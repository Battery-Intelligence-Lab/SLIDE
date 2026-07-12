/**
 * @file StateArena.hpp
 * @brief v4 core state storage: SoA arena with variable-major layout (PLAN.md §3.1, D-01).
 *
 * A CellBatch holds N cells ("lanes") of identical model composition. All state lives in
 * one contiguous, 64-byte-aligned arena, variable-major: the value of state-row r for
 * lane c is data[r * stride + c]. One row = one variable across all lanes (SIMD sweep).
 *
 * Performance contract (PLAN.md §1.1): PC-1 zero per-step allocations (arena allocates
 * exactly once, at build), PC-3 contiguous/64B-aligned/SoA, PC-7 rollback = memcpy.
 *
 * @date 2026-07-07
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace slide::core {

//!< Q1 (PLAN.md §7, DECIDED 2026-07-07): f64 everywhere in v4.0. The alias is the
//!< insurance: f32-storage/f64-accumulate later is a new arena instantiation, not a rewrite.
using real_t = double;

//!< Unit tag carried by every declared state row (CLAUDE.md: dimensional consistency).
//!< Cold-path metadata only — never touched by kernels.
enum class Unit : int {
  none,   //!< dimensionless (lithium fractions, transformed z-modes, SOC)
  A,      //!< Ampere
  V,      //!< Volt
  K,      //!< Kelvin
  s,      //!< second
  Ah,     //!< Ampere-hour (cumulative charge throughput)
  Wh,     //!< Watt-hour (cumulative energy throughput)
  C,      //!< Coulomb / Ampere-second (lost lithium inventory)
  W,      //!< Watt (heat flows)
  J,      //!< Joule (integrated heat)
  Pa,     //!< Pascal (stress-history states used by ageing models)
  m,      //!< metre (layer thicknesses, SEI thickness)
  m2,     //!< square metre (crack surface)
  m2_s,   //!< square metre per second (solid diffusion coefficient)
  inv_m,  //!< inverse metre (specific interfacial area, m2/m3)
  ohm_m2, //!< area-specific resistance
  mol_m3  //!< mol/m^3 (concentrations)
};

//!< Integration role for a state row (PLAN.md §3.12, D-23). The role is cold-path
//!< layout metadata used by steppers to keep algebraic, cumulative and externally
//!< supplied rows out of generic ODE integrators.
enum class StateRole : unsigned char {
  ode,
  algebraic,
  cumulative,
  input
};

//!< A model component's declaration of one state variable (PLAN.md §3.1).
//!< `name` must outlive the builder that receives it — use string literals.
struct StateSpec
{
  std::string_view name{}; //!< unique within a batch, e.g. "zp", "T", "delta"
  int rows{ 1 };           //!< e.g. "zn" has nch rows
  Unit unit{ Unit::none };
  StateRole role{ StateRole::ode };
};

//!< Integer handle a model component holds after layout (PLAN.md §3.1).
struct StateSlice
{
  int row_begin{ 0 };
  int rows{ 0 };
};

/**
 * Owns the doubles. Allocated exactly once (by BatchBuilder::build); geometry is
 * immutable afterwards. snapshot/restore of a row range is a single memcpy because
 * consecutive rows are contiguous in variable-major layout.
 */
class StateArena
{
  struct CheckedShape
  {
    int n_rows{};
    int n_lanes{};
    int stride{};
    std::size_t elements{};
  };

public:
  static constexpr std::size_t alignment = 64; //!< cache line / AVX-512 vector
  static constexpr int lanes_per_block = static_cast<int>(alignment / sizeof(real_t));

  StateArena() = default;

  StateArena(const StateArena &) = delete;
  StateArena &operator=(const StateArena &) = delete;

  StateArena(StateArena &&other) noexcept
    : n_rows_{ std::exchange(other.n_rows_, 0) },
      n_lanes_{ std::exchange(other.n_lanes_, 0) },
      stride_{ std::exchange(other.stride_, 0) },
      data_{ std::move(other.data_) }
  {}

  StateArena &operator=(StateArena &&other) noexcept
  {
    if (this != &other) {
      n_rows_ = std::exchange(other.n_rows_, 0);
      n_lanes_ = std::exchange(other.n_lanes_, 0);
      stride_ = std::exchange(other.stride_, 0);
      data_ = std::move(other.data_);
    }
    return *this;
  }

  StateArena(int n_rows, int n_lanes)
    : StateArena{ checked_shape(n_rows, n_lanes) }
  {}

  //!< One variable across all lanes. Excludes the alignment padding.
  std::span<real_t> row(int r)
  {
    assert(data_ && 0 <= r && r < n_rows_);
    return { data_.get() + static_cast<std::size_t>(r) * stride_, static_cast<std::size_t>(n_lanes_) };
  }
  std::span<const real_t> row(int r) const
  {
    assert(data_ && 0 <= r && r < n_rows_);
    return { data_.get() + static_cast<std::size_t>(r) * stride_, static_cast<std::size_t>(n_lanes_) };
  }

  //!< Full padded row (stride_ elements); padding lanes stay zero. For SIMD sweeps
  //!< that prefer whole vectors — kernels must not depend on padding values.
  std::span<real_t> row_padded(int r)
  {
    assert(data_ && 0 <= r && r < n_rows_);
    return { data_.get() + static_cast<std::size_t>(r) * stride_, static_cast<std::size_t>(stride_) };
  }

  //!< Element access for cold paths and tests: row `r` (relative to the slice) of lane `c`.
  real_t &at(StateSlice s, int r, int c)
  {
    assert(data_ && 0 <= s.row_begin && 0 <= r && r < s.rows
           && s.row_begin + r < n_rows_ && 0 <= c && c < n_lanes_);
    return data_[static_cast<std::size_t>(s.row_begin + r) * stride_ + static_cast<std::size_t>(c)];
  }
  const real_t &at(StateSlice s, int r, int c) const
  {
    assert(data_ && 0 <= s.row_begin && 0 <= r && r < s.rows
           && s.row_begin + r < n_rows_ && 0 <= c && c < n_lanes_);
    return data_[static_cast<std::size_t>(s.row_begin + r) * stride_ + static_cast<std::size_t>(c)];
  }

  //!< Number of real_t a snapshot buffer for `s` must hold (includes padding lanes).
  std::size_t slice_size(StateSlice s) const { return static_cast<std::size_t>(s.rows) * stride_; }

  //!< PC-7: checkpoint = one memcpy of a contiguous row range (never a tree traversal).
  void snapshot(StateSlice s, real_t *dst) const
  {
    assert(0 <= s.row_begin && s.row_begin + s.rows <= n_rows_);
    std::memcpy(dst, data_.get() + static_cast<std::size_t>(s.row_begin) * stride_, slice_size(s) * sizeof(real_t));
  }
  void restore(StateSlice s, const real_t *src)
  {
    assert(0 <= s.row_begin && s.row_begin + s.rows <= n_rows_);
    std::memcpy(data_.get() + static_cast<std::size_t>(s.row_begin) * stride_, src, slice_size(s) * sizeof(real_t));
  }

  //!< Whole-arena view (recording snapshots, whole-batch rollback).
  std::span<real_t> raw() { return { data_.get(), size() }; }
  std::span<const real_t> raw() const { return { data_.get(), size() }; }

  int n_rows() const { return n_rows_; }
  int n_lanes() const { return n_lanes_; }
  int stride() const { return stride_; }
  std::size_t size() const { return static_cast<std::size_t>(n_rows_) * static_cast<std::size_t>(stride_); }

private:
  explicit StateArena(CheckedShape shape)
    : n_rows_{ shape.n_rows }, n_lanes_{ shape.n_lanes }, stride_{ shape.stride },
      data_{ allocate(shape.elements) }
  {
    std::fill_n(data_.get(), shape.elements, real_t{}); // zero-init incl. padding lanes
  }

  static CheckedShape checked_shape(int n_rows, int n_lanes)
  {
    if (n_rows <= 0 || n_lanes <= 0)
      throw std::invalid_argument{ "StateArena requires positive row and lane counts" };

    constexpr auto maximum = std::numeric_limits<std::size_t>::max();
    const auto lanes = static_cast<std::size_t>(n_lanes);
    const auto block = static_cast<std::size_t>(lanes_per_block);
    const auto stride = ((lanes + block - 1) / block) * block;
    if (stride > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::length_error{ "StateArena padded lane count is not representable" };
    const auto rows = static_cast<std::size_t>(n_rows);
    if (rows > maximum / stride)
      throw std::length_error{ "StateArena element count is not representable" };
    const auto elements = rows * stride;
    if (elements > maximum / sizeof(real_t))
      throw std::length_error{ "StateArena byte count is not representable" };
    return { n_rows, n_lanes, static_cast<int>(stride), elements };
  }

  struct AlignedDelete
  {
    void operator()(real_t *p) const { ::operator delete[](p, std::align_val_t{ alignment }); }
  };

  static std::unique_ptr<real_t[], AlignedDelete> allocate(std::size_t n)
  {
    return std::unique_ptr<real_t[], AlignedDelete>{
      static_cast<real_t *>(::operator new[](n * sizeof(real_t), std::align_val_t{ alignment }))
    };
  }

  int n_rows_{ 0 }, n_lanes_{ 0 }, stride_{ 0 };
  std::unique_ptr<real_t[], AlignedDelete> data_{};
};

} // namespace slide::core
