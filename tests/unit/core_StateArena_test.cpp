/**
 * @file core_StateArena_test.cpp
 * @brief Phase-1 tests for slide::core StateArena/BatchBuilder (PLAN.md §3.1, D-01).
 *
 * Registered bands (written BEFORE first run, CLAUDE.md §3):
 *  - Layout: declared offsets exactly {zp:0, zn:5, T:10, SOC:11, q_ext:12}, rows()==13.
 *  - Geometry: n_lanes=10 -> stride==16 (pad to 64 B / 8 doubles); arena base and every
 *    row pointer 64-byte aligned.
 *  - Snapshot/restore: bit-identical round-trip (operator== on every element).
 *  - Hot path (P1-G2 groundwork): row()/at()/snapshot()/restore() perform EXACTLY 0
 *    heap allocations, counted by the replaced global operator new (all forms).
 *
 * @date 2026-07-07
 */

#include "../../src/core/BatchBuilder.hpp"
#include "../../src/core/BatchView.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// Allocation-counting global new (P1-G2 groundwork). Counts every allocation in
// the binary; tests diff the counter around the measured region only.
// ---------------------------------------------------------------------------
static std::atomic<std::size_t> g_alloc_count{ 0 };

#if defined(_WIN32)
  #include <malloc.h>
static void *aligned_alloc_impl(std::size_t n, std::size_t al) { return _aligned_malloc(n, al); }
static void aligned_free_impl(void *p) { _aligned_free(p); }
#else
static void *aligned_alloc_impl(std::size_t n, std::size_t al)
{
  return std::aligned_alloc(al, ((n + al - 1) / al) * al);
}
static void aligned_free_impl(void *p) { std::free(p); }
#endif

void *operator new(std::size_t n)
{
  ++g_alloc_count;
  if (void *p = std::malloc(n ? n : 1)) return p;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t n)
{
  ++g_alloc_count;
  if (void *p = std::malloc(n ? n : 1)) return p;
  throw std::bad_alloc{};
}
void *operator new(std::size_t n, std::align_val_t al)
{
  ++g_alloc_count;
  if (void *p = aligned_alloc_impl(n ? n : 1, static_cast<std::size_t>(al))) return p;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t n, std::align_val_t al)
{
  ++g_alloc_count;
  if (void *p = aligned_alloc_impl(n ? n : 1, static_cast<std::size_t>(al))) return p;
  throw std::bad_alloc{};
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }
void operator delete(void *p, std::align_val_t) noexcept { aligned_free_impl(p); }
void operator delete(void *p, std::size_t, std::align_val_t) noexcept { aligned_free_impl(p); }
void operator delete[](void *p, std::align_val_t) noexcept { aligned_free_impl(p); }
void operator delete[](void *p, std::size_t, std::align_val_t) noexcept { aligned_free_impl(p); }

// ---------------------------------------------------------------------------

using namespace slide::core;

namespace {

//!< SPM-shaped archetype used by all tests: 5-mode diffusion per electrode + thermal +
//!< SOC + the Q9 thermal-flux seam. Registered offsets: zp:0 zn:5 T:10 SOC:11 q_ext:12.
struct Layout {
  BatchBuilder b{};
  StateSlice zp, zn, T, SOC, q_ext;
  Layout()
  {
    zp = b.declare({ "zp", 5, Unit::none });
    zn = b.declare({ "zn", 5, Unit::none });
    T = b.declare({ "T", 1, Unit::K });
    SOC = b.declare({ "SOC", 1, Unit::none });
    q_ext = b.reserve_thermal_flux();
  }
};

} // namespace

TEST_CASE("BatchBuilder layout: appended rows, exact offsets, name lookup", "[core][P1]")
{
  Layout L;

  // Registered band: exact integer offsets.
  REQUIRE(L.zp.row_begin == 0);
  REQUIRE(L.zp.rows == 5);
  REQUIRE(L.zn.row_begin == 5);
  REQUIRE(L.T.row_begin == 10);
  REQUIRE(L.SOC.row_begin == 11);
  REQUIRE(L.q_ext.row_begin == 12);
  REQUIRE(L.b.rows() == 13);

  // Name lookup (cold path, recording).
  auto zn = L.b.find("zn");
  REQUIRE(zn.has_value());
  REQUIRE(zn->row_begin == 5);
  REQUIRE(zn->rows == 5);
  REQUIRE_FALSE(L.b.find("does_not_exist").has_value());
}

TEST_CASE("Q9 seam: reserve_thermal_flux is idempotent", "[core][P1]")
{
  Layout L;
  const auto again = L.b.reserve_thermal_flux();
  REQUIRE(again.row_begin == L.q_ext.row_begin);
  REQUIRE(again.rows == 1);
  REQUIRE(L.b.rows() == 13); // no second row appended
}

TEST_CASE("StateArena geometry: stride padding and 64-byte alignment", "[core][P1]")
{
  Layout L;
  auto arena = L.b.build(10);

  // Registered band: 10 lanes pad to 16 (64 B / sizeof(double) = 8 lanes per block).
  REQUIRE(arena.n_rows() == 13);
  REQUIRE(arena.n_lanes() == 10);
  REQUIRE(arena.stride() == 16);
  REQUIRE(arena.size() == 13u * 16u);

  REQUIRE(reinterpret_cast<std::uintptr_t>(arena.raw().data()) % StateArena::alignment == 0);
  for (int r = 0; r < arena.n_rows(); r++) {
    REQUIRE(reinterpret_cast<std::uintptr_t>(arena.row(r).data()) % StateArena::alignment == 0);
    REQUIRE(arena.row(r).size() == 10u);
    REQUIRE(arena.row_padded(r).size() == 16u);
  }
}

TEST_CASE("StateArena zero-init and lane addressing", "[core][P1]")
{
  Layout L;
  auto arena = L.b.build(7);

  for (int r = 0; r < arena.n_rows(); r++)
    for (const auto x : arena.row_padded(r))
      REQUIRE(x == 0.0);

  // at(slice, r, c) addresses arena[(row_begin + r)*stride + c].
  arena.at(L.zn, 2, 3) = 42.5;
  REQUIRE(arena.row(L.zn.row_begin + 2)[3] == 42.5);
  REQUIRE(arena.raw()[static_cast<std::size_t>(7) * arena.stride() + 3] == 42.5);
}

TEST_CASE("Snapshot/restore: bit-identical round-trip via memcpy (PC-7)", "[core][P1]")
{
  Layout L;
  auto arena = L.b.build(33); // deliberately not a multiple of 8

  // Deterministic non-degenerate fill (no symmetric/uniform oracle, CLAUDE.md §4).
  for (int r = 0; r < arena.n_rows(); r++) {
    auto row = arena.row(r);
    for (int c = 0; c < static_cast<int>(row.size()); c++)
      row[c] = 1.0 + 0.001 * r + 1e-7 * c + ((r * 31 + c * 7) % 13) * 0.01;
  }

  // Snapshot the diffusion block only; mutate everything; restore.
  std::vector<real_t> buf(arena.slice_size(L.zp) + arena.slice_size(L.zn));
  const StateSlice diffusion{ L.zp.row_begin, L.zp.rows + L.zn.rows }; // contiguous by layout
  arena.snapshot(diffusion, buf.data());

  std::vector<real_t> T_before(arena.row(L.T.row_begin).begin(), arena.row(L.T.row_begin).end());
  for (auto &x : arena.raw()) x = -1.0;
  arena.restore(diffusion, buf.data());

  // Registered band: bit-identical (operator==), every diffusion element.
  for (int r = 0; r < diffusion.rows; r++) {
    auto row = arena.row(diffusion.row_begin + r);
    for (int c = 0; c < static_cast<int>(row.size()); c++)
      REQUIRE(row[c] == 1.0 + 0.001 * r + 1e-7 * c + ((r * 31 + c * 7) % 13) * 0.01);
  }
  // Rows OUTSIDE the slice must NOT be restored (T stays mutated).
  REQUIRE(arena.at(L.T, 0, 0) == -1.0);
  REQUIRE(T_before[0] != -1.0);
}

TEST_CASE("Hot path performs zero heap allocations (P1-G2 groundwork)", "[core][P1]")
{
  Layout L;
  auto arena = L.b.build(10'000); // PLAN.md P1-G2 scale: one batch of 1e4 cells
  std::vector<real_t> buf(arena.slice_size(L.zp)); // scratch preallocated (PC-1)

  const auto n0 = g_alloc_count.load();

  // Representative step work: SIMD-style row sweeps, element access, checkpoint cycle.
  for (int r = 0; r < arena.n_rows(); r++)
    for (auto &x : arena.row_padded(r))
      x += 1.0;
  arena.at(L.T, 0, 4242) = 298.15;
  arena.snapshot(L.zp, buf.data());
  arena.restore(L.zp, buf.data());

  const auto n1 = g_alloc_count.load();

  // Registered band: EXACTLY zero.
  REQUIRE(n1 - n0 == 0);
}

TEST_CASE("Builder may build several arenas of one archetype (layout shared)", "[core][P1]")
{
  Layout L;
  auto a1 = L.b.build(4);
  auto a2 = L.b.build(4096);

  REQUIRE(a1.n_rows() == a2.n_rows());
  REQUIRE(a1.stride() == 8);      // 4 -> pad to 8
  REQUIRE(a2.stride() == 4096);   // already a multiple of 8
  REQUIRE(a1.n_lanes() == 4);
  REQUIRE(a2.n_lanes() == 4096);
}

TEST_CASE("BatchView rebinds to integrator trial vectors without copying", "[core][P1][rhs]")
{
  Layout L;
  auto arena = L.b.build(7);
  const auto shape = BatchShape::from(arena);
  std::vector<real_t> trial_a(shape.storage_size(), 1.0);
  std::vector<real_t> trial_b(shape.storage_size(), 2.0);
  std::vector<real_t> derivative(shape.storage_size(), 9.0);

  RhsViews views{ shape };
  views.rebind(trial_a, derivative);
  REQUIRE(views.y.at(L.T, 0, 3) == 1.0);

  views.rebind(trial_b, derivative);
  REQUIRE(views.y.at(L.T, 0, 3) == 2.0);

  views.zero_derivative();
  for (const auto value : derivative)
    REQUIRE(value == 0.0);
}

TEST_CASE("BatchBuilder records ODE-row roles for generic integrators", "[core][P1][rhs]")
{
  BatchBuilder builder;
  const auto z = builder.declare({ "z", 3, Unit::none, StateRole::ode });
  const auto current = builder.declare({ "I", 1, Unit::A, StateRole::algebraic });
  const auto throughput = builder.declare({ "Ah", 1, Unit::Ah, StateRole::cumulative });
  const auto q_ext = builder.reserve_thermal_flux();

  REQUIRE(z.row_begin == 0);
  for (int row = 0; row < z.rows; ++row)
    REQUIRE(builder.is_ode_row(z.row_begin + row));
  REQUIRE_FALSE(builder.is_ode_row(current.row_begin));
  REQUIRE_FALSE(builder.is_ode_row(throughput.row_begin));
  REQUIRE_FALSE(builder.is_ode_row(q_ext.row_begin));
  REQUIRE(builder.roles().size() == 6);
}

TEST_CASE("Moving StateArena leaves a safe empty source", "[core][P1]")
{
  Layout L;
  auto source = L.b.build(4);
  source.at(L.T, 0, 1) = 301.0;

  StateArena destination{ std::move(source) };
  REQUIRE(destination.at(L.T, 0, 1) == 301.0);
  REQUIRE(source.n_rows() == 0);
  REQUIRE(source.n_lanes() == 0);
  REQUIRE(source.stride() == 0);
  REQUIRE(source.raw().empty());
}
