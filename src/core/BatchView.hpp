/**
 * @file BatchView.hpp
 * @brief Rebindable, non-owning views used by v4 RHS kernels (PLAN.md §3.12, D-23).
 *
 * A BatchView carries only the immutable arena geometry and a trial-vector base pointer.
 * Adaptive integrators evaluate the RHS at their own trial vectors, so the pointer is
 * rebound on every evaluation; slices and strides never change. Rebinding allocates
 * nothing and does not copy state.
 */

#pragma once

#include "StateArena.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <span>

namespace slide::core {

struct BatchShape {
  int n_rows{};
  int n_lanes{};
  int stride{};

  constexpr std::size_t storage_size() const
  {
    return static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(stride);
  }

  static BatchShape from(const StateArena &arena)
  {
    return { arena.n_rows(), arena.n_lanes(), arena.stride() };
  }
};

template <class Real>
class BasicBatchView
{
public:
  using value_type = Real;

  explicit BasicBatchView(BatchShape shape) : shape_{ shape }
  {
    assert(shape.n_rows > 0 && shape.n_lanes > 0 && shape.stride >= shape.n_lanes);
  }

  BasicBatchView(BatchShape shape, std::span<Real> values) : BasicBatchView(shape)
  {
    rebind(values);
  }

  void rebind(std::span<Real> values)
  {
    assert(values.size() == shape_.storage_size());
    values_ = values;
  }

  bool bound() const { return values_.data() != nullptr; }

  std::span<Real> raw() const
  {
    assert(bound());
    return values_;
  }

  std::span<Real> row(int row) const
  {
    assert(bound() && 0 <= row && row < shape_.n_rows);
    return { values_.data() + static_cast<std::size_t>(row) * shape_.stride,
             static_cast<std::size_t>(shape_.n_lanes) };
  }

  Real &at(StateSlice slice, int row, int lane) const
  {
    assert(bound() && 0 <= slice.row_begin && 0 <= row && row < slice.rows
           && slice.row_begin + row < shape_.n_rows
           && 0 <= lane && lane < shape_.n_lanes);
    return values_[static_cast<std::size_t>(slice.row_begin + row) * shape_.stride
                   + static_cast<std::size_t>(lane)];
  }

  BatchShape shape() const { return shape_; }
  int n_rows() const { return shape_.n_rows; }
  int n_lanes() const { return shape_.n_lanes; }
  int stride() const { return shape_.stride; }

private:
  BatchShape shape_{};
  std::span<Real> values_{};
};

using BatchView = BasicBatchView<real_t>;
using ConstBatchView = BasicBatchView<const real_t>;

//!< Frozen operating inputs for one RHS evaluation/outer segment. State such as cell
//!< temperature remains in BatchView; externally imposed current density belongs here.
template <class Real = real_t>
struct BasicStepCtx {
  Real time{};
  Real dt{};
  std::span<const Real> i_app{}; //!< per-lane applied current density [A m^-2]

  void assert_valid_for(int n_lanes) const
  {
    assert(static_cast<int>(i_app.size()) == n_lanes);
  }
};

using StepCtx = BasicStepCtx<real_t>;

/**
 * Pair of views rebound by an RHS adapter for every trial evaluation. zero_derivative()
 * is the mandatory first stage because component kernels accumulate into shared rows.
 */
class RhsViews
{
public:
  explicit RhsViews(BatchShape shape) : y{ shape }, ydot{ shape } {}

  void rebind(std::span<const real_t> trial_y, std::span<real_t> trial_ydot)
  {
    y.rebind(trial_y);
    ydot.rebind(trial_ydot);
  }

  void zero_derivative()
  {
    std::fill(ydot.raw().begin(), ydot.raw().end(), real_t{});
  }

  ConstBatchView y;
  BatchView ydot;
};

} // namespace slide::core
