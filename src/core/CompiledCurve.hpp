/**
 * @file CompiledCurve.hpp
 * @brief Build-time curve canonicalisation with allocation-free O(1) evaluation (PLAN.md §3.11, D-16).
 */

#pragma once

#include "../types/Status.hpp"
#include "Numeric.hpp"
#include "SpmScalarKernels.hpp"
#include "StateArena.hpp"

#include <algorithm>
#include <cassert>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace slide::core {

namespace detail {

  inline constexpr real_t default_lut_tolerance = 1e-6;

  inline bool curve_is_finite(const real_t &value) noexcept
  {
    constexpr std::uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
    // A by-value floating-point parameter can acquire LLVM `nofpclass` attributes under
    // `-ffinite-math-only`, making a NaN call undefined before its body can classify the bits.
    // Keep this boundary reference-based, then retain the integer barrier.
    volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
    return (bits & exponent_mask) != exponent_mask;
  }

#if defined(_MSC_VER)
  __declspec(noinline)
#elif defined(__GNUC__)
  __attribute__((noinline))
#endif
  inline real_t curve_nan() noexcept
  {
    static_assert(sizeof(real_t) == sizeof(std::uint64_t));
    // Keep the sentinel opaque to the finite-math optimiser. If this were a constexpr NaN,
    // Clang could delete the guarding branch because `-ffinite-math-only` declares a NaN return
    // unreachable, exposing the float-to-index conversion again.
    volatile std::uint64_t bits = UINT64_C(0x7ff8000000000000);
    const std::uint64_t copied_bits = bits;
    return std::bit_cast<real_t>(copied_bits);
  }

} // namespace detail

/**
 * Exact piecewise-linear table with a uniform segment-index accelerator. This preserves measured
 * knots exactly (needed by parity and parameter round trips) while replacing a binary search with
 * fma/floor + one indexed segment lookup. The bin width is no larger than the minimum knot spacing,
 * so evaluation needs at most one segment correction.
 */
class IndexedPiecewiseLinear
{
public:
  [[nodiscard]] slide::Status build(std::span<const real_t> x,
                                    std::span<const real_t> y)
  {
    x_.clear();
    y_.clear();
    segment_.clear();
    inv_bin_width_ = 0.0;
    // Adaptive BPX functions can legitimately contain a 1/65536-wide segment,
    // which needs 65537 endpoint-covering bins. Keep a bounded 512 KiB index
    // budget while admitting that canonical 4096-knot refinement result.
    constexpr std::size_t max_bins = 131'072;
    if (x.size() != y.size() || x.size() < 2 || x.size() > max_bins)
      return slide::Status::Invalid_parameters;

    if (!detail::curve_is_finite(x[0]) || !detail::curve_is_finite(y[0]))
      return slide::Status::Invalid_parameters;

    real_t min_spacing = 0.0;
    for (std::size_t i = 1; i < x.size(); ++i) {
      if (!detail::curve_is_finite(x[i]) || !detail::curve_is_finite(y[i]))
        return slide::Status::Invalid_parameters;
      const real_t spacing = x[i] - x[i - 1];
      const real_t delta_y = y[i] - y[i - 1];
      if (!detail::curve_is_finite(spacing) || !(spacing > 0.0)
          || !detail::curve_is_finite(delta_y))
        return slide::Status::Invalid_parameters;
      const real_t slope = delta_y / spacing;
      if (!detail::curve_is_finite(slope))
        return slide::Status::Invalid_parameters;
      min_spacing = i == 1 ? spacing : std::min(min_spacing, spacing);
    }

    const real_t range = x.back() - x.front();
    if (!detail::curve_is_finite(range) || !(range > 0.0))
      return slide::Status::Invalid_parameters;
    const real_t bin_ratio = range / min_spacing;
    if (!detail::curve_is_finite(bin_ratio) || !(bin_ratio > 0.0)
        || bin_ratio > static_cast<real_t>(max_bins - 1))
      return slide::Status::Invalid_parameters;
    const auto required_bins = static_cast<std::size_t>(std::ceil(bin_ratio)) + 1;

    const std::size_t bins = std::max<std::size_t>(256, required_bins);
    const real_t inv_bin_width = static_cast<real_t>(bins) / range;
    if (!detail::curve_is_finite(inv_bin_width) || !(inv_bin_width > 0.0))
      return slide::Status::Invalid_parameters;

    std::vector<std::uint32_t> candidate_segment(bins);
    std::size_t segment_index = 0;
    for (std::size_t bin = 0; bin < bins; ++bin) {
      const real_t left = x.front() + static_cast<real_t>(bin) / inv_bin_width;
      // source endpoints, their finite positive range, and a bin strictly
      // inside that range make this a finite convex-domain coordinate.
      assert(detail::curve_is_finite(left));
      while (segment_index + 1 < x.size() - 1 && x[segment_index + 1] <= left)
        ++segment_index;
      candidate_segment[bin] = static_cast<std::uint32_t>(segment_index);
    }

    x_.assign(x.begin(), x.end());
    y_.assign(y.begin(), y.end());
    segment_ = std::move(candidate_segment);
    inv_bin_width_ = inv_bin_width;
    return slide::Status::Success;
  }

  template <class Real>
  Real eval(const Real &x) const
  {
    assert(valid());
    const real_t xp = primal_value(x);
    if (!detail::curve_is_finite(xp))
      return static_cast<Real>(detail::curve_nan());
    if (xp <= x_.front())
      return static_cast<Real>(y_.front());
    if (xp >= x_.back())
      return static_cast<Real>(y_.back());

    const auto raw_bin = static_cast<std::size_t>((xp - x_.front()) * inv_bin_width_);
    const auto bin = std::min(raw_bin, segment_.size() - 1);
    std::size_t i = segment_[bin];
    if (xp >= x_[i + 1] && i + 1 < x_.size() - 1)
      ++i;
    else if (xp < x_[i] && i > 0)
      --i;

    return spm_scalar::linearInterpolate(x,
                                         static_cast<Real>(x_[i]),
                                         static_cast<Real>(x_[i + 1]),
                                         static_cast<Real>(y_[i]),
                                         static_cast<Real>(y_[i + 1]));
  }

  /** Piecewise-constant slope using the same O(1) segment lookup as eval(). */
  real_t derivative(const real_t &query) const
  {
    assert(valid());
    if (!detail::curve_is_finite(query))
      return detail::curve_nan();
    real_t x = query;
    if (x <= x_.front())
      x = x_.front();
    else if (x >= x_.back())
      x = x_.back();
    const auto raw_bin = static_cast<std::size_t>((x - x_.front()) * inv_bin_width_);
    const auto bin = std::min(raw_bin, segment_.size() - 1);
    std::size_t i = segment_[bin];
    if (i + 1 < x_.size() - 1 && x >= x_[i + 1])
      ++i;
    else if (x < x_[i] && i > 0)
      --i;
    return (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
  }

  bool valid() const { return x_.size() >= 2 && segment_.size() >= 2; }
  real_t x_min() const { return x_.front(); }
  real_t x_max() const { return x_.back(); }
  std::size_t knots() const { return x_.size(); }

private:
  std::vector<real_t> x_{};
  std::vector<real_t> y_{};
  std::vector<std::uint32_t> segment_{};
  real_t inv_bin_width_{};
};

/**
 * Uniform value LUT for smooth injected functions/curves. build() samples the exact piecewise
 * source, then independently validates source knots and interval midpoints against the requested
 * relative tolerance. Failed accuracy leaves the LUT invalid.
 */
class UniformLut
{
public:
  [[nodiscard]] slide::Status build(std::span<const real_t> x,
                                    std::span<const real_t> y,
                                    const real_t &relative_tolerance = detail::default_lut_tolerance,
                                    std::size_t points = 4096)
  {
    values_.clear();
    x_min_ = 0.0;
    x_max_ = 0.0;
    inv_dx_ = 0.0;
    max_relative_error_ = 0.0;
    if (!detail::curve_is_finite(relative_tolerance) || !(relative_tolerance > 0.0)
        || points < 2 || points > 4096)
      return slide::Status::Invalid_parameters;

    IndexedPiecewiseLinear source;
    const auto status = source.build(x, y);
    if (status != slide::Status::Success)
      return status;

    std::vector<real_t> candidate(points);
    const real_t xmin = x.front();
    const real_t xmax = x.back();
    const real_t range = xmax - xmin;
    // IndexedPiecewiseLinear accepted this exact source range above.
    assert(detail::curve_is_finite(range) && range > 0.0);
    const real_t dx = range / static_cast<real_t>(points - 1);
    // source.build() also proved 256/range representable; division by at
    // most 4095 therefore cannot underflow this spacing to zero.
    assert(detail::curve_is_finite(dx) && dx > 0.0);
    const real_t inv_dx = real_t{ 1 } / dx;
    if (!detail::curve_is_finite(inv_dx) || !(inv_dx > 0.0))
      return slide::Status::Invalid_parameters;
    for (std::size_t i = 0; i < points; ++i) {
      const real_t query = xmin + static_cast<real_t>(i) * dx;
      // Floating addition can overflow even when both endpoints and their
      // subtraction are finite (rounding need not preserve the real interval).
      if (!detail::curve_is_finite(query))
        return slide::Status::Numerical_failure;
      candidate[i] = source.eval(query);
      if (!detail::curve_is_finite(candidate[i]))
        return slide::Status::Numerical_failure;
    }

    auto candidate_eval = [&](real_t query) {
      if (!detail::curve_is_finite(query)) return detail::curve_nan();
      if (query <= xmin) return candidate.front();
      if (query >= xmax) return candidate.back();
      const real_t scaled = (query - xmin) / dx;
      if (!detail::curve_is_finite(scaled) || !(scaled >= 0.0))
        return detail::curve_nan();
      const auto i = std::min(static_cast<std::size_t>(scaled), points - 2);
      return candidate[i] + (candidate[i + 1] - candidate[i]) * (query - (xmin + static_cast<real_t>(i) * dx)) / dx;
    };

    real_t max_relative_error = 0.0;
    auto observe_error = [&](real_t query) {
      const real_t expected = source.eval(query);
      const real_t actual = candidate_eval(query);
      const real_t difference = actual - expected;
      const real_t denominator = std::max(std::abs(expected), real_t{ 1e-12 });
      const real_t relative_error = std::abs(difference) / denominator;
      if (!detail::curve_is_finite(expected) || !detail::curve_is_finite(actual)
          || !detail::curve_is_finite(difference)
          || !detail::curve_is_finite(denominator)
          || !detail::curve_is_finite(relative_error))
        return false;
      max_relative_error = std::max(max_relative_error, relative_error);
      return true;
    };
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (!observe_error(x[i]))
        return slide::Status::Numerical_failure;
      if (i + 1 < x.size()) {
        const real_t midpoint = x[i] + real_t{ 0.5 } * (x[i + 1] - x[i]);
        if (!observe_error(midpoint))
          return slide::Status::Numerical_failure;
      }
    }
    if (max_relative_error > relative_tolerance)
      return slide::Status::Numerical_failure;

    values_ = std::move(candidate);
    x_min_ = xmin;
    x_max_ = xmax;
    inv_dx_ = inv_dx;
    max_relative_error_ = max_relative_error;
    return slide::Status::Success;
  }

  template <class Real>
  Real eval(const Real &x) const
  {
    assert(valid());
    const real_t xp = primal_value(x);
    if (!detail::curve_is_finite(xp))
      return static_cast<Real>(detail::curve_nan());
    if (xp <= x_min_) return static_cast<Real>(values_.front());
    if (xp >= x_max_) return static_cast<Real>(values_.back());
    const real_t scaled = (xp - x_min_) * inv_dx_;
    if (!detail::curve_is_finite(scaled) || !(scaled >= 0.0))
      return static_cast<Real>(detail::curve_nan());
    const auto i = std::min(static_cast<std::size_t>(scaled), values_.size() - 2);
    const Real fraction = (x - static_cast<Real>(x_min_)) * static_cast<Real>(inv_dx_)
                          - static_cast<Real>(i);
    return static_cast<Real>(values_[i])
           + static_cast<Real>(values_[i + 1] - values_[i]) * fraction;
  }

  bool valid() const { return values_.size() >= 2; }
  real_t max_relative_error() const { return max_relative_error_; }

private:
  std::vector<real_t> values_{};
  real_t x_min_{};
  real_t x_max_{};
  real_t inv_dx_{};
  real_t max_relative_error_{};
};

} // namespace slide::core
