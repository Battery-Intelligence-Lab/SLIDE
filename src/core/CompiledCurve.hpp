/**
 * @file CompiledCurve.hpp
 * @brief Build-time curve canonicalisation with allocation-free O(1) evaluation (PLAN.md §3.11, D-16).
 */

#pragma once

#include "../types/Status.hpp"
#include "Numeric.hpp"
#include "StateArena.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace slide::core {

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
    if (x.size() != y.size() || x.size() < 2)
      return slide::Status::Invalid_parameters;

    real_t min_spacing = x[1] - x[0];
    if (!is_finite(x[0]) || !is_finite(y[0]) || min_spacing <= 0.0)
      return slide::Status::Invalid_parameters;
    for (std::size_t i = 1; i < x.size(); ++i) {
      if (!is_finite(x[i]) || !is_finite(y[i]) || x[i] <= x[i - 1])
        return slide::Status::Invalid_parameters;
      min_spacing = std::min(min_spacing, x[i] - x[i - 1]);
    }

    const real_t range = x.back() - x.front();
    const auto required_bins = static_cast<std::size_t>(std::ceil(range / min_spacing)) + 1;
    constexpr std::size_t max_bins = 65'536;
    if (required_bins > max_bins)
      return slide::Status::Invalid_parameters;

    x_.assign(x.begin(), x.end());
    y_.assign(y.begin(), y.end());
    const std::size_t bins = std::max<std::size_t>(256, required_bins);
    segment_.resize(bins);
    inv_bin_width_ = static_cast<real_t>(bins) / range;

    std::size_t segment = 0;
    for (std::size_t bin = 0; bin < bins; ++bin) {
      const real_t left = x_.front() + static_cast<real_t>(bin) / inv_bin_width_;
      while (segment + 1 < x_.size() - 1 && x_[segment + 1] <= left)
        ++segment;
      segment_[bin] = static_cast<std::uint32_t>(segment);
    }
    return slide::Status::Success;
  }

  template <class Real>
  Real eval(const Real &x) const
  {
    assert(valid());
    const real_t xp = primal_value(x);
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

    return static_cast<Real>(y_[i])
           + static_cast<Real>(y_[i + 1] - y_[i])
               * (x - static_cast<Real>(x_[i]))
               / static_cast<Real>(x_[i + 1] - x_[i]);
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
                                    real_t relative_tolerance = 1e-6,
                                    std::size_t points = 4096)
  {
    values_.clear();
    x_min_ = 0.0;
    x_max_ = 0.0;
    inv_dx_ = 0.0;
    max_relative_error_ = 0.0;
    if (!(relative_tolerance > 0.0) || points < 2 || points > 4096)
      return slide::Status::Invalid_parameters;

    IndexedPiecewiseLinear source;
    const auto status = source.build(x, y);
    if (status != slide::Status::Success)
      return status;

    std::vector<real_t> candidate(points);
    const real_t xmin = x.front();
    const real_t xmax = x.back();
    const real_t dx = (xmax - xmin) / static_cast<real_t>(points - 1);
    for (std::size_t i = 0; i < points; ++i)
      candidate[i] = source.eval(xmin + static_cast<real_t>(i) * dx);

    auto candidate_eval = [&](real_t query) {
      if (query <= xmin) return candidate.front();
      if (query >= xmax) return candidate.back();
      const real_t scaled = (query - xmin) / dx;
      const auto i = std::min(static_cast<std::size_t>(scaled), points - 2);
      return candidate[i] + (candidate[i + 1] - candidate[i]) * (query - (xmin + static_cast<real_t>(i) * dx)) / dx;
    };

    real_t max_relative_error = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      const real_t expected = source.eval(x[i]);
      max_relative_error = std::max(max_relative_error,
                                    std::abs(candidate_eval(x[i]) - expected)
                                      / std::max(std::abs(expected), real_t{ 1e-12 }));
      if (i + 1 < x.size()) {
        const real_t midpoint = 0.5 * (x[i] + x[i + 1]);
        const real_t mid_expected = source.eval(midpoint);
        max_relative_error = std::max(max_relative_error,
                                      std::abs(candidate_eval(midpoint) - mid_expected)
                                        / std::max(std::abs(mid_expected), real_t{ 1e-12 }));
      }
    }
    if (max_relative_error > relative_tolerance)
      return slide::Status::Numerical_failure;

    values_ = std::move(candidate);
    x_min_ = xmin;
    x_max_ = xmax;
    inv_dx_ = real_t{ 1 } / dx;
    max_relative_error_ = max_relative_error;
    return slide::Status::Success;
  }

  template <class Real>
  Real eval(const Real &x) const
  {
    assert(valid());
    const real_t xp = primal_value(x);
    if (xp <= x_min_) return static_cast<Real>(values_.front());
    if (xp >= x_max_) return static_cast<Real>(values_.back());
    const real_t scaled = (xp - x_min_) * inv_dx_;
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
