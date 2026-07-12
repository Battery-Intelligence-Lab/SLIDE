/**
 * @file ParameterCurve.hpp
 * @brief Exact validation and adaptive sampling for cold parameter curves.
 *
 * M0.7 / 9C-3 cold-path contract: this header owns the one reusable curve
 * validation/sampling implementation shared by built-in absorption and BPX.
 * It preserves the D-16 priority refinement order and 4,096-point ceiling.
 */

#pragma once

#include "../CellDesign.hpp"
#include "../Numeric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <queue>
#include <vector>

namespace slide::core::detail {

inline bool validParameterCurve(const OCVCurve &curve)
{
  if (curve.stoichiometry.size() != curve.value.size()
      || curve.stoichiometry.size() < 2)
    return false;
  for (std::size_t i = 0; i < curve.stoichiometry.size(); ++i)
    if (!is_finite(curve.stoichiometry[i]) || !is_finite(curve.value[i])
        || (i > 0 && curve.stoichiometry[i] <= curve.stoichiometry[i - 1]))
      return false;
  return true;
}

template <class Function>
OCVCurve sampleParameterCurve(Function function)
{
  struct Segment
  {
    real_t left;
    real_t right;
    real_t left_value;
    real_t right_value;
    real_t relative_error;
  };
  const auto make_segment = [&function](real_t left, real_t right, real_t left_value, real_t right_value) {
    real_t error{};
    for (const real_t fraction : {
           real_t{ 0.25 }, real_t{ 0.5 }, real_t{ 0.75 } }) {
      const real_t x = left + fraction * (right - left);
      const real_t exact = function(x);
      const real_t linear = left_value + fraction * (right_value - left_value);
      error = std::max(error,
                       std::abs(linear - exact)
                         / std::max(std::abs(exact), real_t{ 1e-12 }));
    }
    return Segment{ left, right, left_value, right_value, error };
  };
  const auto lower_error = [](const Segment &left, const Segment &right) {
    return left.relative_error < right.relative_error;
  };

  constexpr std::size_t maximum_points = 4096;
  // Quarter-point checks plus this safety factor bound interpolation between
  // the probes without spending the table budget in traversal order.
  constexpr real_t relative_tolerance = 1e-7;

  std::priority_queue<Segment, std::vector<Segment>, decltype(lower_error)>
    work{ lower_error };
  work.push(make_segment(0.0, 1.0, function(0.0), function(1.0)));
  while (work.size() + 1 < maximum_points
         && work.top().relative_error > relative_tolerance) {
    const Segment segment = work.top();
    work.pop();
    const real_t middle = 0.5 * (segment.left + segment.right);
    const real_t middle_value = function(middle);
    work.push(make_segment(segment.left, middle, segment.left_value, middle_value));
    work.push(make_segment(middle, segment.right, middle_value, segment.right_value));
  }

  std::vector<Segment> segments;
  segments.reserve(work.size());
  while (!work.empty()) {
    segments.push_back(work.top());
    work.pop();
  }
  std::sort(segments.begin(), segments.end(), [](const Segment &left, const Segment &right) {
    return left.left < right.left;
  });
  OCVCurve curve;
  curve.stoichiometry.reserve(segments.size() + 1);
  curve.value.reserve(segments.size() + 1);
  for (const auto &segment : segments) {
    curve.stoichiometry.push_back(segment.left);
    curve.value.push_back(segment.left_value);
  }
  curve.stoichiometry.push_back(segments.back().right);
  curve.value.push_back(segments.back().right_value);
  return curve;
}

} // namespace slide::core::detail
