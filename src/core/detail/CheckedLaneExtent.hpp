/**
 * @file CheckedLaneExtent.hpp
 * @brief Cold checked sizing for fixed-field, lane-major scratch storage.
 *
 * Contract (PLAN.md §6 MQ.2/O1): construction-only; validated lane counts
 * enter here, and hot kernels never call this helper.
 * @surface internal
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace slide::core::detail {

inline std::size_t checked_lane_extent(int n_lanes,
                                       std::size_t values_per_lane,
                                       const char *extent_error)
{
  assert(n_lanes > 0 && values_per_lane > 0);
  const auto lanes = static_cast<std::size_t>(n_lanes);
  if (lanes > std::numeric_limits<std::size_t>::max() / values_per_lane)
    throw std::length_error{ extent_error };
  return lanes * values_per_lane;
}

} // namespace slide::core::detail
