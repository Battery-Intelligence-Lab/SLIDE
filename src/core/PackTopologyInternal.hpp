/**
 * @file PackTopologyInternal.hpp
 * @brief Internal cold-path finalisation shared by topology importers and solvers.
 * @surface internal
 */

#pragma once

#include "PackTopology.hpp"

#include <algorithm>
#include <functional>
#include <ranges>

namespace slide::core::detail {

template <class Range, class Projection = std::identity>
  requires std::ranges::random_access_range<const Range>
           && std::ranges::sized_range<const Range>
[[nodiscard]] constexpr bool firstOccurrence(
  const Range &range, std::size_t index, Projection projection = {})
{
  if (index >= std::ranges::size(range))
    return false;
  const auto first = std::ranges::begin(range);
  const auto current =
    first + static_cast<std::ranges::range_difference_t<const Range>>(index);
  return std::ranges::find(
           first, current, std::invoke(projection, *current), projection)
         == current;
}

[[nodiscard]] slide::Status validateElectricalNetlist(
  const CompiledElectricalNetlist &netlist,
  std::size_t cell_count);

/** Finalise metadata and an empty thermal graph for an imported candidate. */
[[nodiscard]] slide::Status finalizeImportedPackTopology(
  CompiledPackTopology &candidate,
  std::uint32_t node_count);

} // namespace slide::core::detail
