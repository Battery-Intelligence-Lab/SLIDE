/**
 * @file PackTopologyInternal.hpp
 * @brief Internal cold-path finalisation shared by topology importers and solvers.
 * @surface internal
 */

#pragma once

#include "PackTopology.hpp"

namespace slide::core::detail {

[[nodiscard]] slide::Status validateElectricalNetlist(
  const CompiledElectricalNetlist &netlist,
  std::size_t cell_count);

/** Finalise metadata and an empty thermal graph for an imported candidate. */
[[nodiscard]] slide::Status finalizeImportedPackTopology(
  CompiledPackTopology &candidate,
  std::uint32_t node_count);

} // namespace slide::core::detail
