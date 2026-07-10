/**
 * @file FuzzSupport.hpp
 * @brief Trap-based deterministic and atomicity oracles for parser fuzzing.
 */

#pragma once

#include "core/Experiment.hpp"
#include "core/NetlistCsv.hpp"
#include "core/ParameterSet.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace slide::fuzz {

[[noreturn]] inline void trap() { __builtin_trap(); }

inline void require(bool condition)
{
  if (!condition)
    trap();
}

inline bool finiteBits(const double &value) noexcept
{
  constexpr std::uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
  // Keep an independent opaque reference boundary and integer barrier. The
  // Release build may otherwise fold std::isfinite under -ffinite-math-only.
  volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  return (bits & exponent_mask) != exponent_mask;
}

struct PoisonControl
{
  std::uint64_t marker{};
  double operator()(const core::ExperimentVariables &) const
  {
    return 783.0;
  }
};

struct PoisonIndicator
{
  std::uint64_t marker{};
  double operator()(const core::ExperimentVariables &) const
  {
    return 784.0;
  }
};

inline bool sameSegment(const core::ExperimentSegment &left,
                        const core::ExperimentSegment &right)
{
  if (left.mode != right.mode || left.direction != right.direction
      || left.value != right.value
      || left.value_is_c_rate != right.value_is_c_rate
      || left.duration != right.duration
      || left.voltage_limit != right.voltage_limit
      || left.current_cutoff != right.current_cutoff
      || left.cutoff_is_c_rate != right.cutoff_is_c_rate
      || left.drive_cycle != right.drive_cycle
      || static_cast<bool>(left.custom_control)
           != static_cast<bool>(right.custom_control)
      || left.custom_terminations.size() != right.custom_terminations.size()
      || left.scheduled_start != right.scheduled_start
      || left.sample_period != right.sample_period || left.source != right.source)
    return false;
  const auto *left_control = left.custom_control.target<PoisonControl>();
  const auto *right_control = right.custom_control.target<PoisonControl>();
  if ((left_control != nullptr || right_control != nullptr)
      && (left_control == nullptr || right_control == nullptr
          || left_control->marker != right_control->marker))
    return false;
  for (std::size_t i = 0; i < left.custom_terminations.size(); ++i)
    if (left.custom_terminations[i].name
          != right.custom_terminations[i].name
        || static_cast<bool>(left.custom_terminations[i].indicator)
             != static_cast<bool>(right.custom_terminations[i].indicator))
      return false;
    else {
      const auto *left_indicator =
        left.custom_terminations[i].indicator.target<PoisonIndicator>();
      const auto *right_indicator =
        right.custom_terminations[i].indicator.target<PoisonIndicator>();
      if ((left_indicator != nullptr || right_indicator != nullptr)
          && (left_indicator == nullptr || right_indicator == nullptr
              || left_indicator->marker != right_indicator->marker))
        return false;
    }
  return true;
}

inline bool sameExperiment(const core::Experiment &left,
                           const core::Experiment &right)
{
  if (left.segments.size() != right.segments.size())
    return false;
  for (std::size_t i = 0; i < left.segments.size(); ++i)
    if (!sameSegment(left.segments[i], right.segments[i]))
      return false;
  return true;
}

inline core::Experiment poisonExperiment()
{
  core::Experiment output;
  output.segments.push_back(
    { .mode = core::ControlMode::custom_differential,
      .direction = core::Direction::charge,
      .value = -777.0,
      .value_is_c_rate = true,
      .duration = -778.0,
      .voltage_limit = -779.0,
      .current_cutoff = -780.0,
      .cutoff_is_c_rate = true,
      .drive_cycle = "fuzz-poison-drive",
      .custom_control = PoisonControl{ .marker = UINT64_C(0xa55aa55aa55aa55a) },
      .custom_terminations = {
        { .name = "fuzz-poison-event",
          .indicator = PoisonIndicator{
            .marker = UINT64_C(0x5aa55aa55aa55aa5) } } },
      .scheduled_start = 781.0,
      .sample_period = 782.0,
      .source = "fuzz-poison-experiment" });
  return output;
}

inline bool validExperiment(const core::Experiment &experiment)
{
  if (experiment.segments.empty() || experiment.segments.size() > 10'000)
    return false;
  std::size_t retained{};
  for (const auto &segment : experiment.segments) {
    if (segment.source.empty() || segment.source == "fuzz-poison-experiment"
        || segment.drive_cycle == "fuzz-poison-drive"
        || segment.source.size() > 65'536 || segment.drive_cycle.size() > 1024
        || segment.source.size() > 4U * 1024U * 1024U - retained
        || segment.drive_cycle.size()
             > 4U * 1024U * 1024U - retained - segment.source.size())
      return false;
    retained += segment.source.size() + segment.drive_cycle.size();
    if (!finiteBits(segment.value) || !finiteBits(segment.duration)
        || !finiteBits(segment.voltage_limit)
        || !finiteBits(segment.current_cutoff)
        || segment.duration < 0.0 || segment.voltage_limit < 0.0
        || segment.current_cutoff < 0.0 || segment.scheduled_start != -1.0
        || segment.sample_period != -1.0 || segment.custom_control
        || !segment.custom_terminations.empty()
        || (segment.current_cutoff == 0.0 && segment.cutoff_is_c_rate))
      return false;
    switch (segment.mode) {
    case core::ControlMode::current: {
      if (!(segment.value > 0.0)
          || !(segment.direction == core::Direction::charge
               || segment.direction == core::Direction::discharge)
          || !segment.drive_cycle.empty()
          || !(segment.duration > 0.0 || segment.voltage_limit > 0.0
               || segment.current_cutoff > 0.0))
        return false;
      break;
    }
    case core::ControlMode::power:
      if (!(segment.value > 0.0) || segment.value_is_c_rate
          || !(segment.direction == core::Direction::charge
               || segment.direction == core::Direction::discharge)
          || !segment.drive_cycle.empty()
          || !(segment.duration > 0.0 || segment.voltage_limit > 0.0
               || segment.current_cutoff > 0.0))
        return false;
      break;
    case core::ControlMode::voltage:
      if (!(segment.value > 0.0) || segment.value_is_c_rate
          || segment.direction != core::Direction::none
          || segment.voltage_limit != 0.0 || !segment.drive_cycle.empty()
          || !(segment.duration > 0.0 || segment.current_cutoff > 0.0))
        return false;
      break;
    case core::ControlMode::rest:
      if (segment.direction != core::Direction::none || segment.value != 0.0
          || segment.value_is_c_rate || !(segment.duration > 0.0)
          || segment.voltage_limit != 0.0 || segment.current_cutoff != 0.0
          || segment.cutoff_is_c_rate || !segment.drive_cycle.empty())
        return false;
      break;
    case core::ControlMode::drive_cycle:
      if (segment.direction != core::Direction::none || segment.value != 0.0
          || segment.value_is_c_rate || segment.duration != 0.0
          || segment.voltage_limit != 0.0 || segment.current_cutoff != 0.0
          || segment.cutoff_is_c_rate || segment.drive_cycle.empty())
        return false;
      break;
    default:
      return false;
    }
  }
  return retained <= 4U * 1024U * 1024U;
}

inline bool sameParameters(const core::ParameterSet &left,
                           const core::ParameterSet &right)
{
  const auto left_values = left.describe();
  const auto right_values = right.describe();
  if (left_values.size() != right_values.size())
    return false;
  for (std::size_t i = 0; i < left_values.size(); ++i)
    if (left_values[i].name != right_values[i].name
        || left_values[i].value != right_values[i].value
        || left_values[i].provenance != right_values[i].provenance)
      return false;
  return true;
}

inline core::ParameterSet poisonParameters()
{
  core::ParameterSet output;
  require(output.set("fuzz-poison-scalar", 777.0, "fuzz-poison-origin")
          == Status::Success);
  core::OCVCurve curve{ .stoichiometry = { 0.0, 0.4, 1.0 },
                        .value = { 2.0, 3.0, 4.0 } };
  require(output.set("fuzz-poison-curve", std::move(curve), "fuzz-poison-origin")
          == Status::Success);
  return output;
}

inline bool hasPoisonParameters(const core::ParameterSet &parameters)
{
  return parameters.contains("fuzz-poison-scalar")
         || parameters.contains("fuzz-poison-curve");
}

inline bool sameBranch(const core::CompiledElectricalBranch &left,
                       const core::CompiledElectricalBranch &right)
{
  return left.node_positive == right.node_positive
         && left.node_negative == right.node_negative && left.kind == right.kind
         && left.cell == right.cell && left.resistance == right.resistance;
}

inline bool sameTopology(const core::CompiledPackTopology &left,
                         const core::CompiledPackTopology &right)
{
  if (left.cells.size() != right.cells.size()
      || left.batch_archetypes != right.batch_archetypes
      || left.electrical.node_count != right.electrical.node_count
      || left.electrical.terminal_positive
           != right.electrical.terminal_positive
      || left.electrical.terminal_negative
           != right.electrical.terminal_negative
      || left.electrical.branches.size() != right.electrical.branches.size()
      || left.electrical.nodal_sparsity != right.electrical.nodal_sparsity
      || left.electrical.ladder_offsets != right.electrical.ladder_offsets
      || left.electrical.ladder_cells != right.electrical.ladder_cells
      || left.electrical.ladder_nodes != right.electrical.ladder_nodes
      || left.electrical.connected != right.electrical.connected
      || left.electrical.index1_candidate
           != right.electrical.index1_candidate
      || left.electrical.series_parallel_ladder
           != right.electrical.series_parallel_ladder
      || left.thermal.cell_count != right.thermal.cell_count
      || left.thermal.boundary_count != right.thermal.boundary_count
      || left.thermal.edges.size() != right.thermal.edges.size()
      || left.thermal.offsets != right.thermal.offsets
      || left.thermal.incidents.size() != right.thermal.incidents.size()
      || left.thermal.edge_flux != right.thermal.edge_flux
      || left.thermal.trial_edge_flux != right.thermal.trial_edge_flux
      || left.thermal.trial_endpoint_heat != right.thermal.trial_endpoint_heat
      || left.thermal.trial_edge_incidence
           != right.thermal.trial_edge_incidence)
    return false;
  for (std::size_t i = 0; i < left.cells.size(); ++i) {
    const auto &a = left.cells[i];
    const auto &b = right.cells[i];
    if (a.path != b.path || a.archetype != b.archetype
        || a.location.batch != b.location.batch
        || a.location.lane != b.location.lane || a.thermal != b.thermal)
      return false;
  }
  for (std::size_t i = 0; i < left.electrical.branches.size(); ++i)
    if (!sameBranch(left.electrical.branches[i],
                    right.electrical.branches[i]))
      return false;
  for (std::size_t i = 0; i < left.thermal.edges.size(); ++i) {
    const auto &a = left.thermal.edges[i];
    const auto &b = right.thermal.edges[i];
    if (a.low != b.low || a.high != b.high
        || a.conductance != b.conductance)
      return false;
  }
  for (std::size_t i = 0; i < left.thermal.incidents.size(); ++i) {
    const auto &a = left.thermal.incidents[i];
    const auto &b = right.thermal.incidents[i];
    if (a.edge != b.edge || a.sign != b.sign)
      return false;
  }
  return true;
}

inline core::CompiledPackTopology poisonTopology()
{
  return {
    .cells = { { .path = "fuzz-poison-cell",
                 .archetype = "fuzz-poison-archetype",
                 .location = { .batch = 3, .lane = 4 },
                 .thermal = true } },
    .batch_archetypes = { "fuzz-poison-batch" },
    .electrical = { .node_count = 9,
                    .terminal_positive = 8,
                    .terminal_negative = 7,
                    .branches = { { .node_positive = 8,
                                    .node_negative = 7,
                                    .kind = core::ElectricalBranchKind::resistor,
                                    .cell = 6,
                                    .resistance = 5.0 } },
                    .nodal_sparsity = { { 7, 8 } },
                    .ladder_offsets = { 2, 3 },
                    .ladder_cells = { 4 },
                    .ladder_nodes = { 8, 7 },
                    .connected = true,
                    .index1_candidate = true,
                    .series_parallel_ladder = true },
    .thermal = { .cell_count = 1,
                 .boundary_count = 2,
                 .edges = { { .low = 0, .high = 1, .conductance = 6.0 } },
                 .offsets = { 0, 1 },
                 .incidents = { { .edge = 0, .sign = 1 } },
                 .edge_flux = { 7.0 },
                 .trial_edge_flux = { 8.0 },
                 .trial_endpoint_heat = { 9.0 },
                 .trial_edge_incidence = { 3 } }
  };
}

inline bool validImportedTopology(const core::CompiledPackTopology &topology)
{
  if (topology.cells.empty() || topology.cells.size() > 100'000
      || topology.batch_archetypes != std::vector<std::string>{ "cell" }
      || topology.thermal.cell_count != topology.cells.size()
      || topology.thermal.boundary_count != 0
      || !topology.thermal.edges.empty() || !topology.thermal.incidents.empty()
      || !topology.thermal.edge_flux.empty()
      || !topology.thermal.trial_edge_flux.empty()
      || !topology.thermal.trial_edge_incidence.empty()
      || topology.thermal.offsets.size() != topology.cells.size() + 1
      || std::any_of(topology.thermal.offsets.begin(),
                     topology.thermal.offsets.end(),
                     [](std::uint32_t value) { return value != 0; })
      || topology.thermal.trial_endpoint_heat.size() != topology.cells.size()
      || std::any_of(topology.thermal.trial_endpoint_heat.begin(),
                     topology.thermal.trial_endpoint_heat.end(),
                     [](double value) { return value != 0.0; }))
    return false;

  std::set<std::string> paths;
  for (std::size_t i = 0; i < topology.cells.size(); ++i) {
    const auto &cell = topology.cells[i];
    if (cell.path.empty() || cell.path == "fuzz-poison-cell"
        || !paths.insert(cell.path).second || cell.archetype != "cell"
        || cell.thermal || cell.location.batch != 0
        || cell.location.lane != i)
      return false;
  }

  const auto &netlist = topology.electrical;
  if (!netlist.connected || !netlist.index1_candidate
      || netlist.node_count < 2 || netlist.node_count > 200'002
      || netlist.terminal_positive >= netlist.node_count
      || netlist.terminal_negative >= netlist.node_count
      || netlist.terminal_positive == netlist.terminal_negative
      || netlist.branches.empty()
      || netlist.branches.size() < netlist.node_count - 1)
    return false;

  const auto unused = std::numeric_limits<std::uint32_t>::max();
  std::vector<unsigned char> seen_cell(topology.cells.size());
  std::vector<std::pair<std::uint32_t, std::uint32_t>> cell_endpoints(
    topology.cells.size(), { unused, unused });
  std::vector<std::vector<std::uint32_t>> adjacency(netlist.node_count);
  std::vector<std::pair<std::uint32_t, std::uint32_t>> expected_sparsity;
  expected_sparsity.reserve(netlist.branches.size() * 3);
  bool all_cells = true;
  for (const auto &branch : netlist.branches) {
    if (branch.node_positive >= netlist.node_count
        || branch.node_negative >= netlist.node_count
        || branch.node_positive == branch.node_negative)
      return false;
    if (branch.kind == core::ElectricalBranchKind::cell) {
      if (branch.cell >= topology.cells.size() || seen_cell[branch.cell] != 0
          || branch.resistance != 0.0)
        return false;
      seen_cell[branch.cell] = 1;
      cell_endpoints[branch.cell] = { branch.node_positive,
                                      branch.node_negative };
    } else if (branch.kind == core::ElectricalBranchKind::resistor) {
      if (branch.cell != 0 || !finiteBits(branch.resistance)
          || !(branch.resistance > 0.0))
        return false;
      all_cells = false;
    } else {
      return false;
    }
    adjacency[branch.node_positive].push_back(branch.node_negative);
    adjacency[branch.node_negative].push_back(branch.node_positive);
    expected_sparsity.emplace_back(branch.node_positive, branch.node_positive);
    expected_sparsity.emplace_back(branch.node_negative, branch.node_negative);
    expected_sparsity.emplace_back(std::min(branch.node_positive,
                                            branch.node_negative),
                                   std::max(branch.node_positive,
                                            branch.node_negative));
  }
  if (std::any_of(seen_cell.begin(), seen_cell.end(), [](unsigned char value) { return value != 1; }))
    return false;

  std::sort(expected_sparsity.begin(), expected_sparsity.end());
  expected_sparsity.erase(
    std::unique(expected_sparsity.begin(), expected_sparsity.end()),
    expected_sparsity.end());
  if (netlist.nodal_sparsity != expected_sparsity)
    return false;

  std::vector<unsigned char> visited(netlist.node_count);
  std::queue<std::uint32_t> pending;
  pending.push(netlist.terminal_positive);
  visited[netlist.terminal_positive] = 1;
  while (!pending.empty()) {
    const auto node = pending.front();
    pending.pop();
    for (const auto next : adjacency[node])
      if (visited[next] == 0) {
        visited[next] = 1;
        pending.push(next);
      }
  }
  if (std::any_of(visited.begin(), visited.end(), [](unsigned char value) { return value == 0; }))
    return false;

  bool ladder_expected = all_cells;
  std::map<std::pair<std::uint32_t, std::uint32_t>,
           std::vector<std::uint32_t>>
    layers;
  std::vector<std::vector<std::uint32_t>> simple_adjacency(
    netlist.node_count);
  if (ladder_expected) {
    for (std::uint32_t cell = 0; cell < cell_endpoints.size(); ++cell) {
      const auto endpoints = std::minmax(cell_endpoints[cell].first,
                                         cell_endpoints[cell].second);
      auto &layer = layers[{ endpoints.first, endpoints.second }];
      if (layer.empty()) {
        simple_adjacency[endpoints.first].push_back(endpoints.second);
        simple_adjacency[endpoints.second].push_back(endpoints.first);
      }
      layer.push_back(cell);
    }
    ladder_expected = layers.size() + 1 == netlist.node_count
                      && simple_adjacency[netlist.terminal_positive].size() == 1
                      && simple_adjacency[netlist.terminal_negative].size() == 1;
    for (std::uint32_t node = 0; ladder_expected && node < netlist.node_count;
         ++node)
      if (node != netlist.terminal_positive
          && node != netlist.terminal_negative
          && simple_adjacency[node].size() != 2)
        ladder_expected = false;
  }

  std::vector<std::uint32_t> expected_offsets, expected_cells, expected_nodes;
  if (ladder_expected) {
    expected_offsets.push_back(0);
    expected_nodes.push_back(netlist.terminal_positive);
    auto previous = unused;
    auto current = netlist.terminal_positive;
    std::vector<unsigned char> visited_path(netlist.node_count);
    while (current != netlist.terminal_negative && ladder_expected) {
      if (visited_path[current] != 0) {
        ladder_expected = false;
        break;
      }
      visited_path[current] = 1;
      const auto &neighbors = simple_adjacency[current];
      const auto next_it = std::find_if(neighbors.begin(), neighbors.end(), [previous](std::uint32_t node) {
        return node != previous;
      });
      if (next_it == neighbors.end()) {
        ladder_expected = false;
        break;
      }
      const auto next = *next_it;
      const auto endpoints = std::minmax(current, next);
      auto layer = layers[{ endpoints.first, endpoints.second }];
      std::sort(layer.begin(), layer.end());
      for (const auto cell : layer) {
        if (cell_endpoints[cell]
            != std::pair<std::uint32_t, std::uint32_t>{ current, next }) {
          ladder_expected = false;
          break;
        }
        expected_cells.push_back(cell);
      }
      if (!ladder_expected)
        break;
      expected_offsets.push_back(
        static_cast<std::uint32_t>(expected_cells.size()));
      expected_nodes.push_back(next);
      previous = current;
      current = next;
    }
    ladder_expected = ladder_expected
                      && expected_nodes.size() == netlist.node_count
                      && expected_cells.size() == topology.cells.size();
  }

  if (netlist.series_parallel_ladder != ladder_expected)
    return false;
  if (ladder_expected) {
    if (netlist.ladder_offsets != expected_offsets
        || netlist.ladder_cells != expected_cells
        || netlist.ladder_nodes != expected_nodes)
      return false;
  } else if (!netlist.ladder_offsets.empty() || !netlist.ladder_cells.empty()
             || !netlist.ladder_nodes.empty()) {
    return false;
  }
  return true;
}

} // namespace slide::fuzz
