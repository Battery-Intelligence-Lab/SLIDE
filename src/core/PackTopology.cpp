/**
 * @file PackTopology.cpp
 * @brief Pack description flattening and D-21 thermal graph compilation.
 */

#include "PackTopology.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <tuple>
#include <unordered_map>

namespace slide::core {

PackNode cell(PackCellSpec spec)
{
  return { .kind = PackNodeKind::cell, .cell = std::move(spec) };
}

namespace {

  PackNode group(PackNodeKind kind, std::vector<PackNode> children, PackLink link)
  {
    return { .kind = kind, .link = link, .children = std::move(children) };
  }

} // namespace

PackNode series(int count, const PackNode &child, PackLink link)
{
  return count > 0 ? group(PackNodeKind::series,
                           std::vector<PackNode>(static_cast<std::size_t>(count), child),
                           link)
                   : group(PackNodeKind::series, {}, link);
}

PackNode parallel(int count, const PackNode &child, PackLink link)
{
  return count > 0 ? group(PackNodeKind::parallel,
                           std::vector<PackNode>(static_cast<std::size_t>(count), child),
                           link)
                   : group(PackNodeKind::parallel, {}, link);
}

PackNode series(std::vector<PackNode> children, PackLink link)
{
  return group(PackNodeKind::series, std::move(children), link);
}

PackNode parallel(std::vector<PackNode> children, PackLink link)
{
  return group(PackNodeKind::parallel, std::move(children), link);
}

namespace {

  struct CompileContext
  {
    CompiledPackTopology result{};
    std::uint32_t next_node{ 2 };
    bool valid{ true };

    std::uint32_t allocateNode() { return next_node++; }

    void resistor(std::uint32_t positive, std::uint32_t negative, real_t resistance)
    {
      if (!(is_finite(resistance) && resistance >= 0.0)) {
        valid = false;
        return;
      }
      if (resistance == 0.0)
        return;
      result.electrical.branches.push_back({ .node_positive = positive,
                                             .node_negative = negative,
                                             .kind = ElectricalBranchKind::resistor,
                                             .resistance = resistance });
    }

    void compile(const PackNode &node, std::uint32_t positive, std::uint32_t negative,
                 const std::string &path)
    {
      if (!valid || positive == negative) {
        valid = false;
        return;
      }
      if (node.kind != PackNodeKind::cell && node.kind != PackNodeKind::series
          && node.kind != PackNodeKind::parallel) {
        valid = false;
        return;
      }
      if (node.kind == PackNodeKind::cell) {
        if (node.cell.archetype.empty() || !node.children.empty()) {
          valid = false;
          return;
        }
        const auto index = static_cast<std::uint32_t>(result.cells.size());
        result.cells.push_back({ .path = path.empty() ? "c00" : path,
                                 .archetype = node.cell.archetype,
                                 .thermal = node.cell.thermal });
        result.electrical.branches.push_back({ .node_positive = positive,
                                               .node_negative = negative,
                                               .kind = ElectricalBranchKind::cell,
                                               .cell = index });
        return;
      }
      if (node.children.empty() || !(is_finite(node.link.resistance) && node.link.resistance >= 0.0)) {
        valid = false;
        return;
      }

      const char prefix = node.kind == PackNodeKind::series ? 's' : 'p';
      auto childPath = [&](std::size_t index) {
        std::string token{ prefix };
        if (index < 10)
          token += '0';
        token += std::to_string(index);
        return path.empty() ? token : path + "." + token;
      };

      if (node.kind == PackNodeKind::parallel) {
        for (std::size_t i = 0; i < node.children.size(); ++i) {
          std::uint32_t child_positive = positive;
          if (node.link.resistance > 0.0) {
            child_positive = allocateNode();
            resistor(positive, child_positive, node.link.resistance);
          }
          compile(node.children[i], child_positive, negative, childPath(i));
        }
        return;
      }

      std::uint32_t cursor = positive;
      for (std::size_t i = 0; i < node.children.size(); ++i) {
        const bool last = i + 1 == node.children.size();
        const auto child_negative = last ? negative : allocateNode();
        compile(node.children[i], cursor, child_negative, childPath(i));
        if (!last && node.link.resistance > 0.0) {
          cursor = allocateNode();
          resistor(child_negative, cursor, node.link.resistance);
        } else {
          cursor = child_negative;
        }
      }
    }
  };

  bool assignBatchLocations(CompiledPackTopology &pack)
  {
    std::map<std::string, std::uint32_t> batches;
    std::map<std::string, bool> thermal;
    std::map<std::string, std::uint32_t> next_lane;
    for (const auto &cell : pack.cells) {
      batches.emplace(cell.archetype, 0);
      const auto [it, inserted] = thermal.emplace(cell.archetype, cell.thermal);
      if (!inserted && it->second != cell.thermal)
        return false;
    }
    std::uint32_t batch{};
    for (auto &[name, index] : batches) {
      index = batch++;
      pack.batch_archetypes.push_back(name);
    }
    for (auto &cell : pack.cells)
      cell.location = { .batch = batches.at(cell.archetype),
                        .lane = next_lane[cell.archetype]++ };
    return true;
  }

  void compileElectricalMetadata(CompiledPackTopology &pack, std::uint32_t node_count)
  {
    auto &netlist = pack.electrical;
    netlist.node_count = node_count;
    std::vector<std::vector<std::uint32_t>> adjacency(node_count);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> sparsity;
    for (const auto &branch : netlist.branches) {
      adjacency[branch.node_positive].push_back(branch.node_negative);
      adjacency[branch.node_negative].push_back(branch.node_positive);
      sparsity.emplace_back(branch.node_positive, branch.node_positive);
      sparsity.emplace_back(branch.node_negative, branch.node_negative);
      sparsity.emplace_back(std::min(branch.node_positive, branch.node_negative),
                            std::max(branch.node_positive, branch.node_negative));
    }
    std::sort(sparsity.begin(), sparsity.end());
    sparsity.erase(std::unique(sparsity.begin(), sparsity.end()), sparsity.end());
    netlist.nodal_sparsity = std::move(sparsity);

    std::vector<unsigned char> seen(node_count);
    std::queue<std::uint32_t> pending;
    pending.push(netlist.terminal_positive);
    seen[netlist.terminal_positive] = 1;
    while (!pending.empty()) {
      const auto node = pending.front();
      pending.pop();
      for (const auto next : adjacency[node])
        if (seen[next] == 0) {
          seen[next] = 1;
          pending.push(next);
        }
    }
    netlist.connected = std::all_of(seen.begin(), seen.end(), [](unsigned char value) {
      return value != 0;
    });
    netlist.index1_candidate = netlist.connected;

    netlist.series_parallel_ladder = false;
    if (!netlist.connected
        || std::any_of(netlist.branches.begin(), netlist.branches.end(), [](const auto &branch) {
             return branch.kind != ElectricalBranchKind::cell;
           }))
      return;
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::uint32_t>> layers;
    std::vector<std::vector<std::uint32_t>> simple_adjacency(node_count);
    for (const auto &branch : netlist.branches) {
      const auto endpoints = std::minmax(branch.node_positive, branch.node_negative);
      auto &cells = layers[{ endpoints.first, endpoints.second }];
      if (cells.empty()) {
        simple_adjacency[endpoints.first].push_back(endpoints.second);
        simple_adjacency[endpoints.second].push_back(endpoints.first);
      }
      cells.push_back(branch.cell);
    }
    if (layers.size() + 1 != node_count
        || simple_adjacency[netlist.terminal_positive].size() != 1
        || simple_adjacency[netlist.terminal_negative].size() != 1)
      return;
    for (std::uint32_t node = 0; node < node_count; ++node)
      if (node != netlist.terminal_positive && node != netlist.terminal_negative
          && simple_adjacency[node].size() != 2)
        return;

    netlist.ladder_offsets.push_back(0);
    netlist.ladder_nodes.push_back(netlist.terminal_positive);
    std::uint32_t previous = node_count;
    std::uint32_t current = netlist.terminal_positive;
    while (current != netlist.terminal_negative) {
      const auto &neighbors = simple_adjacency[current];
      const auto next_it = std::find_if(neighbors.begin(), neighbors.end(), [&](std::uint32_t node) {
        return node != previous;
      });
      if (next_it == neighbors.end()) {
        netlist.ladder_offsets.clear();
        netlist.ladder_cells.clear();
        netlist.ladder_nodes.clear();
        return;
      }
      const auto next = *next_it;
      const auto endpoints = std::minmax(current, next);
      auto cells = layers.at({ endpoints.first, endpoints.second });
      std::sort(cells.begin(), cells.end());
      for (const auto cell : cells) {
        const auto branch = std::find_if(netlist.branches.begin(), netlist.branches.end(), [&](const auto &candidate) {
          return candidate.kind == ElectricalBranchKind::cell
                 && candidate.cell == cell;
        });
        if (branch == netlist.branches.end() || branch->node_positive != current
            || branch->node_negative != next) {
          netlist.ladder_offsets.clear();
          netlist.ladder_cells.clear();
          netlist.ladder_nodes.clear();
          return;
        }
        netlist.ladder_cells.push_back(cell);
      }
      netlist.ladder_offsets.push_back(static_cast<std::uint32_t>(netlist.ladder_cells.size()));
      netlist.ladder_nodes.push_back(next);
      previous = current;
      current = next;
    }
    netlist.series_parallel_ladder = netlist.ladder_cells.size() == pack.cells.size();
  }

  slide::Status compileThermal(const PackDescription &description,
                               CompiledPackTopology &pack)
  {
    auto &graph = pack.thermal;
    graph.cell_count = static_cast<std::uint32_t>(pack.cells.size());
    graph.boundary_count = static_cast<std::uint32_t>(description.thermal_boundaries.size());
    std::unordered_map<std::string, std::uint32_t> endpoint;
    endpoint.reserve(pack.cells.size() + description.thermal_boundaries.size());
    for (std::uint32_t i = 0; i < graph.cell_count; ++i)
      if (!endpoint.emplace(pack.cells[i].path, i).second)
        return slide::Status::Invalid_parameters;
    for (std::uint32_t i = 0; i < graph.boundary_count; ++i) {
      const auto &name = description.thermal_boundaries[i].name;
      if (name.empty() || !endpoint.emplace(name, graph.cell_count + i).second)
        return slide::Status::Invalid_parameters;
    }

    for (const auto &source : description.thermal_links) {
      const auto a = endpoint.find(source.endpoint_a);
      const auto b = endpoint.find(source.endpoint_b);
      if (a == endpoint.end() || b == endpoint.end() || a->second == b->second
          || !is_finite(source.conductance) || source.conductance < 0.0)
        return slide::Status::Invalid_parameters;
      if ((a->second < graph.cell_count && !pack.cells[a->second].thermal)
          || (b->second < graph.cell_count && !pack.cells[b->second].thermal)
          || (a->second >= graph.cell_count && b->second >= graph.cell_count))
        return slide::Status::Invalid_parameters;
      if (source.conductance == 0.0)
        continue;
      graph.edges.push_back({ .low = std::min(a->second, b->second),
                              .high = std::max(a->second, b->second),
                              .conductance = source.conductance });
    }
    std::sort(graph.edges.begin(), graph.edges.end(), [](const auto &a, const auto &b) {
      return std::tie(a.low, a.high) < std::tie(b.low, b.high);
    });
    for (std::size_t i = 1; i < graph.edges.size(); ++i)
      if (graph.edges[i - 1].low == graph.edges[i].low
          && graph.edges[i - 1].high == graph.edges[i].high)
        return slide::Status::Invalid_parameters;

    const auto endpoint_count = static_cast<std::size_t>(graph.cell_count + graph.boundary_count);
    graph.offsets.assign(endpoint_count + 1, 0);
    for (const auto &edge : graph.edges) {
      ++graph.offsets[static_cast<std::size_t>(edge.low) + 1];
      ++graph.offsets[static_cast<std::size_t>(edge.high) + 1];
    }
    std::partial_sum(graph.offsets.begin(), graph.offsets.end(), graph.offsets.begin());
    graph.incidents.resize(2 * graph.edges.size());
    auto cursor = graph.offsets;
    for (std::uint32_t edge = 0; edge < graph.edges.size(); ++edge) {
      const auto &pair = graph.edges[edge];
      graph.incidents[cursor[pair.low]++] = { .edge = edge, .sign = 1 };
      graph.incidents[cursor[pair.high]++] = { .edge = edge, .sign = -1 };
    }
    graph.edge_flux.resize(graph.edges.size());
    graph.trial_edge_flux.resize(graph.edges.size());
    graph.trial_endpoint_heat.resize(endpoint_count);
    graph.trial_edge_incidence.resize(graph.edges.size());
    return slide::Status::Success;
  }

} // namespace

slide::Status CompiledThermalGraph::assemble(std::span<const real_t> cell_temperature,
                                             std::span<const real_t>
                                               boundary_temperature,
                                             std::span<real_t>
                                               q_ext,
                                             std::span<real_t>
                                               boundary_heat)
{
  if (cell_temperature.size() != cell_count || q_ext.size() != cell_count
      || boundary_temperature.size() != boundary_count || boundary_heat.size() != boundary_count)
    return slide::Status::Invalid_parameters;
  const auto endpoint_count = static_cast<std::size_t>(cell_count)
                              + static_cast<std::size_t>(boundary_count);
  if (offsets.size() != endpoint_count + 1 || offsets.empty()
      || offsets.front() != 0 || offsets.back() != incidents.size()
      || edges.size() > std::numeric_limits<std::size_t>::max() / 2
      || incidents.size() != 2 * edges.size()
      || edge_flux.size() != edges.size()
      || trial_edge_flux.size() != edges.size()
      || trial_endpoint_heat.size() != endpoint_count
      || trial_edge_incidence.size() != edges.size())
    return slide::Status::Invalid_parameters;
  for (const auto value : cell_temperature)
    if (!is_finite(value))
      return slide::Status::Invalid_states;
  for (const auto value : boundary_temperature)
    if (!is_finite(value))
      return slide::Status::Invalid_states;

  auto temperature = [&](std::uint32_t endpoint) {
    return endpoint < cell_count ? cell_temperature[endpoint]
                                 : boundary_temperature[endpoint - cell_count];
  };
  for (std::size_t i = 0; i < edges.size(); ++i) {
    const auto &edge = edges[i];
    if (edge.low >= endpoint_count || edge.high >= endpoint_count
        || edge.low >= edge.high || !is_finite(edge.conductance)
        || !(edge.conductance > 0.0))
      return slide::Status::Invalid_parameters;
    const real_t difference = temperature(edge.high) - temperature(edge.low);
    if (!is_finite(difference))
      return slide::Status::Invalid_states;
    const real_t flux = edge.conductance * difference;
    if (!is_finite(flux))
      return slide::Status::Invalid_states;
    trial_edge_flux[i] = flux;
  }
  std::fill(trial_edge_incidence.begin(), trial_edge_incidence.end(), 0);
  for (std::size_t endpoint = 0; endpoint < endpoint_count; ++endpoint) {
    if (offsets[endpoint] > offsets[endpoint + 1]
        || offsets[endpoint + 1] > incidents.size())
      return slide::Status::Invalid_parameters;
    real_t total{};
    for (std::uint32_t i = offsets[endpoint]; i < offsets[endpoint + 1]; ++i) {
      const auto &incident = incidents[i];
      if (incident.edge >= edges.size()
          || (incident.sign != 1 && incident.sign != -1))
        return slide::Status::Invalid_parameters;
      const auto &edge = edges[incident.edge];
      unsigned char incidence_bit{};
      if (endpoint == edge.low && incident.sign == 1)
        incidence_bit = 1;
      else if (endpoint == edge.high && incident.sign == -1)
        incidence_bit = 2;
      else
        return slide::Status::Invalid_parameters;
      if ((trial_edge_incidence[incident.edge] & incidence_bit) != 0)
        return slide::Status::Invalid_parameters;
      trial_edge_incidence[incident.edge] |= incidence_bit;
      const real_t contribution = static_cast<real_t>(incident.sign)
                                  * trial_edge_flux[incident.edge];
      const real_t updated = total + contribution;
      if (!is_finite(contribution) || !is_finite(updated))
        return slide::Status::Invalid_states;
      total = updated;
    }
    trial_endpoint_heat[endpoint] = total;
  }
  if (std::any_of(trial_edge_incidence.begin(),
                  trial_edge_incidence.end(),
                  [](unsigned char incidence) { return incidence != 3; }))
    return slide::Status::Invalid_parameters;

  std::copy(trial_edge_flux.begin(), trial_edge_flux.end(), edge_flux.begin());
  for (std::size_t endpoint = 0; endpoint < endpoint_count; ++endpoint) {
    const real_t total = trial_endpoint_heat[endpoint];
    if (endpoint < cell_count)
      q_ext[endpoint] = total;
    else
      boundary_heat[endpoint - cell_count] = total;
  }
  return slide::Status::Success;
}

slide::Status compilePackDescription(const PackDescription &description,
                                     CompiledPackTopology &output)
{
  CompileContext context;
  context.compile(description.root, 0, 1, {});
  if (!context.valid || context.result.cells.empty())
    return slide::Status::Invalid_parameters;
  if (!assignBatchLocations(context.result))
    return slide::Status::Invalid_parameters;
  compileElectricalMetadata(context.result, context.next_node);
  if (!context.result.electrical.connected)
    return slide::Status::Invalid_parameters;
  const auto status = compileThermal(description, context.result);
  if (status != slide::Status::Success)
    return status;
  output = std::move(context.result);
  return slide::Status::Success;
}

} // namespace slide::core
