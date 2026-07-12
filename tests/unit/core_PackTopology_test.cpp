/**
 * @file core_PackTopology_test.cpp
 * @brief Phase-2 pack flattening and D-21 thermal adjacency gates.
 */

#include "../../src/core/PackTopology.hpp"
#include "../../src/core/PackTopologyInternal.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <set>
#include <string>
#include <string_view>

using namespace slide;

TEST_CASE("nested pack combinators compile to one flat netlist", "[core][pack][compile]")
{
  const auto leaf = core::cell({ .archetype = "spm", .thermal = true });
  const auto nested = core::series(2, core::parallel(2, core::parallel(2, leaf)));
  core::CompiledPackTopology pack;
  REQUIRE(core::compilePackDescription({ .root = nested }, pack) == Status::Success);
  REQUIRE(pack.cells.size() == 8);
  REQUIRE(pack.batch_archetypes == std::vector<std::string>{ "spm" });
  REQUIRE(pack.electrical.connected);
  REQUIRE(pack.electrical.index1_candidate);
  REQUIRE(pack.electrical.series_parallel_ladder);
  REQUIRE(pack.electrical.ladder_offsets == std::vector<std::uint32_t>{ 0, 4, 8 });
  REQUIRE(pack.cells.front().path == "s00.p00.p00");
  REQUIRE(pack.cells.back().path == "s01.p01.p01");
  std::set<std::string> paths;
  for (const auto &cell : pack.cells)
    REQUIRE(paths.insert(cell.path).second);
  for (std::size_t lane = 0; lane < pack.cells.size(); ++lane) {
    REQUIRE(pack.cells[lane].location.batch == 0);
    REQUIRE(pack.cells[lane].location.lane == lane);
  }
}

TEST_CASE("pack links flatten to explicit resistor branches", "[core][pack][compile]")
{
  const auto leaf = core::cell({ .archetype = "ecm" });
  core::CompiledPackTopology parallel_pack;
  REQUIRE(core::compilePackDescription(
            { .root = core::parallel(3, leaf, { .resistance = 0.01 }) }, parallel_pack)
          == Status::Success);
  REQUIRE(parallel_pack.cells.size() == 3);
  REQUIRE(parallel_pack.electrical.branches.size() == 6);
  REQUIRE(parallel_pack.electrical.node_count == 5);

  core::CompiledPackTopology series_pack;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(3, leaf, { .resistance = 0.02 }) }, series_pack)
          == Status::Success);
  REQUIRE(series_pack.cells.size() == 3);
  REQUIRE(series_pack.electrical.branches.size() == 5);
}

TEST_CASE("D-21 thermal compile is canonical and assemble conserves pair energy",
          "[core][pack][thermal]")
{
  const auto root = core::parallel(2, core::cell({ .archetype = "thermal", .thermal = true }));
  const core::PackDescription first{
    .root = root,
    .thermal_boundaries = { { "coolant" } },
    .thermal_links = { { "p01", "coolant", 2.0 }, { "p00", "p01", 3.0 } }
  };
  const core::PackDescription reordered{
    .root = root,
    .thermal_boundaries = { { "coolant" } },
    .thermal_links = { { "p01", "p00", 3.0 }, { "coolant", "p01", 2.0 } }
  };
  core::CompiledPackTopology a, b;
  REQUIRE(core::compilePackDescription(first, a) == Status::Success);
  REQUIRE(core::compilePackDescription(reordered, b) == Status::Success);
  REQUIRE(a.thermal.edges.size() == 2);
  REQUIRE(a.thermal.edges.size() == b.thermal.edges.size());
  REQUIRE(std::memcmp(a.thermal.edges.data(), b.thermal.edges.data(), a.thermal.edges.size() * sizeof(core::ThermalEdge))
          == 0);
  REQUIRE(std::memcmp(a.thermal.incidents.data(), b.thermal.incidents.data(), a.thermal.incidents.size() * sizeof(core::ThermalIncident))
          == 0);

  constexpr std::array cell_temperature{ 300.0, 310.0 };
  constexpr std::array boundary_temperature{ 290.0 };
  std::array<double, 2> q_ext{};
  std::array<double, 1> boundary_heat{};
  REQUIRE(a.thermal.assemble(cell_temperature, boundary_temperature, q_ext, boundary_heat)
          == Status::Success);
  REQUIRE(q_ext[0] == 30.0);
  REQUIRE(q_ext[1] == -70.0);
  REQUIRE(boundary_heat[0] == 40.0);
  REQUIRE(q_ext[0] + q_ext[1] + boundary_heat[0] == 0.0);
}

TEST_CASE("thermal assembly rejects finite derived overflow atomically",
          "[core][pack][thermal][validation][P9-G4]")
{
  const auto root = core::cell({ .archetype = "thermal", .thermal = true });

  SECTION("edge product overflow")
  {
    core::CompiledPackTopology pack;
    REQUIRE(core::compilePackDescription(
              { .root = root,
                .thermal_boundaries = { { "hot" } },
                .thermal_links = { { "c00", "hot", 2.0 } } },
              pack)
            == Status::Success);
    constexpr std::array cell_temperature{ 1e-300 };
    constexpr std::array boundary_temperature{ 1e308 };
    std::array q_ext{ 17.0 };
    std::array boundary_heat{ 23.0 };
    const auto old_edge_flux = pack.thermal.edge_flux;

    REQUIRE(pack.thermal.assemble(
              cell_temperature, boundary_temperature, q_ext, boundary_heat)
            == Status::Invalid_states);
    REQUIRE(q_ext == std::array{ 17.0 });
    REQUIRE(boundary_heat == std::array{ 23.0 });
    REQUIRE(pack.thermal.edge_flux == old_edge_flux);
  }

  SECTION("endpoint accumulation overflow")
  {
    core::CompiledPackTopology pack;
    REQUIRE(core::compilePackDescription(
              { .root = root,
                .thermal_boundaries = { { "hot-a" }, { "hot-b" } },
                .thermal_links = { { "c00", "hot-a", 1.0 },
                                   { "c00", "hot-b", 1.0 } } },
              pack)
            == Status::Success);
    constexpr std::array cell_temperature{ 1e-300 };
    constexpr std::array boundary_temperature{ 1e308, 1e308 };
    std::array q_ext{ 17.0 };
    std::array boundary_heat{ 23.0, 29.0 };
    const auto old_edge_flux = pack.thermal.edge_flux;

    REQUIRE(pack.thermal.assemble(
              cell_temperature, boundary_temperature, q_ext, boundary_heat)
            == Status::Invalid_states);
    REQUIRE(q_ext == std::array{ 17.0 });
    REQUIRE(boundary_heat == std::array{ 23.0, 29.0 });
    REQUIRE(pack.thermal.edge_flux == old_edge_flux);
  }

  SECTION("incidence orientation")
  {
    core::CompiledPackTopology pack;
    REQUIRE(core::compilePackDescription(
              { .root = root,
                .thermal_boundaries = { { "coolant" } },
                .thermal_links = { { "c00", "coolant", 2.0 } } },
              pack)
            == Status::Success);
    REQUIRE(pack.thermal.incidents.size() == 2);
    pack.thermal.incidents[1].sign = 1;
    constexpr std::array cell_temperature{ 300.0 };
    constexpr std::array boundary_temperature{ 290.0 };
    std::array q_ext{ 17.0 };
    std::array boundary_heat{ 23.0 };
    const auto old_edge_flux = pack.thermal.edge_flux;

    REQUIRE(pack.thermal.assemble(
              cell_temperature, boundary_temperature, q_ext, boundary_heat)
            == Status::Invalid_parameters);
    REQUIRE(q_ext == std::array{ 17.0 });
    REQUIRE(boundary_heat == std::array{ 23.0 });
    REQUIRE(pack.thermal.edge_flux == old_edge_flux);
  }
}

TEST_CASE("pack cold validation is atomic", "[core][pack][validation]")
{
  const auto root = core::parallel(2, core::cell({ .archetype = "thermal", .thermal = true }));
  core::CompiledPackTopology output;
  REQUIRE(core::compilePackDescription({ .root = root }, output) == Status::Success);
  const auto old_cells = output.cells.size();
  core::PackDescription duplicate{
    .root = root,
    .thermal_links = { { "p00", "p01", 1.0 }, { "p01", "p00", 2.0 } }
  };
  REQUIRE(core::compilePackDescription(duplicate, output) == Status::Invalid_parameters);
  REQUIRE(output.cells.size() == old_cells);
  REQUIRE(core::compilePackDescription({ .root = core::series(0, core::cell()) }, output)
          == Status::Invalid_parameters);
  REQUIRE(output.cells.size() == old_cells);

  core::PackNode invalid_kind{
    .kind = static_cast<core::PackNodeKind>(255),
    .children = { core::cell() }
  };
  REQUIRE(core::compilePackDescription({ .root = invalid_kind }, output)
          == Status::Invalid_parameters);
  REQUIRE(output.cells.size() == old_cells);
}

TEST_CASE("thermal descriptions reject every invalid endpoint class",
          "[core][pack][thermal][validation][coverage]")
{
  const auto thermal_cell = core::cell({ .archetype = "thermal", .thermal = true });
  const auto isothermal_cell = core::cell({ .archetype = "isothermal" });
  const std::array invalid{
    core::PackDescription{ .root = thermal_cell,
                           .thermal_boundaries = { { "" } } },
    core::PackDescription{ .root = thermal_cell,
                           .thermal_boundaries = { { "coolant" } },
                           .thermal_links = { { "missing", "coolant", 1.0 } } },
    core::PackDescription{ .root = isothermal_cell,
                           .thermal_boundaries = { { "coolant" } },
                           .thermal_links = { { "c00", "coolant", 1.0 } } },
  };

  for (std::size_t index = 0; index < invalid.size(); ++index) {
    CAPTURE(index);
    core::CompiledPackTopology output;
    CHECK(core::compilePackDescription(invalid[index], output)
          == Status::Invalid_parameters);
  }
}

TEST_CASE("electrical validator rejects independently corrupted metadata",
          "[core][pack][validation][coverage]")
{
  core::CompiledPackTopology base;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(2, core::cell()) }, base)
          == Status::Success);

  using Mutator = std::function<void(core::CompiledElectricalNetlist &)>;
  const std::array cases{
    std::pair<std::string_view, Mutator>{ "sparsity", [](auto &netlist) {
                                           netlist.nodal_sparsity.clear();
                                         } },
    std::pair<std::string_view, Mutator>{ "ladder shape", [](auto &netlist) {
                                           netlist.ladder_offsets.clear();
                                         } },
    std::pair<std::string_view, Mutator>{ "ladder node", [](auto &netlist) {
                                           netlist.ladder_nodes[1] = netlist.ladder_nodes[0];
                                         } },
    std::pair<std::string_view, Mutator>{ "ladder range", [](auto &netlist) {
                                           netlist.ladder_offsets[1] = 0;
                                         } },
  };
  for (const auto &[name, mutate] : cases) {
    CAPTURE(name);
    auto netlist = base.electrical;
    mutate(netlist);
    CHECK(core::detail::validateElectricalNetlist(netlist, base.cells.size())
          == Status::Invalid_parameters);
  }

  auto resistor = base.electrical;
  resistor.branches.push_back({ .node_positive = 0,
                                .node_negative = 1,
                                .kind = core::ElectricalBranchKind::resistor,
                                .resistance = 0.0 });
  CHECK(core::detail::validateElectricalNetlist(resistor, base.cells.size())
        == Status::Invalid_parameters);

  auto missing_cell = base.electrical;
  missing_cell.branches[1].kind = core::ElectricalBranchKind::resistor;
  missing_cell.branches[1].resistance = 1.0;
  CHECK(core::detail::validateElectricalNetlist(missing_cell, base.cells.size())
        == Status::Invalid_parameters);

  core::CompiledElectricalNetlist disconnected{
    .node_count = 3,
    .terminal_positive = 0,
    .terminal_negative = 1,
    .branches = { { .node_positive = 0,
                    .node_negative = 1,
                    .kind = core::ElectricalBranchKind::cell,
                    .cell = 0 },
                  { .node_positive = 0,
                    .node_negative = 1,
                    .kind = core::ElectricalBranchKind::cell,
                    .cell = 1 } },
    .nodal_sparsity = { { 0, 0 }, { 0, 1 }, { 1, 1 } },
    .connected = true,
  };
  CHECK(core::detail::validateElectricalNetlist(disconnected, 2)
        == Status::Invalid_parameters);
}

TEST_CASE("imported topology finalization validates names, endpoints, and batch physics",
          "[core][pack][validation][coverage]")
{
  core::CompiledPackTopology base;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(std::vector{
                core::cell({ .archetype = "a" }),
                core::cell({ .archetype = "b" }) }) },
            base)
          == Status::Success);
  const auto node_count = base.electrical.node_count;

  auto invalid_shape = base;
  CHECK(core::detail::finalizeImportedPackTopology(invalid_shape, 1)
        == Status::Invalid_parameters);

  using Mutator = std::function<void(core::CompiledPackTopology &)>;
  const std::array cases{
    std::pair<std::string_view, Mutator>{ "empty path", [](auto &pack) {
                                           pack.cells[0].path.clear();
                                         } },
    std::pair<std::string_view, Mutator>{ "branch endpoint", [=](auto &pack) {
                                           pack.electrical.branches[0].node_positive = node_count;
                                         } },
    std::pair<std::string_view, Mutator>{ "mixed batch physics", [](auto &pack) {
                                           pack.cells[1].archetype = pack.cells[0].archetype;
                                           pack.cells[1].thermal = true;
                                         } },
  };
  for (const auto &[name, mutate] : cases) {
    CAPTURE(name);
    auto candidate = base;
    mutate(candidate);
    CHECK(core::detail::finalizeImportedPackTopology(candidate, node_count)
          == Status::Invalid_parameters);
  }
}

TEST_CASE("thermal assembly rejects every independently corrupted graph field",
          "[core][pack][thermal][validation][coverage]")
{
  const auto make_graph = [] {
    core::CompiledPackTopology pack;
    REQUIRE(core::compilePackDescription(
              { .root = core::cell({ .archetype = "thermal", .thermal = true }),
                .thermal_boundaries = { { "coolant" } },
                .thermal_links = { { "c00", "coolant", 2.0 } } },
              pack)
            == Status::Success);
    return pack.thermal;
  };
  constexpr std::array cell_temperature{ 300.0 };
  constexpr std::array boundary_temperature{ 290.0 };

  SECTION("shape and finite inputs")
  {
    auto graph = make_graph();
    std::array q_ext{ 17.0 };
    std::array boundary_heat{ 23.0 };
    CHECK(graph.assemble({}, boundary_temperature, q_ext, boundary_heat)
          == Status::Invalid_parameters);
    CHECK(q_ext == std::array{ 17.0 });
    CHECK(boundary_heat == std::array{ 23.0 });

    constexpr std::array nonfinite_cell{
      std::bit_cast<double>(UINT64_C(0x7ff8000000000000))
    };
    CHECK(graph.assemble(nonfinite_cell, boundary_temperature, q_ext, boundary_heat)
          == Status::Invalid_states);
    constexpr std::array nonfinite_boundary{
      std::bit_cast<double>(UINT64_C(0x7ff0000000000000))
    };
    CHECK(graph.assemble(cell_temperature, nonfinite_boundary, q_ext, boundary_heat)
          == Status::Invalid_states);
  }

  using Mutator = std::function<void(core::CompiledThermalGraph &)>;
  const std::array cases{
    std::pair<std::string_view, Mutator>{ "workspace", [](auto &graph) {
                                           graph.trial_edge_flux.clear();
                                         } },
    std::pair<std::string_view, Mutator>{ "edge", [](auto &graph) {
                                           graph.edges[0].conductance = 0.0;
                                         } },
    std::pair<std::string_view, Mutator>{ "offset", [](auto &graph) {
                                           graph.offsets[1] = 3;
                                         } },
    std::pair<std::string_view, Mutator>{ "incident", [](auto &graph) {
                                           graph.incidents[0].edge = 1;
                                         } },
  };
  for (const auto &[name, mutate] : cases) {
    CAPTURE(name);
    auto graph = make_graph();
    mutate(graph);
    std::array q_ext{ 17.0 };
    std::array boundary_heat{ 23.0 };
    CHECK(graph.assemble(cell_temperature, boundary_temperature, q_ext, boundary_heat)
          == Status::Invalid_parameters);
    CHECK(q_ext == std::array{ 17.0 });
    CHECK(boundary_heat == std::array{ 23.0 });
  }

  SECTION("finite temperature difference overflow")
  {
    auto graph = make_graph();
    constexpr std::array cold{ -1e308 };
    constexpr std::array hot{ 1e308 };
    std::array q_ext{ 17.0 };
    std::array boundary_heat{ 23.0 };
    CHECK(graph.assemble(cold, hot, q_ext, boundary_heat)
          == Status::Invalid_states);
  }

  SECTION("duplicate edge incidence")
  {
    core::CompiledPackTopology pack;
    REQUIRE(core::compilePackDescription(
              { .root = core::cell({ .archetype = "thermal", .thermal = true }),
                .thermal_boundaries = { { "a" }, { "b" } },
                .thermal_links = { { "c00", "a", 1.0 },
                                   { "c00", "b", 1.0 } } },
              pack)
            == Status::Success);
    pack.thermal.incidents[1] = pack.thermal.incidents[0];
    constexpr std::array boundaries{ 290.0, 295.0 };
    std::array q_ext{ 17.0 };
    std::array boundary_heat{ 23.0, 29.0 };
    CHECK(pack.thermal.assemble(
            cell_temperature, boundaries, q_ext, boundary_heat)
          == Status::Invalid_parameters);
  }
}
