/**
 * @file core_PackTopology_test.cpp
 * @brief Phase-2 pack flattening and D-21 thermal adjacency gates.
 */

#include "../../src/core/PackTopology.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstring>

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
}
