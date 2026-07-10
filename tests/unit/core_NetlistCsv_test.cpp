/**
 * @file core_NetlistCsv_test.cpp
 * @brief Phase-9 bounded liionpack CSV parser and atomicity regressions.
 */

#include "../../src/core/NetlistCsv.hpp"
#include "../../src/core/PackSolver.hpp"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <utility>

using namespace slide;

namespace {

constexpr std::string_view parallel_fixture =
  "desc,node1,node2,value,node1_x,node1_y,node2_x,node2_y\n"
  "V0,2,1,4.2,0,1,0,0\n"
  "Ri0,3,2,0.01,0,2,0,1\n"
  "Rc0,4,3,0.01,0,3,0,2\n"
  "V1,6,5,4.2,1,1,1,0\n"
  "Ri1,7,6,0.01,1,2,1,1\n"
  "Rc1,4,7,0.01,1,3,1,2\n"
  "Rbn0,1,5,0.0001,0,0,1,0\n"
  "Rtp1,4,8,0.00001,0,3,-1,3\n"
  "I0,8,0,80,-1,3,-1,0\n"
  "Rtn1,0,1,0.00001,-1,0,0,0\n";

core::CompiledPackTopology sentinelTopology()
{
  core::CompiledPackTopology result;
  REQUIRE(core::compilePackDescription(
            { .root = core::cell({ .archetype = "sentinel" }) }, result)
          == Status::Success);
  return result;
}

void requireSentinel(const core::CompiledPackTopology &topology)
{
  REQUIRE(topology.cells.size() == 1);
  CHECK(topology.cells[0].path == "c00");
  CHECK(topology.cells[0].archetype == "sentinel");
  CHECK(topology.batch_archetypes == std::vector<std::string>{ "sentinel" });
  CHECK(topology.electrical.node_count == 2);
  CHECK(topology.electrical.terminal_positive == 0);
  CHECK(topology.electrical.terminal_negative == 1);
  CHECK(topology.electrical.branches.size() == 1);
  CHECK(topology.electrical.branches[0].kind
        == core::ElectricalBranchKind::cell);
  CHECK(topology.electrical.branches[0].node_positive == 0);
  CHECK(topology.electrical.branches[0].node_negative == 1);
  CHECK(topology.electrical.connected);
  CHECK(topology.electrical.series_parallel_ladder);
  CHECK(topology.electrical.ladder_offsets
        == std::vector<std::uint32_t>{ 0, 1 });
  CHECK(topology.electrical.ladder_cells
        == std::vector<std::uint32_t>{ 0 });
  CHECK(topology.electrical.ladder_nodes
        == std::vector<std::uint32_t>{ 0, 1 });
  CHECK(topology.thermal.cell_count == 1);
  CHECK(topology.thermal.edges.empty());
}

} // namespace

TEST_CASE("liionpack CSV imports arbitrary resistor/cell graphs atomically",
          "[core][pack][netlist][csv][P9-G2]")
{
  core::CompiledPackTopology topology;
  core::NetlistCsvDiagnostic diagnostic;
  REQUIRE(core::parseLiionpackNetlistCsv(
            parallel_fixture,
            topology,
            diagnostic,
            { .archetype = "spm", .thermal = true })
          == Status::Success);
  CHECK(diagnostic.message.empty());
  REQUIRE(topology.cells.size() == 2);
  CHECK(topology.cells[0].path == "V0");
  CHECK(topology.cells[1].path == "V1");
  CHECK(topology.cells[0].archetype == "spm");
  CHECK(topology.cells[0].thermal);
  CHECK(topology.cells[0].location.batch == 0);
  CHECK(topology.cells[0].location.lane == 0);
  CHECK(topology.cells[1].location.lane == 1);
  CHECK(topology.batch_archetypes == std::vector<std::string>{ "spm" });
  CHECK(topology.electrical.node_count == 9);
  CHECK(topology.electrical.terminal_positive == 8);
  CHECK(topology.electrical.terminal_negative == 0);
  CHECK(topology.electrical.branches.size() == 9);
  CHECK(topology.electrical.connected);
  CHECK_FALSE(topology.electrical.series_parallel_ladder);
  CHECK(topology.thermal.cell_count == 2);
  CHECK(topology.thermal.offsets == std::vector<std::uint32_t>{ 0, 0, 0 });
  core::SolverWorkspace workspace;
  CHECK(workspace.configure(topology.electrical, topology.cells.size())
        == Status::Success);
}

TEST_CASE("liionpack CSV contracts ideal wires and remaps sparse node labels",
          "[core][pack][netlist][csv][P9-G2]")
{
  constexpr std::string_view csv =
    "value,node2,desc,node1,ignored\r\n"
    "4.2,4000000000,\"V0\",4294967295,x\r\n"
    "0,4294967294,Rwire,4294967295,x\r\n"
    "5,4000000000,I0,4294967294,x\r\n";
  core::CompiledPackTopology topology;
  core::NetlistCsvDiagnostic diagnostic;
  REQUIRE(core::parseLiionpackNetlistCsv(csv, topology, diagnostic)
          == Status::Success);
  REQUIRE(topology.cells.size() == 1);
  CHECK(topology.electrical.node_count == 2);
  CHECK(topology.electrical.branches.size() == 1);
  CHECK(topology.electrical.terminal_positive
        == topology.electrical.branches[0].node_positive);
  CHECK(topology.electrical.terminal_negative
        == topology.electrical.branches[0].node_negative);
  CHECK(topology.electrical.series_parallel_ladder);
}

TEST_CASE("resistor-free CSV receives the same ladder fast-path metadata",
          "[core][pack][netlist][csv][ladder][P9-G2]")
{
  constexpr std::string_view csv =
    "desc,node1,node2,value\n"
    "V0,3,2,4.2\n"
    "V1,3,2,4.2\n"
    "V2,2,1,4.2\n"
    "V3,2,1,4.2\n"
    "I0,3,1,8\n";
  core::CompiledPackTopology imported;
  core::NetlistCsvDiagnostic diagnostic;
  REQUIRE(core::parseLiionpackNetlistCsv(csv, imported, diagnostic)
          == Status::Success);

  core::CompiledPackTopology composed;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(2, core::parallel(2, core::cell())) },
            composed)
          == Status::Success);
  CHECK(imported.electrical.series_parallel_ladder);
  CHECK(imported.electrical.ladder_offsets
        == composed.electrical.ladder_offsets);
  CHECK(imported.electrical.ladder_cells == composed.electrical.ladder_cells);
  CHECK(imported.electrical.ladder_nodes.size()
        == composed.electrical.ladder_nodes.size());
}

TEST_CASE("liionpack CSV rejects hostile grammar and topology without publication",
          "[core][pack][netlist][csv][parser][P9-G2]")
{
  const auto reject = [](std::string_view csv) {
    auto topology = sentinelTopology();
    core::NetlistCsvDiagnostic diagnostic;
    CHECK(core::parseLiionpackNetlistCsv(csv, topology, diagnostic)
          == Status::Invalid_parameters);
    CHECK_FALSE(diagnostic.message.empty());
    requireSentinel(topology);
  };

  reject("");
  reject("desc,node1,node2,value\n");
  reject("desc,node1,node2\nV0,1,0\nI0,1,0\n");
  reject("desc,node1,node1,value\nV0,1,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\n\"V0,1,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2,extra\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nX0,1,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\nV0,2,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\nR0,2,1,0.1\nI0,1,0,1\nR0,2,1,0.1\n");
  reject("desc,node1,node2,value\n,1,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV 0,1,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,+1,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,-1,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,4294967296,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1x,0,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,1,4.2\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\nI0,1,1,1\n");
  reject("desc,node1,node2,value\nV0,1,0,nan\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,inf\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,1e999\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\nR0,2,1,-0.1\nI0,2,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\n");
  reject("desc,node1,node2,value\nR0,1,0,0.1\nI0,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\nI0,1,0,1\nI1,1,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\nI0,2,0,1\n");
  reject("desc,node1,node2,value\nV0,1,0,4.2\nR0,2,1,0\nI0,2,1,1\n");
  reject("desc,node1,node2,value\rV0,1,0,4.2\rI0,1,0,1\r");
  reject(std::string{ "desc,node1,node2,value\nV0,1,0,4.2\0I0,1,0,1\n",
                      50 });

  const std::string overlong_descriptor(128, 'V');
  reject("desc,node1,node2,value\n" + overlong_descriptor
         + ",1,0,4.2\nI0,1,0,1\n");

  auto invalid_options = sentinelTopology();
  core::NetlistCsvDiagnostic options_diagnostic;
  CHECK(core::parseLiionpackNetlistCsv(
          parallel_fixture,
          invalid_options,
          options_diagnostic,
          { .archetype = "" })
        == Status::Invalid_parameters);
  requireSentinel(invalid_options);

  std::string long_archetype(1025, 'x');
  CHECK(core::parseLiionpackNetlistCsv(
          parallel_fixture,
          invalid_options,
          options_diagnostic,
          { .archetype = std::move(long_archetype) })
        == Status::Invalid_parameters);
  requireSentinel(invalid_options);

  std::string oversized{ parallel_fixture };
  oversized.append(4U * 1024U * 1024U, 'x');
  reject(oversized);
}

TEST_CASE("liionpack CSV file reads are bounded and exact",
          "[core][pack][netlist][csv][file][P9-G2]")
{
  const auto path = std::filesystem::temp_directory_path()
                    / ("slide_netlist_"
                       + std::filesystem::current_path().filename().string()
                       + ".csv");
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.write(parallel_fixture.data(),
                 static_cast<std::streamsize>(parallel_fixture.size()));
    REQUIRE(output.good());
  }
  core::CompiledPackTopology topology;
  core::NetlistCsvDiagnostic diagnostic;
  REQUIRE(core::loadLiionpackNetlistCsv(path, topology, diagnostic)
          == Status::Success);
  CHECK(topology.cells.size() == 2);

  const auto sentinel = sentinelTopology();
  topology = sentinel;
  {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.seekp(static_cast<std::streamoff>(4U * 1024U * 1024U));
    output.put('x');
    REQUIRE(output.good());
  }
  CHECK(core::loadLiionpackNetlistCsv(path, topology, diagnostic)
        == Status::Invalid_parameters);
  requireSentinel(topology);
  std::filesystem::remove(path, ignored);
}
