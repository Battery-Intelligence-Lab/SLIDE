/**
 * @file core_NetlistCsv_test.cpp
 * @brief Phase-9 bounded liionpack CSV parser and atomicity regressions.
 */

#include "../../src/core/NetlistCsv.hpp"
#include "../../src/core/BoundedFileReader.hpp"
#include "../../src/core/PackSolver.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <istream>
#include <sstream>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

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

class ThrowingReadBuffer : public std::streambuf
{
protected:
  std::streamsize xsgetn(char *, std::streamsize) override
  {
    throw std::ios_base::failure("injected bounded-read failure");
  }
};

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

TEST_CASE("liionpack CSV parses scientific-notation values exactly",
          "[core][pack][netlist][csv]")
{
  //!< pandas/numpy write small floats in exponent form by default, so a liionpack
  //!< `netlist.to_csv()` routinely contains values like 1e-05 for the connection
  //!< resistances. `parseValue` has an exponent branch (NetlistCsv.cpp:185-197),
  //!< but no test had ever driven a value through it: the fixtures elsewhere in
  //!< this file all use plain decimals. This pins the accepted spellings AND the
  //!< resulting magnitudes, so a regression cannot silently reinterpret 1e-05.
  constexpr std::string_view csv =
    "desc,node1,node2,value\n"
    "V0,2,1,4.2\n"
    "R0,3,2,1e-05\n"   //!< numpy's default small-float spelling
    "R1,4,3,1.5E+03\n" //!< uppercase E and an explicit + sign
    "R2,5,4,2.5e2\n"   //!< exponent with no sign at all
    "I0,5,1,80\n";
  core::CompiledPackTopology topology;
  core::NetlistCsvDiagnostic diagnostic;
  REQUIRE(core::parseLiionpackNetlistCsv(csv, topology, diagnostic)
          == Status::Success);
  CHECK(diagnostic.message.empty());

  std::vector<core::real_t> resistances;
  for (const auto &branch : topology.electrical.branches)
    if (branch.kind == core::ElectricalBranchKind::resistor)
      resistances.push_back(branch.resistance);
  std::ranges::sort(resistances);

  //!< Exact equality is the right assertion, not a tolerance: `std::from_chars`
  //!< and the C++ literal must both round to the same nearest double.
  REQUIRE(resistances.size() == 3);
  CHECK(resistances[0] == 1e-05);
  CHECK(resistances[1] == 250.0);
  CHECK(resistances[2] == 1500.0);
}

TEST_CASE("liionpack CSV rejects noncanonical values rather than guessing",
          "[core][pack][netlist][csv]")
{
  //!< Pin the complete strict grammar, including exponent digits, the
  //!< minus-only leading sign, no leading zeros, and nonempty integer/fraction.
  const auto rejects = [](std::string_view value) {
    std::string csv = "desc,node1,node2,value\nV0,2,1,4.2\nR0,3,2,";
    csv += value;
    csv += "\nI0,3,1,80\n";
    core::CompiledPackTopology topology;
    core::NetlistCsvDiagnostic diagnostic;
    const Status status =
      core::parseLiionpackNetlistCsv(csv, topology, diagnostic);
    //!< Rejection must also be atomic: nothing published on the way out.
    CHECK(topology.cells.empty());
    CHECK(topology.electrical.branches.empty());
    return status == Status::Invalid_parameters;
  };
  CHECK(rejects("1e"));
  CHECK(rejects("1e+"));
  CHECK(rejects("1e-"));
  CHECK(rejects("1e+x"));
  CHECK(rejects("1.e5e5"));
  CHECK(rejects("+4.2"));
  CHECK(rejects("01"));
  CHECK(rejects(".5"));
  CHECK(rejects("1."));
}

TEST_CASE("liionpack CSV reports exact first-failure diagnostics",
          "[core][pack][netlist][csv][diagnostic][MQ.2][C1]")
{
  struct ExactFailure
  {
    std::string_view name;
    std::string csv;
    std::size_t row;
    std::size_t offset;
    std::string_view message;
  };

  std::vector<ExactFailure> failures;
  failures.reserve(12);
  failures.push_back(
    { "row width",
      "desc,node1,node2,value\nV0,1,0\n",
      2,
      0,
      "liionpack CSV row width differs from header" });
  failures.push_back(
    { "unsupported descriptor",
      "desc,node1,node2,value\nX0,1,0,4.2\n",
      2,
      0,
      "unsupported liionpack descriptor" });
  failures.push_back(
    { "duplicate descriptor",
      "desc,node1,node2,value\nV0,1,0,4.2\nV0,2,0,4.2\n",
      3,
      0,
      "duplicate liionpack descriptor" });
  failures.push_back(
    { "invalid node",
      "desc,node1,node2,value\nV0,+1,0,4.2\n",
      2,
      0,
      "invalid liionpack node label" });
  failures.push_back(
    { "identical endpoints",
      "desc,node1,node2,value\nV0,1,1,4.2\n",
      2,
      0,
      "liionpack element has identical endpoints" });
  failures.push_back(
    { "invalid value",
      "desc,node1,node2,value\nV0,1,0,+4.2\n",
      2,
      0,
      "invalid liionpack element value" });
  failures.push_back(
    { "negative resistance",
      "desc,node1,node2,value\nV0,1,0,4.2\nR0,2,1,-0.1\n",
      3,
      0,
      "liionpack resistance must be non-negative" });
  failures.push_back(
    { "ideal wire short",
      "desc,node1,node2,value\n"
      "V0,1,0,4.2\n"
      "R0,2,1,0\n"
      "I0,2,1,1\n",
      4,
      0,
      "liionpack ideal wire shorts an element or terminal source" });

  std::string wide_header;
  wide_header.reserve(66);
  for (std::size_t column = 0; column < 33; ++column) {
    if (column != 0)
      wide_header.push_back(',');
    wide_header.push_back('x');
  }
  wide_header.push_back('\n');
  failures.push_back(
    { "column limit",
      std::move(wide_header),
      1,
      64,
      "CSV row exceeds 32 columns" });

  std::string unquoted_field(65'537, 'x');
  failures.push_back(
    { "unquoted field limit",
      std::move(unquoted_field),
      1,
      65'537,
      "CSV field exceeds 65536 bytes" });

  std::string quoted_field{ "\"" };
  quoted_field.append(65'537, 'x');
  failures.push_back(
    { "quoted field limit",
      std::move(quoted_field),
      1,
      65'538,
      "CSV field exceeds 65536 bytes" });
  failures.push_back(
    { "quote inside unquoted field",
      "desc,node1,node2,value\nV0,1,0,4\"2\n",
      2,
      31,
      "quote inside unquoted CSV field" });

  for (const auto &failure : failures) {
    CAPTURE(failure.name);
    core::CompiledPackTopology topology;
    core::NetlistCsvDiagnostic diagnostic{
      .row = 777,
      .offset = 888,
      .message = "poison",
    };
    CHECK(core::parseLiionpackNetlistCsv(
            failure.csv, topology, diagnostic)
          == Status::Invalid_parameters);
    CHECK(diagnostic.row == failure.row);
    CHECK(diagnostic.offset == failure.offset);
    CHECK(diagnostic.message == failure.message);
  }
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
  constexpr char embedded_nul_csv[] =
    "desc,node1,node2,value\nV0,1,0,4.2\0I0,1,0,1\n";
  reject(std::string_view{ embedded_nul_csv, sizeof embedded_nul_csv - 1 });

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
  diagnostic = { .row = 777, .offset = 888, .message = "poison" };
  CHECK(core::loadLiionpackNetlistCsv(path, topology, diagnostic)
        == Status::Invalid_parameters);
  CHECK(diagnostic.row == 0);
  CHECK(diagnostic.offset == 0);
  CHECK(diagnostic.message == "liionpack CSV exceeds 4194304 bytes");
  requireSentinel(topology);
  std::filesystem::remove(path, ignored);

  CHECK(core::loadLiionpackNetlistCsv(path, topology, diagnostic)
        == Status::Invalid_parameters);
  CHECK(diagnostic.message.find("open") != std::string::npos);
  requireSentinel(topology);
}

TEST_CASE("bounded parser reads distinguish exact EOF, limits, and I/O faults",
          "[core][parser][file][coverage]")
{
  std::string output{ "sentinel" };
  std::istringstream exact{ "abcd" };
  const auto exact_result = core::detail::readBoundedStream(exact, 4, output);
  CHECK(exact_result == core::detail::BoundedFileRead::success);
  CHECK(core::detail::boundedFileStatus(exact_result) == Status::Success);
  CHECK(output == "abcd");

  output = "sentinel";
  std::istringstream oversized{ "abcde" };
  const auto oversized_result =
    core::detail::readBoundedStream(oversized, 4, output);
  CHECK(oversized_result == core::detail::BoundedFileRead::too_large);
  CHECK(core::detail::boundedFileStatus(oversized_result)
        == Status::Invalid_parameters);
  CHECK(output == "sentinel");

  ThrowingReadBuffer buffer;
  std::istream failing{ &buffer };
  const auto failure_result =
    core::detail::readBoundedStream(failing, 4, output);
  CHECK(failure_result == core::detail::BoundedFileRead::io_failed);
  CHECK(core::detail::boundedFileStatus(failure_result)
        == Status::Numerical_failure);
  CHECK(output == "sentinel");

  CHECK(core::detail::boundedFileStatus(
          core::detail::BoundedFileRead::open_failed)
        == Status::Invalid_parameters);
}

TEST_CASE("liionpack CSV enforces row and retained archetype budgets",
          "[core][pack][netlist][csv][limits][coverage]")
{
  std::string too_many_rows{ "desc,node1,node2,value\n" };
  too_many_rows.reserve(2U * 1024U * 1024U);
  for (std::size_t row = 0; row <= 100'000; ++row)
    too_many_rows += "V" + std::to_string(row) + ",1,0,4.2\n";
  core::CompiledPackTopology topology = sentinelTopology();
  core::NetlistCsvDiagnostic diagnostic{
    .row = 777,
    .offset = 888,
    .message = "poison",
  };
  CHECK(core::parseLiionpackNetlistCsv(
          too_many_rows, topology, diagnostic)
        == Status::Invalid_parameters);
  CHECK(diagnostic.row == 100'002);
  CHECK(diagnostic.offset == 0);
  CHECK(diagnostic.message == "liionpack CSV exceeds 100000 rows");
  requireSentinel(topology);

  std::string retained{ "desc,node1,node2,value\n" };
  for (std::size_t cell = 0; cell < 4097; ++cell)
    retained += "V" + std::to_string(cell) + ",1,0,4.2\n";
  retained += "I0,1,0,1\n";
  topology = sentinelTopology();
  CHECK(core::parseLiionpackNetlistCsv(
          retained,
          topology,
          diagnostic,
          { .archetype = std::string(1024, 'x') })
        == Status::Invalid_parameters);
  CHECK(diagnostic.message.find("retained-text budget") != std::string::npos);
  requireSentinel(topology);
}
