/**
 * @file core_ParserAllocation_test.cpp
 * @brief Deterministic allocation-failure atomicity for Phase-9 parser gates.
 */

#include "../../src/core/Experiment.hpp"
#include "../../src/core/NetlistCsv.hpp"
#include "../../src/core/ParameterSet.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdlib>
#include <new>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <malloc.h>
#endif

namespace {

thread_local bool measure_matching_allocations{};
thread_local std::size_t matching_allocation_size{};
thread_local std::size_t measured_matching_allocations{};
thread_local bool measure_largest_allocation{};
thread_local std::size_t largest_allocation_size{};
thread_local std::size_t largest_allocation_count{};
thread_local bool fail_matching_allocation{};
thread_local std::size_t matching_allocation_cursor{};
thread_local std::size_t failing_matching_allocation{};
thread_local bool allocation_failure_triggered{};
thread_local bool persist_allocation_failure{};

void beforeAllocation(std::size_t bytes)
{
  if (measure_matching_allocations && bytes == matching_allocation_size)
    ++measured_matching_allocations;
  if (measure_largest_allocation) {
    if (bytes > largest_allocation_size) {
      largest_allocation_size = bytes;
      largest_allocation_count = 1;
    } else if (bytes == largest_allocation_size) {
      ++largest_allocation_count;
    }
  }
  if (fail_matching_allocation && allocation_failure_triggered
      && persist_allocation_failure)
    throw std::bad_alloc{};
  if (fail_matching_allocation && bytes == matching_allocation_size) {
    const std::size_t current = matching_allocation_cursor++;
    if (!allocation_failure_triggered
        && current == failing_matching_allocation) {
      allocation_failure_triggered = true;
      throw std::bad_alloc{};
    }
  }
}

#if defined(_WIN32)
void *alignedAllocate(std::size_t bytes, std::size_t alignment)
{
  return _aligned_malloc(bytes, alignment);
}
void alignedRelease(void *pointer) { _aligned_free(pointer); }
#else
void *alignedAllocate(std::size_t bytes, std::size_t alignment)
{
  return std::aligned_alloc(alignment,
                            ((bytes + alignment - 1) / alignment) * alignment);
}
void alignedRelease(void *pointer) { std::free(pointer); }
#endif

class FailAllocationOfSize
{
public:
  FailAllocationOfSize(std::size_t bytes, std::size_t occurrence = 0,
                       bool persistent = false)
  {
    matching_allocation_size = bytes;
    matching_allocation_cursor = 0;
    failing_matching_allocation = occurrence;
    allocation_failure_triggered = false;
    persist_allocation_failure = persistent;
    fail_matching_allocation = true;
  }
  ~FailAllocationOfSize()
  {
    fail_matching_allocation = false;
    persist_allocation_failure = false;
  }
};

class MeasureAllocationsOfSize
{
public:
  explicit MeasureAllocationsOfSize(std::size_t bytes)
  {
    matching_allocation_size = bytes;
    measured_matching_allocations = 0;
    measure_matching_allocations = true;
  }
  ~MeasureAllocationsOfSize() { measure_matching_allocations = false; }
};

class MeasureLargestAllocation
{
public:
  MeasureLargestAllocation()
  {
    largest_allocation_size = 0;
    largest_allocation_count = 0;
    measure_largest_allocation = true;
  }
  ~MeasureLargestAllocation() { measure_largest_allocation = false; }
};

constexpr std::string_view bpx_fixture = R"json({
  "Header": {"BPX": "1.0.0", "Model": "SPM"},
  "Parameterisation": {
    "Cell": {
      "Electrode area [m2]": 0.1027,
      "Nominal cell capacity [A.h]": 5.0,
      "Reference temperature [K]": 298.15,
      "External surface area [m2]": 0.00531
    },
    "Negative electrode": {
      "Thickness [m]": 8.52e-5,
      "Minimum stoichiometry": 0.026,
      "Maximum stoichiometry": 0.91,
      "Maximum concentration [mol.m-3]": 33133.0,
      "Particle radius [m]": 5.86e-6,
      "Surface area per unit volume [m-1]": 383959.0443686007,
      "Diffusivity [m2.s-1]": 3.3e-14,
      "OCP [V]": {"x": [0.0, 0.5, 1.0], "y": [1.98, 0.25, 0.1]},
      "Reaction rate constant [mol.m-2.s-1]": 6.716e-12
    },
    "Positive electrode": {
      "Thickness [m]": 7.56e-5,
      "Minimum stoichiometry": 0.264,
      "Maximum stoichiometry": 0.854,
      "Maximum concentration [mol.m-3]": 63104.0,
      "Particle radius [m]": 5.22e-6,
      "Surface area per unit volume [m-1]": 382183.908045977,
      "Diffusivity [m2.s-1]": 4e-15,
      "OCP [V]": {"x": [0.0, 0.5, 1.0], "y": [4.5, 3.8, 3.0]},
      "Reaction rate constant [mol.m-2.s-1]": 3.544e-11
    }
  },
  "State": {"Initial conditions": {"Initial state-of-charge": 0.75}}
})json";

} // namespace

void *operator new(std::size_t bytes)
{
  beforeAllocation(bytes);
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}

void *operator new[](std::size_t bytes)
{
  beforeAllocation(bytes);
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}

void *operator new(std::size_t bytes, std::align_val_t alignment)
{
  beforeAllocation(bytes);
  if (void *pointer = alignedAllocate(
        bytes != 0 ? bytes : 1, static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}

void *operator new[](std::size_t bytes, std::align_val_t alignment)
{
  beforeAllocation(bytes);
  if (void *pointer = alignedAllocate(
        bytes != 0 ? bytes : 1, static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}

void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void *pointer) noexcept { std::free(pointer); }
void operator delete[](void *pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::align_val_t) noexcept
{
  alignedRelease(pointer);
}
void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept
{
  alignedRelease(pointer);
}
void operator delete[](void *pointer, std::align_val_t) noexcept
{
  alignedRelease(pointer);
}
void operator delete[](void *pointer, std::size_t, std::align_val_t) noexcept
{
  alignedRelease(pointer);
}

using namespace slide;

TEST_CASE("Experiment parser maps a late expansion allocation failure without publication",
          "[core][experiment][parser][allocation][P9]")
{
  const std::vector<std::string> source{
    "Rest for 1 s * 37" + std::string(128, ' '),
  };
  core::Experiment warm;
  core::ParseDiagnostic diagnostic;
  REQUIRE(core::Experiment::parse(source, warm, diagnostic) == Status::Success);
  REQUIRE(warm.segments.size() == 37);

  core::Experiment measured;
  Status measured_status{};
  {
    MeasureLargestAllocation measure;
    measured_status = core::Experiment::parse(source, measured, diagnostic);
  }
  REQUIRE(measured_status == Status::Success);
  REQUIRE(largest_allocation_size > 0);
  REQUIRE(largest_allocation_count == 1);

  core::Experiment sentinel;
  sentinel.segments.push_back(
    { .mode = core::ControlMode::rest, .duration = 17.0, .source = "sentinel" });
  core::Experiment output = sentinel;
  diagnostic = {};
  Status status{};
  {
    FailAllocationOfSize failure{ largest_allocation_size, 0, true };
    status = core::Experiment::parse(source, output, diagnostic);
  }
  REQUIRE(allocation_failure_triggered);
  CHECK(status == Status::Numerical_failure);
  REQUIRE(output.segments.size() == 1);
  CHECK(output.segments[0].duration == 17.0);
  CHECK(output.segments[0].source == "sentinel");
}

TEST_CASE("BPX parser propagates every ParameterSet node allocation failure atomically",
          "[core][parameters][BPX][parser][allocation][P9]")
{
  core::ParameterSet warm;
  std::string diagnostic;
  REQUIRE(core::ParameterSet::fromBpxJson(bpx_fixture, warm, diagnostic)
          == Status::Success);

  // A scalar insertion with SSO-sized strings performs one repeatable map-node
  // allocation.  Targeting that size avoids faulting one-time CRT/locale
  // internals while reaching every candidate.set() call in the full importer.
  core::ParameterSet probe;
  Status probe_status{};
  {
    MeasureLargestAllocation measure;
    probe_status = probe.set("x", 1.0, "p");
  }
  REQUIRE(probe_status == Status::Success);
  const std::size_t parameter_node_bytes = largest_allocation_size;
  REQUIRE(parameter_node_bytes > 0);
  REQUIRE(largest_allocation_count == 1);

  core::ParameterSet failed_node;
  Status failed_node_status{};
  {
    FailAllocationOfSize failure{ parameter_node_bytes, 0, true };
    failed_node_status = failed_node.set("x", 1.0, "p");
  }
  REQUIRE(allocation_failure_triggered);
  CHECK(failed_node_status == Status::Numerical_failure);
  CHECK(failed_node.size() == 0);

  const std::string long_name(513, 'x');
  core::ParameterSet canonical_measured;
  Status canonical_measured_status{};
  auto measured_name = long_name;
  {
    MeasureLargestAllocation measure;
    canonical_measured_status = canonical_measured.set(
      std::move(measured_name), 2.0, "p");
  }
  REQUIRE(canonical_measured_status == Status::Success);
  const std::size_t canonical_name_bytes = largest_allocation_size;
  REQUIRE(canonical_name_bytes > parameter_node_bytes);
  REQUIRE(largest_allocation_count == 1);

  core::ParameterSet failed_canonical_name;
  auto failing_name = long_name;
  Status failed_canonical_status{};
  {
    FailAllocationOfSize failure{ canonical_name_bytes, 0, true };
    failed_canonical_status = failed_canonical_name.set(
      std::move(failing_name), 2.0, "p");
  }
  REQUIRE(allocation_failure_triggered);
  CHECK(failed_canonical_status == Status::Numerical_failure);
  CHECK(failed_canonical_name.size() == 0);

  core::ParameterSet measured;
  Status measured_status{};
  {
    MeasureAllocationsOfSize measure{ parameter_node_bytes };
    measured_status = core::ParameterSet::fromBpxJson(
      bpx_fixture, measured, diagnostic);
  }
  REQUIRE(measured_status == Status::Success);
  const std::size_t matching_allocations = measured_matching_allocations;
  REQUIRE(matching_allocations > 0);

  core::ParameterSet sentinel;
  REQUIRE(sentinel.set("sentinel", 7.0, "allocation-test")
          == Status::Success);
  for (std::size_t occurrence = 0; occurrence < matching_allocations;
       ++occurrence) {
    CAPTURE(occurrence, matching_allocations, parameter_node_bytes);
    core::ParameterSet output = sentinel;
    diagnostic.clear();
    Status status{};
    {
      FailAllocationOfSize failure{
        parameter_node_bytes, occurrence, true
      };
      status = core::ParameterSet::fromBpxJson(
        bpx_fixture, output, diagnostic);
    }
    REQUIRE(allocation_failure_triggered);
    CHECK(status == Status::Numerical_failure);
    REQUIRE(output.size() == 1);
    const auto *value = output.findScalar("sentinel");
    REQUIRE(value != nullptr);
    CHECK(*value == 7.0);
  }
}

TEST_CASE("netlist CSV maps a late cell-vector allocation failure atomically",
          "[core][pack][netlist][parser][allocation][P9]")
{
  std::string csv{ "desc,node1,node2,value\n" };
  for (std::size_t cell = 0; cell < 37; ++cell)
    csv += "V" + std::to_string(cell) + ",1,0,4.2\n";
  csv += "I0,1,0,5\n";

  core::CompiledPackTopology sentinel;
  REQUIRE(core::compilePackDescription(
            { .root = core::cell({ .archetype = "sentinel" }) }, sentinel)
          == Status::Success);
  core::CompiledPackTopology output = sentinel;
  core::NetlistCsvDiagnostic diagnostic;
  Status status{};
  {
    FailAllocationOfSize failure{ 37U * sizeof(core::CompiledCell), 0, true };
    status = core::parseLiionpackNetlistCsv(csv, output, diagnostic);
  }
  REQUIRE(allocation_failure_triggered);
  CHECK(status == Status::Numerical_failure);
  REQUIRE(output.cells.size() == 1);
  CHECK(output.cells[0].path == "c00");
  CHECK(output.cells[0].archetype == "sentinel");
  CHECK(output.electrical.branches.size() == 1);
}
