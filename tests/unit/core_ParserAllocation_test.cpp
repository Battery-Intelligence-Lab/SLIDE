/**
 * @file core_ParserAllocation_test.cpp
 * @brief Deterministic allocation-failure atomicity for Phase-9 cold transactions.
 */

#include "../../src/core/Experiment.hpp"
#include "../../src/core/NetlistCsv.hpp"
#include "../../src/core/PackSolver.hpp"
#include "../../src/core/ParameterSet.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <new>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
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
thread_local bool fail_allocation_at_occurrence{};
thread_local std::size_t allocation_cursor{};
thread_local std::size_t failing_allocation{};

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
  if (fail_allocation_at_occurrence) {
    const std::size_t current = allocation_cursor++;
    if (current == failing_allocation) {
      allocation_failure_triggered = true;
      throw std::bad_alloc{};
    }
  }
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

class FailAllocationAtOccurrence
{
public:
  explicit FailAllocationAtOccurrence(std::size_t occurrence)
  {
    allocation_cursor = 0;
    failing_allocation = occurrence;
    allocation_failure_triggered = false;
    fail_allocation_at_occurrence = true;
  }
  ~FailAllocationAtOccurrence()
  {
    fail_allocation_at_occurrence = false;
  }
};

class ScopedTemporaryDirectory
{
public:
  ScopedTemporaryDirectory()
  {
    const auto nonce = static_cast<std::uint64_t>(
      std::chrono::steady_clock::now().time_since_epoch().count());
    const auto root = std::filesystem::temp_directory_path();
    for (std::uint64_t attempt = 0; attempt < 1024; ++attempt) {
      const auto candidate = root
                             / ("slide_parser_allocation_"
                                + std::to_string(nonce) + "_"
                                + std::to_string(attempt));
      std::error_code error;
      if (std::filesystem::create_directory(candidate, error)) {
        path_ = candidate;
        return;
      }
      if (error)
        throw std::runtime_error{ "cannot create parser test directory: "
                                  + error.message() };
    }
    throw std::runtime_error{ "cannot reserve a parser test directory" };
  }

  ~ScopedTemporaryDirectory()
  {
    std::error_code ignored;
    std::filesystem::remove_all(path_, ignored);
  }

  ScopedTemporaryDirectory(const ScopedTemporaryDirectory &) = delete;
  ScopedTemporaryDirectory &operator=(const ScopedTemporaryDirectory &) = delete;

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  std::filesystem::path path_{};
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

struct AllocationAffineBatch
{
  std::vector<double> ocv{};
  std::vector<double> resistance{};

  slide::Status linearizeThevenin(std::span<const double> current,
                                  std::span<double>
                                    output_ocv,
                                  std::span<double>
                                    output_resistance)
  {
    if (current.size() != ocv.size())
      return slide::Status::Invalid_parameters;
    std::copy(ocv.begin(), ocv.end(), output_ocv.begin());
    std::copy(resistance.begin(), resistance.end(), output_resistance.begin());
    return slide::Status::Success;
  }
};

slide::Status makeSpmInputSentinel(slide::core::SpmFactoryInput &sentinel)
{
  sentinel = {};
  double marker = 101.0;
  for (const slide::core::Domain domain : slide::core::domains) {
    const auto d = slide::core::domain_index(domain);
    auto &electrode = sentinel.design.electrode[d];
    electrode.active_material.ocv = {
      .stoichiometry = { 0.0, 1.0 },
      .value = { marker, marker + 1.0 },
    };
    electrode.active_material.cs_max = marker + 2.0;
    electrode.thickness = marker + 3.0;
    electrode.aging.push_back(
      { .kind = slide::core::AgingMechanismKind::semi_empirical,
        .name = d == 0 ? "negative sentinel" : "positive sentinel",
        .coefficients = { marker + 4.0 } });
    marker += 10.0;
  }
  sentinel.design.separator.thickness = 201.0;
  sentinel.design.electrolyte.concentration = 202.0;
  sentinel.design.thermal.density = 203.0;
  sentinel.design.capacity_Ah = 204.0;
  sentinel.design.electrode_area = 205.0;
  sentinel.total_entropic_coefficient = {
    .stoichiometry = { 0.0, 1.0 }, .value = { 206.0, 207.0 }
  };
  sentinel.negative_entropic_coefficient = {
    .stoichiometry = { 0.0, 1.0 }, .value = { 208.0, 209.0 }
  };
  sentinel.negative_laresgoiti_stress = {
    .stoichiometry = { 0.0, 1.0 }, .value = { 210.0, 211.0 }
  };
  sentinel.initial_soc = 212.0;
  sentinel.initial_temperature = 213.0;
  sentinel.initial_sei_thickness = 214.0;
  sentinel.initial_lost_lithium = 215.0;
  sentinel.initial_crack_surface_fraction = 216.0;
  sentinel.initial_plated_lithium_thickness = 217.0;
  sentinel.initial_specific_resistance = { 218.0, 219.0 };
  sentinel.initial_current_collector_resistance = 220.0;
  sentinel.initial_stress_interval = 221.0;
  sentinel.sei_resistivity_area = 222.0;
  sentinel.sei.F = 223.0;
  sentinel.surface_crack.F = 224.0;
  sentinel.lam.F = 225.0;
  sentinel.lithium_plating.F = 226.0;
  constexpr std::array x{ 0.0, 1.0 };
  constexpr std::array y{ 227.0, 228.0 };
  return sentinel.lam.positive_ocv.build(x, y);
}

bool sameSpmInputSentinel(const slide::core::SpmFactoryInput &actual,
                          const slide::core::SpmFactoryInput &expected)
{
  for (const slide::core::Domain domain : slide::core::domains) {
    const auto d = slide::core::domain_index(domain);
    const auto &a = actual.design.electrode[d];
    const auto &e = expected.design.electrode[d];
    if (a.active_material.ocv != e.active_material.ocv
        || a.active_material.cs_max != e.active_material.cs_max
        || a.thickness != e.thickness || a.aging.size() != 1
        || e.aging.size() != 1 || a.aging[0].kind != e.aging[0].kind
        || a.aging[0].name != e.aging[0].name
        || a.aging[0].coefficients != e.aging[0].coefficients)
      return false;
  }
  if (actual.design.separator.thickness
        != expected.design.separator.thickness
      || actual.design.electrolyte.concentration
           != expected.design.electrolyte.concentration
      || actual.design.thermal.density != expected.design.thermal.density
      || actual.design.capacity_Ah != expected.design.capacity_Ah
      || actual.design.electrode_area != expected.design.electrode_area
      || actual.total_entropic_coefficient
           != expected.total_entropic_coefficient
      || actual.negative_entropic_coefficient
           != expected.negative_entropic_coefficient
      || actual.negative_laresgoiti_stress
           != expected.negative_laresgoiti_stress
      || actual.initial_soc != expected.initial_soc
      || actual.initial_temperature != expected.initial_temperature
      || actual.initial_sei_thickness != expected.initial_sei_thickness
      || actual.initial_lost_lithium != expected.initial_lost_lithium
      || actual.initial_crack_surface_fraction
           != expected.initial_crack_surface_fraction
      || actual.initial_plated_lithium_thickness
           != expected.initial_plated_lithium_thickness
      || actual.initial_specific_resistance
           != expected.initial_specific_resistance
      || actual.initial_current_collector_resistance
           != expected.initial_current_collector_resistance
      || actual.initial_stress_interval != expected.initial_stress_interval
      || actual.sei_resistivity_area != expected.sei_resistivity_area
      || actual.sei.F != expected.sei.F
      || actual.surface_crack.F != expected.surface_crack.F
      || actual.lam.F != expected.lam.F
      || actual.lithium_plating.F != expected.lithium_plating.F
      || actual.lam.positive_ocv.valid()
           != expected.lam.positive_ocv.valid())
    return false;
  return !actual.lam.positive_ocv.valid()
         || (actual.lam.positive_ocv.knots()
               == expected.lam.positive_ocv.knots()
             && actual.lam.positive_ocv.x_min()
                  == expected.lam.positive_ocv.x_min()
             && actual.lam.positive_ocv.x_max()
                  == expected.lam.positive_ocv.x_max()
             && actual.lam.positive_ocv.eval(0.37)
                  == expected.lam.positive_ocv.eval(0.37));
}

bool sameParameterSentinel(const slide::core::ParameterSet &parameters)
{
  if (parameters.size() != 1)
    return false;
  const auto descriptions = parameters.describe();
  return descriptions.size() == 1 && descriptions[0].name == "sentinel"
         && descriptions[0].provenance == "allocation-test"
         && std::get<double>(descriptions[0].value) == 7.0;
}

slide::core::ExperimentSolution makeExperimentSolutionSentinel()
{
  return {
    .time = { 101.0, 102.0 },
    .voltage = { 201.0, 202.0 },
    .current = { 301.0, 302.0 },
    .sample_segment = { 7, 9 },
    .reason = slide::core::TerminationReason::error,
    .status = slide::Status::Unknown_problem,
    .segment = 11,
    .termination_name = "allocation sentinel",
  };
}

bool sameExperimentSolution(
  const slide::core::ExperimentSolution &left,
  const slide::core::ExperimentSolution &right)
{
  return left.time == right.time && left.voltage == right.voltage
         && left.current == right.current
         && left.sample_segment == right.sample_segment
         && left.reason == right.reason && left.status == right.status
         && left.segment == right.segment
         && left.termination_name == right.termination_name;
}

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

TEST_CASE("Chen2020 translates representative allocation failures atomically",
          "[core][parameters][Chen2020][allocation][9C-3]")
{
  core::ParameterSet warm;
  REQUIRE(core::ParameterSet::chen2020(warm) == Status::Success);

  core::ParameterSet measured;
  Status measured_status{};
  {
    MeasureLargestAllocation measure;
    measured_status = core::ParameterSet::chen2020(measured);
  }
  REQUIRE(measured_status == Status::Success);
  const std::size_t large_allocation = largest_allocation_size;
  const std::size_t large_occurrences = largest_allocation_count;
  REQUIRE(large_allocation > 0);
  REQUIRE(large_occurrences > 0);

  core::ParameterSet sentinel;
  REQUIRE(sentinel.set("sentinel", 7.0, "allocation-test")
          == Status::Success);
  const auto check_failure = [&](bool fail_by_size,
                                 std::size_t occurrence) {
    CAPTURE(fail_by_size, occurrence, large_allocation, large_occurrences);
    core::ParameterSet output = sentinel;
    Status status = Status::Unknown_problem;
    bool escaped{};
    try {
      if (fail_by_size) {
        FailAllocationOfSize injection{ large_allocation, occurrence };
        status = core::ParameterSet::chen2020(output);
      } else {
        FailAllocationAtOccurrence injection{ occurrence };
        status = core::ParameterSet::chen2020(output);
      }
    } catch (const std::bad_alloc &) {
      escaped = true;
    } catch (const std::length_error &) {
      escaped = true;
    }
    REQUIRE(allocation_failure_triggered);
    CHECK_FALSE(escaped);
    CHECK(status == Status::Numerical_failure);
    CHECK(sameParameterSentinel(output));
  };

  // The first allocation exercises caller-side construction before set() can
  // translate it. The largest early/late allocations exercise adaptive curve
  // storage without injecting into implementation-defined noexcept internals.
  check_failure(false, 0);
  check_failure(true, 0);
  if (large_occurrences > 1)
    check_failure(true, large_occurrences - 1);
}

TEST_CASE("SPM input compilation translates representative allocation failures atomically",
          "[core][parameters][factory][allocation][9C-3]")
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  core::SpmFactoryInput warm;
  REQUIRE(parameters.toSpmInput(warm) == Status::Success);

  core::SpmFactoryInput measured;
  Status measured_status{};
  {
    MeasureLargestAllocation measure;
    measured_status = parameters.toSpmInput(measured);
  }
  REQUIRE(measured_status == Status::Success);
  const std::size_t large_allocation = largest_allocation_size;
  const std::size_t large_occurrences = largest_allocation_count;
  REQUIRE(large_allocation > 0);
  REQUIRE(large_occurrences > 0);

  core::SpmFactoryInput sentinel;
  REQUIRE(makeSpmInputSentinel(sentinel) == Status::Success);
  const auto check_failure = [&](bool fail_by_size,
                                 std::size_t occurrence) {
    CAPTURE(fail_by_size, occurrence, large_allocation, large_occurrences);
    core::SpmFactoryInput output = sentinel;
    Status status = Status::Unknown_problem;
    bool escaped{};
    try {
      if (fail_by_size) {
        FailAllocationOfSize injection{ large_allocation, occurrence };
        status = parameters.toSpmInput(output);
      } else {
        FailAllocationAtOccurrence injection{ occurrence };
        status = parameters.toSpmInput(output);
      }
    } catch (const std::bad_alloc &) {
      escaped = true;
    } catch (const std::length_error &) {
      escaped = true;
    }
    REQUIRE(allocation_failure_triggered);
    CHECK_FALSE(escaped);
    CHECK(status == Status::Numerical_failure);
    CHECK(sameSpmInputSentinel(output, sentinel));
  };

  check_failure(false, 0);
  check_failure(true, 0);
  if (large_occurrences > 1)
    check_failure(true, large_occurrences - 1);
}

TEST_CASE("drive-cycle registration owns its copy inside the Status transaction",
          "[core][experiment][allocation][drive-cycle][9C-3]")
{
  core::DriveCycle candidate;
  candidate.name.assign(1024, 'c');
  candidate.time.resize(37);
  candidate.current.resize(37);
  for (std::size_t i = 0; i < candidate.time.size(); ++i) {
    candidate.time[i] = static_cast<double>(i);
    candidate.current[i] = 0.25 * static_cast<double>(i);
  }

  core::CyclerV2 name_probe;
  Status name_probe_status{};
  {
    MeasureLargestAllocation measure;
    name_probe_status = name_probe.registerDriveCycle(candidate);
  }
  REQUIRE(name_probe_status == Status::Success);
  const std::size_t name_allocation = largest_allocation_size;
  const std::size_t name_occurrences = largest_allocation_count;
  REQUIRE(name_allocation > candidate.time.size() * sizeof(double));
  REQUIRE(name_occurrences > 0);

  const std::size_t table_allocation = candidate.time.size() * sizeof(double);
  core::CyclerV2 table_probe;
  Status table_probe_status{};
  {
    MeasureAllocationsOfSize measure{ table_allocation };
    table_probe_status = table_probe.registerDriveCycle(candidate);
  }
  REQUIRE(table_probe_status == Status::Success);
  const std::size_t table_occurrences = measured_matching_allocations;
  REQUIRE(table_occurrences >= 2);

  const core::DriveCycle retained{
    .name = "retained", .time = { 0.0, 1.0 }, .current = { 0.0, 1.0 }
  };
  const auto exercise = [&](std::size_t bytes, std::size_t occurrence) {
    CAPTURE(bytes, occurrence, name_occurrences, table_occurrences);
    core::CyclerV2 cycler;
    REQUIRE(cycler.registerDriveCycle(retained) == Status::Success);
    Status status{};
    {
      FailAllocationOfSize failure{ bytes, occurrence };
      status = cycler.registerDriveCycle(candidate);
    }
    REQUIRE(allocation_failure_triggered);
    CHECK(status == Status::Numerical_failure);
    CHECK(cycler.registerDriveCycle(retained)
          == Status::Invalid_parameters);
    CHECK(cycler.registerDriveCycle(candidate) == Status::Success);
  };
  for (std::size_t occurrence = 0; occurrence < name_occurrences;
       ++occurrence)
    exercise(name_allocation, occurrence);
  for (std::size_t occurrence = 0; occurrence < table_occurrences;
       ++occurrence)
    exercise(table_allocation, occurrence);
}

TEST_CASE("Cycler run allocation failure restores both arenas and output",
          "[core][experiment][allocation][run][9C-3]")
{
  const auto input = test_support::make_legacy_kokam_input(
    0.55, settings::T_ENV, 298.0);
  core::Experiment rest;
  rest.segments.push_back(
    { .mode = core::ControlMode::rest, .duration = 1.0 });

  core::SpmBatch measured_batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, measured_batch)
          == Status::Success);
  core::CyclerV2 measured_cycler;
  REQUIRE(measured_cycler.configure(measured_batch) == Status::Success);
  core::ExperimentSolution warm_output;
  REQUIRE(measured_cycler.run(rest, 1.0, warm_output) == Status::Success);
  core::ExperimentSolution measured_output;
  Status measured_status{};
  constexpr std::size_t second_sample_bytes = 2 * sizeof(double);
  {
    MeasureAllocationsOfSize measure{ second_sample_bytes };
    measured_status = measured_cycler.run(rest, 1.0, measured_output);
  }
  REQUIRE(measured_status == Status::Success);
  const std::size_t sample_allocations = measured_matching_allocations;
  REQUIRE(sample_allocations >= 4);
  const std::size_t first_post_step_allocation = sample_allocations - 4;

  // The final four matching allocations are the time/voltage/current/segment
  // growth after advance(). Earlier 16-byte requests are Debug STL proxies.
  for (std::size_t occurrence = first_post_step_allocation;
       occurrence < sample_allocations;
       ++occurrence) {
    CAPTURE(occurrence, sample_allocations);
    core::SpmBatch batch;
    REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);
    core::CyclerV2 cycler;
    REQUIRE(cycler.configure(batch) == Status::Success);
    core::ExperimentSolution target_warm_output;
    REQUIRE(cycler.run(rest, 1.0, target_warm_output) == Status::Success);
    const std::vector<double> state_before(
      batch.state().raw().begin(), batch.state().raw().end());
    const std::vector<double> derivative_before(
      batch.derivative().raw().begin(), batch.derivative().raw().end());
    const auto sentinel = makeExperimentSolutionSentinel();
    auto output = sentinel;
    Status status{};
    {
      FailAllocationOfSize failure{
        second_sample_bytes, occurrence, true
      };
      status = cycler.run(rest, 1.0, output);
    }
    REQUIRE(allocation_failure_triggered);
    CHECK(status == Status::Numerical_failure);
    CHECK(std::equal(state_before.begin(), state_before.end(), batch.state().raw().begin()));
    CHECK(std::equal(derivative_before.begin(), derivative_before.end(), batch.derivative().raw().begin()));
    CHECK(sameExperimentSolution(output, sentinel));
  }
}

TEST_CASE("callback allocation failure restores the pre-run transaction",
          "[core][experiment][allocation][callback][9C-3]")
{
  const auto input = test_support::make_legacy_kokam_input(
    0.55, settings::T_ENV, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  const std::vector<double> state_before(
    batch.state().raw().begin(), batch.state().raw().end());
  const std::vector<double> derivative_before(
    batch.derivative().raw().begin(), batch.derivative().raw().end());

  core::Experiment experiment;
  experiment.segments = {
    { .mode = core::ControlMode::rest, .duration = 1.0 },
    { .mode = core::ControlMode::custom_explicit,
      .duration = 1.0,
      .custom_control = [](const core::ExperimentVariables &) {
        return 0.0;
      },
      .custom_terminations = {
        { .name = "post-advance allocation", .indicator = [](const core::ExperimentVariables &variables) {
           if (variables.local_time > 0.0)
             throw std::bad_alloc{};
           return 1.0;
         } },
      } },
  };
  const auto sentinel = makeExperimentSolutionSentinel();
  auto output = sentinel;
  CHECK(cycler.run(experiment, 1.0, output)
        == Status::Numerical_failure);
  CHECK(std::equal(state_before.begin(), state_before.end(), batch.state().raw().begin()));
  CHECK(std::equal(derivative_before.begin(), derivative_before.end(), batch.derivative().raw().begin()));
  CHECK(sameExperimentSolution(output, sentinel));

  core::Experiment retry;
  retry.segments.push_back(
    { .mode = core::ControlMode::rest, .duration = 1.0 });
  core::ExperimentSolution retry_output;
  CHECK(cycler.run(retry, 1.0, retry_output) == Status::Success);
}

TEST_CASE("missing drive cycles are rejected before any prior segment advances",
          "[core][experiment][drive-cycle][transaction][9C-3]")
{
  const auto input = test_support::make_legacy_kokam_input(
    0.55, settings::T_ENV, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  const std::vector<double> state_before(
    batch.state().raw().begin(), batch.state().raw().end());
  const std::vector<double> derivative_before(
    batch.derivative().raw().begin(), batch.derivative().raw().end());
  core::Experiment experiment;
  experiment.segments = {
    { .mode = core::ControlMode::rest, .duration = 1.0 },
    { .mode = core::ControlMode::drive_cycle,
      .drive_cycle = "missing" },
  };
  const auto sentinel = makeExperimentSolutionSentinel();
  auto output = sentinel;
  CHECK(cycler.run(experiment, 1.0, output)
        == Status::Invalid_parameters);
  CHECK(std::equal(state_before.begin(), state_before.end(), batch.state().raw().begin()));
  CHECK(std::equal(derivative_before.begin(), derivative_before.end(), batch.derivative().raw().begin()));
  CHECK(sameExperimentSolution(output, sentinel));
}

TEST_CASE("active Cycler callbacks cannot invalidate their own transaction",
          "[core][experiment][callback][reentrant][9C-3]")
{
  const auto input = test_support::make_legacy_kokam_input(
    0.55, settings::T_ENV, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);
  core::CyclerV2 cycler;
  REQUIRE(cycler.configure(batch) == Status::Success);
  const core::DriveCycle candidate{
    .name = "callback cycle",
    .time = { 0.0, 1.0 },
    .current = { 0.0, 1.0 },
  };
  Status register_status = Status::Unknown_problem;
  Status configure_status = Status::Unknown_problem;
  Status nested_status = Status::Unknown_problem;
  core::Experiment experiment;
  experiment.segments.push_back(
    { .mode = core::ControlMode::custom_explicit,
      .duration = 1.0,
      .custom_control = [&](const core::ExperimentVariables &) {
        register_status = cycler.registerDriveCycle(candidate);
        configure_status = cycler.configure(batch);
        core::Experiment nested;
        nested.segments.push_back(
          { .mode = core::ControlMode::rest, .duration = 1.0 });
        core::ExperimentSolution nested_output;
        nested_status = cycler.run(nested, 1.0, nested_output);
        return 0.0;
      } });
  core::ExperimentSolution output;
  REQUIRE(cycler.run(experiment, 1.0, output) == Status::Success);
  CHECK(register_status == Status::Invalid_parameters);
  CHECK(configure_status == Status::Invalid_parameters);
  CHECK(nested_status == Status::Invalid_parameters);
  CHECK(cycler.registerDriveCycle(candidate) == Status::Success);
}

TEST_CASE("Cycler reconfiguration publishes allocation-heavy scratch atomically",
          "[core][experiment][allocation][coverage]")
{
  const auto input = test_support::make_legacy_kokam_input(
    0.55, settings::T_ENV, 298.0);
  core::SpmBatch old_batch;
  REQUIRE(core::buildSpmBatch(input, {}, 1, old_batch) == Status::Success);
  core::SpmModelOptions new_options;
  new_options.nch = 8;
  core::SpmBatch new_batch;
  REQUIRE(core::buildSpmBatch(input, new_options, 1, new_batch)
          == Status::Success);
  REQUIRE(new_batch.state().size() != old_batch.state().size());

  const std::size_t state_bytes = new_batch.state().size() * sizeof(double);
  Status probe_status{};
  {
    core::CyclerV2 probe;
    MeasureAllocationsOfSize measure{ state_bytes };
    probe_status =
      probe.configure(new_batch, core::CyclerIntegrator::exponential);
  }
  REQUIRE(probe_status == Status::Success);
  REQUIRE(measured_matching_allocations > 0);

  core::Experiment rest;
  rest.segments.push_back(
    { .mode = core::ControlMode::rest, .duration = 1.0 });
  for (std::size_t occurrence = 0;
       occurrence < measured_matching_allocations;
       ++occurrence) {
    core::CyclerV2 cycler;
    REQUIRE(cycler.configure(old_batch, core::CyclerIntegrator::euler_legacy)
            == Status::Success);
    const std::vector<double> old_before(old_batch.state().raw().begin(),
                                         old_batch.state().raw().end());
    const std::vector<double> new_before(new_batch.state().raw().begin(),
                                         new_batch.state().raw().end());
    Status status{};
    {
      FailAllocationOfSize failure{ state_bytes, occurrence };
      status = cycler.configure(new_batch,
                                core::CyclerIntegrator::exponential);
    }
    CAPTURE(occurrence, measured_matching_allocations);
    REQUIRE(allocation_failure_triggered);
    CHECK(status == Status::Numerical_failure);
    CHECK(std::equal(old_before.begin(), old_before.end(), old_batch.state().raw().begin()));
    CHECK(std::equal(new_before.begin(), new_before.end(), new_batch.state().raw().begin()));

    core::ExperimentSolution output;
    CHECK(cycler.run(rest, 1.0, output) == Status::Success);
    CHECK_FALSE(std::equal(old_before.begin(), old_before.end(), old_batch.state().raw().begin()));
    CHECK(std::equal(new_before.begin(), new_before.end(), new_batch.state().raw().begin()));
  }
}

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

  core::ParameterSet update_target;
  REQUIRE(update_target.set("sentinel", 7.0, "allocation-test")
          == Status::Success);
  const std::array updates{
    core::ParameterDescription{ .name = "new value",
                                .value = 2.0,
                                .provenance = "allocation-test" }
  };
  Status update_status{};
  {
    FailAllocationOfSize failure{ parameter_node_bytes, 0, true };
    update_status = update_target.update(updates);
  }
  REQUIRE(allocation_failure_triggered);
  CHECK(update_status == Status::Numerical_failure);
  CHECK(update_target.size() == 1);
  CHECK(update_target.findScalar("sentinel") != nullptr);
  CHECK_FALSE(update_target.contains("new value"));

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

TEST_CASE("bounded parser files translate reader allocation failures atomically",
          "[core][parser][file][allocation][P9]")
{
  constexpr std::size_t reader_reserve = 8192;
  std::size_t reader_allocation{};
  {
    std::string probe;
    MeasureLargestAllocation measure;
    probe.reserve(reader_reserve);
    reader_allocation = largest_allocation_size;
  }
  REQUIRE(reader_allocation >= reader_reserve);

  const ScopedTemporaryDirectory temporary;
  const auto bpx_path = temporary.path() / "fixture.json";
  const auto csv_path = temporary.path() / "fixture.csv";
  std::error_code ignored;
  std::filesystem::remove(bpx_path, ignored);
  std::filesystem::remove(csv_path, ignored);
  {
    std::ofstream output(bpx_path, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.write(bpx_fixture.data(),
                 static_cast<std::streamsize>(bpx_fixture.size()));
    REQUIRE(output.good());
  }
  const std::string csv =
    "desc,node1,node2,value\nV0,1,0,4.2\nI0,1,0,1\n";
  {
    std::ofstream output(csv_path, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.write(csv.data(), static_cast<std::streamsize>(csv.size()));
    REQUIRE(output.good());
  }

  core::ParameterSet warm_bpx;
  std::string bpx_diagnostic;
  REQUIRE(core::ParameterSet::fromBpxFile(
            bpx_path, warm_bpx, bpx_diagnostic)
          == Status::Success);
  core::CompiledPackTopology warm_topology;
  core::NetlistCsvDiagnostic csv_diagnostic;
  REQUIRE(core::loadLiionpackNetlistCsv(
            csv_path, warm_topology, csv_diagnostic)
          == Status::Success);

  core::ParameterSet bpx_sentinel;
  REQUIRE(bpx_sentinel.set("sentinel", 7.0, "allocation-test")
          == Status::Success);
  Status bpx_status{};
  {
    FailAllocationOfSize failure{ reader_allocation, 0, true };
    bpx_status = core::ParameterSet::fromBpxFile(
      bpx_path, bpx_sentinel, bpx_diagnostic);
  }
  REQUIRE(allocation_failure_triggered);
  CHECK(bpx_status == Status::Numerical_failure);
  REQUIRE(bpx_sentinel.size() == 1);
  CHECK(bpx_sentinel.findScalar("sentinel") != nullptr);

  core::CompiledPackTopology topology_sentinel;
  REQUIRE(core::compilePackDescription(
            { .root = core::cell({ .archetype = "sentinel" }) },
            topology_sentinel)
          == Status::Success);
  Status csv_status{};
  {
    FailAllocationOfSize failure{ reader_allocation, 0, true };
    csv_status = core::loadLiionpackNetlistCsv(
      csv_path, topology_sentinel, csv_diagnostic);
  }
  REQUIRE(allocation_failure_triggered);
  CHECK(csv_status == Status::Numerical_failure);
  REQUIRE(topology_sentinel.cells.size() == 1);
  CHECK(topology_sentinel.cells[0].archetype == "sentinel");

  std::filesystem::remove(bpx_path, ignored);
  std::filesystem::remove(csv_path, ignored);
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

TEST_CASE("PackSolver late allocation failure preserves the prior configuration",
          "[core][pack][solver][allocation][P9-B37]")
{
  core::CompiledPackTopology sentinel_topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::cell({ .archetype = "sentinel" }) },
            sentinel_topology)
          == Status::Success);
  AllocationAffineBatch sentinel_batch{ .ocv = { 4.0 },
                                        .resistance = { 0.1 } };
  const std::array<core::TheveninBatchView, 1> sentinel_view{
    core::TheveninBatchView::bind(sentinel_batch, 1)
  };

  constexpr std::size_t lanes_per_batch = 37;
  std::vector<core::PackNode> cells;
  cells.reserve(2 * lanes_per_batch);
  for (std::size_t index = 0; index < 2 * lanes_per_batch; ++index)
    cells.push_back(core::cell(
      { .archetype = index % 2 == 0 ? "a" : "b" }));
  core::CompiledPackTopology target_topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::parallel(std::move(cells)) }, target_topology)
          == Status::Success);
  AllocationAffineBatch a{
    .ocv = std::vector<double>(lanes_per_batch, 4.0),
    .resistance = std::vector<double>(lanes_per_batch, 0.1)
  };
  AllocationAffineBatch b = a;
  const std::array<core::TheveninBatchView, 2> target_views{
    core::TheveninBatchView::bind(a, static_cast<int>(lanes_per_batch)),
    core::TheveninBatchView::bind(b, static_cast<int>(lanes_per_batch))
  };

  constexpr std::size_t target_cell_bytes =
    2 * lanes_per_batch * sizeof(double);
  // The first cell-sized allocation builds the candidate solution. The second
  // reaches later candidate scratch, after the old implementation had already
  // published topology, executor, and workspace members.
  constexpr std::size_t late_cell_allocation = 1;
  core::PackSolver measured;
  REQUIRE(measured.configure(sentinel_topology, sentinel_view, 1)
          == Status::Success);
  Status measured_status{};
  {
    MeasureAllocationsOfSize measure{ target_cell_bytes };
    measured_status = measured.configure(target_topology, target_views, 2);
  }
  REQUIRE(measured_status == Status::Success);
  REQUIRE(measured_matching_allocations > 1);

  core::PackSolver solver;
  REQUIRE(solver.configure(sentinel_topology, sentinel_view, 1)
          == Status::Success);
  REQUIRE(solver.solve(1.0, core::PackSolveMode::ladder)
          == Status::Success);
  const auto expected_solution = solver.solution();
  const auto expected_diagnostics = solver.diagnostics();
  const auto expected_workspace_valid = solver.workspace().valid();
  const auto expected_workspace_age = solver.workspace().age();
  const auto expected_numeric_factorizations =
    solver.workspace().numericFactorizations();
  const auto expected_symbolic_factorizations =
    solver.workspace().symbolicFactorizations();
  bool threw{};
  Status status = Status::Unknown_problem;
  {
    FailAllocationOfSize failure{ target_cell_bytes, late_cell_allocation };
    try {
      status = solver.configure(target_topology, target_views, 2);
    } catch (const std::bad_alloc &) {
      threw = true;
    }
  }
  REQUIRE(allocation_failure_triggered);
  CHECK_FALSE(threw);
  CHECK(status == Status::Numerical_failure);
  CHECK(solver.batchWorkerCount() == 1);
  CHECK(solver.solution().cell_current == expected_solution.cell_current);
  CHECK(solver.solution().node_voltage == expected_solution.node_voltage);
  CHECK(solver.solution().terminal_voltage == expected_solution.terminal_voltage);
  CHECK(solver.diagnostics().iterations == expected_diagnostics.iterations);
  CHECK(solver.diagnostics().numeric_factorizations
        == expected_diagnostics.numeric_factorizations);
  CHECK(solver.diagnostics().symbolic_factorizations
        == expected_diagnostics.symbolic_factorizations);
  CHECK(solver.diagnostics().jacobian_refreshes
        == expected_diagnostics.jacobian_refreshes);
  CHECK(solver.diagnostics().source_steps
        == expected_diagnostics.source_steps);
  CHECK(solver.diagnostics().residual_norm
        == expected_diagnostics.residual_norm);
  CHECK(solver.diagnostics().constraint_drift
        == expected_diagnostics.constraint_drift);
  CHECK(solver.diagnostics().constraint_bound
        == expected_diagnostics.constraint_bound);
  CHECK(solver.diagnostics().relaxation_gain
        == expected_diagnostics.relaxation_gain);
  CHECK(solver.workspace().valid() == expected_workspace_valid);
  CHECK(solver.workspace().age() == expected_workspace_age);
  CHECK(solver.workspace().numericFactorizations()
        == expected_numeric_factorizations);
  CHECK(solver.workspace().symbolicFactorizations()
        == expected_symbolic_factorizations);

  REQUIRE(solver.solve(1.0, core::PackSolveMode::ladder)
          == Status::Success);
  CHECK(solver.solution().cell_current == expected_solution.cell_current);
  CHECK(solver.solution().node_voltage == expected_solution.node_voltage);
  CHECK(solver.solution().terminal_voltage == expected_solution.terminal_voltage);

  REQUIRE(solver.configure(sentinel_topology, sentinel_view, 1)
          == Status::Success);
  REQUIRE(solver.solve(1.0, core::PackSolveMode::ladder)
          == Status::Success);
}
