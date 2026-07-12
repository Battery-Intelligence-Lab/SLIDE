/**
 * @file core_P1G2_allocation_test.cpp
 * @brief P1-G2: 10,000-cell accepted step performs exactly zero heap allocations.
 */

#include "../../src/core/EulerLegacy.hpp"
#include "../../src/core/ForwardSensitivity.hpp"
#include "../../src/core/Simulation.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <cstdlib>
#include <new>
#include <vector>

static std::atomic<std::size_t> allocation_count{ 0 };
static thread_local bool fail_matching_allocation{};
static thread_local bool matching_failure_triggered{};
static thread_local std::size_t matching_allocation_size{};

static void before_allocation(std::size_t bytes)
{
  ++allocation_count;
  if (fail_matching_allocation && !matching_failure_triggered
      && bytes == matching_allocation_size) {
    matching_failure_triggered = true;
    throw std::bad_alloc{};
  }
}

#if defined(_WIN32)
#include <malloc.h>
static void *aligned_allocate(std::size_t bytes, std::size_t alignment)
{
  return _aligned_malloc(bytes, alignment);
}
static void aligned_release(void *pointer) { _aligned_free(pointer); }
#else
static void *aligned_allocate(std::size_t bytes, std::size_t alignment)
{
  return std::aligned_alloc(alignment,
                            ((bytes + alignment - 1) / alignment) * alignment);
}
static void aligned_release(void *pointer) { std::free(pointer); }
#endif

void *operator new(std::size_t bytes)
{
  before_allocation(bytes);
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes)
{
  before_allocation(bytes);
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new(std::size_t bytes, std::align_val_t alignment)
{
  before_allocation(bytes);
  if (void *pointer = aligned_allocate(bytes != 0 ? bytes : 1,
                                       static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes, std::align_val_t alignment)
{
  before_allocation(bytes);
  if (void *pointer = aligned_allocate(bytes != 0 ? bytes : 1,
                                       static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void *pointer) noexcept { std::free(pointer); }
void operator delete[](void *pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::align_val_t) noexcept { aligned_release(pointer); }
void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept
{
  aligned_release(pointer);
}
void operator delete[](void *pointer, std::align_val_t) noexcept
{
  aligned_release(pointer);
}
void operator delete[](void *pointer, std::size_t, std::align_val_t) noexcept
{
  aligned_release(pointer);
}

using namespace slide;

namespace {

core::SpmFactoryInput make_input()
{
  core::SpmFactoryInput input;
  input.design.capacity_Ah = 1.0;
  input.design.electrode_area = 0.1;
  input.design.electrolyte.concentration = 1000.0;
  input.design.thermal.reference_temperature = 298.15;
  for (const core::Domain domain : core::domains) {
    const auto d = core::domain_index(domain);
    auto &electrode = input.design.electrode[d];
    electrode.thickness = domain == core::Domain::neg ? 75e-6 : 87e-6;
    electrode.porosity = 0.3;
    electrode.active_fraction = 0.5;
    electrode.particle_radius = domain == core::Domain::neg ? 12.5e-6 : 8.5e-6;
    auto &material = electrode.active_material;
    material.ocv.stoichiometry = { 0.0, 1.0 };
    material.ocv.value = domain == core::Domain::neg
                           ? std::vector<double>{ 0.0, 0.1 }
                           : std::vector<double>{ 3.0, 4.0 };
    material.cs_max = domain == core::Domain::neg ? 30'555.0 : 51'385.0;
    material.x_0 = domain == core::Domain::neg ? 0.2 : 0.8;
    material.x_100 = domain == core::Domain::neg ? 0.8 : 0.2;
    material.D_s = { .reference_value = domain == core::Domain::neg ? 7e-14 : 8e-14,
                     .activation_energy = 0.0,
                     .reference_temperature = 298.15 };
    material.k_ct = { .reference_value = 1e-11,
                      .activation_energy = 0.0,
                      .reference_temperature = 298.15 };
  }
  return input;
}

class FailAllocationOfSize
{
public:
  explicit FailAllocationOfSize(std::size_t bytes)
  {
    matching_allocation_size = bytes;
    matching_failure_triggered = false;
    fail_matching_allocation = true;
  }

  ~FailAllocationOfSize() { fail_matching_allocation = false; }
};

} // namespace

TEST_CASE("P1-G2 10000-lane Euler step allocates nothing", "[core][P1-G2]")
{
  constexpr int lanes = 10'000;
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(make_input(), {}, lanes, batch) == Status::Success);
  core::EulerLegacy stepper;
  REQUIRE(stepper.configure(batch) == Status::Success);
  std::vector<double> current_density(static_cast<std::size_t>(lanes), 0.0);

  const double state_bytes_per_cell = static_cast<double>(batch.state().size() * sizeof(double)) / lanes;
  REQUIRE(state_bytes_per_cell <= 300.0);

  // Warm any toolchain/runtime one-time paths before the decisive counter interval.
  REQUIRE(stepper.step(batch, current_density, 0.0, 1.0) == Status::Success);
  const std::size_t before = allocation_count.load(std::memory_order_relaxed);
  const Status status = stepper.step(batch, current_density, 1.0, 1.0);
  const std::size_t after = allocation_count.load(std::memory_order_relaxed);

  REQUIRE(status == Status::Success);
  REQUIRE(after == before);
}

TEST_CASE("solution builders translate late allocation failures atomically",
          "[core][allocation][simulation][sensitivity][coverage]")
{
  const auto input = make_input();
  constexpr std::size_t samples = 17;

  core::Simulation simulation;
  REQUIRE(simulation.build(input, {}, 1) == Status::Success);
  core::SimulationSolution simulation_output;
  simulation_output.time = { 42.0 };
  Status simulation_status{};
  {
    FailAllocationOfSize failure{ samples * sizeof(double) };
    simulation_status = simulation.solve(
      { .current_A = 0.0, .duration = 16.0, .step = 1.0 },
      simulation_output);
  }
  REQUIRE(matching_failure_triggered);
  CHECK(simulation_status == Status::Numerical_failure);
  CHECK(simulation_output.time == std::vector<double>{ 42.0 });

  constexpr std::array sensitivity_parameters{
    core::SensitivityParameter::nominal_capacity,
    core::SensitivityParameter::contact_resistance,
  };
  core::ForwardSensitivitySolution sensitivity_output;
  sensitivity_output.time = { 42.0 };
  Status sensitivity_status{};
  {
    FailAllocationOfSize failure{
      samples * sensitivity_parameters.size() * sizeof(double)
    };
    sensitivity_status = core::solveCcForwardSensitivities(
      input, 5, 1.0, false, core::Direction::discharge, 16.0, 1.0, sensitivity_parameters, sensitivity_output);
  }
  REQUIRE(matching_failure_triggered);
  CHECK(sensitivity_status == Status::Numerical_failure);
  CHECK(sensitivity_output.time == std::vector<double>{ 42.0 });
}
