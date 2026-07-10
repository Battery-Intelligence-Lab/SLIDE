/**
 * @file core_P1G2_allocation_test.cpp
 * @brief P1-G2: 10,000-cell accepted step performs exactly zero heap allocations.
 */

#include "../../src/core/EulerLegacy.hpp"

#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <cstdlib>
#include <new>
#include <vector>

static std::atomic<std::size_t> allocation_count{ 0 };

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
  ++allocation_count;
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes)
{
  ++allocation_count;
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new(std::size_t bytes, std::align_val_t alignment)
{
  ++allocation_count;
  if (void *pointer = aligned_allocate(bytes != 0 ? bytes : 1,
                                       static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes, std::align_val_t alignment)
{
  ++allocation_count;
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
