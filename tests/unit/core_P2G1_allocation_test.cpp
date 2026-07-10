/**
 * @file core_P2G1_allocation_test.cpp
 * @brief D-21/P2-G1 accepted coupled pack step performs zero heap allocations.
 */

#include "../../src/core/PackStepper.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <cstdlib>
#include <new>

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
void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept { aligned_release(pointer); }
void operator delete[](void *pointer, std::align_val_t) noexcept { aligned_release(pointer); }
void operator delete[](void *pointer, std::size_t, std::align_val_t) noexcept { aligned_release(pointer); }

using namespace slide;

TEST_CASE("P2-G1 coupled thermal pack step allocates nothing",
          "[core][pack][thermal][allocation][P2-G1]")
{
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  input.design.thermal.density = 1626.0;
  input.design.thermal.heat_capacity = 750.0;
  input.design.thermal.volume = 1e-4;
  input.design.thermal.surface_area = 0.0;
  input.design.thermal.h_conv = 0.0;
  input.design.thermal.environment_temperature = 298.0;
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, { .nch = 5, .thermal = true }, 2, batch)
          == Status::Success);
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::parallel(2, core::cell({ .archetype = "spm", .thermal = true })),
              .thermal_links = { { "p00", "p01", 2.0 } } },
            topology)
          == Status::Success);
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper stepper;
  REQUIRE(stepper.configure(topology, batches) == Status::Success);
  REQUIRE(stepper.step(20.0, 0.0, 0.1) == Status::Success);

  const auto before = allocation_count.load(std::memory_order_relaxed);
  const auto status = stepper.step(20.0, 0.1, 0.1);
  const auto after = allocation_count.load(std::memory_order_relaxed);
  REQUIRE(status == Status::Success);
  REQUIRE(after == before);
}
