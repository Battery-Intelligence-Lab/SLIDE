/**
 * @file core_RecorderAllocation_test.cpp
 * @brief Phase-6 accepted snapshot allocation counter.
 */

#include "../../src/core/Recorder.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <cstdlib>
#include <new>
#include <vector>

static std::atomic<std::size_t> allocation_count{};
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

TEST_CASE("Phase-6 accepted snapshots allocate nothing",
          "[core][recorder][allocation][P6-G1]")
{
  core::SpmBatch batch;
  const auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);
  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .capacity = 2 }) == Status::Success);
  const std::array current{ 8.0, -4.0 };
  REQUIRE(recorder.record(0, current) == Status::Success);
  const auto before = allocation_count.load(std::memory_order_relaxed);
  const auto status = recorder.record(1, current);
  const auto after = allocation_count.load(std::memory_order_relaxed);
  REQUIRE(status == Status::Success);
  REQUIRE(after == before);
}

TEST_CASE("Recorder reconfiguration translates allocation failure atomically",
          "[core][recorder][allocation][coverage]")
{
  core::SpmBatch batch;
  const auto input = test_support::make_legacy_kokam_input(
    0.55, 298.0, 298.0);
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);

  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .capacity = 3 }) == Status::Success);
  const std::array current{ 8.0, -4.0 };
  REQUIRE(recorder.record(7, current) == Status::Success);
  const auto before = recorder.snapshot(0);
  const std::vector<double> before_current(before.current_density.begin(),
                                           before.current_density.end());
  const std::vector<double> before_state(before.state.begin(),
                                         before.state.end());

  matching_allocation_size = 2 * batch.state().size() * sizeof(double);
  matching_failure_triggered = false;
  fail_matching_allocation = true;
  const auto reconfigure_status =
    recorder.configure(batch, { .capacity = 2 });
  fail_matching_allocation = false;
  REQUIRE(matching_failure_triggered);
  CHECK(reconfigure_status == Status::Numerical_failure);
  REQUIRE(recorder.configured());
  CHECK(recorder.capacity() == 3);
  REQUIRE(recorder.size() == 1);
  const auto after = recorder.snapshot(0);
  CHECK(after.accepted_step == 7);
  const std::vector<double> after_current(after.current_density.begin(),
                                          after.current_density.end());
  const std::vector<double> after_state(after.state.begin(),
                                        after.state.end());
  CHECK(after_current == before_current);
  CHECK(after_state == before_state);
  std::array<double, 2> voltage{};
  CHECK(recorder.terminalVoltage(0, voltage) == Status::Success);

  REQUIRE(recorder.configure(batch, { .capacity = 2 }) == Status::Success);
  CHECK(recorder.capacity() == 2);
  CHECK(recorder.size() == 0);
}
