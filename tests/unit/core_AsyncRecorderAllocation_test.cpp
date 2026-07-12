/**
 * @file core_AsyncRecorderAllocation_test.cpp
 * @brief P8-G3 simulation-thread allocation gate for async enqueue.
 */

#include "../../src/core/AsyncRecorder.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <new>
#include <string>

static std::atomic<std::size_t> producer_allocations{};
static thread_local bool count_producer_allocations{};
static thread_local bool measure_allocation_size{};
static thread_local std::size_t largest_measured_allocation{};
static thread_local bool fail_matching_allocation{};
static thread_local std::size_t matching_allocation_size{};
static thread_local bool matching_failure_triggered{};

static void before_allocation(std::size_t bytes)
{
  if (measure_allocation_size && bytes > largest_measured_allocation)
    largest_measured_allocation = bytes;
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
  if (count_producer_allocations)
    ++producer_allocations;
  before_allocation(bytes);
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes)
{
  if (count_producer_allocations)
    ++producer_allocations;
  before_allocation(bytes);
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new(std::size_t bytes, std::align_val_t alignment)
{
  if (count_producer_allocations)
    ++producer_allocations;
  before_allocation(bytes);
  if (void *pointer = aligned_allocate(bytes != 0 ? bytes : 1,
                                       static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes, std::align_val_t alignment)
{
  if (count_producer_allocations)
    ++producer_allocations;
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

class FailAllocationOfSize
{
public:
  explicit FailAllocationOfSize(std::size_t bytes)
  {
    matching_allocation_size = bytes;
    matching_failure_triggered = false;
    fail_matching_allocation = true;
  }
  ~FailAllocationOfSize()
  {
    fail_matching_allocation = false;
  }
};

} // namespace

TEST_CASE("Async configure rolls back a path-copy allocation failure",
          "[core][async-recorder][allocation][rollback][P9]")
{
  core::SpmBatch batch;
  const auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);
  const auto path = std::filesystem::temp_directory_path()
                    / ("slide_async_fault_" + std::string(128, 'x')
                       + ".slcmp");
  std::error_code ignored;
  std::filesystem::remove(path, ignored);

  // Keep fault injection out of one-time standard-library locale/thread setup;
  // the contract under test starts at AsyncRecorder's repeatable allocations.
  {
    core::AsyncRecorder warmup;
    REQUIRE(warmup.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    REQUIRE(warmup.finish() == Status::Success);
  }
  std::filesystem::remove(path, ignored);

  std::filesystem::path path_copy;
  largest_measured_allocation = 0;
  measure_allocation_size = true;
  path_copy = path;
  measure_allocation_size = false;
  REQUIRE(path_copy == path);
  REQUIRE(largest_measured_allocation > 0);

  core::AsyncRecorder recorder;
  Status status{};
  {
    FailAllocationOfSize failure{ largest_measured_allocation };
    status = recorder.configure(
      batch,
      path,
      { .ring_slots = 3,
        .backpressure = core::AsyncBackpressurePolicy::block,
        .codec = core::CompressionCodec::none });
  }
  REQUIRE(matching_failure_triggered);
  CHECK(status == Status::Numerical_failure);
  CHECK_FALSE(recorder.configured());

  if (recorder.configured())
    (void)recorder.finish();
  else {
    REQUIRE(recorder.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    REQUIRE(recorder.finish() == Status::Success);
  }
  std::filesystem::remove(path, ignored);
}

TEST_CASE("Compressed reader translates allocation failure",
          "[core][async-recorder][reader][allocation][coverage]")
{
  core::SpmBatch batch;
  const auto input = test_support::make_legacy_kokam_input(
    0.55, 298.0, 298.0);
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);
  const auto path = std::filesystem::temp_directory_path()
                    / "slide_async_reader_allocation.slcmp";
  std::error_code ignored;
  std::filesystem::remove(path, ignored);

  core::AsyncRecorder writer;
  REQUIRE(writer.configure(
            batch,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Success);
  const std::array current{ 1.0, -1.0 };
  REQUIRE(writer.enqueue(0, current) == Status::Success);
  REQUIRE(writer.finish() == Status::Success);

  const std::size_t raw_bytes =
    (batch.state().size() + static_cast<std::size_t>(batch.n_lanes()))
    * sizeof(double);
  REQUIRE(raw_bytes > 0);

  core::CompressedRecording allocation_failure;
  Status status{};
  {
    FailAllocationOfSize failure{ raw_bytes };
    status = allocation_failure.open(path);
  }
  REQUIRE(matching_failure_triggered);
  CHECK(status == Status::Numerical_failure);
  CHECK_FALSE(allocation_failure.valid());

  core::CompressedRecording valid;
  REQUIRE(valid.open(path) == Status::Success);
  CHECK(valid.valid());
  std::filesystem::remove(path, ignored);
}

TEST_CASE("P8-G3 async enqueue allocates nothing on the simulation thread",
          "[core][async-recorder][allocation][P8-G3]")
{
  core::SpmBatch batch;
  const auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);
  const auto path = std::filesystem::temp_directory_path()
                    / "slide_async_allocation.slcmp";
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  core::AsyncRecorder recorder;
  REQUIRE(recorder.configure(
            batch,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Success);
  const std::array current{ 8.0, -4.0 };
  producer_allocations.store(0, std::memory_order_relaxed);
  count_producer_allocations = true;
  const auto status = recorder.enqueue(0, current);
  count_producer_allocations = false;
  CHECK(status == Status::Success);
  CHECK(producer_allocations.load(std::memory_order_relaxed) == 0);
  REQUIRE(recorder.finish() == Status::Success);
  std::filesystem::remove(path, ignored);
}
