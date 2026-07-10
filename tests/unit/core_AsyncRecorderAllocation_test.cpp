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

static std::atomic<std::size_t> producer_allocations{};
static thread_local bool count_producer_allocations{};

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
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes)
{
  if (count_producer_allocations)
    ++producer_allocations;
  if (void *pointer = std::malloc(bytes != 0 ? bytes : 1))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new(std::size_t bytes, std::align_val_t alignment)
{
  if (count_producer_allocations)
    ++producer_allocations;
  if (void *pointer = aligned_allocate(bytes != 0 ? bytes : 1,
                                       static_cast<std::size_t>(alignment)))
    return pointer;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes, std::align_val_t alignment)
{
  if (count_producer_allocations)
    ++producer_allocations;
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
