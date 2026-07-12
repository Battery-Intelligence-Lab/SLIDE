/**
 * @file core_ThreadPool_test.cpp
 * @brief P8-G4 persistent pool, failure, determinism, and diagnostic gates.
 */

#include "../../src/core/ThreadPool.hpp"
#include "../../src/utility/parallelisation.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <vector>

using namespace slide;

namespace slide::core::detail {

struct BatchExecutorTestAccess
{
  static void failPoolCreation(BatchExecutor &executor)
  {
    executor.pool_factory_ = [](unsigned) -> std::unique_ptr<ThreadPool> {
      throw std::bad_alloc{};
    };
  }
};

} // namespace slide::core::detail

TEST_CASE("P8-G4 persistent pool executes each task exactly once and propagates failure",
          "[core][thread-pool][P8-G4]")
{
  core::ThreadPool pool{ 4 };
  REQUIRE(pool.workerCount() == 4);
  std::array<std::atomic<int>, 257> visits{};
  const auto status = pool.parallelFor(visits.size(), [&](std::size_t index) {
    visits[index].fetch_add(1, std::memory_order_relaxed);
    return index == 113 ? Status::Invalid_states : Status::Success;
  });
  CHECK(status == Status::Invalid_states);
  for (const auto &visit : visits)
    CHECK(visit.load(std::memory_order_relaxed) == 1);

  std::vector<double> output(257);
  REQUIRE(pool.parallelFor(output.size(), [&](std::size_t index) {
    output[index] = static_cast<double>(index * index);
  }) == Status::Success);
  for (std::size_t i = 0; i < output.size(); ++i)
    CHECK(output[i] == static_cast<double>(i * i));
}

TEST_CASE("P8-G4 fixed-order reductions are worker-count bit-repeatable",
          "[core][thread-pool][P8-G4]")
{
  std::vector<double> values(4099);
  for (std::size_t i = 0; i < values.size(); ++i)
    values[i] = i % 2 == 0 ? 1.0 / static_cast<double>(i + 1)
                           : -1.0 / static_cast<double>(i + 1);
  const double reference = core::fixedOrderSum(values);
  for (const unsigned workers : { 1U, 2U, 7U }) {
    core::ThreadPool pool{ workers };
    std::vector<double> copied(values.size());
    REQUIRE(pool.parallelFor(values.size(), [&](std::size_t i) {
      copied[i] = values[i];
    }) == Status::Success);
    CHECK(std::bit_cast<std::uint64_t>(core::fixedOrderSum(copied))
          == std::bit_cast<std::uint64_t>(reference));
  }

  const std::array cancellation{ 0x1p54, -0x1p54, 1.0 };
  CHECK(std::bit_cast<std::uint64_t>(core::fixedOrderSum(cancellation))
        == std::bit_cast<std::uint64_t>(1.0));

  std::vector<double> cancellation_stress;
  cancellation_stress.reserve(3U * 1024U);
  for (int repetition = 0; repetition < 1024; ++repetition) {
    cancellation_stress.push_back(0x1p54);
    cancellation_stress.push_back(1.0);
    cancellation_stress.push_back(-0x1p54);
  }
  CHECK(std::bit_cast<std::uint64_t>(
          core::fixedOrderSum(cancellation_stress))
        == std::bit_cast<std::uint64_t>(0.0));
}

TEST_CASE("P9-B30 failure propagation selects the lowest failing task index",
          "[core][thread-pool][P9-B30]")
{
  core::ThreadPool pool{ 2 };
  std::atomic<bool> higher_failure_recorded{};
  const auto status = pool.parallelFor(3, [&](std::size_t index) {
    if (index == 0) {
      while (!higher_failure_recorded.load(std::memory_order_acquire))
        std::this_thread::yield();
      return Status::Numerical_failure;
    }
    if (index == 1)
      return Status::Invalid_parameters;
    higher_failure_recorded.store(true, std::memory_order_release);
    return Status::Success;
  });
  CHECK(status == Status::Numerical_failure);
}

TEST_CASE("P9-B31 const and null function callbacks obey the public contract",
          "[core][thread-pool][P9-B31]")
{
  core::ThreadPool pool{ 2 };
  std::array<std::atomic<int>, 3> visits{};
  const auto callback = [&](std::size_t index) {
    visits[index].fetch_add(1, std::memory_order_relaxed);
    return Status::Success;
  };
  REQUIRE(pool.parallelFor(visits.size(), callback) == Status::Success);
  for (const auto &visit : visits)
    CHECK(visit.load(std::memory_order_relaxed) == 1);

  using FunctionPointer = Status (*)(std::size_t);
  FunctionPointer null_callback{};
  CHECK(pool.parallelFor(1, null_callback) == Status::Invalid_parameters);
  CHECK(pool.parallelFor(0, null_callback) == Status::Success);

  bool invoked{};
  CHECK(pool.parallelFor(std::numeric_limits<std::size_t>::max(),
                         [&](std::size_t) { invoked = true; })
        == Status::Invalid_parameters);
  CHECK_FALSE(invoked);
}

TEST_CASE("P9-B30 exceptions finish every task and the pool remains reusable",
          "[core][thread-pool][P9-B30]")
{
  core::ThreadPool pool{ 3 };
  std::array<std::atomic<int>, 17> visits{};
  const auto status = pool.parallelFor(visits.size(), [&](std::size_t index) {
    visits[index].fetch_add(1, std::memory_order_relaxed);
    if (index == 4)
      throw std::runtime_error("controlled callback failure");
  });
  CHECK(status == Status::Unknown_problem);
  for (const auto &visit : visits)
    CHECK(visit.load(std::memory_order_relaxed) == 1);

  std::atomic<std::size_t> sum{};
  REQUIRE(pool.parallelFor(257, [&](std::size_t index) {
    sum.fetch_add(index + 1, std::memory_order_relaxed);
  }) == Status::Success);
  CHECK(sum.load(std::memory_order_relaxed) == 257U * 258U / 2U);
}

TEST_CASE("P9-B30 a concurrent submission is rejected without corrupting the active generation",
          "[core][thread-pool][P9-B30]")
{
  core::ThreadPool pool{ 2 };
  std::atomic<bool> entered{};
  std::atomic<bool> release{};
  Status primary_status = Status::Unknown_problem;
  std::thread submitter{ [&] {
    primary_status = pool.parallelFor(2, [&](std::size_t) {
      entered.store(true, std::memory_order_release);
      while (!release.load(std::memory_order_acquire))
        std::this_thread::yield();
      return Status::Success;
    });
  } };
  while (!entered.load(std::memory_order_acquire))
    std::this_thread::yield();

  std::atomic<bool> secondary_ran{};
  CHECK(pool.parallelFor(1, [&](std::size_t) {
    secondary_ran.store(true, std::memory_order_relaxed);
  }) == Status::Invalid_states);
  CHECK_FALSE(secondary_ran.load(std::memory_order_relaxed));
  release.store(true, std::memory_order_release);
  submitter.join();
  CHECK(primary_status == Status::Success);
}

TEST_CASE("production batch executor is bounded and reports its selected workers",
          "[core][thread-pool][batch-executor][P9-B32]")
{
  core::BatchExecutor executor;
  std::atomic<int> calls{};
  CHECK(executor.parallelFor(1, [&](std::size_t) { ++calls; })
        == Status::Invalid_parameters);
  CHECK(executor.configure(0, 2) == Status::Invalid_parameters);
  CHECK(executor.workerCount() == 0);

  core::BatchExecutor automatic;
  REQUIRE(automatic.configure(2) == Status::Success);
  CHECK(automatic.workerCount() >= 1);
  CHECK(automatic.workerCount() <= 2);

  REQUIRE(executor.configure(3, 7) == Status::Success);
  CHECK(executor.workerCount() == 3);
  CHECK(executor.parallelFor(4, [&](std::size_t) { ++calls; })
        == Status::Invalid_parameters);
  CHECK(calls.load() == 0);
  REQUIRE(executor.parallelFor(3, [&](std::size_t index) {
    calls.fetch_add(static_cast<int>(index + 1), std::memory_order_relaxed);
  }) == Status::Success);
  CHECK(calls.load() == 6);

  REQUIRE(executor.configure(1, 7) == Status::Success);
  CHECK(executor.workerCount() == 1);
  const auto caller = std::this_thread::get_id();
  std::thread::id invoked{};
  REQUIRE(executor.parallelFor(1, [&](std::size_t) {
    invoked = std::this_thread::get_id();
  }) == Status::Success);
  CHECK(invoked == caller);

  using FunctionPointer = Status (*)(std::size_t);
  FunctionPointer null_callback{};
  CHECK(executor.parallelFor(1, null_callback) == Status::Invalid_parameters);

  core::detail::BatchExecutorTestAccess::failPoolCreation(executor);
  CHECK(executor.configure(2, 2) == Status::Numerical_failure);
  CHECK(executor.workerCount() == 1);

  core::BatchExecutor source;
  REQUIRE(source.configure(2, 2) == Status::Success);
  core::BatchExecutor moved{ std::move(source) };
  CHECK(source.workerCount() == 0);
  CHECK(moved.workerCount() == 2);
}

TEST_CASE("legacy parallel facade never silently selects zero workers and transports exceptions",
          "[legacy][thread-pool][P9-B33]")
{
  CHECK(legacy_parallel_detail::workerCount(0, 0, 5) == 1);
  CHECK(legacy_parallel_detail::workerCount(32, 8, 3) == 3);
  CHECK(legacy_parallel_detail::workerCount(2, 8, 0) == 0);

  for (const unsigned workers : { 1U, 2U }) {
    CAPTURE(workers);
    std::array<std::atomic<int>, 6> visits{};
    bool threw{};
    try {
      run(
        [&](int index) {
          visits[static_cast<std::size_t>(index)].fetch_add(
            1, std::memory_order_relaxed);
          if (index == 1 || index == 4)
            throw std::runtime_error(index == 1 ? "lowest" : "higher");
        },
        static_cast<int>(visits.size()),
        workers);
    } catch (const std::runtime_error &error) {
      threw = true;
      CHECK(std::string_view{ error.what() } == "lowest");
    }
    CHECK(threw);
    for (const auto &visit : visits)
      CHECK(visit.load(std::memory_order_relaxed) == 1);
  }
}

TEST_CASE("P8-G4 parallelisation diagnostic reports a validated measured smoke",
          "[core][thread-pool][P8-G4]")
{
  const auto diagnostic = test::parallelisation(2, 32'768);
  INFO(diagnostic.backend);
  INFO(diagnostic.serial_seconds);
  INFO(diagnostic.parallel_seconds);
  INFO(diagnostic.speedup);
  REQUIRE(diagnostic.status == Status::Success);
  CHECK(diagnostic.workers == 2);
  CHECK_FALSE(diagnostic.backend.empty());
  CHECK(diagnostic.serial_seconds > 0.0);
  CHECK(diagnostic.parallel_seconds > 0.0);
  CHECK(std::isfinite(diagnostic.speedup));
  CHECK(std::isfinite(diagnostic.checksum));
}
