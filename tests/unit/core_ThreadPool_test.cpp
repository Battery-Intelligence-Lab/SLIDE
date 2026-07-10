/**
 * @file core_ThreadPool_test.cpp
 * @brief P8-G4 persistent pool, failure, determinism, and diagnostic gates.
 */

#include "../../src/core/ThreadPool.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <vector>

using namespace slide;

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
