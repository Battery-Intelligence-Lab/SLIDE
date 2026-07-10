/**
 * @file core_only_smoke.cpp
 * @brief External-consumer and disabled-capability smoke gate for SLIDE_CORE_ONLY.
 */

#include <core/AsyncRecorder.hpp>
#include <core/BatchBuilder.hpp>
#include <core/CudaSpmBatch.hpp>
#include <core/PackTopology.hpp>
#include <core/ParameterSet.hpp>
#include <core/ThreadPool.hpp>

#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>

namespace {

int fail(const char *message)
{
  std::cerr << "SLIDE core-only smoke failure: " << message << '\n';
  return 1;
}

} // namespace

int main()
{
  using namespace slide::core;

  if (CudaSpmBatch::available())
    return fail("CUDA reports available in an explicitly disabled build");
  if (compressionCodecAvailable(CompressionCodec::zstd))
    return fail("zstd reports available in an explicitly disabled build");
  if (AsyncRecorderConfig{}.codec != CompressionCodec::none)
    return fail("the dependency-light recorder defaults to an unavailable codec");
  Recorder recorder;
  if (recorder.writeParquet("core-only-disabled.parquet") != slide::Status::NotImplementedYet)
    return fail("Arrow/Parquet does not report its disabled capability explicitly");

  BatchBuilder builder;
  const auto state = builder.declare({ "smoke_state", 2, Unit::none, StateRole::ode });
  auto arena = builder.build(3);
  if (arena.n_rows() != 2 || arena.n_lanes() != 3 || arena.stride() != 8)
    return fail("StateArena geometry or alignment padding is incorrect");
  if (reinterpret_cast<std::uintptr_t>(arena.raw().data()) % StateArena::alignment != 0)
    return fail("StateArena storage is not 64-byte aligned");
  for (int lane = 0; lane < arena.n_lanes(); ++lane)
    arena.at(state, 1, lane) = static_cast<real_t>(lane + 1);

  ThreadPool pool(2);
  std::array<std::atomic<int>, 7> visits{};
  const auto status = pool.parallelFor(visits.size(), [&](std::size_t index) {
    visits[index].fetch_add(1, std::memory_order_relaxed);
  });
  if (status != slide::Status::Success)
    return fail("ThreadPool rejected a valid task batch");
  for (const auto &visit : visits)
    if (visit.load(std::memory_order_relaxed) != 1)
      return fail("ThreadPool did not execute every task exactly once");

  const std::array<real_t, 4> values{ 1.0, -2.0, 3.0, 4.0 };
  if (fixedOrderSum(values) != 6.0)
    return fail("fixed-order reduction produced the wrong result");

  return 0;
}
