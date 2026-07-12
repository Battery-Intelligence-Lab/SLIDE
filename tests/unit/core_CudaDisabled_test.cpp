/**
 * @file core_CudaDisabled_test.cpp
 * @brief Optional-off CUDA facade contract and P9-G3 failure-arm coverage.
 */

#include "../../src/core/CudaSpmBatch.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <filesystem>

using namespace slide;

TEST_CASE("CUDA-disabled facades reject every device operation",
          "[core][cuda][optional-off][P9-G3]")
{
  REQUIRE_FALSE(core::CudaSpmBatch::available());

  core::CudaSpmBatch batch;
  const core::SpmFactoryInput input{};
  const core::SpmModelOptions options{};
  const std::array<double, 1> current_density{ 0.0 };
  CHECK(batch.build(input, options, 1) == Status::NotImplementedYet);
  CHECK(batch.uploadState() == Status::NotImplementedYet);
  CHECK(batch.downloadState() == Status::NotImplementedYet);
  CHECK(batch.step(current_density, 0.0, 1.0) == Status::NotImplementedYet);
  CHECK(batch.synchronize() == Status::NotImplementedYet);
  CHECK(batch.checkpoint() == Status::NotImplementedYet);
  CHECK(batch.restore() == Status::NotImplementedYet);

  core::CudaAsyncRecorder recorder;
  CHECK(recorder.configure(batch,
                           std::filesystem::path{ "cuda-disabled.slrec" })
        == Status::NotImplementedYet);
  CHECK(recorder.enqueue(0) == Status::NotImplementedYet);
  CHECK(recorder.finish() == Status::Success);
  CHECK(recorder.thinnedSnapshots() == 0);
  CHECK(recorder.snapshotsWritten() == 0);
  CHECK_FALSE(recorder.usesPinnedMemory());
  CHECK_FALSE(recorder.usesNonDefaultStream());
}
