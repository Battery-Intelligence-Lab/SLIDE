/**
 * @file CudaSpmBatch.hpp
 * @brief Optional one-lane-per-thread base-isothermal SPM backend.
 */

#pragma once

#include "AsyncRecorder.hpp"

#include <cstddef>
#include <memory>
#include <span>

namespace slide::core {

class CudaAsyncRecorder;

/**
 * Owns one CPU metadata/state mirror and one CUDA arena. CUDA implementation
 * details are hidden behind the pimpl so this public CPU header contains no
 * CUDA types or headers.
 */
class CudaSpmBatch
{
public:
  CudaSpmBatch();
  ~CudaSpmBatch();
  CudaSpmBatch(const CudaSpmBatch &) = delete;
  CudaSpmBatch &operator=(const CudaSpmBatch &) = delete;
  CudaSpmBatch(CudaSpmBatch &&) noexcept;
  CudaSpmBatch &operator=(CudaSpmBatch &&) noexcept;

  static bool available() noexcept;
  [[nodiscard]] slide::Status build(const SpmFactoryInput &input,
                                    const SpmModelOptions &options,
                                    int n_lanes);
  bool valid() const noexcept;

  /** Host mirror used for cold variation setup and host-side PackSolver coupling. */
  SpmBatch &hostBatch();
  const SpmBatch &hostBatch() const;
  [[nodiscard]] slide::Status uploadState();
  [[nodiscard]] slide::Status downloadState();

  /** Queue one accepted exact-modal step. This call allocates and synchronizes nothing. */
  [[nodiscard]] slide::Status step(std::span<const real_t> current_density,
                                   real_t time,
                                   real_t dt);
  /** Wait only for this batch's non-default compute stream and fetch voltages. */
  [[nodiscard]] slide::Status synchronize();
  std::span<const real_t> terminalVoltage() const;

  /** Explicit device-to-device checkpoint/restore, both queued without host sync. */
  [[nodiscard]] slide::Status checkpoint();
  [[nodiscard]] slide::Status restore();

  int nLanes() const;
  std::size_t deviceArenaBytes() const noexcept;
  std::size_t deviceAllocationCount() const noexcept;
  std::size_t deviceWideSynchronizationCount() const noexcept;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_{};
  friend class CudaAsyncRecorder;
};

/** Pinned side-stream bridge into AsyncRecorder's shuffled/compressed writer. */
class CudaAsyncRecorder
{
public:
  CudaAsyncRecorder();
  ~CudaAsyncRecorder();
  CudaAsyncRecorder(const CudaAsyncRecorder &) = delete;
  CudaAsyncRecorder &operator=(const CudaAsyncRecorder &) = delete;

  [[nodiscard]] slide::Status configure(
    CudaSpmBatch &batch,
    const std::filesystem::path &path,
    AsyncRecorderConfig config = {});
  /** Queue a pinned snapshot after the most recently launched accepted step. */
  [[nodiscard]] slide::Status enqueue(std::uint64_t accepted_step);
  [[nodiscard]] slide::Status finish();

  std::uint64_t thinnedSnapshots() const;
  std::uint64_t snapshotsWritten() const;
  bool usesPinnedMemory() const;
  bool usesNonDefaultStream() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_{};
};

} // namespace slide::core
