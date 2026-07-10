/**
 * @file AsyncRecorder.hpp
 * @brief Preallocated asynchronous snapshot ring with shuffled block compression.
 */

#pragma once

#include "Recorder.hpp"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <span>
#include <thread>
#include <vector>

namespace slide::core {

enum class AsyncBackpressurePolicy : unsigned char {
  block,
  thin
};

enum class CompressionCodec : unsigned char {
  none,
  zstd
};

struct AsyncRecorderConfig
{
  std::size_t cadence{ 1 };
  std::size_t ring_slots{ 3 };
  AsyncBackpressurePolicy backpressure{ AsyncBackpressurePolicy::thin };
#if defined(SLIDE_WITH_ZSTD)
  CompressionCodec codec{ CompressionCodec::zstd };
#else
  CompressionCodec codec{ CompressionCodec::none };
#endif
  int compression_level{ 3 };
};

/** True when the requested codec is compiled into this build. */
bool compressionCodecAvailable(CompressionCodec codec);

/** Reversible Blosc-style byte transpose for homogeneous fixed-width values. */
void byteShuffle(std::span<const std::byte> input,
                 std::span<std::byte>
                   output,
                 std::size_t element_width);
void byteUnshuffle(std::span<const std::byte> input,
                   std::span<std::byte>
                     output,
                   std::size_t element_width);

class AsyncRecorder
{
public:
  AsyncRecorder() = default;
  ~AsyncRecorder();
  AsyncRecorder(const AsyncRecorder &) = delete;
  AsyncRecorder &operator=(const AsyncRecorder &) = delete;

  [[nodiscard]] slide::Status configure(
    SpmBatch &batch,
    const std::filesystem::path &path,
    AsyncRecorderConfig config = {});
  [[nodiscard]] slide::Status enqueue(
    std::uint64_t accepted_step,
    std::span<const real_t> total_current_A);
  /** Drain, finalize CRC/header metadata, and close. Idempotent. */
  [[nodiscard]] slide::Status finish();

  bool configured() const { return batch_ != nullptr; }
  std::uint64_t thinnedSnapshots() const
  {
    return thinned_.load(std::memory_order_relaxed);
  }
  std::uint64_t snapshotsWritten() const
  {
    return written_.load(std::memory_order_relaxed);
  }
  slide::Status workerStatus() const
  {
    return static_cast<slide::Status>(
      worker_status_.load(std::memory_order_relaxed));
  }

  /** Cold-path test/instrumentation hook, invoked by the drain thread per block. */
  void setDrainHook(std::function<void()> hook) { drain_hook_ = std::move(hook); }

private:
  enum class SlotState : unsigned char { empty,
                                         filling,
                                         ready,
                                         draining };
  struct Slot
  {
    SlotState state{ SlotState::empty };
    std::uint64_t accepted_step{};
    real_t time{};
    std::vector<real_t> values{}; //!< current density then full padded arena
    std::vector<std::byte> shuffled{};
    std::vector<std::byte> compressed{};
  };

  void drainLoop();
  slide::Status drainSlot(Slot &slot);
  slide::Status finalizeFile();
  void setWorkerFailure(slide::Status status);

  SpmBatch *batch_{};
  AsyncRecorderConfig config_{};
  std::filesystem::path path_{};
  int rows_{};
  int lanes_{};
  int stride_{};
  std::size_t state_values_{};
  std::size_t raw_values_{};
  std::vector<Slot> slots_{};
  std::size_t write_slot_{};
  std::size_t read_slot_{};
  std::uint64_t previous_step_{};
  bool has_previous_step_{};
  bool closing_{};
  bool finished_{};
  std::mutex mutex_{};
  std::condition_variable ready_{};
  std::condition_variable space_{};
  std::thread drain_{};
  std::ofstream output_{};
  std::atomic<std::uint64_t> thinned_{};
  std::atomic<std::uint64_t> written_{};
  std::atomic<int> worker_status_{ static_cast<int>(slide::Status::Success) };
  std::function<void()> drain_hook_{};
  void *codec_context_{};
};

/** Hardened eager reader for the async shuffled/compressed block format. */
class CompressedRecording
{
public:
  [[nodiscard]] slide::Status open(const std::filesystem::path &path);
  void close();
  bool valid() const { return valid_; }
  std::size_t size() const { return steps_.size(); }
  int nRows() const { return rows_; }
  int nLanes() const { return lanes_; }
  int stride() const { return stride_; }
  SnapshotView snapshot(std::size_t index) const;

private:
  bool valid_{};
  int rows_{};
  int lanes_{};
  int stride_{};
  std::size_t state_values_{};
  std::vector<std::uint64_t> steps_{};
  std::vector<real_t> times_{};
  std::vector<real_t> currents_{};
  std::vector<real_t> states_{};
};

} // namespace slide::core
