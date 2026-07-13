/**
 * @file AsyncRecorder.cpp
 * @brief The asynchronous recording ring: slots, the drain worker, and the block writer.
 *
 * Owns: every `AsyncRecorder` member -- configuration, the slot ring, enqueue, the drain loop,
 * and file finalisation. Implements PLAN.md §3.7. Hot: `enqueue` runs on the simulation thread
 * and must not block on I/O. The file format it writes, and the reader that consumes it, live in
 * `detail/AsyncRecordingFormat.hpp` and `AsyncRecordingCodec.cpp`.
 * @surface internal
 */

#include "AsyncRecorder.hpp"
#include "detail/AsyncRecordingFormat.hpp"
#include "detail/CheckedArithmetic.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <system_error>

namespace slide::core {

using detail::allocationFailureStatus;
using detail::blockHeaderCrc;
using detail::block_magic;
using detail::checkedAdd;
using detail::checkedMultiply;
using detail::compressed_magic;
using detail::compressionBound;
using detail::CompressedBlockHeader;
using detail::CompressedFileHeader;
using detail::crc32;
using detail::endian_marker;
using detail::fileHeaderCrc;
using detail::format_major;
using detail::format_minor;

AsyncRecorder::~AsyncRecorder()
{
  (void)finish();
}

std::thread AsyncRecorder::makeDrainThread(AsyncRecorder &recorder)
{
  return std::thread{ [&recorder] { recorder.drainLoop(); } };
}

void AsyncRecorder::noOutputHook(std::ofstream &) noexcept
{}

slide::Status AsyncRecorder::configure(SpmBatch &batch,
                                       const std::filesystem::path &path,
                                       AsyncRecorderConfig config)
{
  if (configured() || !batch.valid() || path.empty() || config.cadence == 0
      || config.ring_slots < 3 || config.ring_slots > 1024
      || !(config.backpressure == AsyncBackpressurePolicy::block
           || config.backpressure == AsyncBackpressurePolicy::thin)
      || !compressionCodecAvailable(config.codec)
      || config.ring_slots > std::vector<Slot>{}.max_size()
      || config.compression_level < -20 || config.compression_level > 22)
    return slide::Status::Invalid_parameters;
  const std::size_t state_values = batch.state().size();
  const std::size_t lanes = static_cast<std::size_t>(batch.n_lanes());
  detail::AsyncBufferLayout layout;
  const auto layout_status = detail::asyncBufferLayout(
    lanes, state_values, config.codec, layout);
  if (layout_status != slide::Status::Success)
    return layout_status;

  try {
    std::filesystem::path candidate_path{ path };
    std::vector<Slot> slots(config.ring_slots);
    for (auto &slot : slots) {
      slot.values.resize(layout.values);
      slot.shuffled.resize(layout.raw_bytes);
      slot.compressed.resize(layout.compressed_bound);
    }
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output)
      return slide::Status::Invalid_parameters;
    const CompressedFileHeader placeholder{};
    before_placeholder_write_(output);
    output.write(reinterpret_cast<const char *>(&placeholder), sizeof(placeholder));
    if (!output)
      return slide::Status::Numerical_failure;

#if defined(SLIDE_WITH_ZSTD)
    void *context = config.codec == CompressionCodec::zstd
                      ? static_cast<void *>(ZSTD_createCCtx())
                      : nullptr;
    if (config.codec == CompressionCodec::zstd && context == nullptr)
      return slide::Status::Numerical_failure;
#else
    void *context = nullptr;
#endif
    batch_ = &batch;
    config_ = config;
    path_.swap(candidate_path);
    rows_ = batch.state().n_rows();
    lanes_ = batch.state().n_lanes();
    stride_ = batch.state().stride();
    state_values_ = state_values;
    raw_values_ = layout.values;
    slots_ = std::move(slots);
    output_ = std::move(output);
    codec_context_ = context;
    write_slot_ = 0;
    read_slot_ = 0;
    previous_step_ = 0;
    has_previous_step_ = false;
    closing_ = false;
    finished_ = false;
    thinned_.store(0, std::memory_order_relaxed);
    written_.store(0, std::memory_order_relaxed);
    worker_status_.store(static_cast<int>(slide::Status::Success),
                         std::memory_order_relaxed);
    try {
      drain_ = drain_thread_factory_(*this);
    } catch (...) {
#if defined(SLIDE_WITH_ZSTD)
      if (codec_context_ != nullptr)
        ZSTD_freeCCtx(static_cast<ZSTD_CCtx *>(codec_context_));
#endif
      codec_context_ = nullptr;
      output_.close();
      slots_.clear();
      batch_ = nullptr;
      return slide::Status::Numerical_failure;
    }
  } catch (const std::bad_alloc &) {
    return allocationFailureStatus();
  } catch (const std::length_error &) {
    return allocationFailureStatus();
  }
  return slide::Status::Success;
}

slide::Status AsyncRecorder::enqueue(
  std::uint64_t accepted_step,
  std::span<const real_t>
    total_current_A)
{
  if (!configured())
    return slide::Status::Invalid_parameters;
  return enqueueValues(accepted_step,
                       batch_->state().at(batch_->layout().elapsed_time, 0, 0),
                       total_current_A,
                       batch_->state().raw(),
                       false);
}

slide::Status AsyncRecorder::enqueueSnapshot(
  std::uint64_t accepted_step,
  real_t time,
  std::span<const real_t>
    current_density,
  std::span<const real_t>
    state)
{
  return enqueueValues(accepted_step,
                       time,
                       current_density,
                       state,
                       true);
}

slide::Status AsyncRecorder::enqueueValues(
  std::uint64_t accepted_step,
  real_t time,
  std::span<const real_t>
    currents,
  std::span<const real_t>
    state,
  bool currents_are_density)
{
  if (!configured() || finished_
      || currents.size() != static_cast<std::size_t>(lanes_)
      || state.size() != state_values_ || !is_finite(time))
    return slide::Status::Invalid_parameters;
  if (accepted_step % config_.cadence != 0)
    return slide::Status::Success;
  for (const real_t current : currents)
    if (!is_finite(current)
        || (!currents_are_density
            && !is_finite(current / batch_->electrode_area())))
      return slide::Status::Invalid_parameters;

  std::unique_lock lock{ mutex_ };
  if (has_previous_step_ && accepted_step <= previous_step_)
    return slide::Status::Invalid_parameters;
  previous_step_ = accepted_step;
  has_previous_step_ = true;
  if (workerStatus() != slide::Status::Success)
    return workerStatus();
  auto available = [&] {
    return slots_[write_slot_].state == SlotState::empty
           || workerStatus() != slide::Status::Success || closing_;
  };
  if (!available()) {
    if (config_.backpressure == AsyncBackpressurePolicy::thin) {
      thinned_.fetch_add(1, std::memory_order_relaxed);
      return slide::Status::Success;
    }
    space_.wait(lock, available);
  }
  if (closing_ || workerStatus() != slide::Status::Success)
    return workerStatus() == slide::Status::Success
             ? slide::Status::Invalid_states
             : workerStatus();
  Slot &slot = slots_[write_slot_];
  slot.state = SlotState::filling;
  slot.accepted_step = accepted_step;
  slot.time = time;
  lock.unlock();

  for (int lane = 0; lane < lanes_; ++lane)
    slot.values[static_cast<std::size_t>(lane)] =
      currents[static_cast<std::size_t>(lane)]
      / (currents_are_density ? 1.0 : batch_->electrode_area());
  std::memcpy(slot.values.data() + lanes_,
              state.data(),
              state_values_ * sizeof(real_t));

  lock.lock();
  slot.state = SlotState::ready;
  write_slot_ = (write_slot_ + 1) % slots_.size();
  lock.unlock();
  ready_.notify_one();
  return slide::Status::Success;
}

void AsyncRecorder::setWorkerFailure(slide::Status status)
{
  int expected = static_cast<int>(slide::Status::Success);
  worker_status_.compare_exchange_strong(expected,
                                         static_cast<int>(status),
                                         std::memory_order_relaxed);
}

slide::Status AsyncRecorder::drainSlot(Slot &slot)
{
  const auto raw = std::as_bytes(std::span{ slot.values });
  if (byteShuffle(raw, slot.shuffled, sizeof(real_t))
      != slide::Status::Success)
    return slide::Status::Invalid_states;
  std::span<const std::byte> payload;
#if defined(SLIDE_WITH_ZSTD)
  if (config_.codec == CompressionCodec::none) {
    std::memcpy(slot.compressed.data(), slot.shuffled.data(), slot.shuffled.size());
    payload = std::span<const std::byte>{ slot.compressed }.first(
      slot.shuffled.size());
  } else {
    assert(config_.codec == CompressionCodec::zstd);
    const std::size_t compressed = ZSTD_compressCCtx(
      static_cast<ZSTD_CCtx *>(codec_context_),
      slot.compressed.data(),
      slot.compressed.size(),
      slot.shuffled.data(),
      slot.shuffled.size(),
      config_.compression_level);
    if (ZSTD_isError(compressed))
      return slide::Status::Numerical_failure;
    payload = std::span<const std::byte>{ slot.compressed }.first(compressed);
  }
#else
  // configure() rejects every codec except none in an optional-off build, and
  // config_ is private and immutable while the worker is running.
  assert(config_.codec == CompressionCodec::none);
  std::memcpy(slot.compressed.data(), slot.shuffled.data(), slot.shuffled.size());
  payload = std::span<const std::byte>{ slot.compressed }.first(
    slot.shuffled.size());
#endif

  CompressedBlockHeader header{
    .magic = block_magic,
    .header_size = sizeof(CompressedBlockHeader),
    .codec = static_cast<std::uint32_t>(config_.codec),
    .flags = 1U, // byte shuffled
    .accepted_step = slot.accepted_step,
    .time = slot.time,
    .raw_bytes = raw.size(),
    .payload_bytes = payload.size(),
    .raw_crc32 = crc32(raw),
    .payload_crc32 = crc32(payload),
  };
  header.header_crc32 = blockHeaderCrc(header);
  before_block_write_(output_);
  output_.write(reinterpret_cast<const char *>(&header), sizeof(header));
  output_.write(reinterpret_cast<const char *>(payload.data()),
                static_cast<std::streamsize>(payload.size()));
  return output_ ? slide::Status::Success
                 : slide::Status::Numerical_failure;
}

void AsyncRecorder::drainLoop()
{
  while (true) {
    std::unique_lock lock{ mutex_ };
    ready_.wait(lock, [&] {
      return slots_[read_slot_].state == SlotState::ready || closing_;
    });
    Slot &slot = slots_[read_slot_];
    if (slot.state != SlotState::ready) {
      if (closing_)
        break;
      continue;
    }
    slot.state = SlotState::draining;
    lock.unlock();

    if (drain_hook_)
      drain_hook_();
    const auto status = drainSlot(slot);
    if (status != slide::Status::Success)
      setWorkerFailure(status);

    lock.lock();
    slot.state = SlotState::empty;
    read_slot_ = (read_slot_ + 1) % slots_.size();
    written_.fetch_add(status == slide::Status::Success ? 1U : 0U,
                       std::memory_order_relaxed);
    lock.unlock();
    space_.notify_all();
  }
}

slide::Status AsyncRecorder::finalizeFile()
{
  // A failed block write sets worker_status_ and finish() skips finalization;
  // otherwise configure() owns this still-open stream until this function.
  assert(output_.is_open() && output_);
  before_finalize_tell_(output_);
  const auto position = output_.tellp();
  if (position < 0) {
    output_.close();
    return slide::Status::Numerical_failure;
  }
  const auto file_size = static_cast<std::uint64_t>(position);
  CompressedFileHeader header{ .magic = compressed_magic,
                               .major = format_major,
                               .minor = format_minor,
                               .endian = endian_marker,
                               .header_size = sizeof(CompressedFileHeader),
                               .rows = static_cast<std::uint32_t>(rows_),
                               .lanes = static_cast<std::uint32_t>(lanes_),
                               .stride = static_cast<std::uint32_t>(stride_),
                               .codec = static_cast<std::uint32_t>(config_.codec),
                               .snapshots = snapshotsWritten(),
                               .file_size = file_size };
  header.header_crc32 = fileHeaderCrc(header);
  before_finalize_write_(output_);
  output_.seekp(0);
  output_.write(reinterpret_cast<const char *>(&header), sizeof(header));
  output_.flush();
  output_.close();
  return output_.good() ? slide::Status::Success
                        : slide::Status::Numerical_failure;
}

slide::Status AsyncRecorder::finish()
{
  if (!configured())
    return slide::Status::Success;
  if (finished_)
    return workerStatus();
  {
    const std::lock_guard lock{ mutex_ };
    closing_ = true;
  }
  ready_.notify_all();
  space_.notify_all();
  if (drain_.joinable())
    drain_.join();
  if (workerStatus() == slide::Status::Success) {
    const auto status = finalizeFile();
    if (status != slide::Status::Success)
      setWorkerFailure(status);
  } else {
    output_.close();
  }
#if defined(SLIDE_WITH_ZSTD)
  if (codec_context_ != nullptr)
    ZSTD_freeCCtx(static_cast<ZSTD_CCtx *>(codec_context_));
#endif
  codec_context_ = nullptr;
  finished_ = true;
  return workerStatus();
}


} // namespace slide::core
