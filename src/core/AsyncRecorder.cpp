/**
 * @file AsyncRecorder.cpp
 * @brief Async ring, byte shuffle, zstd block writer, and hardened reader.
 */

#include "AsyncRecorder.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <system_error>

#if defined(SLIDE_WITH_ZSTD)
#include <zstd.h>
#endif

namespace slide::core {
namespace {

  constexpr std::array<char, 8> compressed_magic{ 'S', 'L', 'I', 'D', 'E', 'C', 'M', 'P' };
  constexpr std::uint32_t block_magic = 0x314b4c42U; // BLK1
  constexpr std::uint32_t endian_marker = 0x01020304U;
  constexpr std::uint16_t format_major = 1;
  constexpr std::uint16_t format_minor = 0;

#pragma pack(push, 1)
  struct CompressedFileHeader
  {
    std::array<char, 8> magic{};
    std::uint16_t major{};
    std::uint16_t minor{};
    std::uint32_t endian{};
    std::uint32_t header_size{};
    std::uint32_t header_crc32{};
    std::uint32_t rows{};
    std::uint32_t lanes{};
    std::uint32_t stride{};
    std::uint32_t codec{};
    std::uint64_t snapshots{};
    std::uint64_t file_size{};
    std::uint64_t reserved{};
  };

  struct CompressedBlockHeader
  {
    std::uint32_t magic{};
    std::uint32_t header_size{};
    std::uint32_t codec{};
    std::uint32_t flags{};
    std::uint64_t accepted_step{};
    real_t time{};
    std::uint64_t raw_bytes{};
    std::uint64_t payload_bytes{};
    std::uint32_t raw_crc32{};
    std::uint32_t payload_crc32{};
    std::uint32_t header_crc32{};
    std::uint32_t reserved{};
  };
#pragma pack(pop)

  static_assert(sizeof(CompressedFileHeader) == 64);
  static_assert(sizeof(CompressedBlockHeader) == 64);

  std::uint32_t crc32(std::span<const std::byte> bytes)
  {
    std::uint32_t crc = 0xffffffffU;
    for (const auto byte : bytes) {
      crc ^= std::to_integer<std::uint8_t>(byte);
      for (int bit = 0; bit < 8; ++bit)
        crc = (crc >> 1U) ^ (0xedb88320U & (0U - (crc & 1U)));
    }
    return ~crc;
  }

  std::uint32_t fileHeaderCrc(CompressedFileHeader header)
  {
    header.header_crc32 = 0;
    return crc32(std::as_bytes(std::span{ &header, 1 }));
  }

  std::uint32_t blockHeaderCrc(CompressedBlockHeader header)
  {
    header.header_crc32 = 0;
    return crc32(std::as_bytes(std::span{ &header, 1 }));
  }

  template <class T>
  bool checkedAdd(T left, T right, T &result)
  {
    if (left > std::numeric_limits<T>::max() - right)
      return false;
    result = left + right;
    return true;
  }

  template <class T>
  bool checkedMultiply(T left, T right, T &result)
  {
    if (left != 0 && right > std::numeric_limits<T>::max() / left)
      return false;
    result = left * right;
    return true;
  }

  slide::Status allocationFailureStatus() noexcept
  {
    return slide::Status::Numerical_failure;
  }

  std::size_t compressionBound(CompressionCodec codec, std::size_t raw_bytes)
  {
    if (codec == CompressionCodec::none)
      return raw_bytes;
#if defined(SLIDE_WITH_ZSTD)
    if (codec == CompressionCodec::zstd)
      return ZSTD_compressBound(raw_bytes);
#else
    (void)raw_bytes;
#endif
    return 0;
  }

} // namespace

bool compressionCodecAvailable(CompressionCodec codec)
{
  if (codec == CompressionCodec::none)
    return true;
#if defined(SLIDE_WITH_ZSTD)
  return codec == CompressionCodec::zstd;
#else
  return false;
#endif
}

slide::Status detail::asyncBufferLayout(std::size_t lanes,
                                        std::size_t state_values,
                                        CompressionCodec codec,
                                        AsyncBufferLayout &output)
{
  if (lanes > std::numeric_limits<std::size_t>::max() - state_values)
    return slide::Status::Invalid_parameters;
  const std::size_t raw_values = lanes + state_values;
  std::size_t raw_bytes{};
  if (!checkedMultiply(raw_values, sizeof(real_t), raw_bytes))
    return slide::Status::Invalid_parameters;
  const std::size_t bound = compressionBound(codec, raw_bytes);
  if (bound == 0 || raw_values > std::vector<real_t>{}.max_size()
      || raw_bytes > std::vector<std::byte>{}.max_size()
      || bound > std::vector<std::byte>{}.max_size())
    return slide::Status::Invalid_parameters;
  output = { .state_values = state_values,
             .values = raw_values,
             .raw_bytes = raw_bytes,
             .compressed_bound = bound };
  return slide::Status::Success;
}

slide::Status detail::compressedRecordingBufferLayout(
  std::size_t rows,
  std::size_t lanes,
  std::size_t stride,
  CompressionCodec codec,
  AsyncBufferLayout &output)
{
  std::size_t state_values{};
  if (!checkedMultiply(rows, stride, state_values))
    return slide::Status::Invalid_parameters;
  return asyncBufferLayout(lanes, state_values, codec, output);
}

slide::Status detail::compressedRecordingStorageLayout(
  std::size_t snapshots,
  std::size_t state_values,
  std::size_t lanes,
  std::size_t raw_bytes,
  std::size_t compressed_bound,
  CompressionCodec codec,
  std::size_t retained_bytes,
  CompressedRecordingStorageLayout &output)
{
  std::size_t state_storage{}, current_storage{};
  std::size_t step_bytes{}, time_bytes{}, state_bytes{}, current_bytes{};
  std::size_t raw_scratch_bytes{}, scratch_bytes{}, resident_bytes{};
  const auto real_vector_limit = std::vector<real_t>{}.max_size();
  const bool valid_codec = codec == CompressionCodec::none
                           || codec == CompressionCodec::zstd;
  const std::size_t workspace_bytes = snapshots == 0 ? 0 : raw_bytes;
  const std::size_t payload_bytes = snapshots == 0
                                      ? 0
                                    : codec == CompressionCodec::none
                                      ? raw_bytes
                                      : compressed_bound;
  if (!valid_codec
      || !checkedMultiply(snapshots, state_values, state_storage)
      || !checkedMultiply(snapshots, lanes, current_storage)
      || !checkedMultiply(
        snapshots, sizeof(std::uint64_t), step_bytes)
      || !checkedMultiply(snapshots, sizeof(real_t), time_bytes)
      || !checkedMultiply(state_storage, sizeof(real_t), state_bytes)
      || !checkedMultiply(current_storage, sizeof(real_t), current_bytes)
      || !checkedMultiply(
        workspace_bytes, std::size_t{ 2 }, raw_scratch_bytes)
      || payload_bytes
           > static_cast<std::size_t>(
             std::numeric_limits<std::streamsize>::max())
      || !checkedAdd(raw_scratch_bytes, payload_bytes, scratch_bytes)
      || !checkedAdd(retained_bytes, step_bytes, resident_bytes)
      || !checkedAdd(resident_bytes, time_bytes, resident_bytes)
      || !checkedAdd(resident_bytes, state_bytes, resident_bytes)
      || !checkedAdd(resident_bytes, current_bytes, resident_bytes)
      || !checkedAdd(resident_bytes, scratch_bytes, resident_bytes)
      || snapshots > std::vector<std::uint64_t>{}.max_size()
      || snapshots > real_vector_limit || state_storage > real_vector_limit
      || current_storage > real_vector_limit
      || resident_bytes > max_eager_recording_bytes)
    return slide::Status::Invalid_parameters;
  output = { .state_values = state_storage,
             .current_values = current_storage,
             .payload_bytes = payload_bytes,
             .workspace_bytes = workspace_bytes,
             .resident_bytes = resident_bytes };
  return slide::Status::Success;
}

slide::Status detail::compressedRecordingRetainedBytes(
  std::size_t step_capacity,
  std::size_t time_capacity,
  std::size_t current_capacity,
  std::size_t state_capacity,
  std::size_t &output)
{
  std::size_t step_bytes{}, time_bytes{}, current_bytes{}, state_bytes{};
  std::size_t total{};
  if (!checkedMultiply(step_capacity, sizeof(std::uint64_t), step_bytes)
      || !checkedMultiply(time_capacity, sizeof(real_t), time_bytes)
      || !checkedMultiply(current_capacity, sizeof(real_t), current_bytes)
      || !checkedMultiply(state_capacity, sizeof(real_t), state_bytes)
      || !checkedAdd(step_bytes, time_bytes, total)
      || !checkedAdd(total, current_bytes, total)
      || !checkedAdd(total, state_bytes, total))
    return slide::Status::Invalid_parameters;
  output = total;
  return slide::Status::Success;
}

slide::Status detail::validateCompressedPayloadSize(
  CompressionCodec codec,
  std::uint64_t payload_bytes,
  const AsyncBufferLayout &layout)
{
  const bool valid = codec == CompressionCodec::none
                       ? payload_bytes == layout.raw_bytes
                       : codec == CompressionCodec::zstd
                           && payload_bytes <= layout.compressed_bound;
  return valid ? slide::Status::Success
               : slide::Status::Invalid_parameters;
}

namespace {

  bool spansOverlap(std::span<const std::byte> input,
                    std::span<std::byte>
                      output)
  {
    if (input.empty())
      return false;
    const auto less = std::less<const std::byte *>{};
    const auto *input_begin = input.data();
    const auto *input_end = input_begin + input.size();
    const auto *output_begin = output.data();
    const auto *output_end = output_begin + output.size();
    return less(input_begin, output_end) && less(output_begin, input_end);
  }

  bool validShuffleArguments(std::span<const std::byte> input,
                             std::span<std::byte>
                               output,
                             std::size_t element_width)
  {
    return element_width > 0 && input.size() == output.size()
           && input.size() % element_width == 0
           && !spansOverlap(input, output);
  }

} // namespace

slide::Status byteShuffle(std::span<const std::byte> input,
                          std::span<std::byte>
                            output,
                          std::size_t element_width)
{
  if (!validShuffleArguments(input, output, element_width))
    return slide::Status::Invalid_parameters;
  const std::size_t elements = input.size() / element_width;
  for (std::size_t byte = 0; byte < element_width; ++byte)
    for (std::size_t element = 0; element < elements; ++element)
      output[byte * elements + element] = input[element * element_width + byte];
  return slide::Status::Success;
}

slide::Status byteUnshuffle(std::span<const std::byte> input,
                            std::span<std::byte>
                              output,
                            std::size_t element_width)
{
  if (!validShuffleArguments(input, output, element_width))
    return slide::Status::Invalid_parameters;
  const std::size_t elements = input.size() / element_width;
  for (std::size_t byte = 0; byte < element_width; ++byte)
    for (std::size_t element = 0; element < elements; ++element)
      output[element * element_width + byte] = input[byte * elements + element];
  return slide::Status::Success;
}

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

slide::Status CompressedRecording::readExact(
  std::istream &input,
  std::span<std::byte>
    destination,
  std::uint64_t &remaining)
{
  // All production callers pass fixed 64-byte headers or payloads whose size
  // was checked against streamsize::max before allocation.
  assert(destination.size()
         <= static_cast<std::size_t>(
           std::numeric_limits<std::streamsize>::max()));
  if (destination.size() > remaining)
    return slide::Status::Invalid_parameters;
  const auto bytes = static_cast<std::streamsize>(destination.size());
  input.read(reinterpret_cast<char *>(destination.data()), bytes);
  const bool complete = input && input.gcount() == bytes;
  if (complete)
    remaining -= destination.size();
  return complete ? slide::Status::Success
                  : slide::Status::Invalid_parameters;
}

slide::Status CompressedRecording::open(const std::filesystem::path &path)
try {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    return slide::Status::Invalid_parameters;
  input.seekg(0, std::ios::end);
  const auto length = input.tellg();
  if (length < static_cast<std::streamoff>(sizeof(CompressedFileHeader)))
    return slide::Status::Invalid_parameters;
  input.seekg(0);
  std::uint64_t remaining = static_cast<std::uint64_t>(length);
  CompressedFileHeader header{};
  const auto header_status = readExact(
    input, std::as_writable_bytes(std::span{ &header, 1 }), remaining);
  if (header_status != slide::Status::Success)
    return header_status;
  if (header.magic != compressed_magic || header.major != format_major
      || header.minor != format_minor || header.endian != endian_marker
      || header.header_size != sizeof(CompressedFileHeader)
      || header.header_crc32 != fileHeaderCrc(header)
      || header.file_size != static_cast<std::uint64_t>(length)
      || header.rows == 0 || header.lanes == 0 || header.stride < header.lanes
      || header.rows > static_cast<std::uint32_t>(std::numeric_limits<int>::max())
      || header.lanes > static_cast<std::uint32_t>(std::numeric_limits<int>::max())
      || header.stride > static_cast<std::uint32_t>(std::numeric_limits<int>::max())
      || header.snapshots > std::numeric_limits<std::size_t>::max()
      || (header.codec != static_cast<std::uint32_t>(CompressionCodec::none)
          && header.codec
               != static_cast<std::uint32_t>(CompressionCodec::zstd))
      || !compressionCodecAvailable(
        static_cast<CompressionCodec>(header.codec)))
    return slide::Status::Invalid_parameters;

  if (header.snapshots
      > (header.file_size - sizeof(CompressedFileHeader))
          / sizeof(CompressedBlockHeader))
    return slide::Status::Invalid_parameters;

  detail::AsyncBufferLayout buffer_layout;
  const auto layout_status = detail::compressedRecordingBufferLayout(
    header.rows,
    header.lanes,
    header.stride,
    static_cast<CompressionCodec>(header.codec),
    buffer_layout);
  if (layout_status != slide::Status::Success)
    return layout_status;
  const std::size_t state_values = buffer_layout.state_values;
  const std::size_t raw_bytes = buffer_layout.raw_bytes;
  const auto snapshots = static_cast<std::size_t>(header.snapshots);
  std::size_t retained_bytes{};
  const auto retained_status = detail::compressedRecordingRetainedBytes(
    steps_.capacity(),
    times_.capacity(),
    currents_.capacity(),
    states_.capacity(),
    retained_bytes);
  if (retained_status != slide::Status::Success)
    return retained_status;
  detail::CompressedRecordingStorageLayout storage_layout;
  const auto storage_status = detail::compressedRecordingStorageLayout(
    snapshots,
    state_values,
    header.lanes,
    buffer_layout.raw_bytes,
    buffer_layout.compressed_bound,
    static_cast<CompressionCodec>(header.codec),
    retained_bytes,
    storage_layout);
  if (storage_status != slide::Status::Success)
    return storage_status;

  {
    std::vector<std::uint64_t> steps(static_cast<std::size_t>(header.snapshots));
    std::vector<real_t> times(steps.size());
    std::vector<real_t> currents(storage_layout.current_values);
    std::vector<real_t> states(storage_layout.state_values);
    // Allocate the full validated payload workspace exactly once. Per-block
    // resize growth could otherwise exceed the resident-byte budget.
    std::vector<std::byte> payload(storage_layout.payload_bytes);
    std::vector<std::byte> shuffled(storage_layout.workspace_bytes);
    std::vector<std::byte> raw(storage_layout.workspace_bytes);
    for (std::size_t index = 0; index < steps.size(); ++index) {
      CompressedBlockHeader block{};
      const auto block_status = readExact(
        input, std::as_writable_bytes(std::span{ &block, 1 }), remaining);
      if (block_status != slide::Status::Success)
        return block_status;
      if (block.magic != block_magic
          || block.header_size != sizeof(CompressedBlockHeader)
          || block.header_crc32 != blockHeaderCrc(block)
          || block.codec != header.codec || block.flags != 1U
          || block.raw_bytes != raw_bytes || block.payload_bytes > remaining
          || block.payload_bytes > std::numeric_limits<std::size_t>::max()
          || block.payload_bytes
               > static_cast<std::uint64_t>(
                 std::numeric_limits<std::streamsize>::max())
          || !is_finite(block.time)
          || (index > 0 && block.accepted_step <= steps[index - 1]))
        return slide::Status::Invalid_parameters;
      const auto codec = static_cast<CompressionCodec>(block.codec);
      const auto payload_layout_status = detail::validateCompressedPayloadSize(
        codec, block.payload_bytes, buffer_layout);
      if (payload_layout_status != slide::Status::Success)
        return payload_layout_status;
      auto payload_view = std::span<std::byte>{ payload }.first(
        static_cast<std::size_t>(block.payload_bytes));
      const auto payload_status = readExact(input, payload_view, remaining);
      if (payload_status != slide::Status::Success)
        return payload_status;
      if (block.payload_crc32 != crc32(payload_view))
        return slide::Status::Invalid_parameters;
#if defined(SLIDE_WITH_ZSTD)
      if (codec == CompressionCodec::none) {
        // validateCompressedPayloadSize() fixes the none payload to raw_bytes.
        assert(payload_view.size() == shuffled.size());
        std::memcpy(shuffled.data(), payload_view.data(), payload_view.size());
      } else {
        assert(codec == CompressionCodec::zstd);
        const std::size_t decoded = ZSTD_decompress(shuffled.data(),
                                                    shuffled.size(),
                                                    payload_view.data(),
                                                    payload_view.size());
        if (ZSTD_isError(decoded) || decoded != shuffled.size())
          return slide::Status::Invalid_parameters;
      }
#else
      // The validated file header excludes zstd before allocation when zstd
      // support is absent; payload validation then fixes its size to raw_bytes.
      assert(codec == CompressionCodec::none);
      assert(payload_view.size() == shuffled.size());
      std::memcpy(shuffled.data(), payload_view.data(), payload_view.size());
#endif
      const auto unshuffle_status = byteUnshuffle(shuffled, raw, sizeof(real_t));
      // Both vectors were created at raw_bytes, are distinct allocations, and
      // raw_bytes is an exact multiple of sizeof(real_t).
      assert(unshuffle_status == slide::Status::Success);
      (void)unshuffle_status;
      if (block.raw_crc32 != crc32(raw))
        return slide::Status::Invalid_parameters;
      steps[index] = block.accepted_step;
      times[index] = block.time;
      std::memcpy(currents.data() + index * header.lanes,
                  raw.data(),
                  header.lanes * sizeof(real_t));
      std::memcpy(states.data() + index * state_values,
                  raw.data() + header.lanes * sizeof(real_t),
                  state_values * sizeof(real_t));
    }
    if (remaining != 0)
      return slide::Status::Invalid_parameters;

    close();
    rows_ = static_cast<int>(header.rows);
    lanes_ = static_cast<int>(header.lanes);
    stride_ = static_cast<int>(header.stride);
    state_values_ = state_values;
    steps_ = std::move(steps);
    times_ = std::move(times);
    currents_ = std::move(currents);
    states_ = std::move(states);
    valid_ = true;
  }
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  return allocationFailureStatus();
} catch (const std::length_error &) {
  return allocationFailureStatus();
}

void CompressedRecording::close()
{
  valid_ = false;
  rows_ = 0;
  lanes_ = 0;
  stride_ = 0;
  state_values_ = 0;
  std::vector<std::uint64_t>{}.swap(steps_);
  std::vector<real_t>{}.swap(times_);
  std::vector<real_t>{}.swap(currents_);
  std::vector<real_t>{}.swap(states_);
}

SnapshotView CompressedRecording::snapshot(std::size_t index) const
{
  assert(valid_ && index < steps_.size());
  return { .accepted_step = steps_[index],
           .time = times_[index],
           .current_density = std::span<const real_t>{ currents_ }.subspan(
             index * static_cast<std::size_t>(lanes_),
             static_cast<std::size_t>(lanes_)),
           .state = std::span<const real_t>{ states_ }.subspan(
             index * state_values_, state_values_) };
}

} // namespace slide::core
