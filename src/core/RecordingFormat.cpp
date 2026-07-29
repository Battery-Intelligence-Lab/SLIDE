/**
 * @file RecordingFormat.cpp
 * @brief The CRC-hardened binary recording format: layout, mmap, writer, and reader.
 *
 * Owns: `detail::binaryRecordingLayout`, the synchronous on-disk header, the platform mmap
 * boundary, `Recorder::writeBinary`, and every `BinaryRecording` member. Implements PLAN.md §3.7.
 * Cold: file I/O only. Writer and reader live together because they share one packed-header,
 * layout, and mapping contract. Checksum, byte-order, and allocation-failure facts are shared
 * with compressed recording through `detail/RecordingFormatCommon.hpp`.
 * @surface internal
 */

#include "Recorder.hpp"
#include "detail/CheckedArithmetic.hpp"
#include "detail/RecordingFormatCommon.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <fstream>
#include <limits>
#include <new>
#include <stdexcept>
#include <system_error>
#include <utility>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace slide::core {

using detail::allocationFailureStatus;
using detail::checkedAdd;
using detail::checkedMultiply;
using detail::align64;
using detail::endian_marker;
using detail::headerCrc;

namespace {

  constexpr std::array<char, 8> recording_magic{ 'S', 'L', 'I', 'D', 'E', 'R', 'E', 'C' };
  constexpr std::uint16_t format_major = 1;
  constexpr std::uint16_t format_minor = 0;
  constexpr std::size_t header_bytes = 64;

#pragma pack(push, 1)
  struct RecordingHeader
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
    std::uint32_t snapshots{};
    std::uint64_t offset_table{};
    std::uint64_t data_offset{};
    std::uint64_t file_size{};
  };
#pragma pack(pop)
  static_assert(sizeof(RecordingHeader) == header_bytes);

  template <class T>
  void store(std::byte *destination, const T &value)
  {
    std::memcpy(destination, &value, sizeof(T));
  }

  template <class T>
  T load(const std::byte *source)
  {
    T value{};
    std::memcpy(&value, source, sizeof(T));
    return value;
  }

} // namespace

slide::Status detail::binaryRecordingLayout(
  std::uint64_t lanes,
  std::uint64_t state_values,
  std::uint64_t snapshots,
  BinaryRecordingLayout &output)
{
  std::uint64_t current_bytes{}, state_bytes{}, record_bytes{};
  std::uint64_t table_entries{}, table_bytes{};
  if (!checkedMultiply(lanes, std::uint64_t{ sizeof(real_t) }, current_bytes)
      || !checkedMultiply(
        state_values, std::uint64_t{ sizeof(real_t) }, state_bytes)
      || !checkedAdd(std::uint64_t{ 16 }, current_bytes, record_bytes)
      || !checkedAdd(record_bytes, state_bytes, record_bytes)
      || !checkedAdd(snapshots, std::uint64_t{ 1 }, table_entries)
      || !checkedMultiply(
        table_entries, std::uint64_t{ sizeof(std::uint64_t) }, table_bytes))
    return slide::Status::Invalid_parameters;

  std::uint64_t table_end{}, data_offset{}, records_bytes{}, file_size{};
  if (!checkedAdd(std::uint64_t{ header_bytes }, table_bytes, table_end)
      || !align64(table_end, data_offset)
      || !checkedMultiply(snapshots, record_bytes, records_bytes)
      || !checkedAdd(data_offset, records_bytes, file_size)
      || file_size > std::numeric_limits<std::size_t>::max())
    return slide::Status::Invalid_parameters;

  output = { .current_bytes = current_bytes,
             .state_bytes = state_bytes,
             .record_bytes = record_bytes,
             .table_bytes = table_bytes,
             .data_offset = data_offset,
             .file_size = file_size };
  return slide::Status::Success;
}

struct BinaryRecording::Mapping
{
  std::byte *data{};
  std::size_t size{};
  bool writable{};
#if defined(_WIN32)
  HANDLE file{ INVALID_HANDLE_VALUE };
  HANDLE mapping{};
#else
  int file{ -1 };
#endif

  ~Mapping() { reset(); }

  void reset()
  {
#if defined(_WIN32)
    if (data != nullptr)
      UnmapViewOfFile(data);
    if (mapping != nullptr)
      CloseHandle(mapping);
    if (file != INVALID_HANDLE_VALUE)
      CloseHandle(file);
    file = INVALID_HANDLE_VALUE;
    mapping = nullptr;
#else
    if (data != nullptr)
      munmap(data, size);
    if (file >= 0)
      ::close(file);
    file = -1;
#endif
    data = nullptr;
    size = 0;
    writable = false;
  }

  bool flush()
  {
    if (!writable || data == nullptr)
      return false;
#if defined(_WIN32)
    return FlushViewOfFile(data, size) != 0 && FlushFileBuffers(file) != 0;
#else
    return msync(data, size, MS_SYNC) == 0;
#endif
  }
};

bool Recorder::flushMapping(void *mapping)
{
  assert(mapping != nullptr);
  return static_cast<BinaryRecording::Mapping *>(mapping)->flush();
}

namespace {

  bool mapWritable(const std::filesystem::path &path, std::size_t size,
                   BinaryRecording::Mapping &output)
  {
    output.reset();
    if (size == 0)
      return false;
#if defined(_WIN32)
    if (size > static_cast<std::size_t>(std::numeric_limits<LONGLONG>::max()))
      return false;
    output.file = CreateFileW(path.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (output.file == INVALID_HANDLE_VALUE)
      return false;
    LARGE_INTEGER length{};
    length.QuadPart = static_cast<LONGLONG>(size);
    if (!SetFilePointerEx(output.file, length, nullptr, FILE_BEGIN)
        || !SetEndOfFile(output.file)) {
      output.reset();
      return false;
    }
    output.mapping = CreateFileMappingW(output.file, nullptr, PAGE_READWRITE, 0, 0, nullptr);
    if (output.mapping == nullptr) {
      output.reset();
      return false;
    }
    output.data = static_cast<std::byte *>(
      MapViewOfFile(output.mapping, FILE_MAP_ALL_ACCESS, 0, 0, size));
#else
    output.file = ::open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (output.file < 0 || ftruncate(output.file, static_cast<off_t>(size)) != 0) {
      output.reset();
      return false;
    }
    void *mapped = mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_SHARED, output.file, 0);
    output.data = mapped == MAP_FAILED ? nullptr : static_cast<std::byte *>(mapped);
#endif
    if (output.data == nullptr) {
      output.reset();
      return false;
    }
    output.size = size;
    output.writable = true;
    return true;
  }

  bool mapReadOnly(const std::filesystem::path &path,
                   BinaryRecording::Mapping &output)
  {
    output.reset();
#if defined(_WIN32)
    output.file = CreateFileW(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (output.file == INVALID_HANDLE_VALUE)
      return false;
    LARGE_INTEGER length{};
    if (!GetFileSizeEx(output.file, &length) || length.QuadPart <= 0
        || static_cast<unsigned long long>(length.QuadPart)
             > std::numeric_limits<std::size_t>::max()) {
      output.reset();
      return false;
    }
    output.size = static_cast<std::size_t>(length.QuadPart);
    output.mapping = CreateFileMappingW(output.file, nullptr, PAGE_READONLY, 0, 0, nullptr);
    if (output.mapping == nullptr) {
      output.reset();
      return false;
    }
    output.data = static_cast<std::byte *>(
      MapViewOfFile(output.mapping, FILE_MAP_READ, 0, 0, output.size));
#else
    output.file = ::open(path.c_str(), O_RDONLY);
    struct stat info{};
    if (output.file < 0 || fstat(output.file, &info) != 0 || info.st_size <= 0
        || static_cast<std::uintmax_t>(info.st_size)
             > std::numeric_limits<std::size_t>::max()) {
      output.reset();
      return false;
    }
    output.size = static_cast<std::size_t>(info.st_size);
    void *mapped = mmap(nullptr, output.size, PROT_READ, MAP_PRIVATE, output.file, 0);
    output.data = mapped == MAP_FAILED ? nullptr : static_cast<std::byte *>(mapped);
#endif
    if (output.data == nullptr) {
      output.reset();
      return false;
    }
    return true;
  }

} // namespace

slide::Status Recorder::writeBinary(const std::filesystem::path &path) const
{
  if (!configured() || count_ > std::numeric_limits<std::uint32_t>::max())
    return slide::Status::Invalid_parameters;
  detail::BinaryRecordingLayout layout;
  const auto layout_status = detail::binaryRecordingLayout(
    static_cast<std::uint64_t>(lanes_),
    static_cast<std::uint64_t>(state_values_),
    static_cast<std::uint64_t>(count_),
    layout);
  if (layout_status != slide::Status::Success)
    return layout_status;

  BinaryRecording::Mapping mapping;
  if (!mapWritable(path, static_cast<std::size_t>(layout.file_size), mapping))
    return slide::Status::Invalid_parameters;
  std::memset(mapping.data, 0, mapping.size);
  auto *table = mapping.data + header_bytes;
  for (std::size_t index = 0; index <= count_; ++index) {
    const std::uint64_t offset = layout.data_offset
                                 + static_cast<std::uint64_t>(index)
                                     * layout.record_bytes;
    store(table + index * sizeof(offset), offset);
  }
  for (std::size_t index = 0; index < count_; ++index) {
    const auto recorded = snapshot(index);
    std::byte *cursor = mapping.data + layout.data_offset
                        + static_cast<std::uint64_t>(index)
                            * layout.record_bytes;
    store(cursor, recorded.accepted_step);
    cursor += sizeof(recorded.accepted_step);
    store(cursor, recorded.time);
    cursor += sizeof(recorded.time);
    std::memcpy(cursor, recorded.current_density.data(), recorded.current_density.size_bytes());
    cursor += recorded.current_density.size_bytes();
    std::memcpy(cursor, recorded.state.data(), recorded.state.size_bytes());
  }

  RecordingHeader header{ .magic = recording_magic,
                          .major = format_major,
                          .minor = format_minor,
                          .endian = endian_marker,
                          .header_size = static_cast<std::uint32_t>(header_bytes),
                          .header_crc32 = 0,
                          .rows = static_cast<std::uint32_t>(rows_),
                          .lanes = static_cast<std::uint32_t>(lanes_),
                          .stride = static_cast<std::uint32_t>(stride_),
                          .snapshots = static_cast<std::uint32_t>(count_),
                          .offset_table = header_bytes,
                          .data_offset = layout.data_offset,
                          .file_size = layout.file_size };
  header.header_crc32 = headerCrc(header);
  store(mapping.data, header);
  return mapping_flush_(&mapping) ? slide::Status::Success
                                  : slide::Status::Numerical_failure;
}

BinaryRecording::BinaryRecording() = default;
BinaryRecording::~BinaryRecording() = default;
BinaryRecording::BinaryRecording(BinaryRecording &&) noexcept = default;
BinaryRecording &BinaryRecording::operator=(BinaryRecording &&) noexcept = default;

void BinaryRecording::close()
{
  mapping_.reset();
  snapshots_ = 0;
  rows_ = lanes_ = stride_ = 0;
  table_offset_ = 0;
  state_values_ = 0;
}

bool BinaryRecording::valid() const
{
  return mapping_ != nullptr && mapping_->data != nullptr;
}

slide::Status BinaryRecording::open(const std::filesystem::path &path)
try {
  auto candidate = std::make_unique<Mapping>();
  if (!mapReadOnly(path, *candidate) || candidate->size < header_bytes)
    return slide::Status::Invalid_parameters;
  const RecordingHeader header = load<RecordingHeader>(candidate->data);
  if (header.magic != recording_magic || header.major != format_major
      || header.minor > format_minor || header.endian != endian_marker
      || header.header_size != header_bytes
      || headerCrc(header) != header.header_crc32
      || header.rows == 0 || header.lanes == 0 || header.stride < header.lanes
      || header.rows > static_cast<std::uint32_t>(std::numeric_limits<int>::max())
      || header.lanes > static_cast<std::uint32_t>(std::numeric_limits<int>::max())
      || header.stride > static_cast<std::uint32_t>(std::numeric_limits<int>::max())
      || header.file_size != candidate->size
      || header.offset_table != header_bytes
      || header.data_offset % StateArena::alignment != 0)
    return slide::Status::Invalid_parameters;
  std::uint64_t state_values{}, current_bytes{}, state_bytes{}, record_bytes{}, entries{}, table_bytes{}, table_end{};
  if (!checkedMultiply(static_cast<std::uint64_t>(header.rows),
                       static_cast<std::uint64_t>(header.stride),
                       state_values)
      || !checkedMultiply(static_cast<std::uint64_t>(header.lanes),
                          std::uint64_t{ sizeof(real_t) },
                          current_bytes)
      || !checkedMultiply(state_values, std::uint64_t{ sizeof(real_t) }, state_bytes)
      || !checkedAdd(std::uint64_t{ 16 }, current_bytes, record_bytes)
      || !checkedAdd(record_bytes, state_bytes, record_bytes)
      || !checkedAdd(static_cast<std::uint64_t>(header.snapshots),
                     std::uint64_t{ 1 },
                     entries)
      || !checkedMultiply(entries, std::uint64_t{ sizeof(std::uint64_t) }, table_bytes)
      || !checkedAdd(header.offset_table, table_bytes, table_end)
      || table_end > header.data_offset || header.data_offset > header.file_size
      || state_values > std::numeric_limits<std::size_t>::max())
    return slide::Status::Invalid_parameters;
  const auto *table = candidate->data + header.offset_table;
  std::uint64_t previous = load<std::uint64_t>(table);
  if (previous != header.data_offset)
    return slide::Status::Invalid_parameters;
  for (std::uint64_t index = 1; index <= header.snapshots; ++index) {
    const std::uint64_t next = load<std::uint64_t>(
      table + index * sizeof(std::uint64_t));
    if (next < previous || next > header.file_size
        || next - previous != record_bytes)
      return slide::Status::Invalid_parameters;
    previous = next;
  }
  if (previous != header.file_size)
    return slide::Status::Invalid_parameters;

  close();
  mapping_ = std::move(candidate);
  snapshots_ = header.snapshots;
  rows_ = static_cast<int>(header.rows);
  lanes_ = static_cast<int>(header.lanes);
  stride_ = static_cast<int>(header.stride);
  table_offset_ = header.offset_table;
  state_values_ = static_cast<std::size_t>(state_values);
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  return allocationFailureStatus();
} catch (const std::length_error &) {
  return allocationFailureStatus();
}

SnapshotView BinaryRecording::snapshot(std::size_t index) const
{
  assert(valid() && index < snapshots_);
  const auto *table = mapping_->data + table_offset_;
  const std::uint64_t offset = load<std::uint64_t>(
    table + index * sizeof(std::uint64_t));
  const std::byte *cursor = mapping_->data + offset;
  const std::uint64_t accepted_step = load<std::uint64_t>(cursor);
  cursor += sizeof(accepted_step);
  const real_t time = load<real_t>(cursor);
  cursor += sizeof(time);
  const auto *current = reinterpret_cast<const real_t *>(cursor);
  cursor += static_cast<std::size_t>(lanes_) * sizeof(real_t);
  const auto *state = reinterpret_cast<const real_t *>(cursor);
  return { .accepted_step = accepted_step,
           .time = time,
           .current_density = { current, static_cast<std::size_t>(lanes_) },
           .state = { state, state_values_ } };
}


} // namespace slide::core
