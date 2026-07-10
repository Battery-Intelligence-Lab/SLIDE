/**
 * @file Recorder.cpp
 * @brief Snapshot recorder, CSV output, and CRC-hardened mmap recording format.
 */

#include "Recorder.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <new>
#include <string>
#include <system_error>
#include <utility>

#if defined(SLIDE_WITH_ARROW)
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

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
namespace {

  constexpr std::array<char, 8> recording_magic{ 'S', 'L', 'I', 'D', 'E', 'R', 'E', 'C' };
  constexpr std::uint32_t endian_marker = 0x01020304U;
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

  template <class T>
  bool checkedAdd(T a, T b, T &result)
  {
    if (a > std::numeric_limits<T>::max() - b)
      return false;
    result = a + b;
    return true;
  }

  template <class T>
  bool checkedMultiply(T a, T b, T &result)
  {
    if (a != 0 && b > std::numeric_limits<T>::max() / a)
      return false;
    result = a * b;
    return true;
  }

  bool align64(std::uint64_t value, std::uint64_t &result)
  {
    std::uint64_t enlarged{};
    if (!checkedAdd(value, std::uint64_t{ 63 }, enlarged))
      return false;
    result = enlarged & ~std::uint64_t{ 63 };
    return true;
  }

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

slide::Status Recorder::configure(SpmBatch &batch, RecorderConfig config)
{
  if (!batch.valid() || config.cadence == 0 || config.capacity == 0)
    return slide::Status::Invalid_parameters;
  const std::size_t state_values = batch.state().size();
  const auto lanes = static_cast<std::size_t>(batch.n_lanes());
  std::size_t state_storage{}, current_storage{};
  if (!checkedMultiply(config.capacity, state_values, state_storage)
      || !checkedMultiply(config.capacity, lanes, current_storage))
    return slide::Status::Invalid_parameters;
  try {
    std::vector<std::uint64_t> steps(config.capacity);
    std::vector<real_t> times(config.capacity);
    std::vector<real_t> currents(current_storage);
    std::vector<real_t> states(state_storage);
    batch_ = &batch;
    config_ = config;
    rows_ = batch.state().n_rows();
    lanes_ = batch.state().n_lanes();
    stride_ = batch.state().stride();
    state_values_ = state_values;
    count_ = 0;
    thinned_ = 0;
    accepted_steps_ = std::move(steps);
    times_ = std::move(times);
    current_density_ = std::move(currents);
    states_ = std::move(states);
  } catch (const std::bad_alloc &) {
    return slide::Status::Numerical_failure;
  }
  return slide::Status::Success;
}

slide::Status Recorder::record(std::uint64_t accepted_step,
                               std::span<const real_t>
                                 total_current_A)
{
  if (!configured() || total_current_A.size() != static_cast<std::size_t>(lanes_))
    return slide::Status::Invalid_parameters;
  if (accepted_step % config_.cadence != 0)
    return slide::Status::Success;
  if (count_ > 0 && accepted_step <= accepted_steps_[count_ - 1])
    return slide::Status::Invalid_parameters;
  for (const real_t current : total_current_A)
    if (!is_finite(current))
      return slide::Status::Invalid_parameters;
  if (count_ == config_.capacity) {
    if (config_.backpressure == BackpressurePolicy::thin) {
      ++thinned_;
      return slide::Status::Success;
    }
    return slide::Status::Numerical_failure;
  }
  accepted_steps_[count_] = accepted_step;
  times_[count_] = batch_->state().at(batch_->layout().elapsed_time, 0, 0);
  const std::size_t current_offset = count_ * static_cast<std::size_t>(lanes_);
  for (int lane = 0; lane < lanes_; ++lane)
    current_density_[current_offset + static_cast<std::size_t>(lane)] = total_current_A[static_cast<std::size_t>(lane)] / batch_->electrode_area();
  std::memcpy(states_.data() + count_ * state_values_,
              batch_->state().raw().data(),
              state_values_ * sizeof(real_t));
  ++count_;
  return slide::Status::Success;
}

SnapshotView Recorder::snapshot(std::size_t index) const
{
  assert(index < count_);
  return { .accepted_step = accepted_steps_[index],
           .time = times_[index],
           .current_density = std::span<const real_t>{ current_density_ }.subspan(
             index * static_cast<std::size_t>(lanes_),
             static_cast<std::size_t>(lanes_)),
           .state = std::span<const real_t>{ states_ }.subspan(
             index * state_values_, state_values_) };
}

void Recorder::clear()
{
  count_ = 0;
  thinned_ = 0;
}

slide::Status Recorder::terminalVoltage(std::size_t index,
                                        std::span<real_t>
                                          output)
{
  if (!configured() || index >= count_
      || output.size() != static_cast<std::size_t>(lanes_))
    return slide::Status::Invalid_parameters;
  const auto recorded = snapshot(index);
  return batch_->terminalVoltageAt(recorded.state, recorded.current_density, output);
}

slide::Status Recorder::writeCsv(const std::filesystem::path &path)
{
  if (!configured())
    return slide::Status::Invalid_parameters;
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output)
    return slide::Status::Invalid_parameters;
  output << "accepted_step,time_s";
  for (int lane = 0; lane < lanes_; ++lane)
    output << ",current_density_lane" << lane << "_A_m2";
  for (int row = 0; row < rows_; ++row)
    for (int lane = 0; lane < lanes_; ++lane)
      output << ",state_r" << row << "_lane" << lane;
  for (int lane = 0; lane < lanes_; ++lane)
    output << ",terminal_voltage_lane" << lane << "_V";
  output << '\n'
         << std::setprecision(std::numeric_limits<real_t>::max_digits10)
         << std::scientific;
  std::vector<real_t> voltage(static_cast<std::size_t>(lanes_));
  for (std::size_t index = 0; index < count_; ++index) {
    const auto recorded = snapshot(index);
    const auto status = terminalVoltage(index, voltage);
    if (status != slide::Status::Success)
      return status;
    output << recorded.accepted_step << ',' << recorded.time;
    for (const real_t current : recorded.current_density)
      output << ',' << current;
    for (int row = 0; row < rows_; ++row)
      for (int lane = 0; lane < lanes_; ++lane)
        output << ',' << recorded.state[static_cast<std::size_t>(row * stride_ + lane)];
    for (const real_t value : voltage)
      output << ',' << value;
    output << '\n';
  }
  return output.good() ? slide::Status::Success
                       : slide::Status::Numerical_failure;
}

slide::Status Recorder::writeBinary(const std::filesystem::path &path) const
{
  if (!configured() || count_ > std::numeric_limits<std::uint32_t>::max())
    return slide::Status::Invalid_parameters;
  std::uint64_t current_bytes{}, state_bytes{}, record_bytes{}, table_entries{}, table_bytes{};
  if (!checkedMultiply(static_cast<std::uint64_t>(lanes_),
                       std::uint64_t{ sizeof(real_t) },
                       current_bytes)
      || !checkedMultiply(static_cast<std::uint64_t>(state_values_),
                          std::uint64_t{ sizeof(real_t) },
                          state_bytes)
      || !checkedAdd(std::uint64_t{ 16 }, current_bytes, record_bytes)
      || !checkedAdd(record_bytes, state_bytes, record_bytes)
      || !checkedAdd(static_cast<std::uint64_t>(count_), std::uint64_t{ 1 }, table_entries)
      || !checkedMultiply(table_entries, std::uint64_t{ sizeof(std::uint64_t) }, table_bytes))
    return slide::Status::Invalid_parameters;
  std::uint64_t table_end{}, data_offset{}, records_bytes{}, file_size{};
  if (!checkedAdd(std::uint64_t{ header_bytes }, table_bytes, table_end)
      || !align64(table_end, data_offset)
      || !checkedMultiply(static_cast<std::uint64_t>(count_), record_bytes, records_bytes)
      || !checkedAdd(data_offset, records_bytes, file_size)
      || file_size > std::numeric_limits<std::size_t>::max())
    return slide::Status::Invalid_parameters;

  BinaryRecording::Mapping mapping;
  if (!mapWritable(path, static_cast<std::size_t>(file_size), mapping))
    return slide::Status::Invalid_parameters;
  std::memset(mapping.data, 0, mapping.size);
  auto *table = mapping.data + header_bytes;
  for (std::size_t index = 0; index <= count_; ++index) {
    const std::uint64_t offset = data_offset
                                 + static_cast<std::uint64_t>(index) * record_bytes;
    store(table + index * sizeof(offset), offset);
  }
  for (std::size_t index = 0; index < count_; ++index) {
    const auto recorded = snapshot(index);
    std::byte *cursor = mapping.data + data_offset
                        + static_cast<std::uint64_t>(index) * record_bytes;
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
                          .data_offset = data_offset,
                          .file_size = file_size };
  header.header_crc32 = crc32(std::as_bytes(std::span{ &header, 1 }));
  store(mapping.data, header);
  return mapping.flush() ? slide::Status::Success
                         : slide::Status::Numerical_failure;
}

slide::Status Recorder::writeParquet(const std::filesystem::path &path)
{
#if !defined(SLIDE_WITH_ARROW)
  (void)path;
  return slide::Status::NotImplementedYet;
#else
  if (!configured())
    return slide::Status::Invalid_parameters;
  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> columns;
  auto append_uint64 = [&](const std::string &name, auto value) {
    arrow::UInt64Builder builder;
    for (std::size_t index = 0; index < count_; ++index)
      if (!builder.Append(value(index)).ok())
        return false;
    std::shared_ptr<arrow::Array> array;
    if (!builder.Finish(&array).ok())
      return false;
    fields.push_back(arrow::field(name, arrow::uint64()));
    columns.push_back(std::move(array));
    return true;
  };
  auto append_double = [&](const std::string &name, auto value) {
    arrow::DoubleBuilder builder;
    for (std::size_t index = 0; index < count_; ++index)
      if (!builder.Append(value(index)).ok())
        return false;
    std::shared_ptr<arrow::Array> array;
    if (!builder.Finish(&array).ok())
      return false;
    fields.push_back(arrow::field(name, arrow::float64()));
    columns.push_back(std::move(array));
    return true;
  };
  if (!append_uint64("accepted_step", [&](std::size_t i) { return accepted_steps_[i]; })
      || !append_double("time_s", [&](std::size_t i) { return times_[i]; }))
    return slide::Status::Numerical_failure;
  for (int lane = 0; lane < lanes_; ++lane) {
    if (!append_double("current_density_lane" + std::to_string(lane) + "_A_m2",
                       [&](std::size_t i) {
                         return snapshot(i).current_density[static_cast<std::size_t>(lane)];
                       }))
      return slide::Status::Numerical_failure;
  }
  for (int row = 0; row < rows_; ++row)
    for (int lane = 0; lane < lanes_; ++lane)
      if (!append_double("state_r" + std::to_string(row) + "_lane"
                           + std::to_string(lane),
                         [&](std::size_t i) {
                           return snapshot(i).state[static_cast<std::size_t>(
                             row * stride_ + lane)];
                         }))
        return slide::Status::Numerical_failure;
  std::vector<std::vector<real_t>> voltage(
    static_cast<std::size_t>(lanes_), std::vector<real_t>(count_));
  std::vector<real_t> one_voltage(static_cast<std::size_t>(lanes_));
  for (std::size_t index = 0; index < count_; ++index) {
    if (terminalVoltage(index, one_voltage) != slide::Status::Success)
      return slide::Status::Numerical_failure;
    for (int lane = 0; lane < lanes_; ++lane)
      voltage[static_cast<std::size_t>(lane)][index] = one_voltage[static_cast<std::size_t>(lane)];
  }
  for (int lane = 0; lane < lanes_; ++lane)
    if (!append_double("terminal_voltage_lane" + std::to_string(lane) + "_V",
                       [&](std::size_t i) {
                         return voltage[static_cast<std::size_t>(lane)][i];
                       }))
      return slide::Status::Numerical_failure;

  const auto table = arrow::Table::Make(arrow::schema(std::move(fields)),
                                        std::move(columns));
  auto output_result = arrow::io::FileOutputStream::Open(path.string());
  if (!output_result.ok())
    return slide::Status::Invalid_parameters;
  auto output = *output_result;
  const auto status = parquet::arrow::WriteTable(
    *table, arrow::default_memory_pool(), output, std::max<std::int64_t>(1, static_cast<std::int64_t>(count_)));
  if (!status.ok() || !output->Close().ok())
    return slide::Status::Numerical_failure;
  return slide::Status::Success;
#endif
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
{
  auto candidate = std::make_unique<Mapping>();
  if (!mapReadOnly(path, *candidate) || candidate->size < header_bytes)
    return slide::Status::Invalid_parameters;
  const RecordingHeader header = load<RecordingHeader>(candidate->data);
  RecordingHeader crc_header = header;
  const std::uint32_t expected_crc = std::exchange(crc_header.header_crc32, 0U);
  if (header.magic != recording_magic || header.major != format_major
      || header.minor > format_minor || header.endian != endian_marker
      || header.header_size != header_bytes || expected_crc == 0
      || crc32(std::as_bytes(std::span{ &crc_header, 1 })) != expected_crc
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
