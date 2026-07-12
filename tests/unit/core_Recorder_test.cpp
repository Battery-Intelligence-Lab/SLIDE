/**
 * @file core_Recorder_test.cpp
 * @brief Phase-6 lazy-observable and hardened recording gates.
 */

#include "../../src/core/EulerLegacy.hpp"
#include "../../src/core/Recorder.hpp"
#include "../support/CoreSpmTestHarness.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

using namespace slide;

namespace slide::core::detail {

struct RecorderTestAccess
{
  static Status writeCsvStream(Recorder &recorder, std::ostream &output)
  {
    return recorder.writeCsvStream(output);
  }

  static void failMappingFlush(Recorder &recorder)
  {
    recorder.mapping_flush_ = [](void *) { return false; };
  }
};

} // namespace slide::core::detail

namespace {

std::filesystem::path temporary(const std::string &suffix)
{
  return std::filesystem::temp_directory_path()
         / ("slide_core_recorder_" + suffix);
}

void flipByte(const std::filesystem::path &path, std::uint64_t offset)
{
  std::fstream stream(path, std::ios::binary | std::ios::in | std::ios::out);
  REQUIRE(stream.good());
  stream.seekg(static_cast<std::streamoff>(offset));
  char value{};
  stream.read(&value, 1);
  REQUIRE(stream.gcount() == 1);
  value ^= static_cast<char>(0x5a);
  stream.seekp(static_cast<std::streamoff>(offset));
  stream.write(&value, 1);
  REQUIRE(stream.good());
}

std::uint32_t testCrc32(std::span<const std::byte> bytes)
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
void writeScalar(std::span<std::byte> bytes, std::size_t offset, T value)
{
  static_assert(std::is_trivially_copyable_v<T>);
  REQUIRE(offset <= bytes.size());
  REQUIRE(sizeof(T) <= bytes.size() - offset);
  std::memcpy(bytes.data() + offset, &value, sizeof(value));
}

void sealHeader(std::span<std::byte> bytes)
{
  REQUIRE(bytes.size() >= 64);
  auto header = bytes.first(64);
  writeScalar<std::uint32_t>(header, 20, 0U);
  writeScalar<std::uint32_t>(header, 20, testCrc32(header));
}

std::vector<std::byte> readBytes(const std::filesystem::path &path)
{
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  REQUIRE(input.good());
  const auto length = input.tellg();
  REQUIRE(length >= 0);
  std::vector<std::byte> bytes(static_cast<std::size_t>(length));
  input.seekg(0);
  input.read(reinterpret_cast<char *>(bytes.data()), length);
  REQUIRE(input.good());
  return bytes;
}

void writeBytes(const std::filesystem::path &path,
                std::span<const std::byte>
                  bytes)
{
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  REQUIRE(output.good());
  output.write(reinterpret_cast<const char *>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  REQUIRE(output.good());
}

} // namespace

TEST_CASE("Recorder rejects invalid derived metadata and preserves thin ordering",
          "[core][recorder][validation][P9]")
{
  SECTION("invalid backpressure policy is rejected before configuration")
  {
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::Recorder recorder;
    CHECK(recorder.configure(
            batch,
            { .capacity = 1,
              .backpressure = static_cast<core::BackpressurePolicy>(255) })
          == Status::Invalid_parameters);
    CHECK_FALSE(recorder.configured());
  }

  SECTION("overflowing capacity is rejected before allocation")
  {
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::Recorder recorder;
    CHECK(recorder.configure(
            batch,
            { .capacity = std::numeric_limits<std::size_t>::max() })
          == Status::Invalid_parameters);
    CHECK_FALSE(recorder.configured());
  }

  SECTION("unconfigured calls and invalid sinks are rejected")
  {
    core::Recorder recorder;
    const std::array current{ 1.0, -1.0 };
    std::array<double, 2> voltage{};
    CHECK(recorder.record(0, current) == Status::Invalid_parameters);
    CHECK(recorder.terminalVoltage(0, voltage)
          == Status::Invalid_parameters);
    CHECK(recorder.writeCsv(temporary("unconfigured.csv"))
          == Status::Invalid_parameters);
    CHECK(recorder.writeBinary(temporary("unconfigured.slrec"))
          == Status::Invalid_parameters);

    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    REQUIRE(recorder.configure(batch, { .capacity = 1 })
            == Status::Success);
    CHECK(recorder.terminalVoltage(0, voltage)
          == Status::Invalid_parameters);
    CHECK(recorder.writeCsv(std::filesystem::temp_directory_path())
          == Status::Invalid_parameters);
    CHECK(recorder.writeBinary(std::filesystem::temp_directory_path())
          == Status::Invalid_parameters);
  }

  SECTION("thin mode tracks the last omitted cadence point")
  {
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::Recorder recorder;
    REQUIRE(recorder.configure(
              batch,
              { .capacity = 1,
                .backpressure = core::BackpressurePolicy::thin })
            == Status::Success);
    const std::array current{ 1.0, -1.0 };
    REQUIRE(recorder.record(0, current) == Status::Success);
    REQUIRE(recorder.record(2, current) == Status::Success);
    REQUIRE(recorder.thinnedSnapshots() == 1);
    CHECK(recorder.record(2, current) == Status::Invalid_parameters);
    CHECK(recorder.record(1, current) == Status::Invalid_parameters);
    CHECK(recorder.thinnedSnapshots() == 1);

    recorder.clear();
    CHECK(recorder.record(0, current) == Status::Success);
  }

  SECTION("non-finite elapsed time is rejected atomically")
  {
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::Recorder recorder;
    REQUIRE(recorder.configure(batch, { .capacity = 1 })
            == Status::Success);
    batch.state().at(batch.layout().elapsed_time, 0, 0) =
      std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
    const std::array current{ 1.0, -1.0 };
    CHECK(recorder.record(0, current) == Status::Invalid_parameters);
    CHECK(recorder.size() == 0);
  }

  SECTION("finite current and area cannot store infinite current density")
  {
    auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    input.design.electrode_area = 1e-300;
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::Recorder recorder;
    REQUIRE(recorder.configure(batch, { .capacity = 1 })
            == Status::Success);
    const std::array current{ 1e300, -1e300 };
    CHECK(recorder.record(0, current) == Status::Invalid_parameters);
    CHECK(recorder.size() == 0);
  }
}

TEST_CASE("Recorder checks binary layouts and propagates deterministic sink faults",
          "[core][recorder][fault-injection][coverage]")
{
  SECTION("the pure mmap layout rejects each overflow phase atomically")
  {
    core::detail::BinaryRecordingLayout layout{
      .current_bytes = 1,
      .state_bytes = 2,
      .record_bytes = 3,
      .table_bytes = 4,
      .data_offset = 5,
      .file_size = 6,
    };
    CHECK(core::detail::binaryRecordingLayout(
            1,
            std::numeric_limits<std::uint64_t>::max(),
            1,
            layout)
          == Status::Invalid_parameters);
    CHECK(layout.file_size == 6);

    CHECK(core::detail::binaryRecordingLayout(
            1,
            std::numeric_limits<std::uint64_t>::max() / 16,
            3,
            layout)
          == Status::Invalid_parameters);
    CHECK(layout.file_size == 6);

    REQUIRE(core::detail::binaryRecordingLayout(2, 10, 1, layout)
            == Status::Success);
    CHECK(layout.current_bytes == 2 * sizeof(double));
    CHECK(layout.state_bytes == 10 * sizeof(double));
    CHECK(layout.record_bytes == 16 + 12 * sizeof(double));
    CHECK(layout.data_offset % 64 == 0);
  }

  SECTION("a failed CSV stream reports the write failure")
  {
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::Recorder recorder;
    REQUIRE(recorder.configure(batch, { .capacity = 1 })
            == Status::Success);
    std::ostream failed{ nullptr };
    CHECK(core::detail::RecorderTestAccess::writeCsvStream(recorder, failed)
          == Status::Numerical_failure);
  }

  SECTION("an mmap synchronization failure is not reported as success")
  {
    const auto path = temporary("flush_failure.slrec");
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::Recorder recorder;
    REQUIRE(recorder.configure(batch, { .capacity = 1 })
            == Status::Success);
    const std::array current{ 1.0, -1.0 };
    REQUIRE(recorder.record(0, current) == Status::Success);
    core::detail::RecorderTestAccess::failMappingFlush(recorder);
    CHECK(recorder.writeBinary(path) == Status::Numerical_failure);
    std::filesystem::remove(path, ignored);
  }
}

TEST_CASE("P6-G1 recorded states lazily reproduce live voltage",
          "[core][recorder][lazy][P6-G1]")
{
  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .cadence = 2, .capacity = 3 })
          == Status::Success);
  const std::array current{ 8.0, -4.0 };
  std::array<double, 2> voltage_density_scratch{};
  std::array<double, 2> voltage0{};
  test_support::requireTerminalVoltage(
    batch,
    test_support::CurrentA{ std::span<const core::real_t>{ current } },
    0.0,
    voltage_density_scratch,
    voltage0);
  REQUIRE(recorder.record(0, current) == Status::Success);

  core::EulerLegacy stepper{ batch };
  const std::array density{ current[0] / batch.electrode_area(),
                            current[1] / batch.electrode_area() };
  REQUIRE(stepper.step(batch, density, 0.0, 1.0) == Status::Success);
  REQUIRE(recorder.record(1, current) == Status::Success); // cadence skip
  REQUIRE(stepper.step(batch, density, 1.0, 1.0) == Status::Success);
  std::array<double, 2> voltage2{};
  test_support::requireTerminalVoltage(
    batch,
    test_support::CurrentA{ std::span<const core::real_t>{ current } },
    0.0,
    voltage_density_scratch,
    voltage2);
  REQUIRE(recorder.record(2, current) == Status::Success);
  REQUIRE(recorder.size() == 2);

  // Mutate the subscribed batch; lazy reads must remain snapshot-pure.
  REQUIRE(stepper.step(batch, density, 2.0, 10.0) == Status::Success);
  std::array<double, 2> derived{};
  REQUIRE(recorder.terminalVoltage(0, derived) == Status::Success);
  CHECK(derived[0] == voltage0[0]);
  CHECK(derived[1] == voltage0[1]);
  REQUIRE(recorder.terminalVoltage(1, derived) == Status::Success);
  CHECK(derived[0] == voltage2[0]);
  CHECK(derived[1] == voltage2[1]);
  CHECK(recorder.snapshot(0).time == 0.0);
  CHECK(recorder.snapshot(1).time == 2.0);

  core::Recorder thin;
  REQUIRE(thin.configure(batch, { .cadence = 2, .capacity = 1, .backpressure = core::BackpressurePolicy::thin })
          == Status::Success);
  REQUIRE(thin.record(0, current) == Status::Success);
  REQUIRE(thin.record(1, current) == Status::Success);
  REQUIRE(thin.record(2, current) == Status::Success);
  CHECK(thin.size() == 1);
  CHECK(thin.thinnedSnapshots() == 1);

  core::Recorder stop;
  REQUIRE(stop.configure(batch, { .cadence = 1, .capacity = 1 })
          == Status::Success);
  REQUIRE(stop.record(0, current) == Status::Success);
  CHECK(stop.record(1, current) == Status::Numerical_failure);
}

TEST_CASE("P6-G1 CSV and mmap recordings preserve snapshots",
          "[core][recorder][csv][mmap][P6-G1]")
{
  const auto csv_path = temporary("valid.csv");
  const auto binary_path = temporary("valid.slrec");
  std::error_code ignored;
  std::filesystem::remove(csv_path, ignored);
  std::filesystem::remove(binary_path, ignored);

  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .capacity = 3 }) == Status::Success);
  core::EulerLegacy stepper{ batch };
  const std::array current{ 8.0, -4.0 };
  const std::array density{ current[0] / batch.electrode_area(),
                            current[1] / batch.electrode_area() };
  for (std::uint64_t step = 0; step < 3; ++step) {
    REQUIRE(recorder.record(step, current) == Status::Success);
    if (step + 1 < 3)
      REQUIRE(stepper.step(batch, density, static_cast<double>(step), 1.0)
              == Status::Success);
  }

  REQUIRE(recorder.writeCsv(csv_path) == Status::Success);
#if !defined(SLIDE_WITH_ARROW)
  CHECK(recorder.writeParquet(temporary("disabled.parquet"))
        == Status::NotImplementedYet);
#endif
  std::ifstream csv(csv_path, std::ios::binary);
  const std::string contents{ std::istreambuf_iterator<char>{ csv },
                              std::istreambuf_iterator<char>{} };
  CHECK(contents.starts_with("accepted_step,time_s,current_density_lane0_A_m2"));
  CHECK(contents.find("terminal_voltage_lane1_V") != std::string::npos);
  CHECK(static_cast<std::size_t>(std::count(contents.begin(), contents.end(), '\n'))
        == recorder.size() + 1);

  REQUIRE(recorder.writeBinary(binary_path) == Status::Success);
  core::BinaryRecording mapped;
  REQUIRE(mapped.open(binary_path) == Status::Success);
  REQUIRE(mapped.valid());
  REQUIRE(mapped.size() == recorder.size());
  CHECK(mapped.nRows() == recorder.nRows());
  CHECK(mapped.nLanes() == recorder.nLanes());
  CHECK(mapped.stride() == recorder.stride());
  for (std::size_t index = 0; index < recorder.size(); ++index) {
    const auto expected = recorder.snapshot(index);
    const auto actual = mapped.snapshot(index);
    CHECK(actual.accepted_step == expected.accepted_step);
    CHECK(actual.time == expected.time);
    CHECK(std::equal(actual.current_density.begin(), actual.current_density.end(), expected.current_density.begin(), expected.current_density.end()));
    CHECK(std::equal(actual.state.begin(), actual.state.end(), expected.state.begin(), expected.state.end()));
  }
  mapped.close();
  std::filesystem::remove(csv_path, ignored);
  std::filesystem::remove(binary_path, ignored);
}

TEST_CASE("P6-G1 mmap open rejects truncated CRC and offset corruption",
          "[core][recorder][mmap][hardened][P6-G1]")
{
  const auto valid = temporary("hardened_valid.slrec");
  const auto header_corrupt = temporary("header_corrupt.slrec");
  const auto offset_corrupt = temporary("offset_corrupt.slrec");
  const auto record_offset_corrupt = temporary("record_offset_corrupt.slrec");
  const auto truncated = temporary("truncated.slrec");
  std::error_code ignored;
  for (const auto &path : { valid,
                            header_corrupt,
                            offset_corrupt,
                            record_offset_corrupt,
                            truncated })
    std::filesystem::remove(path, ignored);

  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .capacity = 1 }) == Status::Success);
  const std::array current{ 8.0, -4.0 };
  REQUIRE(recorder.record(0, current) == Status::Success);
  REQUIRE(recorder.writeBinary(valid) == Status::Success);
  REQUIRE(std::filesystem::copy_file(valid, header_corrupt));
  REQUIRE(std::filesystem::copy_file(valid, offset_corrupt));
  REQUIRE(std::filesystem::copy_file(valid, record_offset_corrupt));
  REQUIRE(std::filesystem::copy_file(valid, truncated));
  flipByte(header_corrupt, 24);        // first dimension byte; CRC must catch it
  flipByte(offset_corrupt, 64);        // first absolute offset
  flipByte(record_offset_corrupt, 72); // second absolute offset
  const auto truncated_size = std::filesystem::file_size(truncated);
  REQUIRE(truncated_size > 64);
  std::filesystem::resize_file(truncated, truncated_size - 1);

  core::BinaryRecording recording;
  REQUIRE(recording.open(valid) == Status::Success);
  CHECK(recording.open(header_corrupt) == Status::Invalid_parameters);
  CHECK(recording.valid()); // failed open is atomic
  CHECK(recording.open(offset_corrupt) == Status::Invalid_parameters);
  CHECK(recording.open(record_offset_corrupt) == Status::Invalid_parameters);
  CHECK(recording.open(truncated) == Status::Invalid_parameters);
  recording.close();

  for (const auto &path : { valid,
                            header_corrupt,
                            offset_corrupt,
                            record_offset_corrupt,
                            truncated })
    std::filesystem::remove(path, ignored);
}

TEST_CASE("Binary recording rejects absent and header-short files",
          "[core][recorder][mmap][coverage]")
{
  const auto absent = temporary("absent.slrec");
  const auto short_file = temporary("short.slrec");
  std::error_code ignored;
  std::filesystem::remove(absent, ignored);
  std::filesystem::remove(short_file, ignored);
  {
    std::ofstream output(short_file, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.put('x');
  }
  core::BinaryRecording recording;
  CHECK(recording.open(absent) == Status::Invalid_parameters);
  CHECK(recording.open(short_file) == Status::Invalid_parameters);
  CHECK_FALSE(recording.valid());
  std::filesystem::remove(short_file, ignored);
}

TEST_CASE("Binary recording rejects overflowing layouts and trailing bytes",
          "[core][recorder][mmap][coverage]")
{
  const auto valid = temporary("layout_valid.slrec");
  const auto overflow = temporary("layout_overflow.slrec");
  const auto trailing = temporary("layout_trailing.slrec");
  std::error_code ignored;
  for (const auto &path : { valid, overflow, trailing })
    std::filesystem::remove(path, ignored);

  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::Recorder writer;
  REQUIRE(writer.configure(batch, { .capacity = 1 }) == Status::Success);
  const std::array current{ 1.0, -1.0 };
  REQUIRE(writer.record(0, current) == Status::Success);
  REQUIRE(writer.writeBinary(valid) == Status::Success);
  const auto original = readBytes(valid);
  REQUIRE(original.size() > 72);

  auto overflow_bytes = original;
  constexpr auto maximum_dimension =
    static_cast<std::uint32_t>(std::numeric_limits<int>::max());
  writeScalar<std::uint32_t>(overflow_bytes, 24, maximum_dimension);
  writeScalar<std::uint32_t>(overflow_bytes, 32, maximum_dimension);
  sealHeader(overflow_bytes);
  writeBytes(overflow, overflow_bytes);

  auto trailing_bytes = original;
  trailing_bytes.push_back(std::byte{});
  writeScalar<std::uint64_t>(trailing_bytes, 56, trailing_bytes.size());
  sealHeader(trailing_bytes);
  writeBytes(trailing, trailing_bytes);

  core::BinaryRecording recording;
  REQUIRE(recording.open(valid) == Status::Success);
  CHECK(recording.open(overflow) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(trailing) == Status::Invalid_parameters);
  CHECK(recording.valid());

  for (const auto &path : { valid, overflow, trailing })
    std::filesystem::remove(path, ignored);
}
