/**
 * @file core_AsyncRecorder_test.cpp
 * @brief P8-G3 async ring, backpressure, compression, and corruption gates.
 */

#include "../../src/core/AsyncRecorder.hpp"
#include "../../src/core/EulerLegacy.hpp"
#include "../support/CoreSpmTestHarness.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <iterator>
#include <limits>
#include <mutex>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

using namespace slide;

namespace slide::core::detail {

struct AsyncRecorderTestAccess
{
  static void failThreadCreation(AsyncRecorder &recorder)
  {
    recorder.drain_thread_factory_ = [](AsyncRecorder &) -> std::thread {
      throw std::runtime_error{ "injected thread construction failure" };
    };
  }

  static void restoreThreadCreation(AsyncRecorder &recorder)
  {
    recorder.drain_thread_factory_ = &AsyncRecorder::makeDrainThread;
  }

  static void failPlaceholderWrite(AsyncRecorder &recorder)
  {
    recorder.before_placeholder_write_ = &failOutput;
  }

  static void failFinalizeTell(AsyncRecorder &recorder)
  {
    recorder.before_finalize_tell_ = &failOutput;
  }

  static void failFinalizeWrite(AsyncRecorder &recorder)
  {
    recorder.before_finalize_write_ = &failOutput;
  }

  static void failBlockWrite(AsyncRecorder &recorder)
  {
    recorder.before_block_write_ = &failOutput;
  }

  static void truncateShuffleBuffer(AsyncRecorder &recorder)
  {
    // The worker cannot touch this slot until enqueue() publishes ready under
    // mutex_; this mutation is sequenced before that release.
    recorder.slots_.front().shuffled.pop_back();
  }

  static Status readExact(std::istream &input,
                          std::span<std::byte>
                            destination,
                          std::uint64_t &remaining)
  {
    return CompressedRecording::readExact(input, destination, remaining);
  }

  static std::size_t retainedBytes(const CompressedRecording &recording)
  {
    std::size_t bytes{};
    const auto status = compressedRecordingRetainedBytes(
      recording.steps_.capacity(),
      recording.times_.capacity(),
      recording.currents_.capacity(),
      recording.states_.capacity(),
      bytes);
    assert(status == Status::Success);
    return bytes;
  }

private:
  static void failOutput(std::ofstream &output)
  {
    output.setstate(std::ios::badbit);
  }
};

} // namespace slide::core::detail

namespace {

std::filesystem::path temporary(std::string_view suffix)
{
  return std::filesystem::temp_directory_path()
         / ("slide_async_" + std::string{ suffix });
}

core::CompressionCodec testCodec()
{
  return core::compressionCodecAvailable(core::CompressionCodec::zstd)
           ? core::CompressionCodec::zstd
           : core::CompressionCodec::none;
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

template <class T>
T readScalar(std::span<const std::byte> bytes, std::size_t offset)
{
  static_assert(std::is_trivially_copyable_v<T>);
  REQUIRE(offset <= bytes.size());
  REQUIRE(sizeof(T) <= bytes.size() - offset);
  T value{};
  std::memcpy(&value, bytes.data() + offset, sizeof(value));
  return value;
}

void sealHeader(std::span<std::byte> bytes,
                std::size_t offset,
                std::size_t crc_offset)
{
  REQUIRE(offset <= bytes.size());
  REQUIRE(64 <= bytes.size() - offset);
  auto header = bytes.subspan(offset, 64);
  writeScalar<std::uint32_t>(header, crc_offset, 0U);
  writeScalar<std::uint32_t>(header, crc_offset, testCrc32(header));
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

std::vector<std::byte> readBytes(const std::filesystem::path &path)
{
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  REQUIRE(input.good());
  const auto length = input.tellg();
  REQUIRE(length >= 0);
  std::vector<std::byte> bytes(static_cast<std::size_t>(length));
  input.seekg(0);
  input.read(reinterpret_cast<char *>(bytes.data()),
             static_cast<std::streamsize>(length));
  REQUIRE(input.good());
  return bytes;
}

template <class Function>
Status checkedShuffleCall(Function &&function)
{
  using Result = std::invoke_result_t<Function &>;
  if constexpr (std::is_same_v<Result, Status>)
    return function();
  else {
    function();
    return Status::Success;
  }
}

} // namespace

TEST_CASE("Async recorder rejects adversarial metadata without exceptions",
          "[core][async-recorder][validation][P9]")
{
  std::error_code ignored;

  SECTION("impossible snapshot count is rejected before allocation")
  {
    const auto path = temporary("p9_count.slcmp");
    std::filesystem::remove(path, ignored);
    std::array<std::byte, 64> header{};
    constexpr std::array magic{ 'S', 'L', 'I', 'D', 'E', 'C', 'M', 'P' };
    std::memcpy(header.data(), magic.data(), magic.size());
    writeScalar<std::uint16_t>(header, 8, 1U);
    writeScalar<std::uint16_t>(header, 10, 0U);
    writeScalar<std::uint32_t>(header, 12, 0x01020304U);
    writeScalar<std::uint32_t>(header, 16, 64U);
    writeScalar<std::uint32_t>(header, 24, 1U);
    writeScalar<std::uint32_t>(header, 28, 1U);
    writeScalar<std::uint32_t>(header, 32, 1U);
    writeScalar<std::uint32_t>(header, 36, 0U);
    writeScalar<std::uint64_t>(
      header, 40, std::numeric_limits<std::uint64_t>::max());
    writeScalar<std::uint64_t>(header, 48, header.size());
    sealHeader(header, 0, 20);
    writeBytes(path, header);

    core::CompressedRecording recording;
    Status status = Status::Success;
    CHECK_NOTHROW(status = recording.open(path));
    CHECK(status == Status::Invalid_parameters);
    CHECK_FALSE(recording.valid());
    std::filesystem::remove(path, ignored);
  }

  SECTION("checked buffer layout rejects each overflowing intermediate")
  {
    core::detail::AsyncBufferLayout layout;
    CHECK(core::detail::asyncBufferLayout(
            std::numeric_limits<std::size_t>::max(),
            1,
            core::CompressionCodec::none,
            layout)
          == Status::Invalid_parameters);
    CHECK(core::detail::asyncBufferLayout(
            std::numeric_limits<std::size_t>::max() / sizeof(double) + 1,
            0,
            core::CompressionCodec::none,
            layout)
          == Status::Invalid_parameters);
    CHECK(core::detail::asyncBufferLayout(
            1,
            1,
            static_cast<core::CompressionCodec>(255),
            layout)
          == Status::Invalid_parameters);
    REQUIRE(core::detail::asyncBufferLayout(
              2, 8, core::CompressionCodec::none, layout)
            == Status::Success);
    CHECK(layout.values == 10);
    CHECK(layout.state_values == 8);
    CHECK(layout.raw_bytes == 10 * sizeof(double));

    layout = { .state_values = 7,
               .values = 8,
               .raw_bytes = 9,
               .compressed_bound = 10 };
    CHECK(core::detail::compressedRecordingBufferLayout(
            std::numeric_limits<std::size_t>::max(),
            1,
            2,
            core::CompressionCodec::none,
            layout)
          == Status::Invalid_parameters);
    CHECK(layout.state_values == 7);

    core::detail::CompressedRecordingStorageLayout storage{
      .state_values = 11,
      .current_values = 12,
    };
    CHECK(core::detail::compressedRecordingStorageLayout(
            std::numeric_limits<std::size_t>::max(),
            2,
            1,
            24,
            24,
            core::CompressionCodec::none,
            0,
            storage)
          == Status::Invalid_parameters);
    CHECK(storage.state_values == 11);
    REQUIRE(core::detail::compressedRecordingStorageLayout(
              3, 8, 2, 80, 80, core::CompressionCodec::none, 0, storage)
            == Status::Success);
    CHECK(storage.state_values == 24);
    CHECK(storage.current_values == 6);
    CHECK(storage.payload_bytes == 80);
    CHECK(storage.workspace_bytes == 80);
    CHECK(storage.resident_bytes == 528);

    REQUIRE(core::detail::compressedRecordingStorageLayout(
              0,
              12'000'000,
              1,
              96'000'008,
              96'000'008,
              core::CompressionCodec::none,
              0,
              storage)
            == Status::Success);
    CHECK(storage.state_values == 0);
    CHECK(storage.current_values == 0);
    CHECK(storage.payload_bytes == 0);
    CHECK(storage.workspace_bytes == 0);
    CHECK(storage.resident_bytes == 0);

    CHECK(core::detail::compressedRecordingStorageLayout(
            1,
            12'000'000,
            1,
            96'000'008,
            96'000'008,
            core::CompressionCodec::none,
            0,
            storage)
          == Status::Invalid_parameters);
    CHECK(core::detail::max_eager_recording_bytes == 256U * 1024U * 1024U);

    std::size_t retained = 123;
    CHECK(core::detail::compressedRecordingRetainedBytes(
            std::numeric_limits<std::size_t>::max(), 1, 1, 1, retained)
          == Status::Invalid_parameters);
    CHECK(retained == 123);
    REQUIRE(core::detail::compressedRecordingRetainedBytes(1, 2, 3, 4, retained)
            == Status::Success);
    CHECK(retained == 80);
    CHECK(core::detail::compressedRecordingStorageLayout(
            1,
            1,
            1,
            16,
            16,
            core::CompressionCodec::none,
            core::detail::max_eager_recording_bytes - 79,
            storage)
          == Status::Invalid_parameters);

    core::detail::AsyncBufferLayout payload_layout{
      .raw_bytes = 80,
      .compressed_bound = 96,
    };
    CHECK(core::detail::validateCompressedPayloadSize(
            core::CompressionCodec::none, 79, payload_layout)
          == Status::Invalid_parameters);
    CHECK(core::detail::validateCompressedPayloadSize(
            core::CompressionCodec::none, 80, payload_layout)
          == Status::Success);
    CHECK(core::detail::validateCompressedPayloadSize(
            core::CompressionCodec::zstd, 97, payload_layout)
          == Status::Invalid_parameters);
    CHECK(core::detail::validateCompressedPayloadSize(
            core::CompressionCodec::zstd, 96, payload_layout)
          == Status::Success);
  }

  SECTION("unconfigured and finished producers reject snapshots")
  {
    const auto path = temporary("p9_snapshot_guards.slcmp");
    std::filesystem::remove(path, ignored);
    core::AsyncRecorder recorder;
    const std::array current{ 1.0, -1.0 };
    CHECK(recorder.enqueue(0, current) == Status::Invalid_parameters);

    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    REQUIRE(recorder.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    CHECK(recorder.enqueueSnapshot(
            0, 0.0, std::span<const double>{}, batch.state().raw())
          == Status::Invalid_parameters);
    REQUIRE(recorder.enqueue(1, current) == Status::Success);
    CHECK(recorder.enqueue(1, current) == Status::Invalid_parameters);
    REQUIRE(recorder.finish() == Status::Success);
    CHECK(recorder.enqueue(2, current) == Status::Invalid_parameters);
    std::filesystem::remove(path, ignored);
  }

  SECTION("a directory cannot be configured as an output file")
  {
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder recorder;
    CHECK(recorder.configure(
            batch,
            std::filesystem::temp_directory_path(),
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Invalid_parameters);
    CHECK_FALSE(recorder.configured());
  }

  SECTION("wide codec values cannot alias a supported byte-sized enum")
  {
    const auto path = temporary("p9_codec.slcmp");
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder writer;
    REQUIRE(writer.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    const std::array current{ 1.0, -1.0 };
    REQUIRE(writer.enqueue(0, current) == Status::Success);
    REQUIRE(writer.finish() == Status::Success);
    auto bytes = readBytes(path);
    REQUIRE(bytes.size() >= 128);
    writeScalar<std::uint32_t>(bytes, 36, 256U);
    writeScalar<std::uint32_t>(bytes, 64 + 8, 256U);
    sealHeader(bytes, 64, 56);
    sealHeader(bytes, 0, 20);
    writeBytes(path, bytes);

    core::CompressedRecording recording;
    CHECK(recording.open(path) == Status::Invalid_parameters);
    CHECK_FALSE(recording.valid());
    std::filesystem::remove(path, ignored);
  }

  SECTION("invalid policy and overflowing derived density are rejected")
  {
    const auto invalid_path = temporary("p9_policy.slcmp");
    const auto density_path = temporary("p9_density.slcmp");
    std::filesystem::remove(invalid_path, ignored);
    std::filesystem::remove(density_path, ignored);

    const auto invalid_input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      invalid_input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder invalid;
    CHECK(invalid.configure(
            batch,
            invalid_path,
            { .ring_slots = 3,
              .backpressure = static_cast<core::AsyncBackpressurePolicy>(255),
              .codec = core::CompressionCodec::none })
          == Status::Invalid_parameters);
    CHECK_FALSE(invalid.configured());

    auto density_input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    density_input.design.electrode_area = 1e-300;
    auto tiny_area_batch = test_support::requireSpmBatch(
      density_input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder density;
    REQUIRE(density.configure(
              tiny_area_batch,
              density_path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    const std::array current{ 1e300, -1e300 };
    CHECK(density.enqueue(0, current) == Status::Invalid_parameters);
    CHECK(density.finish() == Status::Success);
    CHECK(density.snapshotsWritten() == 0);
    std::filesystem::remove(invalid_path, ignored);
    std::filesystem::remove(density_path, ignored);
  }
}

TEST_CASE("Async recorder translates deterministic worker and stream faults",
          "[core][async-recorder][fault-injection][coverage]")
{
  std::error_code ignored;
  const std::array current{ 1.0, -1.0 };

  SECTION("placeholder write failure is atomic")
  {
    const auto path = temporary("fault_placeholder.slcmp");
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder recorder;
    core::detail::AsyncRecorderTestAccess::failPlaceholderWrite(recorder);
    CHECK(recorder.configure(
            batch,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Numerical_failure);
    CHECK_FALSE(recorder.configured());
    std::filesystem::remove(path, ignored);
  }

  SECTION("thread construction failure rolls configuration back")
  {
    const auto path = temporary("fault_thread.slcmp");
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder recorder;
    core::detail::AsyncRecorderTestAccess::failThreadCreation(recorder);
    CHECK(recorder.configure(
            batch,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Numerical_failure);
    CHECK_FALSE(recorder.configured());

    core::detail::AsyncRecorderTestAccess::restoreThreadCreation(recorder);
    REQUIRE(recorder.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    REQUIRE(recorder.enqueue(0, current) == Status::Success);
    REQUIRE(recorder.finish() == Status::Success);
    core::CompressedRecording recording;
    REQUIRE(recording.open(path) == Status::Success);
    REQUIRE(recording.size() == 1);
    CHECK(recording.snapshot(0).accepted_step == 0);
    std::filesystem::remove(path, ignored);
  }

  SECTION("corrupt private shuffle storage is detected by the worker")
  {
    const auto path = temporary("fault_shuffle.slcmp");
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder recorder;
    REQUIRE(recorder.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    core::detail::AsyncRecorderTestAccess::truncateShuffleBuffer(recorder);
    REQUIRE(recorder.enqueue(0, current) == Status::Success);
    CHECK(recorder.finish() == Status::Invalid_states);
    CHECK(recorder.snapshotsWritten() == 0);
    std::filesystem::remove(path, ignored);
  }

  SECTION("block sink failure propagates from the drain thread")
  {
    const auto path = temporary("fault_block.slcmp");
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder recorder;
    REQUIRE(recorder.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    core::detail::AsyncRecorderTestAccess::failBlockWrite(recorder);
    REQUIRE(recorder.enqueue(0, current) == Status::Success);
    CHECK(recorder.finish() == Status::Numerical_failure);
    CHECK(recorder.snapshotsWritten() == 0);
    std::filesystem::remove(path, ignored);
  }

  SECTION("final tell failure propagates")
  {
    const auto path = temporary("fault_finalize_tell.slcmp");
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder recorder;
    core::detail::AsyncRecorderTestAccess::failFinalizeTell(recorder);
    REQUIRE(recorder.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    CHECK(recorder.finish() == Status::Numerical_failure);
    std::filesystem::remove(path, ignored);
  }

  SECTION("final header write failure propagates")
  {
    const auto path = temporary("fault_finalize_write.slcmp");
    std::filesystem::remove(path, ignored);
    const auto input =
      test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
    auto batch = test_support::requireSpmBatch(
      input, core::SpmModelOptions{}, 2);
    core::AsyncRecorder recorder;
    core::detail::AsyncRecorderTestAccess::failFinalizeWrite(recorder);
    REQUIRE(recorder.configure(
              batch,
              path,
              { .ring_slots = 3,
                .backpressure = core::AsyncBackpressurePolicy::block,
                .codec = core::CompressionCodec::none })
            == Status::Success);
    CHECK(recorder.finish() == Status::Numerical_failure);
    std::filesystem::remove(path, ignored);
  }

  SECTION("short exact read is rejected")
  {
    std::istringstream input{ "x" };
    std::array<std::byte, 2> destination{};
    std::uint64_t remaining = destination.size();
    CHECK(core::detail::AsyncRecorderTestAccess::readExact(
            input, destination, remaining)
          == Status::Invalid_parameters);
    CHECK(remaining == destination.size());

    std::istringstream complete{ "xx" };
    remaining = 1;
    CHECK(core::detail::AsyncRecorderTestAccess::readExact(
            complete, destination, remaining)
          == Status::Invalid_parameters);
    CHECK(remaining == 1);
  }
}

TEST_CASE("P8-G3 byte shuffle is a bitwise involution",
          "[core][async-recorder][P8-G3]")
{
  const std::array<double, 7> values{
    0.0, -0.0, 1.0, -2.5, 1e-300, 1e300, 3.141592653589793
  };
  const auto bytes = std::as_bytes(std::span{ values });
  std::vector<std::byte> shuffled(bytes.size()), restored(bytes.size());
  REQUIRE(core::byteShuffle(bytes, shuffled, sizeof(double))
          == Status::Success);
  REQUIRE(core::byteUnshuffle(shuffled, restored, sizeof(double))
          == Status::Success);
  CHECK(std::equal(bytes.begin(), bytes.end(), restored.begin(), restored.end()));
}

TEST_CASE("byte shuffle rejects degenerate public spans",
          "[core][async-recorder][shuffle][P9]")
{
  std::array<std::byte, 8> storage{};
  CHECK(checkedShuffleCall([&] {
          return core::byteShuffle({}, {}, 0);
        })
        == Status::Invalid_parameters);
  CHECK(checkedShuffleCall([&] {
          return core::byteShuffle(
            std::span<const std::byte>{ storage },
            std::span<std::byte>{ storage }.first(7),
            sizeof(double));
        })
        == Status::Invalid_parameters);
  CHECK(checkedShuffleCall([&] {
          return core::byteUnshuffle(
            std::span<const std::byte>{ storage },
            std::span<std::byte>{ storage },
            sizeof(double));
        })
        == Status::Invalid_parameters);
}

TEST_CASE("P8-G3 async blocks round-trip accepted snapshots bitwise",
          "[core][async-recorder][P8-G3]")
{
  const auto path = temporary("roundtrip.slcmp");
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::Recorder reference;
  REQUIRE(reference.configure(batch, { .capacity = 8 }) == Status::Success);
  core::AsyncRecorder async;
  REQUIRE(async.configure(batch,
                          path,
                          { .ring_slots = 3,
                            .backpressure = core::AsyncBackpressurePolicy::block,
                            .codec = testCodec() })
          == Status::Success);
  core::EulerLegacy stepper{ batch };
  const std::array current{ 8.0, -4.0 };
  const std::array density{ current[0] / batch.electrode_area(),
                            current[1] / batch.electrode_area() };
  for (std::uint64_t step = 0; step < 8; ++step) {
    REQUIRE(reference.record(step, current) == Status::Success);
    REQUIRE(async.enqueue(step, current) == Status::Success);
    if (step + 1 < 8)
      REQUIRE(stepper.step(batch, density, static_cast<double>(step), 1.0)
              == Status::Success);
  }
  REQUIRE(async.finish() == Status::Success);
  CHECK(async.snapshotsWritten() == 8);
  CHECK(async.thinnedSnapshots() == 0);

  core::CompressedRecording decoded;
  REQUIRE(decoded.open(path) == Status::Success);
  REQUIRE(decoded.size() == reference.size());
  CHECK(decoded.nRows() == reference.nRows());
  CHECK(decoded.nLanes() == reference.nLanes());
  CHECK(decoded.stride() == reference.stride());
  for (std::size_t index = 0; index < reference.size(); ++index) {
    const auto expected = reference.snapshot(index);
    const auto actual = decoded.snapshot(index);
    CHECK(actual.accepted_step == expected.accepted_step);
    CHECK(actual.time == expected.time);
    CHECK(std::equal(actual.current_density.begin(),
                     actual.current_density.end(),
                     expected.current_density.begin(),
                     expected.current_density.end()));
    CHECK(std::equal(actual.state.begin(),
                     actual.state.end(),
                     expected.state.begin(),
                     expected.state.end()));
  }
  CHECK(core::detail::AsyncRecorderTestAccess::retainedBytes(decoded) > 0);
  decoded.close();
  CHECK(core::detail::AsyncRecorderTestAccess::retainedBytes(decoded) == 0);
  REQUIRE(decoded.open(path) == Status::Success);
  CHECK(decoded.size() == reference.size());
  std::filesystem::remove(path, ignored);
}

TEST_CASE("empty async recording opens without decode workspace",
          "[core][async-recorder][empty][coverage]")
{
  const auto path = temporary("empty.slcmp");
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::AsyncRecorder writer;
  REQUIRE(writer.configure(
            batch,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Success);
  REQUIRE(writer.finish() == Status::Success);
  CHECK(std::filesystem::file_size(path) == 64);

  core::CompressedRecording recording;
  REQUIRE(recording.open(path) == Status::Success);
  CHECK(recording.valid());
  CHECK(recording.size() == 0);
  CHECK(core::detail::AsyncRecorderTestAccess::retainedBytes(recording) == 0);
  std::filesystem::remove(path, ignored);
}

TEST_CASE("P8-G3 thin never waits silently and block stays lossless",
          "[core][async-recorder][backpressure][P8-G3]")
{
  using namespace std::chrono_literals;
  const auto thin_path = temporary("thin.slcmp");
  const auto block_path = temporary("block.slcmp");
  std::error_code ignored;
  std::filesystem::remove(thin_path, ignored);
  std::filesystem::remove(block_path, ignored);
  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  const std::array current{ 1.0, 1.0 };

  core::AsyncRecorder thin;
  thin.setDrainHook([] { std::this_thread::sleep_for(20ms); });
  REQUIRE(thin.configure(batch,
                         thin_path,
                         { .ring_slots = 3,
                           .backpressure = core::AsyncBackpressurePolicy::thin,
                           .codec = testCodec() })
          == Status::Success);
  const auto start = std::chrono::steady_clock::now();
  for (std::uint64_t step = 0; step < 50; ++step)
    REQUIRE(thin.enqueue(step, current) == Status::Success);
  const auto producer_time = std::chrono::steady_clock::now() - start;
  REQUIRE(thin.finish() == Status::Success);
  CHECK(producer_time < 100ms);
  CHECK(thin.thinnedSnapshots() > 0);
  CHECK(thin.snapshotsWritten() + thin.thinnedSnapshots() == 50);

  core::AsyncRecorder block;
  block.setDrainHook([] { std::this_thread::sleep_for(2ms); });
  REQUIRE(block.configure(batch,
                          block_path,
                          { .ring_slots = 3,
                            .backpressure = core::AsyncBackpressurePolicy::block,
                            .codec = testCodec() })
          == Status::Success);
  for (std::uint64_t step = 0; step < 8; ++step)
    REQUIRE(block.enqueue(step, current) == Status::Success);
  REQUIRE(block.finish() == Status::Success);
  CHECK(block.snapshotsWritten() == 8);
  CHECK(block.thinnedSnapshots() == 0);

  core::CompressedRecording decoded;
  REQUIRE(decoded.open(block_path) == Status::Success);
  REQUIRE(decoded.size() == 8);
  for (std::size_t i = 0; i < decoded.size(); ++i)
    CHECK(decoded.snapshot(i).accepted_step == i);
  std::filesystem::remove(thin_path, ignored);
  std::filesystem::remove(block_path, ignored);
}

TEST_CASE("closing wakes a producer blocked behind a full async ring",
          "[core][async-recorder][backpressure][coverage]")
{
  const auto path = temporary("closing_wakeup.slcmp");
  std::error_code ignored;
  std::filesystem::remove(path, ignored);
  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);

  std::promise<void> drain_entered;
  auto entered = drain_entered.get_future();
  std::promise<void> release_drain;
  auto release = release_drain.get_future().share();
  std::once_flag first_drain;
  core::AsyncRecorder recorder;
  struct DrainRelease
  {
    std::promise<void> &promise;
    bool released{};

    void release() noexcept
    {
      if (!released) {
        promise.set_value();
        released = true;
      }
    }

    ~DrainRelease() { release(); }
  } unblock{ release_drain };
  recorder.setDrainHook([&] {
    std::call_once(first_drain, [&] {
      drain_entered.set_value();
      release.wait();
    });
  });
  REQUIRE(recorder.configure(
            batch,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Success);
  const std::array current{ 1.0, -1.0 };
  REQUIRE(recorder.enqueue(0, current) == Status::Success);
  REQUIRE(entered.wait_for(std::chrono::seconds{ 2 })
          == std::future_status::ready);
  REQUIRE(recorder.enqueue(1, current) == Status::Success);
  REQUIRE(recorder.enqueue(2, current) == Status::Success);

  Status producer_status{};
  std::promise<void> producer_done;
  auto producer_finished = producer_done.get_future();
  std::thread producer([&] {
    producer_status = recorder.enqueue(3, current);
    producer_done.set_value();
  });
  Status finish_status{};
  std::thread finisher([&] { finish_status = recorder.finish(); });
  const bool producer_woke =
    producer_finished.wait_for(std::chrono::seconds{ 2 })
    == std::future_status::ready;
  unblock.release();
  producer.join();
  finisher.join();
  CHECK(producer_woke);
  CHECK(producer_status == Status::Invalid_states);
  CHECK(finish_status == Status::Success);
  std::filesystem::remove(path, ignored);
}

TEST_CASE("P8-G3 compressed reader atomically rejects corruption and truncation",
          "[core][async-recorder][hardened][P8-G3]")
{
  const auto valid = temporary("valid.slcmp");
  const auto header_corrupt = temporary("header.slcmp");
  const auto payload_corrupt = temporary("payload.slcmp");
  const auto truncated = temporary("truncated.slcmp");
  std::error_code ignored;
  for (const auto &path : { valid, header_corrupt, payload_corrupt, truncated })
    std::filesystem::remove(path, ignored);
  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::AsyncRecorder recorder;
  REQUIRE(recorder.configure(batch,
                             valid,
                             { .ring_slots = 3,
                               .backpressure = core::AsyncBackpressurePolicy::block,
                               .codec = testCodec() })
          == Status::Success);
  const std::array current{ 1.0, 1.0 };
  REQUIRE(recorder.enqueue(0, current) == Status::Success);
  REQUIRE(recorder.finish() == Status::Success);
  REQUIRE(std::filesystem::copy_file(valid, header_corrupt));
  REQUIRE(std::filesystem::copy_file(valid, payload_corrupt));
  REQUIRE(std::filesystem::copy_file(valid, truncated));
  flipByte(header_corrupt, 24);
  flipByte(payload_corrupt, 128);
  const auto size = std::filesystem::file_size(truncated);
  REQUIRE(size > 128);
  std::filesystem::resize_file(truncated, size - 1);

  core::CompressedRecording recording;
  REQUIRE(recording.open(valid) == Status::Success);
  CHECK(recording.open(header_corrupt) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(payload_corrupt) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(truncated) == Status::Invalid_parameters);
  CHECK(recording.valid());
  for (const auto &path : { valid, header_corrupt, payload_corrupt, truncated })
    std::filesystem::remove(path, ignored);
}

TEST_CASE("compressed reader rejects exact block and layout corruptions",
          "[core][async-recorder][hardened][coverage]")
{
  const auto valid = temporary("coverage_valid.slcmp");
  const auto invalid_block = temporary("coverage_block.slcmp");
  const auto short_payload = temporary("coverage_payload_size.slcmp");
  const auto raw_corrupt = temporary("coverage_raw_crc.slcmp");
  const auto trailing = temporary("coverage_trailing.slcmp");
  const auto overflow = temporary("coverage_overflow.slcmp");
  const auto budget = temporary("coverage_budget.slcmp");
  const auto absent = temporary("coverage_absent.slcmp");
  const auto short_file = temporary("coverage_short.slcmp");
  std::error_code ignored;
  for (const auto &path : { valid,
                            invalid_block,
                            short_payload,
                            raw_corrupt,
                            trailing,
                            overflow,
                            budget,
                            absent,
                            short_file })
    std::filesystem::remove(path, ignored);

  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  auto batch = test_support::requireSpmBatch(
    input, core::SpmModelOptions{}, 2);
  core::AsyncRecorder recorder;
  REQUIRE(recorder.configure(
            batch,
            valid,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Success);
  const std::array current{ 1.0, -1.0 };
  REQUIRE(recorder.enqueue(0, current) == Status::Success);
  REQUIRE(recorder.finish() == Status::Success);
  const auto original = readBytes(valid);
  REQUIRE(original.size() > 128);

  auto block_bytes = original;
  writeScalar<double>(block_bytes, 64 + 24, std::bit_cast<double>(UINT64_C(0x7ff8000000000000)));
  sealHeader(block_bytes, 64, 56);
  writeBytes(invalid_block, block_bytes);

  auto short_payload_bytes = original;
  const auto raw_bytes = readScalar<std::uint64_t>(short_payload_bytes, 64 + 32);
  REQUIRE(raw_bytes > sizeof(double));
  const auto payload_bytes = raw_bytes - sizeof(double);
  writeScalar<std::uint64_t>(short_payload_bytes, 64 + 40, payload_bytes);
  const auto payload = std::span<const std::byte>{ short_payload_bytes }.subspan(
    128, static_cast<std::size_t>(payload_bytes));
  writeScalar<std::uint32_t>(short_payload_bytes, 64 + 52, testCrc32(payload));
  sealHeader(short_payload_bytes, 64, 56);
  writeBytes(short_payload, short_payload_bytes);

  auto raw_corrupt_bytes = original;
  raw_corrupt_bytes[128] ^= std::byte{ 0x5a };
  const auto complete_payload = std::span<const std::byte>{ raw_corrupt_bytes }.subspan(
    128, static_cast<std::size_t>(raw_bytes));
  writeScalar<std::uint32_t>(raw_corrupt_bytes,
                             64 + 52,
                             testCrc32(complete_payload));
  sealHeader(raw_corrupt_bytes, 64, 56);
  writeBytes(raw_corrupt, raw_corrupt_bytes);

  auto trailing_bytes = original;
  trailing_bytes.push_back(std::byte{});
  writeScalar<std::uint64_t>(trailing_bytes, 48, trailing_bytes.size());
  sealHeader(trailing_bytes, 0, 20);
  writeBytes(trailing, trailing_bytes);

  std::array<std::byte, 64> overflow_header{};
  constexpr std::array magic{ 'S', 'L', 'I', 'D', 'E', 'C', 'M', 'P' };
  std::memcpy(overflow_header.data(), magic.data(), magic.size());
  writeScalar<std::uint16_t>(overflow_header, 8, 1U);
  writeScalar<std::uint16_t>(overflow_header, 10, 0U);
  writeScalar<std::uint32_t>(overflow_header, 12, 0x01020304U);
  writeScalar<std::uint32_t>(overflow_header, 16, 64U);
  constexpr auto maximum_dimension =
    static_cast<std::uint32_t>(std::numeric_limits<int>::max());
  writeScalar<std::uint32_t>(overflow_header, 24, maximum_dimension);
  writeScalar<std::uint32_t>(overflow_header, 28, maximum_dimension);
  writeScalar<std::uint32_t>(overflow_header, 32, maximum_dimension);
  writeScalar<std::uint32_t>(overflow_header, 36, 0U);
  writeScalar<std::uint64_t>(overflow_header, 40, 0U);
  writeScalar<std::uint64_t>(overflow_header, 48, 64U);
  sealHeader(overflow_header, 0, 20);
  writeBytes(overflow, overflow_header);

  std::array<std::byte, 128> budget_bytes{};
  std::memcpy(budget_bytes.data(), magic.data(), magic.size());
  writeScalar<std::uint16_t>(budget_bytes, 8, 1U);
  writeScalar<std::uint16_t>(budget_bytes, 10, 0U);
  writeScalar<std::uint32_t>(budget_bytes, 12, 0x01020304U);
  writeScalar<std::uint32_t>(budget_bytes, 16, 64U);
  writeScalar<std::uint32_t>(budget_bytes, 24, 1U);
  writeScalar<std::uint32_t>(budget_bytes, 28, 1U);
  writeScalar<std::uint32_t>(budget_bytes, 32, 12'000'000U);
  writeScalar<std::uint32_t>(budget_bytes, 36, 0U);
  writeScalar<std::uint64_t>(budget_bytes, 40, 1U);
  writeScalar<std::uint64_t>(budget_bytes, 48, budget_bytes.size());
  sealHeader(budget_bytes, 0, 20);
  writeBytes(budget, budget_bytes);

  {
    std::ofstream output(short_file, std::ios::binary | std::ios::trunc);
    REQUIRE(output.good());
    output.put('x');
  }

  core::CompressedRecording recording;
  REQUIRE(recording.open(valid) == Status::Success);
  CHECK(recording.open(absent) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(short_file) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(overflow) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(budget) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(invalid_block) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(short_payload) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(raw_corrupt) == Status::Invalid_parameters);
  CHECK(recording.valid());
  CHECK(recording.open(trailing) == Status::Invalid_parameters);
  CHECK(recording.valid());

  for (const auto &path : { valid,
                            invalid_block,
                            short_payload,
                            raw_corrupt,
                            trailing,
                            overflow,
                            budget,
                            short_file })
    std::filesystem::remove(path, ignored);
}
