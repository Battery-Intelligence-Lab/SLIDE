/**
 * @file core_AsyncRecorder_test.cpp
 * @brief P8-G3 async ring, backpressure, compression, and corruption gates.
 */

#include "../../src/core/AsyncRecorder.hpp"
#include "../../src/core/EulerLegacy.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <span>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

using namespace slide;

namespace {

core::SpmBatch makeBatch(double electrode_area = -1.0)
{
  core::SpmBatch batch;
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  if (electrode_area > 0.0)
    input.design.electrode_area = electrode_area;
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);
  return batch;
}

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

  SECTION("wide codec values cannot alias a supported byte-sized enum")
  {
    const auto path = temporary("p9_codec.slcmp");
    std::filesystem::remove(path, ignored);
    auto batch = makeBatch();
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

    auto batch = makeBatch();
    core::AsyncRecorder invalid;
    CHECK(invalid.configure(
            batch,
            invalid_path,
            { .ring_slots = 3,
              .backpressure = static_cast<core::AsyncBackpressurePolicy>(255),
              .codec = core::CompressionCodec::none })
          == Status::Invalid_parameters);
    CHECK_FALSE(invalid.configured());

    auto tiny_area_batch = makeBatch(1e-300);
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
  auto batch = makeBatch();
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
  auto batch = makeBatch();
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
  auto batch = makeBatch();
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
