/**
 * @file core_RecordingFormatCommon_test.cpp
 * @brief Independent byte and eager-index oracles for recording-common ownership.
 */

#include "../../src/core/AsyncRecorder.hpp"
#include "../../src/core/Recorder.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using namespace slide;

namespace {

constexpr std::size_t header_bytes = 64;
constexpr std::uint32_t expected_endian_marker = 0x01020304U;

struct ExpectedSnapshot
{
  std::uint64_t accepted_step{};
  core::real_t time{};
  std::array<core::real_t, 2> current_density{};
  std::vector<core::real_t> state{};
};

struct ByteFingerprint
{
  std::size_t bytes{};
  std::uint64_t fnv1a{};
  std::uint64_t mixed{};
};

struct BlockCrcPair
{
  std::uint32_t raw{};
  std::uint32_t payload{};
};

struct RecordingOracle
{
  ByteFingerprint csv{};
  ByteFingerprint binary{};
  ByteFingerprint async{};
  BlockCrcPair first_block{};
  BlockCrcPair second_block{};
};

#if defined(SLIDE_TEST_IPO)
constexpr RecordingOracle expected_oracle{
  .csv = { 4072, UINT64_C(113020753288656565), UINT64_C(1220644351662519581) },
  .binary = { 3904, UINT64_C(4980822913020495434), UINT64_C(5810136181434859582) },
  .async = { 3936, UINT64_C(15459230203036584454), UINT64_C(6543818133775003568) },
  .first_block = { 1528093825U, 3638640079U },
  .second_block = { 936358080U, 170786859U },
};
#elif defined(SLIDE_TEST_RELEASE)
constexpr RecordingOracle expected_oracle{
  .csv = { 4072, UINT64_C(9836986554487851297), UINT64_C(10195804735363157979) },
  .binary = { 3904, UINT64_C(7644568223896275882), UINT64_C(9507010167800238048) },
  .async = { 3936, UINT64_C(12852073895150042970), UINT64_C(2881486469091071908) },
  .first_block = { 1129188246U, 4052193948U },
  .second_block = { 798320599U, 591966072U },
};
#else
constexpr RecordingOracle expected_oracle{
  .csv = { 4072, UINT64_C(9687440020754251757), UINT64_C(4395482851448939133) },
  .binary = { 3904, UINT64_C(5616759404010358714), UINT64_C(2046236715092195141) },
  .async = { 3936, UINT64_C(674564965144750450), UINT64_C(6807551341330539846) },
  .first_block = { 4136208346U, 118461792U },
  .second_block = { 2589125531U, 3586173060U },
};
#endif

static_assert(expected_oracle.first_block.raw
              != expected_oracle.first_block.payload);
static_assert(expected_oracle.second_block.raw
              != expected_oracle.second_block.payload);

std::filesystem::path temporary(std::string_view suffix)
{
  return std::filesystem::temp_directory_path()
         / ("slide_recording_common_" + std::string{ suffix });
}

core::SpmBatch makeBatch()
{
  const auto input =
    test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  core::SpmBatch batch;
  if (core::buildSpmBatch(input, core::SpmModelOptions{}, 2, batch)
      != Status::Success)
    throw std::runtime_error{ "deterministic recording fixture failed to build" };
  return batch;
}

ExpectedSnapshot prepareSnapshot(core::SpmBatch &batch, std::size_t index)
{
  auto raw = batch.state().raw();
  const auto base = static_cast<core::real_t>((index + 1) * 1000);
  for (int row = 0; row < batch.state().n_rows(); ++row)
    for (int lane = batch.state().n_lanes();
         lane < batch.state().stride();
         ++lane) {
      const auto value =
        static_cast<std::size_t>(row)
          * static_cast<std::size_t>(batch.state().stride())
        + static_cast<std::size_t>(lane);
      raw[value] = base + static_cast<core::real_t>(value);
    }

  const core::real_t time = 10.25 + static_cast<core::real_t>(index);
  batch.state().at(batch.layout().elapsed_time, 0, 0) = time;
  const std::array densities =
    index == 0 ? std::array<core::real_t, 2>{ 1.0, -2.0 }
               : std::array<core::real_t, 2>{ 3.0, -4.0 };
  return { .accepted_step = index == 0 ? 7U : 11U,
           .time = time,
           .current_density = densities,
           .state = { raw.begin(), raw.end() } };
}

std::array<core::real_t, 2> totalCurrent(
  const core::SpmBatch &batch,
  const ExpectedSnapshot &snapshot)
{
  return { batch.electrode_area() * snapshot.current_density[0],
           batch.electrode_area() * snapshot.current_density[1] };
}

std::string expectedCsvHeader(int rows, int lanes)
{
  std::ostringstream output;
  output << "accepted_step,time_s";
  for (int lane = 0; lane < lanes; ++lane)
    output << ",current_density_lane" << lane << "_A_m2";
  for (int row = 0; row < rows; ++row)
    for (int lane = 0; lane < lanes; ++lane)
      output << ",state_r" << row << "_lane" << lane;
  for (int lane = 0; lane < lanes; ++lane)
    output << ",terminal_voltage_lane" << lane << "_V";
  output << '\n';
  return output.str();
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
             static_cast<std::streamsize>(bytes.size()));
  REQUIRE(input.good());
  return bytes;
}

ByteFingerprint fingerprint(std::span<const std::byte> bytes)
{
  std::uint64_t fnv1a = UINT64_C(14695981039346656037);
  std::uint64_t mixed = UINT64_C(0x6a09e667f3bcc909);
  for (const auto byte : bytes) {
    const auto value =
      static_cast<std::uint64_t>(std::to_integer<std::uint8_t>(byte));
    fnv1a = (fnv1a ^ value) * UINT64_C(1099511628211);
    mixed = std::rotl(
      mixed ^ (value + UINT64_C(0x9e3779b97f4a7c15)), 17);
    mixed *= UINT64_C(0xbf58476d1ce4e5b9);
  }
  return { .bytes = bytes.size(), .fnv1a = fnv1a, .mixed = mixed };
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
std::optional<T> scalarAt(std::span<const std::byte> bytes,
                          std::size_t offset)
{
  if (offset > bytes.size() || sizeof(T) > bytes.size() - offset)
    return std::nullopt;
  T value{};
  std::memcpy(&value, bytes.data() + offset, sizeof(value));
  return value;
}

bool headerCrcMatches(std::span<const std::byte> bytes,
                      std::size_t header_offset,
                      std::size_t crc_offset)
{
  if (header_offset > bytes.size()
      || header_bytes > bytes.size() - header_offset
      || crc_offset > header_bytes
      || sizeof(std::uint32_t) > header_bytes - crc_offset)
    return false;
  const auto stored =
    scalarAt<std::uint32_t>(bytes, header_offset + crc_offset);
  if (!stored.has_value() || *stored == 0)
    return false;
  std::array<std::byte, header_bytes> header{};
  std::ranges::copy(
    bytes.subspan(header_offset, header_bytes), header.begin());
  std::ranges::fill(
    std::span{ header }.subspan(crc_offset, sizeof(std::uint32_t)),
    std::byte{});
  return testCrc32(header) == *stored;
}

bool asyncBlockMatches(std::span<const std::byte> bytes,
                       std::size_t block_offset,
                       const ExpectedSnapshot &expected,
                       BlockCrcPair &observed)
{
  std::vector<core::real_t> raw_values;
  raw_values.reserve(
    expected.current_density.size() + expected.state.size());
  raw_values.insert(raw_values.end(),
                    expected.current_density.begin(),
                    expected.current_density.end());
  raw_values.insert(
    raw_values.end(), expected.state.begin(), expected.state.end());
  const auto raw = std::as_bytes(std::span{ raw_values });
  const auto stored_raw_bytes =
    scalarAt<std::uint64_t>(bytes, block_offset + 32);
  const auto stored_payload_bytes =
    scalarAt<std::uint64_t>(bytes, block_offset + 40);
  const auto stored_raw_crc =
    scalarAt<std::uint32_t>(bytes, block_offset + 48);
  const auto stored_payload_crc =
    scalarAt<std::uint32_t>(bytes, block_offset + 52);
  const std::size_t payload_offset = block_offset + header_bytes;
  if (!stored_raw_bytes.has_value() || !stored_payload_bytes.has_value()
      || !stored_raw_crc.has_value() || !stored_payload_crc.has_value()
      || *stored_raw_bytes != raw.size()
      || *stored_payload_bytes != raw.size()
      || payload_offset > bytes.size()
      || raw.size() > bytes.size() - payload_offset)
    return false;
  observed = { .raw = *stored_raw_crc, .payload = *stored_payload_crc };
  const auto payload = bytes.subspan(payload_offset, raw.size());
  return headerCrcMatches(bytes, block_offset, 56)
         && *stored_raw_crc == testCrc32(raw)
         && *stored_payload_crc == testCrc32(payload)
         && *stored_raw_crc != *stored_payload_crc;
}

void checkSnapshot(const core::SnapshotView &actual,
                   const ExpectedSnapshot &expected)
{
  CHECK(actual.accepted_step == expected.accepted_step);
  CHECK(actual.time == expected.time);
  CHECK(std::ranges::equal(
    actual.current_density, expected.current_density));
  CHECK(std::ranges::equal(actual.state, expected.state));
}

} // namespace

TEST_CASE("Recording common format preserves synchronous bytes and eager indexing",
          "[core][recorder][recording-common][MQ2-R1]")
{
  const auto csv_path = temporary("oracle.csv");
  const auto binary_path = temporary("oracle.slrec");
  std::error_code ignored;
  std::filesystem::remove(csv_path, ignored);
  std::filesystem::remove(binary_path, ignored);

  auto batch = makeBatch();
  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .capacity = 2 })
          == Status::Success);

  std::array<ExpectedSnapshot, 2> expected;
  for (std::size_t index = 0; index < expected.size(); ++index) {
    expected[index] = prepareSnapshot(batch, index);
    const auto current = totalCurrent(batch, expected[index]);
    REQUIRE(recorder.record(expected[index].accepted_step, current)
            == Status::Success);
  }
  for (std::size_t index = 0; index < expected.size(); ++index)
    checkSnapshot(recorder.snapshot(index), expected[index]);

  REQUIRE(recorder.writeCsv(csv_path) == Status::Success);
  const auto csv_bytes = readBytes(csv_path);
  const std::string csv_text{
    reinterpret_cast<const char *>(csv_bytes.data()), csv_bytes.size()
  };
  const auto csv_fingerprint = fingerprint(csv_bytes);
  CAPTURE(csv_fingerprint.bytes,
          csv_fingerprint.fnv1a,
          csv_fingerprint.mixed);
  CHECK((csv_text.starts_with(
           expectedCsvHeader(recorder.n_rows(), recorder.n_lanes()))
         && csv_fingerprint.bytes == expected_oracle.csv.bytes
         && csv_fingerprint.fnv1a == expected_oracle.csv.fnv1a
         && csv_fingerprint.mixed == expected_oracle.csv.mixed));

  REQUIRE(recorder.writeBinary(binary_path) == Status::Success);
  const auto binary_bytes = readBytes(binary_path);
  const auto binary_fingerprint = fingerprint(binary_bytes);
  CAPTURE(binary_fingerprint.bytes,
          binary_fingerprint.fnv1a,
          binary_fingerprint.mixed);
  CHECK(binary_fingerprint.bytes == expected_oracle.binary.bytes);
  CHECK(binary_fingerprint.fnv1a == expected_oracle.binary.fnv1a);
  CHECK(binary_fingerprint.mixed == expected_oracle.binary.mixed);
  const auto endian = scalarAt<std::uint32_t>(binary_bytes, 12);
  CHECK((endian.has_value() && *endian == expected_endian_marker));
  CHECK(headerCrcMatches(binary_bytes, 0, 20));

  std::filesystem::remove(csv_path, ignored);
  std::filesystem::remove(binary_path, ignored);
}

TEST_CASE("Recording common format preserves async bytes and eager indexing",
          "[core][async-recorder][recording-common][MQ2-R1]")
{
  const auto path = temporary("oracle.slcmp");
  std::error_code ignored;
  std::filesystem::remove(path, ignored);

  auto batch = makeBatch();
  core::AsyncRecorder recorder;
  REQUIRE(recorder.configure(
            batch,
            path,
            { .ring_slots = 3,
              .backpressure = core::AsyncBackpressurePolicy::block,
              .codec = core::CompressionCodec::none })
          == Status::Success);
  std::array<ExpectedSnapshot, 2> expected;
  for (std::size_t index = 0; index < expected.size(); ++index) {
    expected[index] = prepareSnapshot(batch, index);
    const auto current = totalCurrent(batch, expected[index]);
    REQUIRE(recorder.enqueue(expected[index].accepted_step, current)
            == Status::Success);
  }
  REQUIRE(recorder.finish() == Status::Success);
  CHECK(recorder.snapshotsWritten() == expected.size());

  const auto bytes = readBytes(path);
  const std::size_t raw_bytes =
    (expected.front().current_density.size()
     + expected.front().state.size())
    * sizeof(core::real_t);
  const std::size_t first_block = header_bytes;
  const std::size_t second_block =
    first_block + header_bytes + raw_bytes;
  BlockCrcPair first_crc, second_crc;
  const bool first_block_matches =
    asyncBlockMatches(bytes, first_block, expected[0], first_crc);
  const bool second_block_matches =
    asyncBlockMatches(bytes, second_block, expected[1], second_crc);
  const auto async_fingerprint = fingerprint(bytes);
  CAPTURE(async_fingerprint.bytes,
          async_fingerprint.fnv1a,
          async_fingerprint.mixed,
          first_crc.raw,
          first_crc.payload,
          second_crc.raw,
          second_crc.payload);
  CHECK(async_fingerprint.bytes == expected_oracle.async.bytes);
  CHECK(async_fingerprint.fnv1a == expected_oracle.async.fnv1a);
  CHECK(async_fingerprint.mixed == expected_oracle.async.mixed);
  const auto endian = scalarAt<std::uint32_t>(bytes, 12);
  CHECK((endian.has_value() && *endian == expected_endian_marker));
  CHECK(headerCrcMatches(bytes, 0, 20));
  CHECK((first_block_matches
         && first_crc.raw == expected_oracle.first_block.raw
         && first_crc.payload == expected_oracle.first_block.payload));
  CHECK((second_block_matches
         && second_crc.raw == expected_oracle.second_block.raw
         && second_crc.payload == expected_oracle.second_block.payload));

  core::CompressedRecording decoded;
  REQUIRE(decoded.open(path) == Status::Success);
  REQUIRE(decoded.size() == expected.size());
  CHECK(decoded.n_rows() == batch.state().n_rows());
  CHECK(decoded.n_lanes() == batch.state().n_lanes());
  CHECK(decoded.stride() == batch.state().stride());
  for (std::size_t index = 0; index < expected.size(); ++index)
    checkSnapshot(decoded.snapshot(index), expected[index]);

  std::filesystem::remove(path, ignored);
}
