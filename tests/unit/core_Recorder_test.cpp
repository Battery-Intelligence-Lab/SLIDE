/**
 * @file core_Recorder_test.cpp
 * @brief Phase-6 lazy-observable and hardened recording gates.
 */

#include "../../src/core/EulerLegacy.hpp"
#include "../../src/core/Recorder.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
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

std::array<double, 2> liveVoltage(core::SpmBatch &batch,
                                  const std::array<double, 2> &current)
{
  const std::array density{ current[0] / batch.electrode_area(),
                            current[1] / batch.electrode_area() };
  std::array<double, 2> voltage{};
  REQUIRE(batch.terminalVoltage({ .i_app = density }, voltage) == Status::Success);
  return voltage;
}

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

} // namespace

TEST_CASE("Recorder rejects invalid derived metadata and preserves thin ordering",
          "[core][recorder][validation][P9]")
{
  SECTION("invalid backpressure policy is rejected before configuration")
  {
    auto batch = makeBatch();
    core::Recorder recorder;
    CHECK(recorder.configure(
            batch,
            { .capacity = 1,
              .backpressure = static_cast<core::BackpressurePolicy>(255) })
          == Status::Invalid_parameters);
    CHECK_FALSE(recorder.configured());
  }

  SECTION("thin mode tracks the last omitted cadence point")
  {
    auto batch = makeBatch();
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
    auto batch = makeBatch();
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
    auto batch = makeBatch(1e-300);
    core::Recorder recorder;
    REQUIRE(recorder.configure(batch, { .capacity = 1 })
            == Status::Success);
    const std::array current{ 1e300, -1e300 };
    CHECK(recorder.record(0, current) == Status::Invalid_parameters);
    CHECK(recorder.size() == 0);
  }
}

TEST_CASE("P6-G1 recorded states lazily reproduce live voltage",
          "[core][recorder][lazy][P6-G1]")
{
  auto batch = makeBatch();
  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .cadence = 2, .capacity = 3 })
          == Status::Success);
  const std::array current{ 8.0, -4.0 };
  const auto voltage0 = liveVoltage(batch, current);
  REQUIRE(recorder.record(0, current) == Status::Success);

  core::EulerLegacy stepper{ batch };
  const std::array density{ current[0] / batch.electrode_area(),
                            current[1] / batch.electrode_area() };
  REQUIRE(stepper.step(batch, density, 0.0, 1.0) == Status::Success);
  REQUIRE(recorder.record(1, current) == Status::Success); // cadence skip
  REQUIRE(stepper.step(batch, density, 1.0, 1.0) == Status::Success);
  const auto voltage2 = liveVoltage(batch, current);
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

  auto batch = makeBatch();
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
  const auto truncated = temporary("truncated.slrec");
  std::error_code ignored;
  for (const auto &path : { valid, header_corrupt, offset_corrupt, truncated })
    std::filesystem::remove(path, ignored);

  auto batch = makeBatch();
  core::Recorder recorder;
  REQUIRE(recorder.configure(batch, { .capacity = 1 }) == Status::Success);
  const std::array current{ 8.0, -4.0 };
  REQUIRE(recorder.record(0, current) == Status::Success);
  REQUIRE(recorder.writeBinary(valid) == Status::Success);
  REQUIRE(std::filesystem::copy_file(valid, header_corrupt));
  REQUIRE(std::filesystem::copy_file(valid, offset_corrupt));
  REQUIRE(std::filesystem::copy_file(valid, truncated));
  flipByte(header_corrupt, 24); // first dimension byte; CRC must catch it
  flipByte(offset_corrupt, 64); // first absolute offset
  const auto truncated_size = std::filesystem::file_size(truncated);
  REQUIRE(truncated_size > 64);
  std::filesystem::resize_file(truncated, truncated_size - 1);

  core::BinaryRecording recording;
  REQUIRE(recording.open(valid) == Status::Success);
  CHECK(recording.open(header_corrupt) == Status::Invalid_parameters);
  CHECK(recording.valid()); // failed open is atomic
  CHECK(recording.open(offset_corrupt) == Status::Invalid_parameters);
  CHECK(recording.open(truncated) == Status::Invalid_parameters);
  recording.close();

  for (const auto &path : { valid, header_corrupt, offset_corrupt, truncated })
    std::filesystem::remove(path, ignored);
}
