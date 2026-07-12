/**
 * @file core_P7G3_PyBaMM_test.cpp
 * @brief Hermetic Chen2020 behaviour gates against committed PyBaMM traces.
 */

#include "../../src/core/ExponentialModal.hpp"
#include "../../src/core/ParameterSet.hpp"
#include "../support/CoreSpmTestHarness.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <span>
#include <sstream>
#include <string>
#include <vector>

using namespace slide;

namespace {

struct ReferenceTrace
{
  std::vector<double> time;
  std::vector<double> voltage;
};

ReferenceTrace loadTrace(std::string_view filename)
{
  const auto path = std::filesystem::path{ SLIDE_SOURCE_DIR } / "tests"
                    / "reference" / filename;
  std::ifstream stream(path);
  REQUIRE(stream.good());
  std::string line;
  REQUIRE(std::getline(stream, line));
  REQUIRE(line == "segment,time_s,local_time_s,current_A,voltage_V");
  ReferenceTrace result;
  while (std::getline(stream, line)) {
    std::stringstream row{ line };
    std::array<std::string, 5> fields;
    for (std::size_t i = 0; i < fields.size(); ++i) {
      REQUIRE(std::getline(row, fields[i], ','));
    }
    result.time.push_back(std::stod(fields[1]));
    result.voltage.push_back(std::stod(fields[4]));
  }
  REQUIRE(result.time.size() > 2);
  return result;
}

std::pair<double, double> compareTrace(const ReferenceTrace &reference,
                                       double c_rate)
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  REQUIRE(parameters.set("Initial state-of-charge", 1.0, "P7-G3 fixture")
          == Status::Success);
  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);
  core::SpmModelOptions options;
  options.nch = 12;
  auto batch = test_support::requireSpmBatch(input, options, 1);

  const double current = c_rate * batch.capacity_Ah();
  const std::array current_density{ current / batch.electrode_area() };
  std::vector<double> actual_voltage(reference.time.size());
  core::ExponentialModal stepper;
  test_support::requireConstantCurrentTrace(
    batch,
    stepper,
    test_support::CurrentDensityApm2{
      std::span<const double>{ current_density } },
    reference.time,
    actual_voltage);

  test_support::VoltageError error{};
  REQUIRE(test_support::computeVoltageError(
    actual_voltage, reference.voltage, error));
  return { error.maximum_absolute_V, error.rms_V };
}

} // namespace

TEST_CASE("P7-G3 Chen2020 C/50 matches committed PyBaMM 26.6.2.0 fixture",
          "[core][parameters][PyBaMM][P7-G3]")
{
  const auto [maximum, rms] = compareTrace(
    loadTrace("pybamm_chen2020_c050.csv"), 1.0 / 50.0);
  CAPTURE(maximum, rms);
  CHECK(maximum <= 2e-3);
  CHECK(rms <= 1e-3);
}

TEST_CASE("P7-G3 Chen2020 1C matches committed PyBaMM 26.6.2.0 fixture",
          "[core][parameters][PyBaMM][P7-G3]")
{
  const auto [maximum, rms] = compareTrace(
    loadTrace("pybamm_chen2020_c1.csv"), 1.0);
  CAPTURE(maximum, rms);
  CHECK(maximum <= 15e-3);
  CHECK(rms <= 8e-3);
}
