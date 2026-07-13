/**
 * @file core_ForwardSensitivity_test.cpp
 * @brief P7-G2 dual propagation, production primal, and centered-FD arbiter.
 */

#include "../../src/core/ForwardSensitivity.hpp"
#include "../../src/core/ParameterSet.hpp"
#include "../support/CoreSpmTestHarness.hpp"
#include "../support/RecordedBits.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <vector>

using namespace slide;

namespace {

std::vector<double> productionTrace(const core::SpmFactoryInput &input,
                                    double c_rate,
                                    double duration,
                                    double sample_step,
                                    int nch = 12)
{
  core::SpmModelOptions options;
  options.nch = nch;
  auto batch = test_support::requireSpmBatch(input, options, 1);
  const double current = c_rate * batch.capacity_Ah();
  const std::array density{ current / batch.electrode_area() };
  const auto samples = static_cast<std::size_t>(std::ceil(duration / sample_step)) + 1;
  std::vector<double> sample_time(samples);
  double time{};
  for (std::size_t sample = 1; sample < samples; ++sample) {
    const double dt = std::min(sample_step, duration - time);
    time += dt;
    sample_time[sample] = time;
    // The harness re-derives dt by differencing this grid; on a grid where the
    // accumulation rounds, that would silently integrate a different problem.
    REQUIRE(sample_time[sample] - sample_time[sample - 1] == dt);
  }
  std::vector<double> voltage(samples);
  core::ExponentialModal stepper;
  test_support::requireConstantCurrentTrace(
    batch,
    stepper,
    test_support::CurrentDensityApm2{
      std::span<const double>{ density } },
    sample_time,
    voltage);
  return voltage;
}

void setParameter(core::SpmFactoryInput &input,
                  core::SensitivityParameter parameter,
                  double value)
{
  const auto neg = core::domain_index(core::Domain::neg);
  const auto pos = core::domain_index(core::Domain::pos);
  switch (parameter) {
  case core::SensitivityParameter::negative_diffusivity:
    input.design.electrode[neg].active_material.D_s.reference_value = value;
    break;
  case core::SensitivityParameter::positive_diffusivity:
    input.design.electrode[pos].active_material.D_s.reference_value = value;
    break;
  case core::SensitivityParameter::negative_reaction_rate:
    input.design.electrode[neg].active_material.k_ct.reference_value = value;
    break;
  case core::SensitivityParameter::positive_reaction_rate:
    input.design.electrode[pos].active_material.k_ct.reference_value = value;
    break;
  case core::SensitivityParameter::contact_resistance:
    input.initial_current_collector_resistance = value * input.design.electrode_area;
    break;
  case core::SensitivityParameter::nominal_capacity:
    input.design.capacity_Ah = value;
    break;
  case core::SensitivityParameter::negative_minimum_stoichiometry:
    input.design.electrode[neg].active_material.x_0 = value;
    break;
  case core::SensitivityParameter::negative_maximum_stoichiometry:
    input.design.electrode[neg].active_material.x_100 = value;
    break;
  case core::SensitivityParameter::positive_minimum_stoichiometry:
    input.design.electrode[pos].active_material.x_100 = value;
    break;
  case core::SensitivityParameter::positive_maximum_stoichiometry:
    input.design.electrode[pos].active_material.x_0 = value;
    break;
  }
}

double characteristicScale(core::SensitivityParameter parameter)
{
  switch (parameter) {
  case core::SensitivityParameter::negative_diffusivity:
  case core::SensitivityParameter::positive_diffusivity:
    return 1e-14;
  case core::SensitivityParameter::negative_reaction_rate:
  case core::SensitivityParameter::positive_reaction_rate:
    return 1e-12;
  case core::SensitivityParameter::contact_resistance:
    return 1e-3;
  case core::SensitivityParameter::nominal_capacity:
    return 1.0;
  default:
    return 0.1;
  }
}

} // namespace

TEST_CASE("Dual observation includes entropic OCV value and tangent away from reference temperature",
          "[core][sensitivity][PC-10][entropic]")
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  REQUIRE(parameters.set("Initial state-of-charge", 0.8, "PC-10 entropic")
          == Status::Success);
  core::SpmFactoryInput zero;
  REQUIRE(parameters.toSpmInput(zero) == Status::Success);
  zero.initial_temperature = zero.design.thermal.reference_temperature + 10.0;
  core::SpmFactoryInput entropic = zero;
  entropic.total_entropic_coefficient = {
    .stoichiometry = { 0.0, 1.0 },
    .value = { 1e-4, 1.1e-3 },
  };

  constexpr auto parameter =
    core::SensitivityParameter::positive_minimum_stoichiometry;
  constexpr std::array selected{ parameter };
  core::ForwardSensitivitySolution zero_dual;
  core::ForwardSensitivitySolution entropic_dual;
  REQUIRE(core::solveCcForwardSensitivities(
            zero, 5, 1.0, true, core::Direction::discharge, 0.0, 1.0, selected, zero_dual)
          == Status::Success);
  REQUIRE(core::solveCcForwardSensitivities(
            entropic, 5, 1.0, true, core::Direction::discharge, 0.0, 1.0, selected, entropic_dual)
          == Status::Success);

  const double dual_primal_increment = entropic_dual.terminal_voltage[0]
                                       - zero_dual.terminal_voltage[0];
  const double production_primal_increment =
    productionTrace(entropic, 1.0, 0.0, 1.0, 5)[0]
    - productionTrace(zero, 1.0, 0.0, 1.0, 5)[0];
  CAPTURE(dual_primal_increment, production_primal_increment);
  CHECK(std::abs(dual_primal_increment - production_primal_increment) <= 2e-12);

  constexpr double h = 1e-5;
  auto lower_zero = zero;
  auto upper_zero = zero;
  auto lower_entropic = entropic;
  auto upper_entropic = entropic;
  const double value = core::sensitivityParameterValue(entropic, parameter);
  setParameter(lower_zero, parameter, value - h);
  setParameter(lower_entropic, parameter, value - h);
  setParameter(upper_zero, parameter, value + h);
  setParameter(upper_entropic, parameter, value + h);
  const double lower_increment = productionTrace(lower_entropic, 1.0, 0.0, 1.0, 5)[0]
                                 - productionTrace(lower_zero, 1.0, 0.0, 1.0, 5)[0];
  const double upper_increment = productionTrace(upper_entropic, 1.0, 0.0, 1.0, 5)[0]
                                 - productionTrace(upper_zero, 1.0, 0.0, 1.0, 5)[0];
  const double finite_difference = (upper_increment - lower_increment) / (2.0 * h);
  const double dual_tangent_increment = entropic_dual.derivative[0]
                                        - zero_dual.derivative[0];
  CAPTURE(finite_difference, dual_tangent_increment);
  CHECK(std::abs(finite_difference - 0.008) <= 2e-9);
  CHECK(std::abs(dual_tangent_increment - finite_difference) <= 2e-9);
}

TEST_CASE("P7-G2 ten dual sensitivities match centered-FD arbiters",
          "[core][sensitivity][P7-G2]")
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  REQUIRE(parameters.set("Initial state-of-charge", 0.8, "P7-G2")
          == Status::Success);
  REQUIRE(parameters.set("Contact resistance [Ohm]", 1e-3, "P7-G2")
          == Status::Success);
  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);

  constexpr double c_rate = 1.0;
  constexpr double duration = 600.0;
  constexpr double sample_step = 10.0;
  core::ForwardSensitivitySolution sensitivity;
  REQUIRE(core::solveCcForwardSensitivities(
            input, 12, c_rate, true, core::Direction::discharge, duration, sample_step, core::supported_sensitivity_parameters, sensitivity)
          == Status::Success);

  test_support::RecordedBits recorded;
  recorded.append(sensitivity.time);
  recorded.append(sensitivity.terminal_voltage);
  recorded.append(sensitivity.derivative);
  CAPTURE(recorded.values, recorded.fnv1a, recorded.mixed);
  REQUIRE(recorded.values == 732);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#if defined(SLIDE_TEST_IPO) && defined(__FAST_MATH__)
  constexpr auto expected_fnv = UINT64_C(0x56b17a9abb12f827);
  constexpr auto expected_mixed = UINT64_C(0x0ba390f9803cf187);
#elif defined(__FAST_MATH__)
  constexpr auto expected_fnv = UINT64_C(0x06dc78ab74d5e7a4);
  constexpr auto expected_mixed = UINT64_C(0xd10596df0d2c2c5a);
#else
  constexpr auto expected_fnv = UINT64_C(0x6dd17b952b7e315a);
  constexpr auto expected_mixed = UINT64_C(0xc286fa543130b91b);
#endif
  REQUIRE(recorded.fnv1a == expected_fnv);
  CHECK(recorded.mixed == expected_mixed);
#endif

  const auto production = productionTrace(input, c_rate, duration, sample_step);
  REQUIRE(production.size() == sensitivity.time.size());
  double maximum_primal_error{};
  for (std::size_t sample = 0; sample < production.size(); ++sample)
    maximum_primal_error = std::max(
      maximum_primal_error,
      std::abs(production[sample] - sensitivity.terminal_voltage[sample]));
  CAPTURE(maximum_primal_error);
  REQUIRE(maximum_primal_error <= 2e-12);

  const double root_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
  for (std::size_t p = 0; p < sensitivity.parameters.size(); ++p) {
    const auto parameter = sensitivity.parameters[p];
    const double value = core::sensitivityParameterValue(input, parameter);
    const double h = root_epsilon
                     * std::max(std::abs(value), characteristicScale(parameter));
    double worst_excess{};
    double maximum_difference{};
    for (const double factor : { 0.5, 1.0, 2.0 }) {
      const double fd_step = factor * h;
      core::SpmFactoryInput lower = input;
      core::SpmFactoryInput upper = input;
      setParameter(lower, parameter, value - fd_step);
      setParameter(upper, parameter, value + fd_step);
      const auto lower_voltage = productionTrace(lower, c_rate, duration, sample_step);
      const auto upper_voltage = productionTrace(upper, c_rate, duration, sample_step);
      for (std::size_t sample = 0; sample < production.size(); ++sample) {
        const double fd = value
                          * (upper_voltage[sample] - lower_voltage[sample])
                          / (2.0 * fd_step);
        const double dual = value * sensitivity.derivativeAt(sample)[p];
        const double difference = std::abs(dual - fd);
        const double tolerance = 2e-6 + 2e-4 * std::abs(fd);
        maximum_difference = std::max(maximum_difference, difference);
        worst_excess = std::max(worst_excess, difference - tolerance);
      }
    }
    CAPTURE(core::sensitivityParameterName(parameter), value, h, maximum_difference, worst_excess);
    CHECK(worst_excess <= 0.0);
  }
}

TEST_CASE("Forward sensitivities reject malformed requests before publication",
          "[core][sensitivity][validation]")
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);

  core::SensitivityParameter parsed{};
  REQUIRE(core::parseSensitivityParameter("not a parameter", parsed)
          == Status::Invalid_parameters);

  core::ForwardSensitivitySolution unchanged;
  unchanged.time = { 42.0 };
  REQUIRE(core::solveCcForwardSensitivities(input, 5, 1.0, true, core::Direction::discharge, 0.0, 1.0, {}, unchanged)
          == Status::Invalid_parameters);
  REQUIRE(unchanged.time == std::vector<double>{ 42.0 });

  constexpr std::array duplicate{
    core::SensitivityParameter::nominal_capacity,
    core::SensitivityParameter::nominal_capacity,
  };
  REQUIRE(core::solveCcForwardSensitivities(input, 5, 1.0, true, core::Direction::discharge, 0.0, 1.0, duplicate, unchanged)
          == Status::Invalid_parameters);

  auto invalid_input = input;
  invalid_input.design.capacity_Ah = 0.0;
  constexpr std::array one_parameter{
    core::SensitivityParameter::nominal_capacity,
  };
  REQUIRE(core::solveCcForwardSensitivities(invalid_input, 5, 1.0, true, core::Direction::discharge, 0.0, 1.0, one_parameter, unchanged)
          == Status::Invalid_parameters);

  REQUIRE(core::solveCcForwardSensitivities(
            input, 5, 1.0, true, core::Direction::discharge, std::numeric_limits<double>::max(), 1.0, one_parameter, unchanged)
          == Status::Invalid_parameters);
  REQUIRE(core::solveCcForwardSensitivities(
            input, 5, 1.0, true, core::Direction::discharge, std::ldexp(1.0, 64), 1.0, one_parameter, unchanged)
          == Status::Invalid_parameters);
  REQUIRE(core::solveCcForwardSensitivities(
            input, 5, 1.0, true, core::Direction::discharge, std::ldexp(1.0, 61), 1.0, core::supported_sensitivity_parameters, unchanged)
          == Status::Invalid_parameters);
  REQUIRE(core::solveCcForwardSensitivities(input, 7, 1.0, true, core::Direction::discharge, 0.0, 1.0, one_parameter, unchanged)
          == Status::Invalid_parameters);
  REQUIRE(unchanged.time == std::vector<double>{ 42.0 });

  const double tiny_duration = 1e-301;
  const double rounded_seventh = 1.4285714285714284e-302;
  core::ForwardSensitivitySolution rounded_schedule;
  REQUIRE(core::solveCcForwardSensitivities(
            input, 5, 1.0, false, core::Direction::discharge, tiny_duration, rounded_seventh, one_parameter, rounded_schedule)
          == Status::Success);
  REQUIRE(rounded_schedule.time.size() == 8);
  REQUIRE(rounded_schedule.time.back() == tiny_duration);

  core::ForwardSensitivitySolution underflowed_ratio;
  REQUIRE(core::solveCcForwardSensitivities(
            input, 5, 1.0, false, core::Direction::discharge, std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), one_parameter, underflowed_ratio)
          == Status::Success);
  REQUIRE(underflowed_ratio.time.size() == 2);
  REQUIRE(underflowed_ratio.time.back()
          == std::numeric_limits<double>::min());

  for (const double rounded_duration : {
         std::nextafter(3.0, 0.0),
         std::nextafter(3.0, 4.0) }) {
    core::ForwardSensitivitySolution snapped_horizon;
    REQUIRE(core::solveCcForwardSensitivities(
              input, 5, 1.0, false, core::Direction::discharge, rounded_duration, 1.0, one_parameter, snapped_horizon)
            == Status::Success);
    REQUIRE(snapped_horizon.time.size() == 4);
    REQUIRE(snapped_horizon.time.back() == rounded_duration);
  }
}

TEST_CASE("Forward sensitivities propagate invalid electrochemical states",
          "[core][sensitivity][validation]")
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  REQUIRE(parameters.set("Initial state-of-charge", 0.8, "coverage")
          == Status::Success);
  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);
  constexpr std::array selected{
    core::SensitivityParameter::nominal_capacity,
  };
  core::ForwardSensitivitySolution solution;

  REQUIRE(core::solveCcForwardSensitivities(input, 5, 1e100, true, core::Direction::discharge, 0.0, 1.0, selected, solution)
          == Status::Invalid_states);

  core::ForwardSensitivitySolution initial;
  REQUIRE(core::solveCcForwardSensitivities(input, 5, 1.0, true, core::Direction::discharge, 0.0, 3600.0, selected, initial)
          == Status::Success);
  REQUIRE(core::solveCcForwardSensitivities(input, 5, 1.0, true, core::Direction::discharge, 3600.0, 3600.0, selected, solution)
          == Status::Invalid_states);
}
