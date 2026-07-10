/**
 * @file core_ForwardSensitivity_test.cpp
 * @brief P7-G2 dual propagation, production primal, and centered-FD arbiter.
 */

#include "../../src/core/ForwardSensitivity.hpp"
#include "../../src/core/ParameterSet.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace slide;

namespace {

std::vector<double> productionTrace(const core::SpmFactoryInput &input,
                                    double c_rate,
                                    double duration,
                                    double sample_step)
{
  core::SpmModelOptions options;
  options.nch = 12;
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, options, 1, batch) == Status::Success);
  const double current = c_rate * batch.capacity_Ah();
  const std::array density{ current / batch.electrode_area() };
  const auto samples = static_cast<std::size_t>(std::ceil(duration / sample_step)) + 1;
  std::vector<double> voltage(samples);
  const core::StepCtx initial{ .time = 0.0, .dt = 0.0, .i_app = density };
  REQUIRE(batch.terminalVoltage(initial, std::span<double>{ voltage }.first(1))
          == Status::Success);
  core::ExponentialModal stepper;
  REQUIRE(stepper.configure(batch) == Status::Success);
  double time{};
  for (std::size_t sample = 1; sample < samples; ++sample) {
    const double dt = std::min(sample_step, duration - time);
    REQUIRE(stepper.step(batch, density, time, dt) == Status::Success);
    time += dt;
    voltage[sample] = stepper.terminalVoltage()[0];
  }
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
