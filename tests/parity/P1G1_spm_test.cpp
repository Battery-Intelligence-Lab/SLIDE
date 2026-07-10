/**
 * @file P1G1_spm_test.cpp
 * @brief P1-G1 full-trajectory parity for a Kokam SPM cell in legacy-Euler mode.
 *
 * Registered before first run (PLAN.md sections 5.2 and 6):
 *   - max absolute terminal-voltage error <= 1e-12 V;
 *   - every mapped state satisfies |core-legacy| <= 1e-15 + 1e-12*max(|core|,|legacy|).
 *
 * Scenarios cover a 1200 s mid-SOC 1C discharge and the steep low-SOC OCV tail.
 */

#include "../../src/core/EulerLegacy.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <span>
#include <string_view>

using namespace slide;

namespace {

struct Drift
{
  double voltage_abs{};
  double state_abs{};
  double state_scaled{};
  bool state_within_band{ true };

  void state(double legacy, double current)
  {
    const double error = std::abs(current - legacy);
    const double scale = std::max(std::abs(current), std::abs(legacy));
    const double tolerance = 1e-15 + 1e-12 * scale;
    state_abs = std::max(state_abs, error);
    state_scaled = std::max(state_scaled, error / tolerance);
    state_within_band = state_within_band && error <= tolerance;
  }
};

void compare_state(State_SPM &legacy, const core::SpmBatch &batch, Drift &drift)
{
  const auto &state = batch.state();
  const auto &layout = batch.layout();
  for (const auto legacy_domain : { neg, pos }) {
    const auto domain = legacy_domain == neg ? core::Domain::neg : core::Domain::pos;
    const auto d = core::domain_index(domain);
    for (int mode = 0; mode < static_cast<int>(settings::nch); ++mode)
      drift.state(legacy.z(static_cast<std::size_t>(mode), legacy_domain),
                  state.at(layout.spm.z[d], mode, 0));
    drift.state(legacy.e(legacy_domain), state.at(layout.spm.active_fraction[d], 0, 0));
    drift.state(legacy.D(legacy_domain), state.at(layout.spm.diffusion_coefficient[d], 0, 0));
    drift.state(legacy.thick(legacy_domain), state.at(layout.spm.electrode_thickness[d], 0, 0));
    drift.state(legacy.a(legacy_domain), state.at(layout.spm.specific_surface_area[d], 0, 0));
    drift.state(legacy.rDC(legacy_domain), state.at(layout.spm.specific_resistance[d], 0, 0));
  }
  drift.state(legacy.T(), state.at(layout.spm.temperature, 0, 0));
  drift.state(legacy.delta(), state.at(layout.spm.sei_thickness, 0, 0));
  drift.state(legacy.LLI(), state.at(layout.spm.lost_lithium, 0, 0));
  drift.state(legacy.CS(), state.at(layout.spm.crack_surface, 0, 0));
  drift.state(legacy.delta_pl(), state.at(layout.spm.plated_lithium_thickness, 0, 0));
  drift.state(legacy.rDCcc(), state.at(layout.spm.current_collector_resistance, 0, 0));
  drift.state(legacy.time(), state.at(layout.elapsed_time, 0, 0));
  drift.state(legacy.Ah(), state.at(layout.charge_throughput, 0, 0));
  drift.state(legacy.Wh(), state.at(layout.energy_throughput, 0, 0));
}

struct Scenario
{
  std::string_view name;
  double initial_soc;
  double current_A;
  int steps;
};

} // namespace

TEST_CASE("P1-G1 Kokam SPM trajectories match legacy Euler", "[parity][core][P1-G1]")
{
  constexpr double dt = 1.0;
  constexpr std::array scenarios{
    Scenario{ "mid-SOC 1C discharge", 0.5, 16.0, 1200 },
    Scenario{ "low-SOC steep-tail discharge", 0.15, 16.0, 300 },
  };

  for (const auto &scenario : scenarios) {
    Cell_SPM legacy;
    REQUIRE(test_support::initialize_legacy_kokam(
              legacy, scenario.initial_soc, scenario.current_A)
            == Status::Success);

    double environment_temperature{}, reference_temperature{};
    legacy.getTemperatures(&environment_temperature, &reference_temperature);
    core::SpmBatch batch;
    const auto input = test_support::make_legacy_kokam_input(
      scenario.initial_soc, legacy.T(), reference_temperature);
    REQUIRE(core::buildSpmBatch(input, {}, 1, batch) == Status::Success);
    core::EulerLegacy stepper;
    REQUIRE(stepper.configure(batch) == Status::Success);

    Drift drift;
    compare_state(legacy.getStateObj(), batch, drift);
    const std::array current_density{ scenario.current_A / batch.electrode_area() };
    std::array<double, 1> initial_voltage{};
    const core::StepCtx initial_ctx{ .time = 0.0, .dt = 0.0, .i_app = current_density };
    REQUIRE(batch.terminalVoltage(initial_ctx, initial_voltage) == Status::Success);
    drift.voltage_abs = std::abs(initial_voltage[0] - legacy.V());

    for (int step = 0; step < scenario.steps; ++step) {
      legacy.timeStep_CC(dt);
      REQUIRE(stepper.step(batch, current_density, step * dt, dt) == Status::Success);
      compare_state(legacy.getStateObj(), batch, drift);
      drift.voltage_abs = std::max(
        drift.voltage_abs,
        std::abs(stepper.terminalVoltage()[0] - legacy.V()));
    }

    std::printf("P1-G1 %.*s: dV=%.17g dstate=%.17g scaled=%.17g\n",
                static_cast<int>(scenario.name.size()),
                scenario.name.data(),
                drift.voltage_abs,
                drift.state_abs,
                drift.state_scaled);
    CAPTURE(scenario.name, drift.voltage_abs, drift.state_abs, drift.state_scaled);
    REQUIRE(drift.voltage_abs <= 1e-12);
    REQUIRE(drift.state_within_band);
  }
}
