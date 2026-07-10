/**
 * @file core_Simulation_test.cpp
 * @brief EulerLegacy row-role, cumulative, rollback, and Simulation façade tests.
 */

#include "../../src/core/Simulation.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

using namespace slide;

namespace {

core::SpmFactoryInput make_simulation_input()
{
  core::SpmFactoryInput input;
  input.design.capacity_Ah = 1.0;
  input.design.electrode_area = 0.1;
  input.design.electrolyte.concentration = 1000.0;
  input.design.thermal = { .density = 2.0,
                           .heat_capacity = 3.0,
                           .volume = 4.0,
                           .surface_area = 5.0,
                           .h_conv = 6.0,
                           .reference_temperature = 298.15,
                           .environment_temperature = 300.0 };
  for (const core::Domain domain : core::domains) {
    const auto d = core::domain_index(domain);
    auto &electrode = input.design.electrode[d];
    electrode.thickness = domain == core::Domain::neg ? 75e-6 : 87e-6;
    electrode.porosity = 0.3;
    electrode.active_fraction = 0.5;
    electrode.particle_radius = domain == core::Domain::neg ? 12.5e-6 : 8.5e-6;
    auto &material = electrode.active_material;
    material.ocv.stoichiometry = { 0.0, 1.0 };
    material.ocv.value = domain == core::Domain::neg
                           ? std::vector<double>{ 0.0, 0.1 }
                           : std::vector<double>{ 3.0, 4.0 };
    material.cs_max = domain == core::Domain::neg ? 30'555.0 : 51'385.0;
    material.x_0 = domain == core::Domain::neg ? 0.2 : 0.8;
    material.x_100 = domain == core::Domain::neg ? 0.8 : 0.2;
    material.D_s = { .reference_value = domain == core::Domain::neg ? 7e-14 : 8e-14,
                     .activation_energy = domain == core::Domain::neg ? 7000.0 : 29'000.0,
                     .reference_temperature = 298.15 };
    material.k_ct = { .reference_value = domain == core::Domain::neg ? 1.764e-11 : 5e-11,
                      .activation_energy = domain == core::Domain::neg ? 20'000.0 : 58'000.0,
                      .reference_temperature = 298.15 };
    input.initial_specific_resistance[d] = 2.8e-3;
  }
  input.initial_current_collector_resistance = 0.2325e-3;
  return input;
}

} // namespace

TEST_CASE("EulerLegacy advances only ODE rows and owns cumulative updates",
          "[core][stepper]")
{
  const auto input = make_simulation_input();
  core::SpmModelOptions options{ .thermal = true };
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, options, 2, batch) == Status::Success);
  core::EulerLegacy stepper{ batch };
  const auto &layout = batch.layout();
  REQUIRE(batch.roles()[static_cast<std::size_t>(layout.elapsed_time.row_begin)]
          == core::StateRole::cumulative);
  REQUIRE(batch.roles()[static_cast<std::size_t>(layout.charge_throughput.row_begin)]
          == core::StateRole::cumulative);
  REQUIRE(batch.roles()[static_cast<std::size_t>(layout.energy_throughput.row_begin)]
          == core::StateRole::cumulative);

  batch.state().at(layout.thermal.external_heat_flow, 0, 0) = 2.0;
  batch.state().at(layout.thermal.external_heat_flow, 0, 1) = -1.0;
  const std::array zero_density{ 0.0, 0.0 };
  REQUIRE(stepper.step(batch, zero_density, 0.0, 1.0) == Status::Success);
  REQUIRE(batch.state().at(layout.thermal.external_heat_flow, 0, 0) == 2.0);
  REQUIRE(batch.state().at(layout.thermal.external_heat_flow, 0, 1) == -1.0);
  for (int lane = 0; lane < 2; ++lane) {
    REQUIRE(batch.state().at(layout.elapsed_time, 0, lane) == 1.0);
    REQUIRE(batch.state().at(layout.charge_throughput, 0, lane) == 0.0);
    REQUIRE(batch.state().at(layout.energy_throughput, 0, lane) == 0.0);
  }
  REQUIRE(batch.state().at(layout.spm.temperature, 0, 0) != 298.15);

  const std::array discharge_density{ 1.0, 1.0 }; // 0.1 A because area=0.1 m2
  REQUIRE(stepper.step(batch, discharge_density, 1.0, 2.0) == Status::Success);
  const double expected_Ah = 0.1 * 2.0 / 3600.0;
  for (int lane = 0; lane < 2; ++lane) {
    CAPTURE(lane);
    REQUIRE(batch.state().at(layout.elapsed_time, 0, lane) == 3.0);
    REQUIRE(std::abs(batch.state().at(layout.charge_throughput, 0, lane)
                     - expected_Ah)
            <= 1e-16);
    REQUIRE(std::abs(batch.state().at(layout.energy_throughput, 0, lane)
                     - expected_Ah
                         * stepper.terminalVoltage()[static_cast<std::size_t>(lane)])
            <= 1e-16);
  }
}

TEST_CASE("Simulation façade solves a partial final CC step", "[core][simulation]")
{
  core::Simulation simulation;
  const auto input = make_simulation_input();
  REQUIRE(simulation.build(input, {}, 2) == Status::Success);
  core::SimulationSolution solution;
  const core::ConstantCurrentExperiment experiment{ .current_A = 0.0,
                                                    .duration = 2.5,
                                                    .step = 1.0 };
  REQUIRE(simulation.solve(experiment, solution) == Status::Success);
  REQUIRE(solution.termination == Status::Success);
  REQUIRE(solution.n_lanes == 2);
  REQUIRE(solution.time == std::vector<double>{ 0.0, 1.0, 2.0, 2.5 });
  REQUIRE(solution.terminal_voltage.size() == 8);
  for (const double voltage : solution.terminal_voltage)
    REQUIRE(std::isfinite(voltage));
  REQUIRE(simulation.batch().state().at(simulation.batch().layout().elapsed_time, 0, 0)
          == 2.5);
  REQUIRE(solution.voltageAt(3).size() == 2);

  core::SimulationSolution unchanged;
  unchanged.time = { 42.0 };
  REQUIRE(simulation.solve({ .current_A = 0.0, .duration = 1.0, .step = 0.0 },
                           unchanged)
          == Status::Invalid_parameters);
  REQUIRE(unchanged.time == std::vector<double>{ 42.0 });

  core::Simulation decimal_steps;
  REQUIRE(decimal_steps.build(input, {}, 1) == Status::Success);
  core::SimulationSolution decimal_solution;
  REQUIRE(decimal_steps.solve({ .current_A = 0.0, .duration = 0.3, .step = 0.1 },
                              decimal_solution)
          == Status::Success);
  REQUIRE(decimal_solution.time.size() == 4);
  REQUIRE(std::abs(decimal_solution.time.back() - 0.3) <= 1e-15);
}
