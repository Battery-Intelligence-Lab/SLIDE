/**
 * @file core_P1G4_restart_test.cpp
 * @brief P1-G4 bitwise restart across thermal and all ageing state classes.
 */

#include "../../src/core/EulerLegacy.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstring>
#include <vector>

using namespace slide;

namespace {

core::SpmFactoryInput make_restart_input()
{
  core::SpmFactoryInput input;
  input.design.capacity_Ah = 1.0;
  input.design.electrode_area = 0.1;
  input.design.electrolyte.concentration = 1000.0;
  input.design.thermal = { .density = 1626.0,
                           .heat_capacity = 750.0,
                           .volume = 1e-4,
                           .surface_area = 0.02,
                           .h_conv = 30.0,
                           .reference_temperature = 298.15,
                           .environment_temperature = 298.15 };
  for (const core::Domain domain : core::domains) {
    const auto d = core::domain_index(domain);
    auto &electrode = input.design.electrode[d];
    electrode.thickness = domain == core::Domain::neg ? 75e-6 : 87e-6;
    electrode.porosity = 0.3;
    electrode.active_fraction = 0.5;
    electrode.particle_radius = domain == core::Domain::neg ? 12.5e-6 : 8.5e-6;
    electrode.stress = { .youngs_modulus = domain == core::Domain::neg ? 15e9 : 10e9,
                         .poisson_ratio = 0.3,
                         .partial_molar_volume = 1e-6 };
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
    input.lam.model1_stress_coefficient[d] = 1e-24;
  }
  input.initial_current_collector_resistance = 0.2325e-3;
  input.sei_resistivity_area = 2037.4;
  input.negative_laresgoiti_stress.stoichiometry = { 0.0, 1.0 };
  input.negative_laresgoiti_stress.value = { 0.0, 1e6 };
  return input;
}

core::SpmModelOptions restart_options()
{
  return { .nch = 5,
           .thermal = true,
           .sei_model_mask = core::sei_model_bit(1),
           .surface_crack_model_mask = core::surface_crack_model_bit(1),
           .lam_model_mask = core::lam_model_bit(1),
           .lithium_plating = true };
}

void advance(core::SpmBatch &batch,
             core::EulerLegacy &stepper,
             int first_step,
             int step_count,
             std::span<const double>
               current_density,
             double dt)
{
  for (int step = first_step; step < first_step + step_count; ++step)
    REQUIRE(stepper.step(batch,
                         current_density,
                         static_cast<double>(step) * dt,
                         dt)
            == Status::Success);
}

} // namespace

TEST_CASE("P1-G4 arena-only restart is bitwise identical", "[core][P1-G4]")
{
  constexpr int lanes = 3;
  constexpr int total_steps = 20;
  constexpr int split_step = 10;
  constexpr double dt = 0.5;
  const auto input = make_restart_input();
  const auto options = restart_options();
  const std::array current_density{ 1.0, 1.5, 2.0 };

  core::SpmBatch continuous;
  core::SpmBatch split;
  REQUIRE(core::buildSpmBatch(input, options, lanes, continuous) == Status::Success);
  REQUIRE(core::buildSpmBatch(input, options, lanes, split) == Status::Success);
  core::EulerLegacy continuous_stepper{ continuous };
  core::EulerLegacy split_stepper{ split };

  advance(continuous, continuous_stepper, 0, total_steps, current_density, dt);
  advance(split, split_stepper, 0, split_step, current_density, dt);

  const core::StateSlice whole_arena{ 0, split.state().n_rows() };
  std::vector<double> checkpoint(split.state().slice_size(whole_arena));
  split.state().snapshot(whole_arena, checkpoint.data());

  // Rebuild every non-state object, then restore only the arena bytes.
  core::SpmBatch restarted;
  REQUIRE(core::buildSpmBatch(input, options, lanes, restarted) == Status::Success);
  restarted.state().restore(whole_arena, checkpoint.data());
  core::EulerLegacy restarted_stepper{ restarted };
  advance(restarted,
          restarted_stepper,
          split_step,
          total_steps - split_step,
          current_density,
          dt);

  REQUIRE(continuous.state().raw().size() == restarted.state().raw().size());
  REQUIRE(std::memcmp(continuous.state().raw().data(),
                      restarted.state().raw().data(),
                      continuous.state().raw().size() * sizeof(double))
          == 0);

  std::array<double, lanes> continuous_voltage{};
  std::array<double, lanes> restarted_voltage{};
  const core::StepCtx final_ctx{ .time = total_steps * dt,
                                 .dt = 0.0,
                                 .i_app = current_density };
  REQUIRE(continuous.terminalVoltage(final_ctx, continuous_voltage) == Status::Success);
  REQUIRE(restarted.terminalVoltage(final_ctx, restarted_voltage) == Status::Success);
  REQUIRE(std::memcmp(continuous_voltage.data(),
                      restarted_voltage.data(),
                      sizeof(continuous_voltage))
          == 0);
}
