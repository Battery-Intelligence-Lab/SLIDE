/**
 * @file core_SpmFactory_test.cpp
 * @brief D-02 registry selection, cold validation, initialization, and batch dispatch tests.
 */

#include "../../src/core/SpmFactory.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <span>
#include <vector>

using namespace slide;

namespace {

core::SpmFactoryInput make_input()
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
    electrode.stress = { .youngs_modulus = 10e9,
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
  }
  input.initial_current_collector_resistance = 0.2325e-3;
  input.negative_laresgoiti_stress.stoichiometry = { 0.0, 1.0 };
  input.negative_laresgoiti_stress.value = { 0.0, 1e6 };
  return input;
}

} // namespace

TEST_CASE("SPM factory selects registered compositions and initializes batches",
          "[core][factory][registry]")
{
  const auto input = make_input();
  core::SpmBatch isothermal;
  core::SpmModelOptions options;
  REQUIRE(core::buildSpmBatch(input, options, 3, isothermal) == Status::Success);
  REQUIRE(isothermal.valid());
  REQUIRE(isothermal.nch() == 5);
  REQUIRE(isothermal.n_lanes() == 3);
  REQUIRE(isothermal.composition() == core::SpmComposition::isothermal);
  REQUIRE(static_cast<int>(isothermal.roles().size()) == isothermal.state().n_rows());
  REQUIRE(isothermal.layout().spm.z[core::domain_index(core::Domain::neg)].rows == 5);

  constexpr std::array zero_current{ 0.0, 0.0, 0.0 };
  const core::StepCtx ctx{ .time = 0.0, .dt = 1.0, .i_app = zero_current };
  REQUIRE(isothermal.evaluate(ctx) == Status::Success);

  core::SpmBatch thermal;
  options.thermal = true;
  REQUIRE(core::buildSpmBatch(input, options, 3, thermal) == Status::Success);
  REQUIRE(thermal.composition() == core::SpmComposition::thermal);
  REQUIRE(thermal.state().n_rows() == isothermal.state().n_rows() + 3);
  REQUIRE(thermal.roles()[static_cast<std::size_t>(thermal.layout().thermal.external_heat_flow.row_begin)]
          == core::StateRole::input);

  core::SpmBatch ageing;
  options.sei_model_mask = core::sei_model_bit(1);
  options.surface_crack_model_mask = core::surface_crack_model_bit(1);
  REQUIRE(core::buildSpmBatch(input, options, 3, ageing) == Status::Success);
  REQUIRE(ageing.composition() == core::SpmComposition::thermal_ageing);
  REQUIRE(ageing.state().n_rows() == thermal.state().n_rows() + 4);
  REQUIRE(ageing.storeStressHistory(2.0) == Status::Success);
  for (int lane = 0; lane < ageing.n_lanes(); ++lane)
    REQUIRE(ageing.state().at(ageing.layout().stress_history.interval, 0, lane) == 2.0);
}

TEST_CASE("SPM registry exposes all validated spatial orders", "[core][factory][registry]")
{
  const auto input = make_input();
  for (const int nch : core::registered_spm_nch) {
    core::SpmModelOptions options{ .nch = nch };
    core::SpmBatch batch;
    REQUIRE(core::buildSpmBatch(input, options, 1, batch) == Status::Success);
    CAPTURE(nch);
    REQUIRE(batch.nch() == nch);
    REQUIRE(batch.layout().spm.z[core::domain_index(core::Domain::pos)].rows == nch);
  }
}

TEST_CASE("SPM factory rejects malformed cold inputs without replacing output",
          "[core][factory][validation]")
{
  auto input = make_input();
  core::SpmModelOptions options;
  core::SpmBatch output;
  REQUIRE(core::buildSpmBatch(input, options, 2, output) == Status::Success);
  REQUIRE(output.n_lanes() == 2);

  options.nch = 7;
  REQUIRE(core::buildSpmBatch(input, options, 4, output) == Status::Invalid_parameters);
  REQUIRE(output.valid());
  REQUIRE(output.n_lanes() == 2);
  REQUIRE(output.nch() == 5);

  options.nch = 5;
  options.sei_model_mask = 0x80;
  REQUIRE(core::buildSpmBatch(input, options, 4, output) == Status::Invalid_parameters);
  options = {};
  options.sei_porosity = true;
  REQUIRE(core::buildSpmBatch(input, options, 4, output) == Status::Invalid_parameters);
  options = {};
  options.surface_crack_diffusivity = true;
  REQUIRE(core::buildSpmBatch(input, options, 4, output) == Status::Invalid_parameters);
  REQUIRE(core::buildSpmBatch(input, {}, 0, output) == Status::Invalid_parameters);

  options = {};
  input.design.electrode[core::domain_index(core::Domain::neg)]
    .active_material.ocv.value.pop_back();
  REQUIRE(core::buildSpmBatch(input, options, 4, output) == Status::Invalid_parameters);
  REQUIRE(output.n_lanes() == 2);
}

TEST_CASE("SPM batch validates its generic RHS boundary", "[core][factory][rhs]")
{
  const auto input = make_input();
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);
  constexpr std::array current{ 0.0, 0.0 };
  const core::StepCtx ctx{ .dt = 1.0, .i_app = current };
  std::vector<double> too_short(batch.state().size() - 1);
  REQUIRE(batch.rhs(too_short, batch.derivative().raw(), ctx)
          == Status::Invalid_parameters);
  REQUIRE(batch.storeStressHistory(0.0) == Status::Invalid_parameters);
}
