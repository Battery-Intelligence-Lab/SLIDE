/**
 * @file core_AgeingKernel_test.cpp
 * @brief M0.6 pre-refactor whole-pipeline ageing trace.
 *
 * The recorded case deliberately enables every ageing model over seven
 * heterogeneous lanes. It freezes model-major/lane-major accumulation, the
 * SEI-to-crack dependency, additive RHS publication, and scratch reset before
 * the 9C-2 scaffolding refactor.
 */

#include "../support/RecordedBits.hpp"
#include "../../src/core/SpmFactory.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <span>
#include <vector>

using namespace slide;

namespace {

constexpr int nch = 5;
constexpr int lanes = 7;

core::SpmFactoryInput make_ageing_input()
{
  core::SpmFactoryInput input;
  input.design.capacity_Ah = 1.0;
  input.design.electrode_area = 0.1;
  input.design.electrolyte.concentration = 1000.0;
  input.design.thermal.reference_temperature = 298.15;
  input.design.thermal.environment_temperature = 298.15;
  input.initial_soc = 0.5;
  input.initial_temperature = 298.15;
  input.initial_sei_thickness = 1.2e-9;
  input.initial_crack_surface_fraction = 0.012;
  input.initial_current_collector_resistance = 0.2325e-3;
  input.initial_stress_interval = 30.0;
  input.sei_resistivity_area = 2037.4;
  input.negative_laresgoiti_stress = {
    .stoichiometry = { 0.0, 1.0 }, .value = { 0.0, 1.0e6 }
  };

  for (const core::Domain domain : core::domains) {
    const auto d = core::domain_index(domain);
    auto &electrode = input.design.electrode[d];
    electrode.thickness = domain == core::Domain::neg ? 75e-6 : 87e-6;
    electrode.porosity = 0.3;
    electrode.active_fraction = domain == core::Domain::neg ? 0.55 : 0.52;
    electrode.particle_radius = domain == core::Domain::neg ? 12.5e-6 : 8.5e-6;
    electrode.stress = { .youngs_modulus = domain == core::Domain::neg ? 15e9 : 10e9,
                         .poisson_ratio = 0.3,
                         .partial_molar_volume = domain == core::Domain::neg ? 1.1e-6 : 0.9e-6 };
    auto &material = electrode.active_material;
    material.ocv.stoichiometry = { 0.0, 1.0 };
    material.ocv.value = domain == core::Domain::neg
                           ? std::vector<double>{ 0.05, 0.15 }
                           : std::vector<double>{ 3.2, 4.2 };
    material.cs_max = domain == core::Domain::neg ? 30'555.0 : 51'385.0;
    material.x_0 = domain == core::Domain::neg ? 0.2 : 0.8;
    material.x_100 = domain == core::Domain::neg ? 0.8 : 0.2;
    material.D_s = { .reference_value = domain == core::Domain::neg ? 7e-14 : 8e-14,
                     .activation_energy = domain == core::Domain::neg ? 7000.0 : 29'000.0,
                     .reference_temperature = 298.15 };
    material.k_ct = { .reference_value = domain == core::Domain::neg ? 1.764e-11 : 5e-11,
                      .activation_energy = domain == core::Domain::neg ? 20'000.0 : 58'000.0,
                      .reference_temperature = 298.15 };
    input.initial_specific_resistance[d] = domain == core::Domain::neg ? 2.8e-3 : 3.1e-3;

    input.lam.model1_stress_coefficient[d] = domain == core::Domain::neg ? 1.1e-24 : 0.8e-24;
    input.lam.model2_linear_flux[d] = domain == core::Domain::neg ? -2.0e-7 : -1.5e-7;
    input.lam.model2_sqrt_flux[d] = domain == core::Domain::neg ? -3.0e-11 : -2.0e-11;
    input.lam.model4_area_coefficient[d] = domain == core::Domain::neg ? 2.0e-10 : 1.0e-10;
  }
  input.lam.model2_activation = 9000.0;
  input.lam.model3_k = 2.5e-12;
  input.lam.model3_k_activation = 5000.0;
  input.lam.model3_equilibrium_potential = 4.1;
  return input;
}

core::SpmModelOptions all_ageing_options()
{
  return { .nch = nch,
           .sei_model_mask = 0x0f,
           .sei_porosity = true,
           .surface_crack_model_mask = 0x1f,
           .surface_crack_diffusivity = true,
           .lam_model_mask = 0x0f,
           .lithium_plating = true };
}

void make_lanes_heterogeneous(core::SpmBatch &batch)
{
  constexpr std::array temperatures{ 286.0, 292.0, 298.15, 304.0, 311.0, 319.0, 327.0 };
  constexpr std::array negative_scales{ 0.42, 0.58, 0.82, 1.0, 1.18, 1.38, 1.56 };
  constexpr std::array positive_scales{ 1.08, 1.04, 1.0, 0.97, 0.94, 0.91, 0.88 };
  constexpr std::array crack_scales{ 0.6, 0.8, 1.0, 1.25, 1.5, 1.8, 2.1 };
  constexpr std::array intervals{ 11.0, 17.0, 23.0, 31.0, 43.0, 59.0, 71.0 };

  auto &state = batch.state();
  const auto &layout = batch.layout();
  const auto neg = core::domain_index(core::Domain::neg);
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    state.at(layout.spm.temperature, 0, lane) = temperatures[i];
    state.at(layout.spm.sei_thickness, 0, lane) *= 0.75 + 0.1 * lane;
    state.at(layout.spm.crack_surface, 0, lane) *= crack_scales[i];
    state.at(layout.stress_history.interval, 0, lane) = intervals[i];
    state.at(layout.stress_history.previous_laresgoiti_negative, 0, lane) +=
      (static_cast<double>(lane) - 3.0) * 2.5e4;
    state.at(layout.stress_history.previous_dai[neg], 0, lane) +=
      (3.0 - static_cast<double>(lane)) * 1.75e4;
    for (const core::Domain domain : core::domains) {
      const auto d = core::domain_index(domain);
      const double scale = domain == core::Domain::neg ? negative_scales[i]
                                                        : positive_scales[i];
      for (int mode = 0; mode < nch; ++mode)
        state.at(layout.spm.z[d], mode, lane) *= scale;
      state.at(layout.spm.specific_surface_area[d], 0, lane) *= 0.94 + 0.02 * lane;
      state.at(layout.spm.electrode_thickness[d], 0, lane) *= 0.97 + 0.01 * lane;
    }
  }
}

slide::test_support::RecordedBits recorded_ageing_trace()
{
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(make_ageing_input(), all_ageing_options(), lanes, batch)
          == Status::Success);
  make_lanes_heterogeneous(batch);

  constexpr std::array current_a{ -2.1, -1.2, -0.35, 0.0, 0.45, 1.15, 2.0 };
  constexpr std::array current_b{ 1.6, -0.8, 0.25, -1.45, 0.9, -2.2, 0.55 };
  const std::array<std::span<const double>, 3> currents{
    std::span<const double>{ current_a }, std::span<const double>{ current_b },
    std::span<const double>{ current_a }
  };

  slide::test_support::RecordedBits recorded;
  recorded.append(std::span<const double>{ batch.state().raw() });
  std::vector<double> first_derivative;
  std::array<double, lanes> voltage{};
  for (std::size_t pass = 0; pass < currents.size(); ++pass) {
    const core::StepCtx ctx{ .time = 2.5 * static_cast<double>(pass),
                             .dt = 0.25,
                             .i_app = currents[pass] };
    REQUIRE(batch.evaluate(ctx) == Status::Success);
    REQUIRE(batch.terminalVoltage(ctx, voltage) == Status::Success);
    REQUIRE(std::all_of(batch.derivative().raw().begin(),
                        batch.derivative().raw().end(),
                        [](double value) { return core::is_finite(value); }));
    recorded.append(std::span<const double>{ batch.derivative().raw() });
    recorded.append(voltage);
    if (pass == 0)
      first_derivative.assign(batch.derivative().raw().begin(),
                              batch.derivative().raw().end());
    if (pass == 2) {
      REQUIRE(first_derivative.size() == batch.derivative().raw().size());
      REQUIRE(std::memcmp(first_derivative.data(),
                          batch.derivative().raw().data(),
                          first_derivative.size() * sizeof(double))
              == 0);
    }
  }
  return recorded;
}

} // namespace

TEST_CASE("9C-2 all-mask ageing trace retains its pre-refactor bits",
          "[core][ageing][9C-2][recorded]")
{
  const auto recorded = recorded_ageing_trace();
  CAPTURE(recorded.values, recorded.fnv1a, recorded.mixed);
  REQUIRE(recorded.values == 1077);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#if defined(SLIDE_TEST_IPO) && defined(__FAST_MATH__)
  constexpr auto expected_fnv = UINT64_C(0x0be580cf849e57e1);
  constexpr auto expected_mixed = UINT64_C(0x4c6451971b74f789);
#elif defined(__FAST_MATH__)
  constexpr auto expected_fnv = UINT64_C(0xdc8f92a59dd8d67f);
  constexpr auto expected_mixed = UINT64_C(0x59085638e06725d3);
#else
  constexpr auto expected_fnv = UINT64_C(0x87119b1b6fa83b81);
  constexpr auto expected_mixed = UINT64_C(0x55296fa07292792f);
#endif
  CHECK(recorded.fnv1a == expected_fnv);
  CHECK(recorded.mixed == expected_mixed);
#endif
}
