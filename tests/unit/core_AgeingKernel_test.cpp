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
#include "../../src/core/AgeingKernel.hpp"
#include "../../src/core/Dual.hpp"
#include "../../src/core/SpmFactory.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
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

template <class Real>
std::array<Real, 2> crack_diffusivity_at(Real crack_surface)
{
  constexpr int local_nch = 1;
  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<local_nch>(builder);
  const auto history = core::declareStressHistory(builder);
  auto geometry = builder.build(1);
  const core::BatchShape shape = core::BatchShape::from(geometry);
  std::vector<Real> storage(shape.storage_size());
  core::BasicBatchView<Real> mutable_state{ shape, storage };
  mutable_state.at(layout.temperature, 0, 0) = Real{ 298.15 };
  mutable_state.at(layout.sei_thickness, 0, 0) = Real{ 1e-9 };
  mutable_state.at(layout.crack_surface, 0, 0) = crack_surface;
  const auto neg = core::domain_index(core::Domain::neg);
  mutable_state.at(layout.specific_surface_area[neg], 0, 0) = Real{ 4.0e4 };
  mutable_state.at(layout.electrode_thickness[neg], 0, 0) = Real{ 75e-6 };
  mutable_state.at(layout.diffusion_coefficient[neg], 0, 0) = Real{ 7e-14 };
  mutable_state.at(history.interval, 0, 0) = Real{ 30.0 };
  const core::BasicBatchView<const Real> state{ shape, storage };

  core::SpmObservableScratch<local_nch, Real> observable_scratch{ 1 };
  auto observables = observable_scratch.view();
  core::SpmStressScratch<Real> stress_scratch{ 1 };
  auto stress = stress_scratch.view();
  core::SurfaceCrackScratch<Real> output_scratch{ 1 };
  auto output = output_scratch.view();
  const std::array<Real, 1> current_density{ Real{ 2.0 } };
  const core::BasicStepCtx<Real> ctx{ .i_app = current_density };

  core::SurfaceCrackParams params;
  params.model_mask = core::surface_crack_model_bit(4);
  params.reduce_negative_diffusivity = true;
  params.electrode_area = 0.1;
  params.negative_cs_max = 30'555.0;
  params.model4_alpha = 4.0e-8;
  params.model4_max_surface = 0.03;
  params.diffusion_exponent = 2.3;
  REQUIRE(core::computeSurfaceCrack(params,
                                    state,
                                    layout,
                                    history,
                                    ctx,
                                    observables,
                                    stress,
                                    output)
          == Status::Success);
  return { output.crack_surface_rate[0], output.negative_diffusivity_rate[0] };
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

TEST_CASE("9C-2 surface-crack diffusivity is scalar-generic for Dual",
          "[core][ageing][9C-2][dual]")
{
  constexpr double crack_surface = 0.01;
  const auto dual = crack_diffusivity_at(core::Dual{ crack_surface, 1.0 });
  const auto primal = crack_diffusivity_at(crack_surface);
  CAPTURE(std::bit_cast<std::uint64_t>(dual[0].value),
          std::bit_cast<std::uint64_t>(primal[0]),
          std::bit_cast<std::uint64_t>(dual[1].value),
          std::bit_cast<std::uint64_t>(primal[1]));
  REQUIRE(dual[0].value == primal[0]);
  REQUIRE(std::abs(dual[1].value - primal[1])
          <= 2.0 * std::numeric_limits<double>::epsilon()
               * std::max(std::abs(primal[1]), 1e-300));

  constexpr double h = 1e-6;
  const auto below = crack_diffusivity_at(crack_surface - h);
  const auto above = crack_diffusivity_at(crack_surface + h);
  const double crack_fd = (above[0] - below[0]) / (2.0 * h);
  const double diffusion_fd = (above[1] - below[1]) / (2.0 * h);
  CAPTURE(dual[0].derivative, crack_fd, dual[1].derivative, diffusion_fd);
  REQUIRE(std::abs(dual[0].derivative - crack_fd)
          <= 1e-8 * std::max(std::abs(crack_fd), 1e-30));
  REQUIRE(std::abs(dual[1].derivative - diffusion_fd)
          <= 1e-8 * std::max(std::abs(diffusion_fd), 1e-30));
}

TEST_CASE("9C-2 common ageing scaffold fixes mask, scratch, and traversal order",
          "[core][ageing][9C-2][scaffold]")
{
  STATIC_REQUIRE(core::ageing_model_bit<4>(0) == 0);
  STATIC_REQUIRE(core::ageing_model_bit<4>(1) == 0x01);
  STATIC_REQUIRE(core::ageing_model_bit<4>(4) == 0x08);
  STATIC_REQUIRE(core::ageing_model_bit<4>(5) == 0);
  STATIC_REQUIRE(core::ageing_model_mask<5>() == 0x1f);
  STATIC_REQUIRE(core::valid_ageing_model_mask<4>(0x0f));
  STATIC_REQUIRE_FALSE(core::valid_ageing_model_mask<4>(0));
  STATIC_REQUIRE_FALSE(core::valid_ageing_model_mask<4>(0x10));
  STATIC_REQUIRE(core::valid_optional_ageing_model_mask<4>(0));

  REQUIRE_THROWS_AS((core::detail::AgeingScratchStorage<double, 3>{ 0 }),
                    std::invalid_argument);
  REQUIRE_THROWS_AS((core::detail::AgeingScratchStorage<double, 3>{ -1 }),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(core::SeiScratch<>{ 0 }, std::invalid_argument);
  REQUIRE_THROWS_AS(core::SeiScratch<>{ -1 }, std::invalid_argument);
  REQUIRE_THROWS_AS(core::SurfaceCrackScratch<>{ 0 }, std::invalid_argument);
  REQUIRE_THROWS_AS(core::SurfaceCrackScratch<>{ -1 }, std::invalid_argument);
  REQUIRE_THROWS_AS(core::LamScratch<>{ 0 }, std::invalid_argument);
  REQUIRE_THROWS_AS(core::LamScratch<>{ -1 }, std::invalid_argument);
  REQUIRE_THROWS_AS(core::LithiumPlatingScratch<>{ 0 }, std::invalid_argument);
  REQUIRE_THROWS_AS(core::LithiumPlatingScratch<>{ -1 }, std::invalid_argument);
  core::detail::AgeingScratchStorage<double, 3> scratch{ 4 };
  REQUIRE(scratch.n_lanes() == 4);
  for (std::size_t field = 0; field < scratch.field_count; ++field) {
    REQUIRE(scratch.field(field).size() == 4);
    std::fill(scratch.field(field).begin(), scratch.field(field).end(),
              static_cast<double>(field + 1));
  }
  core::detail::clear_ageing_fields<double, 3>(
    4, { scratch.field(0), scratch.field(1), scratch.field(2) });
  for (std::size_t field = 0; field < scratch.field_count; ++field)
    REQUIRE(std::all_of(scratch.field(field).begin(),
                        scratch.field(field).end(),
                        [](double value) { return value == 0.0; }));

  std::vector<std::array<unsigned, 2>> visits;
  const auto status = core::detail::for_each_enabled_ageing_model_lane<4>(
    0x05, 3, [&](unsigned model, int lane) {
      visits.push_back({ model, static_cast<unsigned>(lane) });
      return Status::Success;
    });
  REQUIRE(status == Status::Success);
  const std::vector<std::array<unsigned, 2>> expected{
    { 1, 0 }, { 1, 1 }, { 1, 2 }, { 3, 0 }, { 3, 1 }, { 3, 2 }
  };
  REQUIRE(visits == expected);

  int last_lane = -1;
  const auto failure = core::detail::for_each_ageing_lane_while_success(
    5, [&](int lane) {
      last_lane = lane;
      return lane == 2 ? Status::Invalid_states : Status::Success;
    });
  REQUIRE(failure == Status::Invalid_states);
  REQUIRE(last_lane == 2);
}
