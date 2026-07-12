/**
 * @file core_ExponentialModal_test.cpp
 * @brief Phase-3 exact modal, conservation, and stability gates.
 */

#include "../../src/core/ExponentialModal.hpp"
#include "../../src/core/EulerLegacy.hpp"
#include "../../src/core/ParameterSet.hpp"
#include "../../src/core/SpectralModel.hpp"
#include "../support/CoreSpmTestHarness.hpp"
#include "../support/KokamSpmFixture.hpp"
#include "../support/RecordedBits.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace slide;

namespace {

template <int NCH>
void exactModalGate()
{
  constexpr double dt = 25.0;
  constexpr double current = 2.0;
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, { .nch = NCH }, 1, batch) == Status::Success);

  core::PerDomain<double> radius{};
  for (const auto domain : core::domains)
    radius[core::domain_index(domain)] = input.design.electrode[core::domain_index(domain)].particle_radius;
  core::CompiledSpectralModel<NCH> model;
  REQUIRE(core::compileSpectralModel<NCH>(radius, model) == Status::Success);

  core::PerDomain<std::array<double, NCH>> initial{};
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    for (int mode = 0; mode < NCH; ++mode) {
      auto &value = batch.state().at(batch.layout().spm.z[d], mode, 0);
      value += 0.01 * (mode + 1);
      initial[d][static_cast<std::size_t>(mode)] = value;
    }
  }

  const std::array density{ current / batch.electrode_area() };
  core::ExponentialModal stepper{ batch };
  REQUIRE(stepper.step(batch, density, 0.0, dt) == Status::Success);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    const auto &electrode = input.design.electrode[d];
    const double area = 3.0 * electrode.active_fraction / electrode.particle_radius;
    const double flux = static_cast<double>(core::molar_flux_sign(domain)) * density[0]
                        / (area * 96487.0 * electrode.thickness);
    const double diffusivity = electrode.active_material.D_s.reference_value;
    for (int mode = 0; mode < NCH; ++mode) {
      // Evaluate the diagonal closed form through a deliberately independent
      // path: long-double libm, with no copy of the production small-x Taylor
      // branch. The shared compiled A/B coefficients are the ODE definition;
      // this gate arbitrates the time propagator, not spectral compilation.
      const long double rate =
        static_cast<long double>(diffusivity)
        * static_cast<long double>(model.A[d][static_cast<std::size_t>(mode)]);
      const long double forcing =
        static_cast<long double>(model.B[d][static_cast<std::size_t>(mode)])
        * static_cast<long double>(flux);
      const long double h = static_cast<long double>(dt);
      const long double z0 =
        static_cast<long double>(initial[d][static_cast<std::size_t>(mode)]);
      const long double expected_long = rate == 0.0L
                                          ? z0 + h * forcing
                                          : std::exp(rate * h) * z0
                                              + std::expm1(rate * h) / rate * forcing;
      const double expected = static_cast<double>(expected_long);
      const double actual = batch.state().at(batch.layout().spm.z[d], mode, 0);
      CAPTURE(NCH, d, mode, rate, expected, actual);
      REQUIRE(std::abs(actual - expected)
              <= 2e-12 * std::max(1.0, std::abs(expected)));
    }
  }
}

core::SpmFactoryInput thermalInput()
{
  auto input = test_support::make_legacy_kokam_input(0.55, 310.0, 298.0);
  input.design.thermal = { .density = 1000.0,
                           .heat_capacity = 1000.0,
                           .volume = 1e-3,
                           .surface_area = 1.0,
                           .h_conv = 10.0,
                           .reference_temperature = 298.0,
                           .environment_temperature = 300.0 };
  return input;
}

template <int NCH>
test_support::RecordedBits recordedCpuKernelTrace()
{
  core::ParameterSet parameters;
  REQUIRE(core::ParameterSet::chen2020(parameters) == Status::Success);
  core::SpmFactoryInput input;
  REQUIRE(parameters.toSpmInput(input) == Status::Success);
  input.initial_soc = 0.55;

  constexpr int lanes = 4;
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, { .nch = NCH }, lanes, batch)
          == Status::Success);
  auto &state = batch.state();
  const auto &layout = batch.layout().spm;
  constexpr std::array target_soc{ 0.31, 0.47, 0.66, 0.79 };
  constexpr std::array temperature{ 291.25, 298.15, 304.75, 313.5 };
  constexpr std::array diffusion_scale{ 0.91, 1.03, 0.97, 1.08 };
  constexpr std::array resistance_scale{ 1.07, 0.94, 1.02, 0.89 };
  for (int lane = 0; lane < lanes; ++lane) {
    state.at(layout.temperature, 0, lane) = temperature[static_cast<std::size_t>(lane)];
    state.at(layout.current_collector_resistance, 0, lane) *=
      resistance_scale[static_cast<std::size_t>(lane)];
    for (const core::Domain domain : core::domains) {
      const auto d = core::domain_index(domain);
      const auto &material = input.design.electrode[d].active_material;
      const double base = material.x_0
                          + input.initial_soc * (material.x_100 - material.x_0);
      const double varied = material.x_0
                            + target_soc[static_cast<std::size_t>(lane)]
                                * (material.x_100 - material.x_0);
      for (int mode = 0; mode < NCH; ++mode)
        state.at(layout.z[d], mode, lane) *= varied / base;
      state.at(layout.diffusion_coefficient[d], 0, lane) *=
        diffusion_scale[static_cast<std::size_t>(lane)];
    }
  }

  const double one_c_density = input.design.capacity_Ah
                               / input.design.electrode_area;
  const std::array density{ -0.45 * one_c_density,
                            0.0,
                            0.35 * one_c_density,
                            0.85 * one_c_density };
  std::array<double, lanes> voltage{};
  const core::StepCtx initial{ .time = 0.0, .dt = 0.0, .i_app = density };
  REQUIRE(batch.terminalVoltage(initial, voltage) == Status::Success);

  test_support::RecordedBits bits;
  bits.append(std::span<const double>{ state.raw() });
  bits.append(voltage);

  core::ExponentialModal stepper{ batch };
  double time{};
  for (const double dt : { 1e-5, 7.25, 19.0 }) {
    REQUIRE(stepper.step(batch, density, time, dt) == Status::Success);
    time += dt;
    bits.append(std::span<const double>{ state.raw() });
    bits.append(stepper.terminalVoltage());
  }
  return bits;
}

} // namespace

TEST_CASE("PC-10 CPU traces retain portable framing and capture-host bits",
          "[core][integrator][PC-10][recorded]")
{
  const std::array traces{ recordedCpuKernelTrace<5>(),
                           recordedCpuKernelTrace<8>(),
                           recordedCpuKernelTrace<12>() };
  constexpr std::array nch{ 5, 8, 12 };
  constexpr std::array<std::size_t, 3> expected_values{ 944, 1136, 1392 };
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#if defined(SLIDE_TEST_IPO) && defined(__FAST_MATH__)
  constexpr std::array expected_fnv{ UINT64_C(0x00a78a61e6bd84d6),
                                     UINT64_C(0xf4d3b8eb337172a6),
                                     UINT64_C(0x8dca38c943fc058a) };
  constexpr std::array expected_mixed{ UINT64_C(0x196a34e9b1d00094),
                                       UINT64_C(0x4e06b19d061124f7),
                                       UINT64_C(0xac4b6059db321d83) };
#elif defined(__FAST_MATH__)
  constexpr std::array expected_fnv{ UINT64_C(0xcd87b2262663a983),
                                     UINT64_C(0x5c766159af330d3e),
                                     UINT64_C(0x21b0119f4f1ec474) };
  constexpr std::array expected_mixed{ UINT64_C(0x7966eedf67f341b5),
                                       UINT64_C(0xad91f8f42bf550ce),
                                       UINT64_C(0x1feb3cf7a32e5b35) };
#else
  constexpr std::array expected_fnv{ UINT64_C(0xf461aa18b8f29e9d),
                                     UINT64_C(0x81d095bfc708d897),
                                     UINT64_C(0x983fa84ba68fc867) };
  constexpr std::array expected_mixed{ UINT64_C(0xe99152e8e9b27262),
                                       UINT64_C(0xde4eb9a14a355bc4),
                                       UINT64_C(0xdfc61b3e371e80f1) };
#endif
#endif
  for (std::size_t i = 0; i < traces.size(); ++i) {
    CAPTURE(nch[i], traces[i].values, traces[i].fnv1a, traces[i].mixed);
    REQUIRE(traces[i].values == expected_values[i]);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
    CHECK(traces[i].fnv1a == expected_fnv[i]);
    CHECK(traces[i].mixed == expected_mixed[i]);
#endif
  }
}

TEST_CASE("P3-G1 exponential modal update matches its independent closed form",
          "[core][integrator][exponential][P3-G1]")
{
  exactModalGate<5>();
  exactModalGate<8>();
  exactModalGate<12>();
}

TEST_CASE("P3-G2 a full signed-current cycle conserves the modal inventory",
          "[core][integrator][conservation][P3-G2]")
{
  constexpr int nch = 8;
  constexpr double current = 1.0;
  constexpr double dt = 10.0;
  constexpr int half_steps = 360;
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, { .nch = nch }, 1, batch) == Status::Success);
  core::PerDomain<double> radius{};
  for (const auto domain : core::domains)
    radius[core::domain_index(domain)] = input.design.electrode[core::domain_index(domain)].particle_radius;
  core::CompiledSpectralModel<nch> model;
  REQUIRE(core::compileSpectralModel<nch>(radius, model) == Status::Success);
  core::PerDomain<double> initial_mass_mode{};
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    initial_mass_mode[d] = batch.state().at(batch.layout().spm.z[d], model.zero_mode[d], 0);
  }

  core::ExponentialModal stepper{ batch };
  for (int step = 0; step < 2 * half_steps; ++step) {
    const double signed_current = step < half_steps ? current : -current;
    const std::array density{ signed_current / batch.electrode_area() };
    REQUIRE(stepper.step(batch, density, step * dt, dt) == Status::Success);
  }
  const double charge_in = current * half_steps * dt / 3600.0;
  const double charge_out = current * half_steps * dt / 3600.0;
  const auto neg = core::domain_index(core::Domain::neg);
  const auto &negative = input.design.electrode[neg];
  const double specific_area = 3.0 * negative.active_fraction
                               / negative.particle_radius;
  const double final_negative_mode = batch.state().at(
    batch.layout().spm.z[neg], model.zero_mode[neg], 0);
  const double stored_charge = (final_negative_mode - initial_mass_mode[neg])
                               * input.design.electrode_area * specific_area
                               * 96487.0 * negative.thickness
                               / (3600.0 * model.B[neg][static_cast<std::size_t>(model.zero_mode[neg])]
                                  * static_cast<double>(core::molar_flux_sign(
                                    core::Domain::neg)));
  REQUIRE(std::abs(charge_in - charge_out - stored_charge) < 1e-9);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    const double final_mode = batch.state().at(batch.layout().spm.z[d], model.zero_mode[d], 0);
    REQUIRE(std::abs(final_mode - initial_mass_mode[d]) <= 1e-9);
  }
}

TEST_CASE("P3-G3 nch12 exponential stepping survives an Euler-unstable step",
          "[core][integrator][stability][P3-G3]")
{
  constexpr int nch = 12;
  auto input = test_support::make_legacy_kokam_input(0.55, 298.0, 298.0);
  const core::SpmModelOptions options{ .nch = nch };
  auto exponential_batch =
    test_support::requireSpmBatch(input, options, 1);
  auto euler_batch = test_support::requireSpmBatch(input, options, 1);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    exponential_batch.state().at(exponential_batch.layout().spm.z[d], nch - 1, 0) += 1.0;
    euler_batch.state().at(euler_batch.layout().spm.z[d], nch - 1, 0) += 1.0;
  }
  const std::array density{ 0.0 };
  core::ExponentialModal exponential{ exponential_batch };
  core::EulerLegacy euler{ euler_batch };
  constexpr double unstable_dt = 1000.0;
  REQUIRE(exponential.step(exponential_batch, density, 0.0, unstable_dt)
          == Status::Success);
  const auto euler_status = euler.step(euler_batch, density, 0.0, unstable_dt);
  double exponential_norm{};
  double euler_norm{};
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    exponential_norm = std::max(exponential_norm,
                                std::abs(exponential_batch.state().at(
                                  exponential_batch.layout().spm.z[d], nch - 1, 0)));
    euler_norm = std::max(euler_norm,
                          std::abs(euler_batch.state().at(
                            euler_batch.layout().spm.z[d], nch - 1, 0)));
  }
  CAPTURE(euler_status, exponential_norm, euler_norm);
  REQUIRE(std::isfinite(exponential_norm));
  REQUIRE((euler_status != Status::Success || euler_norm > 1e6 * exponential_norm));
}

TEST_CASE("Strang slow split is second order and adaptive steps align to events",
          "[core][integrator][strang][adaptive][event]")
{
  const auto input = thermalInput();
  const core::SpmModelOptions options{ .nch = 5, .thermal = true };
  auto coarse = test_support::requireSpmBatch(input, options, 1);
  auto fine = test_support::requireSpmBatch(input, options, 1);
  auto adaptive = test_support::requireSpmBatch(input, options, 1);
  const std::array density{ 0.0 };
  core::ExponentialModal coarse_stepper{ coarse };
  core::ExponentialModal fine_stepper{ fine };
  constexpr double duration = 10.0;
  REQUIRE(coarse_stepper.step(coarse, density, 0.0, duration) == Status::Success);
  REQUIRE(fine_stepper.step(fine, density, 0.0, 0.5 * duration) == Status::Success);
  REQUIRE(fine_stepper.step(fine, density, 0.5 * duration, 0.5 * duration)
          == Status::Success);
  const double rate = input.design.thermal.h_conv * input.design.thermal.surface_area
                      / (input.design.thermal.density
                         * input.design.thermal.heat_capacity
                         * input.design.thermal.volume);
  const double exact = input.design.thermal.environment_temperature
                       + (input.initial_temperature
                          - input.design.thermal.environment_temperature)
                           * std::exp(-rate * duration);
  const double coarse_error = std::abs(
    coarse.state().at(coarse.layout().spm.temperature, 0, 0) - exact);
  const double fine_error = std::abs(
    fine.state().at(fine.layout().spm.temperature, 0, 0) - exact);
  CAPTURE(coarse_error, fine_error);
  REQUIRE(fine_error < 0.35 * coarse_error);

  core::ExponentialModal adaptive_stepper{ adaptive };
  double accepted{}, next{};
  REQUIRE(adaptive_stepper.stepAdaptive(adaptive, density, 0.0, 100.0, 1e-10, 1e-8, accepted, next)
          == Status::Success);
  REQUIRE(accepted > 0.0);
  REQUIRE(accepted <= 100.0);
  REQUIRE(next > 0.0);
  const double adaptive_exact = input.design.thermal.environment_temperature
                                + (input.initial_temperature
                                   - input.design.thermal.environment_temperature)
                                    * std::exp(-rate * accepted);
  REQUIRE(std::abs(adaptive.state().at(adaptive.layout().spm.temperature, 0, 0)
                   - adaptive_exact)
          <= 1e-6);
  REQUIRE(core::ExponentialModal::alignToEvent(3.0, 5.0, 4.25) == 1.25);
  REQUIRE(core::ExponentialModal::alignToEvent(3.0, 0.5, 4.25) == 0.5);
}

TEST_CASE("ExponentialModal validates inputs and restores state after rejected trials",
          "[core][integrator][exponential][validation]")
{
  core::ExponentialModal stepper;
  const core::SpmBatch empty;
  REQUIRE(stepper.configure(empty) == Status::Invalid_parameters);

  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(test_support::make_legacy_kokam_input(0.55, 298.0, 298.0),
                              {},
                              1,
                              batch)
          == Status::Success);
  REQUIRE(stepper.configure(batch) == Status::Success);
  constexpr std::array zero_current{ 0.0 };
  REQUIRE(stepper.step(batch, zero_current, 0.0, 0.0)
          == Status::Invalid_parameters);

  const double quiet_nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  const std::array invalid_current{ quiet_nan };
  REQUIRE(stepper.step(batch, invalid_current, 0.0, 1.0)
          == Status::Invalid_parameters);

  double accepted = 7.0;
  double next = 11.0;
  REQUIRE(stepper.stepAdaptive(batch, zero_current, 0.0, 1.0, 0.0, 1e-6, accepted, next)
          == Status::Invalid_parameters);
  REQUIRE(accepted == 7.0);
  REQUIRE(next == 11.0);

  const std::vector<double> before(batch.state().raw().begin(),
                                   batch.state().raw().end());
  constexpr std::array extreme_current{ 1e200 };
  REQUIRE(stepper.stepAdaptive(batch, extreme_current, 0.0, 1.0, 1e-12, 1e-12, accepted, next)
          == Status::Numerical_failure);
  REQUIRE(std::equal(before.begin(), before.end(), batch.state().raw().begin()));
  REQUIRE(accepted == 7.0);
  REQUIRE(next == 11.0);
}
