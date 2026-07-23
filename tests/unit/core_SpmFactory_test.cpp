/**
 * @file core_SpmFactory_test.cpp
 * @brief D-02 registry selection, cold validation, initialization, and batch dispatch tests.
 */

#include "../../src/core/SpmFactory.hpp"

#include <catch2/catch_test_macros.hpp>

#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#include <boost/hash2/sha2.hpp>
#endif

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <span>
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#include <string>
#endif
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

core::SpmFactoryInput make_factory_refactor_oracle_input()
{
  auto input = make_input();
  input.initial_soc = 0.53;
  input.initial_temperature = 305.25;
  input.initial_sei_thickness = 1.2e-9;
  input.initial_lost_lithium = 2.3e-8;
  input.initial_crack_surface_fraction = 0.012;
  input.initial_plated_lithium_thickness = 3.4e-10;
  input.initial_stress_interval = 30.0;
  input.sei_resistivity_area = 2037.4;

  // These valid sentinels deliberately disagree with the electrical/design
  // owners. The factory must replace every one while compiling each mechanism.
  input.sei.F = 11'111.0;
  input.sei.Rg = 7.11;
  input.sei.n = 0.71;
  input.sei.reference_temperature = 289.0;
  input.sei.electrode_area = 0.81;
  input.sei.negative_particle_radius = 8.1e-6;
  input.sei.sei_resistivity_area = 17.0;

  input.surface_crack.F = 22'222.0;
  input.surface_crack.Rg = 7.22;
  input.surface_crack.reference_temperature = 290.0;
  input.surface_crack.electrode_area = 0.92;
  input.surface_crack.negative_cs_max = 22'222.0;
  input.surface_crack.sei_resistivity_area = 18.0;

  input.lam.F = 33'333.0;
  input.lam.Rg = 7.33;
  input.lam.n = 1.33;
  input.lam.reference_temperature = 291.0;
  input.lam.particle_radius = { 9.1e-6, 9.2e-6 };

  input.lithium_plating.F = 44'444.0;
  input.lithium_plating.Rg = 7.44;
  input.lithium_plating.n = 1.44;
  input.lithium_plating.reference_temperature = 292.0;
  input.lithium_plating.electrode_area = 1.04;
  input.lithium_plating.sei_resistivity_area = 19.0;

  // Proven nonzero coefficients from the all-mask AgeingKernel fixture keep
  // every LAM alternative observable in the post-step digest.
  for (const core::Domain domain : core::domains) {
    const auto d = core::domain_index(domain);
    input.lam.model1_stress_coefficient[d] =
      domain == core::Domain::neg ? 1.1e-24 : 0.8e-24;
    input.lam.model2_linear_flux[d] =
      domain == core::Domain::neg ? -2.0e-7 : -1.5e-7;
    input.lam.model2_sqrt_flux[d] =
      domain == core::Domain::neg ? -3.0e-11 : -2.0e-11;
    input.lam.model4_area_coefficient[d] =
      domain == core::Domain::neg ? 2.0e-10 : 1.0e-10;
  }
  input.lam.model2_activation = 9000.0;
  input.lam.model3_k = 2.5e-12;
  input.lam.model3_k_activation = 5000.0;
  input.lam.model3_equilibrium_potential = 4.1;
  return input;
}

core::SpmModelOptions factory_refactor_oracle_options()
{
  return { .nch = 12,
           .thermal = true,
           .sei_model_mask = core::ageing_model_mask<4>(),
           .sei_porosity = true,
           .surface_crack_model_mask = core::ageing_model_mask<5>(),
           .surface_crack_diffusivity = true,
           .lam_model_mask = core::ageing_model_mask<4>(),
           .lithium_plating = true };
}

void require_positive_zero_padding(core::StateArena &arena)
{
  constexpr auto positive_zero_bits = std::bit_cast<std::uint64_t>(0.0);
  for (int row = 0; row < arena.n_rows(); ++row) {
    const auto padded = arena.row_padded(row);
    for (int lane = arena.n_lanes(); lane < arena.stride(); ++lane) {
      const auto bits =
        std::bit_cast<std::uint64_t>(padded[static_cast<std::size_t>(lane)]);
      CAPTURE(row, lane, bits);
      REQUIRE(bits == positive_zero_bits);
    }
  }
}

#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
std::string sha256(std::span<const double> values)
{
  boost::hash2::sha2_256 hash;
  hash.update(values.data(), values.size_bytes());
  return boost::hash2::to_string(hash.result());
}
#endif

} // namespace

TEST_CASE("MQ.2 factory refactor retains the full padded all-mechanism arenas",
          "[core][factory][MQ.2][recorded]")
{
  constexpr int lanes = 9;
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(make_factory_refactor_oracle_input(),
                              factory_refactor_oracle_options(),
                              lanes,
                              batch)
          == Status::Success);

#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
  // Capture before any explicit evaluation: these are the exact arenas returned
  // by the factory, including every padding lane.
  const auto initial_state_sha256 = sha256(batch.state().raw());
  const auto initial_derivative_sha256 = sha256(batch.derivative().raw());
#endif

  REQUIRE(batch.state().n_rows() == 50);
  REQUIRE(batch.state().n_lanes() == lanes);
  REQUIRE(batch.state().stride() == 16);
  REQUIRE(batch.state().raw().size() == 800);
  REQUIRE(batch.derivative().n_rows() == 50);
  REQUIRE(batch.derivative().n_lanes() == lanes);
  REQUIRE(batch.derivative().stride() == 16);
  REQUIRE(batch.derivative().raw().size() == 800);
  require_positive_zero_padding(batch.state());
  require_positive_zero_padding(batch.derivative());

  constexpr std::array currents{
    -2.25, -1.60, -0.95, -0.35, 0.15, 0.55, 1.05, 1.55, 2.15
  };
  const core::StepCtx ctx{ .time = 2.5, .dt = 1e-3, .i_app = currents };
  std::array<double, lanes> voltage{};
  REQUIRE(batch.exponentialStep(ctx, 1e-3, voltage) == Status::Success);
  REQUIRE(std::all_of(voltage.begin(), voltage.end(), [](double value) {
    return core::is_finite(value);
  }));

#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
  const auto stepped_state_sha256 = sha256(batch.state().raw());
  CAPTURE(initial_state_sha256, initial_derivative_sha256, stepped_state_sha256);

  // Each branch freezes pre/post identity within one compiler/build
  // configuration. Different branches are not a claim of cross-mode equality.
#if defined(SLIDE_TEST_RELEASE) && defined(SLIDE_TEST_IPO)
  constexpr auto expected_initial_state_sha256 =
    "a1c372842363e6109d609a217067b51e6af24a0a77480236f58f99c8d76d7955";
  constexpr auto expected_initial_derivative_sha256 =
    "3f5de3d1146e51b5d8d5755130753f329944a7c9cada50e5909aa348e76ed178";
  constexpr auto expected_stepped_state_sha256 =
    "5f72e3a5fbc7ff60d7c63d21f7ee7b11a9299b1aaba7930c3c0ac64b1af53705";
#elif defined(SLIDE_TEST_RELEASE)
  constexpr auto expected_initial_state_sha256 =
    "fec245fb6f6751b570f02edbcb00c2f84394b7c28e5652e5d151c464f9461f18";
  constexpr auto expected_initial_derivative_sha256 =
    "ad6f936b1087ff2a8647c6f1f69b5ef2e03ac374f7ce9c540f03f1b4095f7a08";
  constexpr auto expected_stepped_state_sha256 =
    "119011a365f6bb438e63d13aafa3faafcaa3f8bcf65cae0df153a161182b94f4";
#else
  constexpr auto expected_initial_state_sha256 =
    "f9c0df074580cad37eab525a428eb75641626e1398d3f3fb71eeffe7c10278fa";
  constexpr auto expected_initial_derivative_sha256 =
    "d7746f27ecd3de2efccc65458bc429421b80d45f3984bb2fe4f28d3d2c01d476";
  constexpr auto expected_stepped_state_sha256 =
    "d5de6e0fb2fc2b550728922066c349a3678032421cc3a7962ecca2d412896110";
#endif
  CHECK(initial_state_sha256 == expected_initial_state_sha256);
  CHECK(initial_derivative_sha256 == expected_initial_derivative_sha256);
  CHECK(stepped_state_sha256 == expected_stepped_state_sha256);
#endif

  require_positive_zero_padding(batch.state());
  require_positive_zero_padding(batch.derivative());
}

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

TEST_CASE("SPM batch validates every public dispatch boundary",
          "[core][factory][validation]")
{
  const auto input = make_input();
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, {}, 2, batch) == Status::Success);
  constexpr std::array current{ 0.0, 0.0 };
  const core::StepCtx ctx{ .dt = 1.0, .i_app = current };
  std::array<double, 2> voltage{};

  REQUIRE(batch.fusedEuler(ctx, 0.0, voltage) == Status::Invalid_parameters);
  REQUIRE(batch.exponentialStep(ctx, 0.0, voltage)
          == Status::Invalid_parameters);
  REQUIRE(batch.terminalVoltage(ctx, std::span<double>{ voltage }.first(1))
          == Status::Invalid_parameters);
  REQUIRE(batch.terminalVoltageAt(std::span<const double>{ batch.state().raw() }.first(
                                    batch.state().size() - 1),
                                  current,
                                  voltage)
          == Status::Invalid_parameters);
  REQUIRE(batch.linearizeThevenin(current, std::span<double>{ voltage }.first(1), voltage)
          == Status::Invalid_parameters);
  REQUIRE(batch.setTrustedLanePeriod(0) == Status::Invalid_parameters);
}

TEST_CASE("SPM factory classifies each cold-input validation layer",
          "[core][factory][validation]")
{
  core::SpmBatch output;

  auto input = make_input();
  input.design.capacity_Ah = 0.0;
  REQUIRE(core::buildSpmBatch(input, {}, 1, output)
          == Status::Invalid_parameters);

  input = make_input();
  input.design.electrode[core::domain_index(core::Domain::neg)].thickness = 0.0;
  REQUIRE(core::buildSpmBatch(input, {}, 1, output)
          == Status::Invalid_parameters);

  input = make_input();
  input.design.electrode[core::domain_index(core::Domain::neg)]
    .stress.poisson_ratio = 1.0;
  core::SpmModelOptions stress_options;
  stress_options.lam_model_mask = core::lam_model_bit(1);
  REQUIRE(core::buildSpmBatch(input, stress_options, 1, output)
          == Status::Invalid_parameters);

  input = make_input();
  core::SpmModelOptions unregistered;
  unregistered.nch = 7;
  REQUIRE(core::buildSpmBatch(input, unregistered, 1, output)
          == Status::Invalid_parameters);
}

TEST_CASE("SPM batch rejects false lane periods and non-passive linearizations",
          "[core][factory][validation]")
{
  auto input = make_input();
  core::SpmBatch periodic;
  REQUIRE(core::buildSpmBatch(input, {}, 4, periodic) == Status::Success);
  periodic.state().at(periodic.layout().spm.temperature, 0, 2) += 1.0;
  REQUIRE(periodic.setTrustedLanePeriod(2) == Status::Invalid_states);

  for (const core::Domain domain : core::domains) {
    auto &curve = input.design.electrode[core::domain_index(domain)]
                    .active_material.ocv;
    curve.value = { 0.0, 1e100 };
  }
  core::SpmBatch non_passive;
  REQUIRE(core::buildSpmBatch(input, {}, 1, non_passive) == Status::Success);
  constexpr std::array current{ 0.0 };
  std::array<double, 1> intercept{}, resistance{};
  REQUIRE(non_passive.linearizeThevenin(current, intercept, resistance)
          == Status::Invalid_states);
}
