/**
 * @file core_CoreSpmTestHarness_test.cpp
 * @brief Independent forwarding, unit, framing, ownership, and metric gates for
 *        the M0.8 shared SPM test harness.
 */

#include "../support/CoreSpmTestHarness.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

using namespace slide;

namespace {

core::SpmFactoryInput makeHarnessInput()
{
  core::SpmFactoryInput input;
  input.design.capacity_Ah = 1.75;
  input.design.electrode_area = 0.25;
  input.design.electrolyte.concentration = 1000.0;
  input.design.thermal = { .density = 1626.0,
                           .heat_capacity = 750.0,
                           .volume = 1e-4,
                           .surface_area = 0.02,
                           .h_conv = 25.0,
                           .reference_temperature = 298.15,
                           .environment_temperature = 301.0 };
  for (const core::Domain domain : core::domains) {
    const auto d = core::domain_index(domain);
    auto &electrode = input.design.electrode[d];
    electrode.thickness = domain == core::Domain::neg ? 75e-6 : 87e-6;
    electrode.porosity = 0.3;
    electrode.active_fraction = 0.5;
    electrode.particle_radius =
      domain == core::Domain::neg ? 12.5e-6 : 8.5e-6;
    auto &material = electrode.active_material;
    material.ocv.stoichiometry = { 0.0, 1.0 };
    material.ocv.value = domain == core::Domain::neg
                           ? std::vector<double>{ 0.0, 0.1 }
                           : std::vector<double>{ 3.0, 4.0 };
    material.cs_max = domain == core::Domain::neg ? 30'555.0 : 51'385.0;
    material.x_0 = domain == core::Domain::neg ? 0.2 : 0.8;
    material.x_100 = domain == core::Domain::neg ? 0.8 : 0.2;
    material.D_s = {
      .reference_value = domain == core::Domain::neg ? 7e-14 : 8e-14,
      .activation_energy = domain == core::Domain::neg ? 7000.0 : 29'000.0,
      .reference_temperature = 298.15,
    };
    material.k_ct = {
      .reference_value = domain == core::Domain::neg ? 1.764e-11 : 5e-11,
      .activation_energy = domain == core::Domain::neg ? 20'000.0 : 58'000.0,
      .reference_temperature = 298.15,
    };
    input.initial_specific_resistance[d] = 2.8e-3;
  }
  input.initial_soc = 0.61;
  input.initial_temperature = 304.0;
  input.initial_current_collector_resistance = 0.2325e-3;
  return input;
}

constexpr core::SpmModelOptions harness_options{ .nch = 8, .thermal = true };
constexpr int harness_lanes = 3;
constexpr core::SpmModelOptions alternate_harness_options{
  .nch = 12, .thermal = false
};
constexpr int alternate_harness_lanes = 2;

core::real_t runtimeRealFromBits(std::uint64_t bits) noexcept
{
  volatile std::uint64_t runtime_bits = bits;
  const std::uint64_t observed_bits = runtime_bits;
  return std::bit_cast<core::real_t>(observed_bits);
}

} // namespace

TEST_CASE("Core SPM harness forwards explicit build and distinguishes current units",
          "[core][test-harness][9C-4][units]")
{
  const auto input = makeHarnessInput();
  auto batch =
    test_support::requireSpmBatch(input, harness_options, harness_lanes);

  REQUIRE(batch.valid());
  CHECK(batch.nch() == 8);
  CHECK(batch.n_lanes() == harness_lanes);
  CHECK(batch.composition() == core::SpmComposition::thermal);
  CHECK(batch.electrode_area() == 0.25);

  const std::vector state_before(
    batch.state().raw().begin(), batch.state().raw().end());
  constexpr std::array<core::real_t, harness_lanes> current_A{
    0.25, -0.5, 0.125
  };
  constexpr std::array<core::real_t, harness_lanes> expected_density_Apm2{
    1.0, -2.0, 0.5
  };
  std::array<core::real_t, harness_lanes> density_scratch_Apm2{
    -7.0, -7.0, -7.0
  };
  std::array<core::real_t, harness_lanes> voltage_from_A_V{};
  std::array<core::real_t, harness_lanes> voltage_from_density_V{};

  test_support::requireTerminalVoltage(
    batch,
    test_support::CurrentA{ std::span<const core::real_t>{ current_A } },
    0.0,
    density_scratch_Apm2,
    voltage_from_A_V);
  CHECK(density_scratch_Apm2 == expected_density_Apm2);

  test_support::requireTerminalVoltage(
    batch,
    test_support::CurrentDensityApm2{
      std::span<const core::real_t>{ expected_density_Apm2 } },
    0.0,
    voltage_from_density_V);
  CHECK(voltage_from_A_V == voltage_from_density_V);
  CHECK(std::ranges::equal(batch.state().raw(), state_before));
}

TEST_CASE("Core SPM harness forwards a second independent build and area tuple",
          "[core][test-harness][9C-4][units][alternate]")
{
  auto input = makeHarnessInput();
  input.design.electrode_area = 0.4;
  auto batch = test_support::requireSpmBatch(
    input, alternate_harness_options, alternate_harness_lanes);

  REQUIRE(batch.valid());
  CHECK(batch.nch() == 12);
  CHECK(batch.n_lanes() == alternate_harness_lanes);
  CHECK(batch.composition() == core::SpmComposition::isothermal);
  CHECK(batch.electrode_area() == 0.4);

  constexpr std::array<core::real_t, alternate_harness_lanes> current_A{
    0.4, -0.2
  };
  constexpr std::array<core::real_t, alternate_harness_lanes>
    expected_density_Apm2{ 1.0, -0.5 };
  std::array<core::real_t, alternate_harness_lanes> density_scratch_Apm2{};
  std::array<core::real_t, alternate_harness_lanes> terminal_voltage_V{};
  test_support::requireTerminalVoltage(
    batch,
    test_support::CurrentA{ std::span<const core::real_t>{ current_A } },
    1.25,
    density_scratch_Apm2,
    terminal_voltage_V);
  CHECK(density_scratch_Apm2 == expected_density_Apm2);
}

TEST_CASE("Core SPM harness uses total pointer ordering for disjoint buffers",
          "[core][test-harness][9C-4][ownership]")
{
  std::array<core::real_t, 6> storage{};
  std::array<core::real_t, 2> separate_storage{};
  const std::span<const core::real_t> whole{ storage };
  const auto first = whole.first(3);
  const auto overlapping = whole.subspan(2, 3);
  const auto adjacent = whole.last(3);

  CHECK(test_support::detail::spansOverlap(first, overlapping));
  CHECK(test_support::detail::spansOverlap(overlapping, first));
  CHECK_FALSE(test_support::detail::spansOverlap(first, adjacent));
  CHECK_FALSE(test_support::detail::spansOverlap(first, separate_storage));
  CHECK_FALSE(test_support::detail::spansOverlap(
    std::span<const core::real_t>{}, first));

  const core::real_t maximum =
    runtimeRealFromBits(UINT64_C(0x7fefffffffffffff));
  const core::real_t minimum_subnormal =
    runtimeRealFromBits(UINT64_C(0x0000000000000001));
  CHECK_FALSE(test_support::detail::canDivideByPositiveNormal(maximum, 0.5));
  CHECK(test_support::detail::canDivideByPositiveNormal(0.25, 0.25));
  CHECK_FALSE(
    test_support::detail::canAddSquare(maximum, minimum_subnormal));
}

TEST_CASE("Core SPM harness owns every nonuniform trace frame and interval",
          "[core][test-harness][9C-4][trace]")
{
  const auto input = makeHarnessInput();
  auto batch =
    test_support::requireSpmBatch(input, harness_options, harness_lanes);
  core::SpmBatch manual_batch;
  REQUIRE(core::buildSpmBatch(
            input, harness_options, harness_lanes, manual_batch)
          == Status::Success);

  constexpr std::array<core::real_t, harness_lanes> current_density_Apm2{
    0.4, -0.2, 0.1
  };
  constexpr std::array<core::real_t, 5> sample_time_s{
    2.0, 2.125, 2.5, 2.75, 2.875
  };
  std::array<core::real_t, sample_time_s.size() * harness_lanes>
    expected_voltage_V{};

  const core::StepCtx initial_observation{
    .time = sample_time_s.front(),
    .dt = 0.0,
    .i_app = current_density_Apm2,
  };
  REQUIRE(manual_batch.terminalVoltage(
            initial_observation,
            std::span<core::real_t>{ expected_voltage_V }.first(
              harness_lanes))
          == Status::Success);
  core::ExponentialModal manual_stepper;
  REQUIRE(manual_stepper.configure(manual_batch) == Status::Success);
  for (std::size_t sample = 1; sample < sample_time_s.size(); ++sample) {
    const core::real_t start_time_s = sample_time_s[sample - 1];
    const core::real_t dt_s = sample_time_s[sample] - start_time_s;
    REQUIRE(dt_s > 0.0);
    REQUIRE(manual_stepper.step(
              manual_batch, current_density_Apm2, start_time_s, dt_s)
            == Status::Success);
    std::ranges::copy(
      manual_stepper.terminalVoltage(),
      expected_voltage_V.begin()
        + static_cast<std::ptrdiff_t>(sample * harness_lanes));
  }

  std::array<core::real_t, sample_time_s.size() * harness_lanes>
    actual_voltage_V{};
  std::ranges::fill(actual_voltage_V, -777.0);
  const auto *const state_storage = batch.state().raw().data();
  const auto *const output_storage = actual_voltage_V.data();
  core::ExponentialModal stepper;
  test_support::requireConstantCurrentTrace(
    batch,
    stepper,
    test_support::CurrentDensityApm2{
      std::span<const core::real_t>{ current_density_Apm2 } },
    sample_time_s,
    actual_voltage_V);

  CHECK(batch.state().raw().data() == state_storage);
  CHECK(actual_voltage_V.data() == output_storage);
  CHECK(actual_voltage_V == expected_voltage_V);
  CHECK(std::ranges::equal(batch.state().raw(), manual_batch.state().raw()));
  CHECK(std::ranges::equal(batch.derivative().raw(),
                           manual_batch.derivative().raw()));
  for (int lane = 0; lane < harness_lanes; ++lane)
    CHECK(batch.state().at(batch.layout().elapsed_time, 0, lane) == 0.875);
}

TEST_CASE("Core SPM harness voltage metrics are literal finite and end-sensitive",
          "[core][test-harness][9C-4][metric]")
{
  constexpr std::array<core::real_t, 3> zero_voltage_V{ 0.0, 0.0, 0.0 };
  constexpr std::array<core::real_t, 3> literal_error_V{ 1.0, 2.0, -1.0 };
  test_support::VoltageError error{ -7.0, -11.0 };
  REQUIRE(test_support::computeVoltageError(
    literal_error_V, zero_voltage_V, error));
  CHECK(error.maximum_absolute_V == 2.0);
  CHECK(error.rms_V == std::sqrt(2.0));

  constexpr std::array<core::real_t, 3> final_maximum_error_V{
    1.0, -1.0, 2.0
  };
  REQUIRE(test_support::computeVoltageError(
    final_maximum_error_V, zero_voltage_V, error));
  CHECK(error.maximum_absolute_V == 2.0);
  CHECK(error.rms_V == std::sqrt(2.0));

  const core::real_t nan =
    runtimeRealFromBits(UINT64_C(0x7ff8000000000000));
  const core::real_t infinity =
    runtimeRealFromBits(UINT64_C(0x7ff0000000000000));
  const std::array<core::real_t, 3> nan_voltage_V{ 1.0, nan, 2.0 };
  const std::array<core::real_t, 3> infinite_voltage_V{
    1.0, 2.0, infinity
  };
  const test_support::VoltageError sentinel{ 7.0, 11.0 };
  error = sentinel;
  CHECK_FALSE(test_support::computeVoltageError(
    nan_voltage_V, zero_voltage_V, error));
  CHECK(error.maximum_absolute_V == sentinel.maximum_absolute_V);
  CHECK(error.rms_V == sentinel.rms_V);
  CHECK_FALSE(test_support::computeVoltageError(
    zero_voltage_V, infinite_voltage_V, error));
  CHECK(error.maximum_absolute_V == sentinel.maximum_absolute_V);
  CHECK(error.rms_V == sentinel.rms_V);

  CHECK_FALSE(test_support::computeVoltageError(
    std::span<const core::real_t>{ literal_error_V }.first(2),
    zero_voltage_V,
    error));
  CHECK_FALSE(test_support::computeVoltageError(
    std::span<const core::real_t>{},
    std::span<const core::real_t>{},
    error));
}

TEST_CASE("Core SPM metric rejects subtraction overflow before evaluation",
          "[core][test-harness][9C-4][metric][fast-math][subtraction]")
{
  const std::array<core::real_t, 1> positive_V{
    runtimeRealFromBits(UINT64_C(0x7fefffffffffffff))
  };
  const std::array<core::real_t, 1> negative_V{
    runtimeRealFromBits(UINT64_C(0xffefffffffffffff))
  };
  const test_support::VoltageError sentinel{ 7.0, 11.0 };
  auto output = sentinel;
  CHECK_FALSE(
    test_support::computeVoltageError(positive_V, negative_V, output));
  CHECK(output.maximum_absolute_V == sentinel.maximum_absolute_V);
  CHECK(output.rms_V == sentinel.rms_V);
}

TEST_CASE("Core SPM metric rejects a finite difference whose square overflows",
          "[core][test-harness][9C-4][metric][fast-math][square]")
{
  const core::real_t two_to_512 =
    runtimeRealFromBits(UINT64_C(0x5ff0000000000000));
  REQUIRE(core::is_finite(two_to_512));
  CHECK_FALSE(test_support::detail::canSquare(two_to_512));
  const std::array<core::real_t, 1> actual_V{ two_to_512 };
  const std::array<core::real_t, 1> expected_V{ 0.0 };
  const test_support::VoltageError sentinel{ 7.0, 11.0 };
  auto output = sentinel;
  CHECK_FALSE(
    test_support::computeVoltageError(actual_V, expected_V, output));
  CHECK(output.maximum_absolute_V == sentinel.maximum_absolute_V);
  CHECK(output.rms_V == sentinel.rms_V);
}

TEST_CASE("Core SPM metric rejects overflow of individually finite squares",
          "[core][test-harness][9C-4][metric][fast-math][sum]")
{
  const core::real_t two_to_511 =
    runtimeRealFromBits(UINT64_C(0x5fe0000000000000));
  core::real_t witness_square{};
  REQUIRE(core::try_multiply_nonnegative(
    two_to_511, two_to_511, witness_square));
  REQUIRE(core::is_finite(witness_square));
  const std::array<core::real_t, 4> actual_V{
    two_to_511, two_to_511, two_to_511, two_to_511
  };
  const std::array<core::real_t, 4> expected_V{};
  const test_support::VoltageError sentinel{ 7.0, 11.0 };
  auto output = sentinel;
  CHECK_FALSE(
    test_support::computeVoltageError(actual_V, expected_V, output));
  CHECK(output.maximum_absolute_V == sentinel.maximum_absolute_V);
  CHECK(output.rms_V == sentinel.rms_V);
}
