/**
 * @file core_SpmPipeline_test.cpp
 * @brief Fixed-pipeline contract: mandatory zeroing, shared diffusion inputs, and rebind safety.
 */

#include "../../src/core/SpmPipeline.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

using namespace slide;

// Force every optional branch through compilation; runtime cases below focus on the pipeline
// contract independently of individual mechanism parity tests.
template class slide::core::SpmPipeline<1, true, true, true, true, true>;

namespace {

constexpr int direct_euler_nch = 1;
constexpr int direct_euler_lanes = 2;
using DirectEulerPipeline =
  core::SpmPipeline<direct_euler_nch, false, false, false, false, false>;

struct DirectEulerFixture
{
  core::SpmBatchLayout layout;
  core::StateArena state;
  DirectEulerPipeline pipeline;
};

double runtimeRealFromBits(std::uint64_t bits)
{
  volatile std::uint64_t opaque_bits = bits;
  const std::uint64_t copied_bits = opaque_bits;
  return std::bit_cast<double>(copied_bits);
}

DirectEulerFixture makeDirectEulerFixture(int lanes)
{
  core::BatchBuilder builder;
  const auto layout = DirectEulerPipeline::declareLayout(builder);
  auto state = builder.build(lanes);

  core::SpmPipelineParams<direct_euler_nch> params;
  params.electrical.F = 1.0;
  params.electrical.Rg = 1.0;
  params.electrical.n = 1.0;
  params.electrical.electrolyte_concentration = 1.0;
  params.electrical.reference_temperature = 1.0;
  params.electrical.electrode_area = 1.0;
  params.electrical.concentration.F = 1.0;
  params.electrical.concentration.Rg = 1.0;
  params.electrical.concentration.n = 1.0;
  params.electrical.concentration.T_ref = 1.0;

  constexpr std::array curve_x{ 0.0, 1.0 };
  constexpr std::array zero_curve{ 0.0, 0.0 };
  constexpr std::array positive_ocv{ 3.0, 4.0 };
  REQUIRE(params.electrical.total_entropic_coefficient.build(curve_x, zero_curve)
          == Status::Success);
  REQUIRE(params.electrical.negative_entropic_coefficient.build(curve_x, zero_curve)
          == Status::Success);

  for (const core::Domain domain : core::domains) {
    const auto d = core::domain_index(domain);
    params.electrical.concentration.C[d][0][0] = 1.0;
    params.electrical.concentration.C[d][1][0] = 1.0;
    params.electrical.concentration.R[d] = 0.0;
    params.electrical.electrode[d].cs_max = 1.0;
    params.electrical.electrode[d].reaction_rate_ref = 1.0;
    const auto &ocv = domain == core::Domain::neg ? zero_curve : positive_ocv;
    REQUIRE(params.electrical.electrode_ocv[d].build(curve_x, ocv)
            == Status::Success);

    for (int lane = 0; lane < lanes; ++lane) {
      state.at(layout.spm.z[d], 0, lane) = 0.5;
      state.at(layout.spm.diffusion_coefficient[d], 0, lane) = 2.0;
      state.at(layout.spm.specific_surface_area[d], 0, lane) = 4.0;
      state.at(layout.spm.electrode_thickness[d], 0, lane) = 5.0;
    }
  }
  const auto neg = core::domain_index(core::Domain::neg);
  const auto pos = core::domain_index(core::Domain::pos);
  params.diffusion.A[neg][0] = -3.0;
  params.diffusion.A[pos][0] = -4.0;
  params.diffusion.B[neg][0] = 2.0;
  params.diffusion.B[pos][0] = 3.0;
  for (int lane = 0; lane < lanes; ++lane)
    state.at(layout.spm.temperature, 0, lane) = 1.0;

  return { layout,
           std::move(state),
           DirectEulerPipeline{ std::move(params), layout, lanes } };
}

} // namespace

TEST_CASE("SPM pipeline zeroes and evaluates the rebound trial vector", "[core][pipeline]")
{
  constexpr int lanes = 1;
  auto fixture = makeDirectEulerFixture(lanes);
  const auto &layout = fixture.layout;
  auto &arena = fixture.state;
  auto &pipeline = fixture.pipeline;
  core::StateArena derivative{ arena.n_rows(), lanes };
  const auto neg = core::domain_index(core::Domain::neg);
  const auto pos = core::domain_index(core::Domain::pos);
  core::RhsViews views{ core::BatchShape::from(arena) };
  std::fill(derivative.raw().begin(), derivative.raw().end(), 9.0);
  views.rebind(arena.raw(), derivative.raw());
  constexpr std::array current_density{ 2.0 };
  const core::StepCtx ctx{ .i_app = current_density };

  REQUIRE(pipeline.evaluate(views, ctx) == Status::Success);
  REQUIRE(std::abs(derivative.at(layout.spm.z[neg], 0, 0) - (-2.8)) <= 1e-13);
  REQUIRE(std::abs(derivative.at(layout.spm.z[pos], 0, 0) - (-4.3)) <= 1e-13);
  REQUIRE(derivative.at(layout.spm.temperature, 0, 0) == 0.0);
  REQUIRE(std::count(derivative.raw().begin(), derivative.raw().end(), 9.0) == 0);

  // Adaptive steppers own their trial vectors. Rebinding must affect every RHS stage.
  std::vector<core::real_t> trial(arena.raw().begin(), arena.raw().end());
  trial[static_cast<std::size_t>(layout.spm.z[neg].row_begin) * arena.stride()] = 0.25;
  trial[static_cast<std::size_t>(layout.spm.specific_surface_area[neg].row_begin)
        * arena.stride()] = 8.0;
  std::fill(derivative.raw().begin(), derivative.raw().end(), -7.0);
  views.rebind(std::span<const core::real_t>{ trial }, derivative.raw());

  REQUIRE(pipeline.evaluate(views, ctx) == Status::Success);
  REQUIRE(std::abs(derivative.at(layout.spm.z[neg], 0, 0) - (-1.4)) <= 1e-13);
  REQUIRE(derivative.at(layout.spm.temperature, 0, 0) == 0.0);
  REQUIRE(std::count(derivative.raw().begin(), derivative.raw().end(), -7.0) == 0);

  core::BatchView mutable_state{ core::BatchShape::from(arena), arena.raw() };
  core::BatchView mutable_derivative{ core::BatchShape::from(derivative),
                                      derivative.raw() };
  const core::ConstBatchView const_state{ core::BatchShape::from(arena),
                                          arena.raw() };
  std::array<double, 1> resistance{}, intercept{};
  REQUIRE(pipeline.advanceExponential(mutable_state, mutable_derivative, ctx, 1.0, {})
          == Status::Invalid_parameters);
  REQUIRE(pipeline.observeTerminalVoltage(const_state, ctx, {})
          == Status::Invalid_parameters);
  REQUIRE(pipeline.linearizeThevenin(const_state, current_density, {}, resistance)
          == Status::Invalid_parameters);

  const double quiet_nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  const std::array invalid_current{ quiet_nan };
  REQUIRE(pipeline.linearizeThevenin(const_state, invalid_current, intercept, resistance)
          == Status::Invalid_parameters);
  REQUIRE(pipeline.setTrustedLanePeriod(const_state, 0)
          == Status::Invalid_parameters);
}

TEST_CASE("SPM direct Euler validates every input before mutation",
          "[core][pipeline][euler][validation]")
{
  auto fixture = makeDirectEulerFixture(direct_euler_lanes);
  const int expected_rows = fixture.layout.energy_throughput.row_begin
                            + fixture.layout.energy_throughput.rows;
  REQUIRE(fixture.state.n_rows() == expected_rows);
  const core::BatchShape valid_shape = core::BatchShape::from(fixture.state);
  const std::vector<double> state_before(fixture.state.raw().begin(),
                                         fixture.state.raw().end());
  constexpr std::array terminal_before{ 101.0, -202.0, 303.0 };
  std::array terminal_voltage = terminal_before;
  constexpr std::array valid_current{ 0.0, 0.0 };
  constexpr double valid_time = 0.5;
  constexpr double valid_context_dt = 0.25;
  constexpr double valid_explicit_dt = 0.001;
  const double quiet_nan = runtimeRealFromBits(UINT64_C(0x7ff8000000000000));
  const double positive_infinity =
    runtimeRealFromBits(UINT64_C(0x7ff0000000000000));
  const double negative_infinity =
    runtimeRealFromBits(UINT64_C(0xfff0000000000000));

  // Preregistered exact matrix for a two-lane otherwise-valid direct pipeline:
  // - shape: final energy-throughput end minus/plus one row, or one/three lanes;
  // - current/output sizes: empty (0), short (1), or long (3), versus exact size 2;
  // - current values: opaque qNaN 0x7ff8000000000000, +Inf 0x7ff0000000000000,
  //   or -Inf 0xfff0000000000000 in an otherwise-zero exact-size span;
  // - explicit dt: 0, -0.25, opaque qNaN, or +Inf;
  // - ctx.dt: opaque qNaN or +Inf; ctx.time: opaque qNaN or +Inf.
  // Every invalid call must return Invalid_parameters with the entire padded state arena
  // and all three terminal sentinel values byte-identical. The exact-shape/size, finite
  // control uses current {0,0}, ctx {time=0.5, dt=0.25}, and explicit dt=0.001.
  const auto require_rejected =
    [&](std::string_view label,
        core::BatchShape shape,
        std::span<const double>
          current,
        double time,
        double context_dt,
        double explicit_dt,
        std::size_t output_size) {
      INFO(label);
      REQUIRE(output_size <= terminal_voltage.size());
      std::vector<double> candidate_state(shape.storage_size(), -404.0);
      std::copy_n(state_before.begin(),
                  std::min(state_before.size(), candidate_state.size()),
                  candidate_state.begin());
      const std::vector<double> candidate_before = candidate_state;
      terminal_voltage = terminal_before;
      core::BatchView state{ shape, candidate_state };
      const core::StepCtx ctx{ .time = time,
                               .dt = context_dt,
                               .i_app = current };

      CHECK(fixture.pipeline.advanceEuler(
              state,
              ctx,
              explicit_dt,
              std::span{ terminal_voltage }.first(output_size))
            == Status::Invalid_parameters);
      CHECK(std::memcmp(candidate_state.data(),
                        candidate_before.data(),
                        candidate_before.size() * sizeof(double))
            == 0);
      CHECK(std::memcmp(terminal_voltage.data(),
                        terminal_before.data(),
                        sizeof(terminal_voltage))
            == 0);
    };

  struct RejectionCase
  {
    std::string_view label;
    int row_delta{};
    int lane_delta{};
    std::array<double, 3> current{};
    std::size_t current_size{ direct_euler_lanes };
    std::size_t output_size{ direct_euler_lanes };
    double time{ 0.5 };
    double context_dt{ 0.25 };
    double explicit_dt{ 0.001 };
  };
  const std::array rejection_cases{
    RejectionCase{ .label = "wrong state row count", .row_delta = -1 },
    RejectionCase{ .label = "oversized state row count", .row_delta = 1 },
    RejectionCase{ .label = "wrong state lane count", .lane_delta = -1 },
    RejectionCase{ .label = "oversized state lane count", .lane_delta = 1 },
    RejectionCase{ .label = "empty current span", .current_size = 0 },
    RejectionCase{ .label = "short current span", .current_size = 1 },
    RejectionCase{ .label = "long current span", .current_size = 3 },
    RejectionCase{ .label = "empty terminal output span", .output_size = 0 },
    RejectionCase{ .label = "short terminal output span", .output_size = 1 },
    RejectionCase{ .label = "long terminal output span", .output_size = 3 },
    RejectionCase{ .label = "qNaN current", .current = { 0.0, quiet_nan, 0.0 } },
    RejectionCase{ .label = "+Inf current",
                   .current = { positive_infinity, 0.0, 0.0 } },
    RejectionCase{ .label = "-Inf current",
                   .current = { 0.0, negative_infinity, 0.0 } },
    RejectionCase{ .label = "zero explicit dt", .explicit_dt = 0.0 },
    RejectionCase{ .label = "negative explicit dt", .explicit_dt = -0.25 },
    RejectionCase{ .label = "qNaN explicit dt", .explicit_dt = quiet_nan },
    RejectionCase{ .label = "+Inf explicit dt",
                   .explicit_dt = positive_infinity },
    RejectionCase{ .label = "qNaN ctx.dt", .context_dt = quiet_nan },
    RejectionCase{ .label = "+Inf ctx.dt", .context_dt = positive_infinity },
    RejectionCase{ .label = "qNaN ctx.time", .time = quiet_nan },
    RejectionCase{ .label = "+Inf ctx.time", .time = positive_infinity }
  };
  for (const auto &test : rejection_cases) {
    const core::BatchShape shape{
      expected_rows + test.row_delta,
      direct_euler_lanes + test.lane_delta,
      fixture.state.stride()
    };
    const auto current =
      test.current_size == 0
        ? std::span<const double>{}
        : std::span<const double>{ test.current }.first(test.current_size);
    require_rejected(test.label,
                     shape,
                     current,
                     test.time,
                     test.context_dt,
                     test.explicit_dt,
                     test.output_size);
  }

  std::copy(state_before.begin(), state_before.end(), fixture.state.raw().begin());
  terminal_voltage = terminal_before;
  core::BatchView state{ valid_shape, fixture.state.raw() };
  const core::StepCtx valid_ctx{ .time = valid_time,
                                 .dt = valid_context_dt,
                                 .i_app = valid_current };
  REQUIRE(fixture.pipeline.advanceEuler(
            state,
            valid_ctx,
            valid_explicit_dt,
            std::span{ terminal_voltage }.first<direct_euler_lanes>())
          == Status::Success);
  CHECK(std::memcmp(fixture.state.raw().data(),
                    state_before.data(),
                    state_before.size() * sizeof(double))
        != 0);
  CHECK(core::is_finite(terminal_voltage[0]));
  CHECK(core::is_finite(terminal_voltage[1]));
  CHECK(terminal_voltage[0] != terminal_before[0]);
  CHECK(terminal_voltage[1] != terminal_before[1]);
  CHECK(terminal_voltage[2] == terminal_before[2]);
}
