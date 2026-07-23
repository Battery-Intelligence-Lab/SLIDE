/**
 * @file core_ThermalLumped_test.cpp
 * @brief Phase-1 analytic and state-contract tests for the lumped thermal RHS.
 *
 * REGISTERED BEFORE FIRST RUN:
 *  - dT/dt matches (Qgen + qext + hA*(Tenv-T))/(rho*Cp*V) to 1e-14 relative.
 *  - internal generated energy and thermal elapsed time are arena derivatives, not hidden state.
 *  - malformed physical descriptions return Invalid_parameters and clear compiled output.
 *  - overflowing h*A returns Invalid_parameters and atomically clears all compiled output.
 *  - with freshly compiled parameters and positive temperature, opaque bit-built qNaN in
 *    either internal or external heat returns Invalid_states; an all-finite control succeeds.
 *  - in the adiabatic zero-heat limit, dT/dt is exactly 0 and d(t_thermal)/dt is exactly 1
 *    after the derivative arena is zeroed.
 */

#include "../../src/core/ThermalLumped.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>

using namespace slide;

TEST_CASE("ThermalLumped RHS matches the analytic energy balance", "[core][thermal]")
{
  constexpr int lanes = 3;
  const core::ThermalDesign design{ .density = 1626.0,
                                    .heat_capacity = 750.0,
                                    .volume = 1.6850e-4 * 0.62,
                                    .surface_area = 0.02,
                                    .h_conv = 45.0,
                                    .environment_temperature = 298.15 };
  core::ThermalLumpedParams params;
  REQUIRE(core::compileThermalLumped(design, params) == Status::Success);

  core::BatchBuilder builder;
  const auto temperature = builder.declare({ "T", 1, core::Unit::K });
  const auto layout = core::declareThermalLumped(builder, temperature);
  auto arena = builder.build(lanes);
  core::StateArena ydot{ arena.n_rows(), lanes };

  constexpr std::array temperatures{ 298.15, 308.15, 288.15 };
  constexpr std::array internal_heat{ 4.0, 2.0, -1.0 };
  constexpr std::array external_heat{ 1.0, -3.0, 0.5 };
  std::array current_density{ 0.0, 0.0, 0.0 };
  for (int lane = 0; lane < lanes; ++lane) {
    arena.at(layout.temperature, 0, lane) = temperatures[static_cast<std::size_t>(lane)];
    arena.at(layout.external_heat_flow, 0, lane) = external_heat[static_cast<std::size_t>(lane)];
  }

  core::RhsViews views{ core::BatchShape::from(arena) };
  views.rebind(arena.raw(), ydot.raw());
  views.zero_derivative();
  core::SpmObservableScratch<1> scratch{ lanes };
  auto observables = scratch.view();
  for (int lane = 0; lane < lanes; ++lane)
    observables.total_heat[static_cast<std::size_t>(lane)] = internal_heat[static_cast<std::size_t>(lane)];
  const core::StepCtx ctx{ .i_app = current_density };

  REQUIRE(core::addThermalLumpedRhs(params, views.y, views.ydot, layout, observables, ctx)
          == Status::Success);
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    const double expected = (internal_heat[i] + external_heat[i]
                             + params.environment_conductance
                                 * (params.environment_temperature - temperatures[i]))
                            / params.thermal_capacity;
    const double actual = views.ydot.at(layout.temperature, 0, lane);
    CAPTURE(lane, expected, actual);
    REQUIRE(std::abs(actual - expected) <= 1e-14 * std::max(1.0, std::abs(expected)));
    REQUIRE(views.ydot.at(layout.generated_heat_energy, 0, lane) == internal_heat[i]);
    REQUIRE(views.ydot.at(layout.thermal_elapsed_time, 0, lane) == 1.0);
  }

  REQUIRE(builder.is_ode_row(layout.temperature.row_begin));
  REQUIRE(builder.is_ode_row(layout.generated_heat_energy.row_begin));
  REQUIRE(builder.is_ode_row(layout.thermal_elapsed_time.row_begin));
  REQUIRE_FALSE(builder.is_ode_row(layout.external_heat_flow.row_begin));

  // Component kernels accumulate: a second load adds to, rather than overwrites, shared rows.
  REQUIRE(core::addThermalLumpedRhs(params, views.y, views.ydot, layout, observables, ctx)
          == Status::Success);
  const double first_expected = (internal_heat[0] + external_heat[0])
                                / params.thermal_capacity;
  REQUIRE(views.ydot.at(layout.temperature, 0, 0) == 2.0 * first_expected);
  REQUIRE(views.ydot.at(layout.generated_heat_energy, 0, 0) == 2.0 * internal_heat[0]);
  REQUIRE(views.ydot.at(layout.thermal_elapsed_time, 0, 0) == 2.0);
}

TEST_CASE("ThermalLumped validates cold parameters and hot state", "[core][thermal]")
{
  core::ThermalDesign design{ .density = 1.0,
                              .heat_capacity = 1.0,
                              .volume = 1.0,
                              .surface_area = 1.0,
                              .h_conv = 1.0,
                              .environment_temperature = 300.0 };
  core::ThermalLumpedParams params;
  REQUIRE(core::compileThermalLumped(design, params) == Status::Success);
  design.volume = 0.0;
  REQUIRE(core::compileThermalLumped(design, params) == Status::Invalid_parameters);
  REQUIRE(params.thermal_capacity == 0.0);

  design.volume = 1.0;
  // Construct the IEEE payload through its bits: under -ffast-math clang may fold
  // numeric_limits<double>::quiet_NaN() itself to zero before validation sees it.
  design.h_conv = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  REQUIRE(core::compileThermalLumped(design, params) == Status::Invalid_parameters);

  design.h_conv = 1.0;
  REQUIRE(core::compileThermalLumped(design, params) == Status::Success);

  design.density = std::numeric_limits<double>::max();
  design.heat_capacity = 2.0;
  REQUIRE(core::compileThermalLumped(design, params)
          == Status::Invalid_parameters);
  REQUIRE(params.thermal_capacity == 0.0);
  design.density = 1.0;
  design.heat_capacity = 1.0;

  core::BatchBuilder builder;
  const auto temperature = builder.declare({ "T", 1, core::Unit::K });
  const auto layout = core::declareThermalLumped(builder, temperature);
  auto arena = builder.build(1);
  core::StateArena ydot{ arena.n_rows(), 1 };
  arena.at(layout.temperature, 0, 0) = -1.0;
  std::array<double, 1> current{};
  core::SpmObservableScratch<1> scratch{ 1 };
  auto observables = scratch.view();
  observables.total_heat[0] = 0.0;
  core::RhsViews views{ core::BatchShape::from(arena) };
  views.rebind(arena.raw(), ydot.raw());
  views.zero_derivative();
  REQUIRE(core::addThermalLumpedRhs(params, views.y, views.ydot, layout, observables, core::StepCtx{ .i_app = current })
          == Status::Invalid_states);

  design.h_conv = 1.0;
  design.surface_area = 1.0;
  REQUIRE(core::compileThermalLumped(design, params) == Status::Success);
  design.h_conv = std::numeric_limits<double>::max();
  design.surface_area = 2.0;
  REQUIRE(core::compileThermalLumped(design, params)
          == Status::Invalid_parameters);
  REQUIRE(params.thermal_capacity == 0.0);
  REQUIRE(params.environment_conductance == 0.0);
  REQUIRE(params.environment_temperature == 0.0);

  design.h_conv = 1.0;
  design.surface_area = 1.0;
  REQUIRE(core::compileThermalLumped(design, params) == Status::Success);
  arena.at(layout.temperature, 0, 0) = 300.0;
  arena.at(layout.external_heat_flow, 0, 0) = 0.0;
  observables.total_heat[0] = 0.0;
  views.zero_derivative();
  REQUIRE(core::addThermalLumpedRhs(params, views.y, views.ydot, layout, observables, core::StepCtx{ .i_app = current })
          == Status::Success);
  REQUIRE(std::isfinite(views.ydot.at(layout.temperature, 0, 0)));
  REQUIRE(std::isfinite(views.ydot.at(layout.generated_heat_energy, 0, 0)));
  REQUIRE(std::isfinite(views.ydot.at(layout.thermal_elapsed_time, 0, 0)));

  const double opaque_qnan =
    std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  observables.total_heat[0] = opaque_qnan;
  REQUIRE(core::addThermalLumpedRhs(params, views.y, views.ydot, layout, observables, core::StepCtx{ .i_app = current })
          == Status::Invalid_states);
  observables.total_heat[0] = 0.0;
  arena.at(layout.external_heat_flow, 0, 0) = opaque_qnan;
  REQUIRE(core::addThermalLumpedRhs(params, views.y, views.ydot, layout, observables, core::StepCtx{ .i_app = current })
          == Status::Invalid_states);
  arena.at(layout.external_heat_flow, 0, 0) = 0.0;

  design.h_conv = 0.0;
  REQUIRE(core::compileThermalLumped(design, params) == Status::Success);
  REQUIRE(params.environment_conductance == 0.0);
  arena.at(layout.temperature, 0, 0) = 310.0;
  views.zero_derivative();
  REQUIRE(core::addThermalLumpedRhs(params, views.y, views.ydot, layout, observables, core::StepCtx{ .i_app = current })
          == Status::Success);
  REQUIRE(views.ydot.at(layout.temperature, 0, 0) == 0.0);
  REQUIRE(views.ydot.at(layout.thermal_elapsed_time, 0, 0) == 1.0);
}
