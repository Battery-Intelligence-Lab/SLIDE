/**
 * @file core_SpmPipeline_test.cpp
 * @brief Fixed-pipeline contract: mandatory zeroing, shared diffusion inputs, and rebind safety.
 */

#include "../../src/core/SpmPipeline.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <utility>
#include <vector>

using namespace slide;

// Force every optional branch through compilation; runtime cases below focus on the pipeline
// contract independently of individual mechanism parity tests.
template class slide::core::SpmPipeline<1, true, true, true, true, true>;

TEST_CASE("SPM pipeline zeroes and evaluates the rebound trial vector", "[core][pipeline]")
{
  constexpr int NCH = 1;
  constexpr int lanes = 1;
  using Pipeline = core::SpmPipeline<NCH, false, false, false, false, false>;

  core::BatchBuilder builder;
  const auto layout = Pipeline::declareLayout(builder);
  auto arena = builder.build(lanes);
  core::StateArena derivative{ arena.n_rows(), lanes };

  core::SpmPipelineParams<NCH> params;
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

    arena.at(layout.spm.z[d], 0, 0) = 0.5;
    arena.at(layout.spm.diffusion_coefficient[d], 0, 0) = 2.0;
    arena.at(layout.spm.specific_surface_area[d], 0, 0) = 4.0;
    arena.at(layout.spm.electrode_thickness[d], 0, 0) = 5.0;
  }
  const auto neg = core::domain_index(core::Domain::neg);
  const auto pos = core::domain_index(core::Domain::pos);
  params.diffusion.A[neg][0] = -3.0;
  params.diffusion.A[pos][0] = -4.0;
  params.diffusion.B[neg][0] = 2.0;
  params.diffusion.B[pos][0] = 3.0;
  arena.at(layout.spm.temperature, 0, 0) = 1.0;

  Pipeline pipeline{ std::move(params), layout, lanes };
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
}
