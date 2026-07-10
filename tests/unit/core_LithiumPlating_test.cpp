/**
 * @file core_LithiumPlating_test.cpp
 * @brief Direct legacy parity for Yang lithium plating.
 *
 * REGISTERED BEFORE FIRST RUN: plating side-current matches Cell_SPM::LiPlating
 * within 1e-13 relative at charge and discharge operating points.
 */

#include "../../src/slide.hpp"
#include "../../src/core/LithiumPlating.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>

using namespace slide;

namespace {

class PlatingReferenceCell : public Cell_SPM
{
public:
  core::LithiumPlatingParams params() const
  {
    return { .F = PhyConst::F,
             .Rg = PhyConst::Rg,
             .n = PhyConst::n,
             .n_plating = npl,
             .alpha_plating = alphapl,
             .reference_temperature = T_ref,
             .electrode_area = geo.elec_surf,
             .sei_resistivity_area = rsei,
             .equilibrium_potential = OCVpl,
             .plated_lithium_molar_density = rhopl,
             .reaction_rate_ref = pl_p.pl1k,
             .reaction_rate_activation = pl_p.pl1k_T };
  }

  double evaluate(double ocv_negative, double eta_negative)
  {
    deg_id.pl_id = 1;
    return LiPlating(ocv_negative, eta_negative);
  }
};

} // namespace

TEST_CASE("Lithium plating matches legacy Cell_SPM", "[core][ageing][plating]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr double ocv_negative = 0.12;
  constexpr double eta_negative = -0.045;
  for (const double current : { -10.0, 15.0 }) {
    PlatingReferenceCell cell;
    cell.setT(311.0);
    REQUIRE(cell.setCurrent(current, false, false) == Status::Success);
    const auto params = cell.params();
    REQUIRE(core::validateLithiumPlatingParams(params) == Status::Success);
    const double expected = cell.evaluate(ocv_negative, eta_negative);

    core::BatchBuilder builder;
    const auto layout = core::declareSpmState<NCH>(builder);
    auto arena = builder.build(1);
    auto &legacy = cell.getStateObj();
    arena.at(layout.temperature, 0, 0) = legacy.T();
    arena.at(layout.sei_thickness, 0, 0) = legacy.delta();
    const auto neg = core::domain_index(core::Domain::neg);
    arena.at(layout.electrode_thickness[neg], 0, 0) = legacy.thickn();
    arena.at(layout.specific_surface_area[neg], 0, 0) = legacy.an();
    core::SpmObservableScratch<NCH> scratch{ 1 };
    auto observables = scratch.view();
    observables.electrode_ocv[neg][0] = ocv_negative;
    observables.negative_entropic_coefficient[0] = 0.0;
    observables.overpotential[neg][0] = eta_negative;
    const std::array current_density{ current / params.electrode_area };
    const core::StepCtx ctx{ .i_app = current_density };
    std::array<double, 1> actual{};
    const core::ConstBatchView state{ core::BatchShape::from(arena), arena.raw() };
    REQUIRE(core::computeLithiumPlating(params, state, layout, ctx, observables, std::span<double>{ actual })
            == Status::Success);
    CAPTURE(current, expected, actual[0]);
    REQUIRE(std::abs(actual[0] - expected)
              / std::max(std::abs(expected), 1e-30)
            <= 1e-13);
  }
}

TEST_CASE("Lithium-plating RHS maps current into all affected states",
          "[core][ageing][plating]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  PlatingReferenceCell cell;
  core::LithiumPlatingRhsParams<NCH> params;
  params.mechanism = cell.params();
  auto *model = Model_SPM<>::makeModel();
  for (int mode = 0; mode < NCH; ++mode)
    params.negative_input_map[static_cast<std::size_t>(mode)] = model->B[neg](mode);
  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(1);
  core::StateArena ydot{ arena.n_rows(), 1 };
  const auto neg = core::domain_index(core::Domain::neg);
  arena.at(layout.electrode_thickness[neg], 0, 0) = 2.0;
  arena.at(layout.specific_surface_area[neg], 0, 0) = 3.0;
  core::RhsViews views{ core::BatchShape::from(arena) };
  views.rebind(arena.raw(), ydot.raw());
  views.zero_derivative();
  const std::array<double, 1> side_current{ 3e-8 };
  core::addLithiumPlatingRhs(params, views.y, views.ydot, layout, std::span<const double>{ side_current });
  const double charge_rate = side_current[0] / (params.mechanism.n_plating * params.mechanism.F);
  REQUIRE(views.ydot.at(layout.plated_lithium_thickness, 0, 0)
          == charge_rate / params.mechanism.plated_lithium_molar_density);
  const double expected_lost_lithium = side_current[0] * params.mechanism.electrode_area * 2.0 * 3.0;
  REQUIRE(std::abs(views.ydot.at(layout.lost_lithium, 0, 0)
                   - expected_lost_lithium)
          <= 1e-13 * expected_lost_lithium);
}
