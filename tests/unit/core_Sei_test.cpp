/**
 * @file core_Sei_test.cpp
 * @brief Direct legacy parity and RHS mapping tests for SEI mechanisms 1--4.
 *
 * REGISTERED BEFORE FIRST RUN: side-reaction current and active-fraction rate match
 * Cell_SPM::SEI to relative error <=1e-13 for every model, with and without porosity loss.
 */

#include "../../src/slide.hpp"
#include "../../src/core/Sei.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <span>

using namespace slide;

namespace {

class SeiReferenceCell : public Cell_SPM
{
public:
  core::SeiParams params(unsigned model, bool porosity) const
  {
    core::SeiParams p;
    p.model_mask = core::sei_model_bit(model);
    p.reduce_active_fraction = porosity;
    p.reference_temperature = T_ref;
    p.electrode_area = geo.elec_surf;
    p.negative_particle_radius = geo.Rn;
    p.sei_resistivity_area = rsei;
    p.sei_equilibrium_potential = OCVsei;
    p.sei_molar_volume = rhosei;
    p.electrolyte_reactant_concentration = c_elec0;
    p.main_molar_volume = Vmain;
    p.side_molar_volume = Vsei;
    p.porosity_coefficient = sei_p.sei_porosity;
    p.model1_k = sei_p.sei1k;
    p.model1_k_activation = sei_p.sei1k_T;
    p.model2_k = sei_p.sei2k;
    p.model2_k_activation = sei_p.sei2k_T;
    p.model2_D = sei_p.sei2D;
    p.model2_D_activation = sei_p.sei2D_T;
    p.model3_k = sei_p.sei3k;
    p.model3_k_activation = sei_p.sei3k_T;
    p.model3_D = sei_p.sei3D;
    p.model3_D_activation = sei_p.sei3D_T;
    p.model4_k = sei_p.sei4k;
    p.model4_k_activation = sei_p.sei4k_T;
    p.model4_D = sei_p.sei4D;
    p.model4_D_activation = sei_p.sei4D_T;
    return p;
  }

  std::pair<double, double> evaluate(unsigned model, bool porosity,
                                     double ocv_negative, double eta_negative)
  {
    deg_id.SEI_id.add_model(static_cast<DEG_ID::data_t>(model));
    deg_id.SEI_porosity = static_cast<DEG_ID::data_t>(porosity);
    double side_current{}, active_fraction_rate{};
    SEI(ocv_negative, eta_negative, &side_current, &active_fraction_rate);
    return { side_current, active_fraction_rate };
  }
};

void copy_required_state(SeiReferenceCell &cell, const core::SpmStateLayout &layout,
                         core::StateArena &arena)
{
  auto &legacy = cell.getStateObj();
  const auto neg = core::domain_index(core::Domain::neg);
  arena.at(layout.temperature, 0, 0) = legacy.T();
  arena.at(layout.sei_thickness, 0, 0) = legacy.delta();
  arena.at(layout.electrode_thickness[neg], 0, 0) = legacy.thickn();
  arena.at(layout.specific_surface_area[neg], 0, 0) = legacy.an();
  arena.at(layout.active_fraction[neg], 0, 0) = legacy.en();
}

} // namespace

TEST_CASE("SEI mechanisms and porosity coupling match legacy Cell_SPM", "[core][ageing][SEI]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr double ocv_negative = 0.18;
  constexpr double eta_negative = -0.035;

  for (unsigned model = 1; model <= 4; ++model) {
    for (const bool porosity : { false, true }) {
      SeiReferenceCell cell;
      cell.setT(310.0);
      REQUIRE(cell.setCurrent(-12.0, false, false) == Status::Success);
      const auto params = cell.params(model, porosity);
      REQUIRE(core::validateSeiParams(params) == Status::Success);
      const auto expected = cell.evaluate(model, porosity, ocv_negative, eta_negative);

      core::BatchBuilder builder;
      const auto layout = core::declareSpmState<NCH>(builder);
      auto arena = builder.build(1);
      copy_required_state(cell, layout, arena);
      const std::array current_density{ cell.I() / params.electrode_area };
      const core::StepCtx ctx{ .i_app = current_density };
      core::SpmObservableScratch<NCH> observable_scratch{ 1 };
      auto observables = observable_scratch.view();
      const auto neg = core::domain_index(core::Domain::neg);
      observables.electrode_ocv[neg][0] = ocv_negative;
      observables.negative_entropic_coefficient[0] = 0.0;
      observables.overpotential[neg][0] = eta_negative;
      core::SeiScratch scratch{ 1 };
      const auto output = scratch.view();
      const core::ConstBatchView state{ core::BatchShape::from(arena), arena.raw() };
      REQUIRE(core::computeSei(params, state, layout, ctx, observables, output)
              == Status::Success);

      const double current_scale = std::max(std::abs(expected.first), 1e-30);
      const double fraction_scale = std::max(std::abs(expected.second), 1e-30);
      CAPTURE(model, porosity, expected.first, output.side_reaction_current[0], expected.second, output.active_fraction_rate[0]);
      REQUIRE(std::abs(output.side_reaction_current[0] - expected.first) / current_scale
              <= 1e-13);
      REQUIRE(std::abs(output.active_fraction_rate[0] - expected.second) / fraction_scale
              <= 1e-13);
    }
  }
}

TEST_CASE("SEI maps side current additively into every affected arena row",
          "[core][ageing][SEI]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  SeiReferenceCell cell;
  REQUIRE(cell.setCurrent(-8.0, false, false) == Status::Success);
  core::SeiRhsParams<NCH> params;
  params.mechanism = cell.params(1, true);
  auto *model = Model_SPM<>::makeModel();
  for (int mode = 0; mode < NCH; ++mode)
    params.negative_input_map[static_cast<std::size_t>(mode)] = model->B[neg](mode);

  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(1);
  copy_required_state(cell, layout, arena);
  core::StateArena ydot{ arena.n_rows(), 1 };
  core::RhsViews views{ core::BatchShape::from(arena) };
  views.rebind(arena.raw(), ydot.raw());
  views.zero_derivative();

  std::array<double, 1> current{ 2.5e-8 };
  std::array<double, 1> fraction_rate{ -3.0e-10 };
  const core::BasicSeiOutput output{ std::span<double>{ current },
                                     std::span<double>{ fraction_rate } };
  core::addSeiRhs(params, views.y, views.ydot, layout, output);

  const auto neg_index = core::domain_index(core::Domain::neg);
  const double charge_rate = current[0] / (params.mechanism.n_sei * params.mechanism.F);
  for (int mode = 0; mode < NCH; ++mode) {
    const double expected = params.negative_input_map[static_cast<std::size_t>(mode)]
                            * charge_rate;
    const double actual = views.ydot.at(layout.z[neg_index], mode, 0);
    REQUIRE(std::abs(actual - expected)
            <= 1e-13 * std::max(std::abs(expected), 1e-30));
  }
  REQUIRE(views.ydot.at(layout.sei_thickness, 0, 0)
          == charge_rate / params.mechanism.sei_molar_volume);
  REQUIRE(views.ydot.at(layout.active_fraction[neg_index], 0, 0) == fraction_rate[0]);
  REQUIRE(views.ydot.at(layout.specific_surface_area[neg_index], 0, 0)
          == 3.0 / params.mechanism.negative_particle_radius * fraction_rate[0]);
  REQUIRE(views.ydot.at(layout.lost_lithium, 0, 0)
          == current[0] * params.mechanism.electrode_area
               * arena.at(layout.electrode_thickness[neg_index], 0, 0)
               * arena.at(layout.specific_surface_area[neg_index], 0, 0));
}

TEST_CASE("SEI rejects invalid model masks", "[core][ageing][SEI]")
{
  core::SeiParams params;
  params.model_mask = 0;
  REQUIRE(core::validateSeiParams(params) == Status::Invalid_parameters);
  params.model_mask = 0x80;
  REQUIRE(core::validateSeiParams(params) == Status::Invalid_parameters);
}
