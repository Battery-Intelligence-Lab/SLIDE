/**
 * @file core_Lam_test.cpp
 * @brief Direct legacy parity for LAM mechanisms 1--4.
 *
 * REGISTERED BEFORE FIRST RUN: all six raw geometry rates match Cell_SPM::LAM
 * within 1e-12 relative for every mechanism.
 */

#include "../../src/slide.hpp"
#include "../../src/core/Lam.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>

using namespace slide;

namespace {

class LamReferenceCell : public Cell_SPM
{
public:
  core::LamParams params() const
  {
    core::LamParams p;
    p.reference_temperature = T_ref;
    core::domain_value(p.particle_radius, core::Domain::pos) = geo.Rp;
    core::domain_value(p.particle_radius, core::Domain::neg) = geo.Rn;
    core::domain_value(p.model1_stress_coefficient, core::Domain::pos) = lam_p.lam1p;
    core::domain_value(p.model1_stress_coefficient, core::Domain::neg) = lam_p.lam1n;
    core::domain_value(p.model2_linear_flux, core::Domain::pos) = lam_p.lam2ap;
    core::domain_value(p.model2_linear_flux, core::Domain::neg) = lam_p.lam2an;
    core::domain_value(p.model2_sqrt_flux, core::Domain::pos) = lam_p.lam2bp;
    core::domain_value(p.model2_sqrt_flux, core::Domain::neg) = lam_p.lam2bn;
    p.model2_activation = lam_p.lam2t;
    p.model3_k = lam_p.lam3k;
    p.model3_k_activation = lam_p.lam3k_T;
    p.model3_equilibrium_potential = OCVnmc;
    core::domain_value(p.model4_area_coefficient, core::Domain::pos) = lam_p.lam4p;
    core::domain_value(p.model4_area_coefficient, core::Domain::neg) = lam_p.lam4n;
    REQUIRE(p.positive_ocv.build(OCV_curves.OCV_pos.x, OCV_curves.OCV_pos.y)
            == Status::Success);
    return p;
  }

  std::array<double, 6> evaluate(unsigned model, double z_positive,
                                 double eta_positive, double current_stress,
                                 double previous_stress, double interval)
  {
    deg_id.LAM_id.add_model(static_cast<DEG_ID::data_t>(model));
    sparam.s_dt = interval;
    sparam.s_dai_update = true;
    sparam.s_dai_p = current_stress;
    sparam.s_dai_n = -current_stress;
    sparam.s_dai_p_prev = previous_stress;
    sparam.s_dai_n_prev = -previous_stress;
    std::array<double, 6> result{};
    LAM(false, z_positive, eta_positive, &result[0], &result[1], &result[2], &result[3], &result[4], &result[5]);
    return result;
  }
};

} // namespace

TEST_CASE("LAM mechanisms match legacy Cell_SPM", "[core][ageing][LAM]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr double z_positive = 0.58;
  constexpr double eta_positive = -0.025;
  constexpr double current_stress = 11.0;
  constexpr double previous_stress = 2.0;
  constexpr double interval = 45.0;

  for (unsigned model_id = 1; model_id <= 4; ++model_id) {
    LamReferenceCell cell;
    cell.setT(309.0);
    REQUIRE(cell.setCurrent(14.0, false, false) == Status::Success);
    auto params = cell.params();
    params.model_mask = core::lam_model_bit(model_id);
    REQUIRE(core::validateLamParams(params) == Status::Success);
    const auto expected = cell.evaluate(model_id, z_positive, eta_positive, current_stress, previous_stress, interval);

    core::BatchBuilder builder;
    const auto layout = core::declareSpmState<NCH>(builder);
    const auto history_layout = core::declareStressHistory(builder);
    auto arena = builder.build(1);
    auto &legacy = cell.getStateObj();
    for (const auto domain : core::domains) {
      const auto d = core::domain_index(domain);
      const auto legacy_domain = domain == core::Domain::pos ? pos : neg;
      arena.at(layout.electrode_thickness[d], 0, 0) = legacy.thick(legacy_domain);
      arena.at(layout.specific_surface_area[d], 0, 0) = legacy.a(legacy_domain);
    }
    arena.at(layout.temperature, 0, 0) = legacy.T();
    arena.at(history_layout.previous_dai[core::domain_index(core::Domain::pos)], 0, 0) = previous_stress;
    arena.at(history_layout.previous_dai[core::domain_index(core::Domain::neg)], 0, 0) = -previous_stress;
    arena.at(history_layout.interval, 0, 0) = interval;

    core::SpmObservableScratch<NCH> observable_scratch{ 1 };
    auto observables = observable_scratch.view();
    const auto pos_index = core::domain_index(core::Domain::pos);
    observables.surface_stoichiometry[pos_index][0] = z_positive;
    observables.overpotential[pos_index][0] = eta_positive;
    core::SpmStressScratch stress_scratch{ 1 };
    auto stress = stress_scratch.view();
    stress.dai_maximum_hydrostatic[pos_index][0] = current_stress;
    stress.dai_maximum_hydrostatic[core::domain_index(core::Domain::neg)][0] = -current_stress;
    constexpr double electrode_area = 0.1 * 0.2 * 31;
    const std::array current_density{ cell.I() / electrode_area };
    const core::StepCtx ctx{ .i_app = current_density };
    core::LamScratch scratch{ 1 };
    const auto output = scratch.view();
    const core::ConstBatchView state{ core::BatchShape::from(arena), arena.raw() };
    REQUIRE(core::computeLam(params, state, layout, history_layout, ctx, observables, stress, output)
            == Status::Success);

    const auto neg_index = core::domain_index(core::Domain::neg);
    const std::array actual{ output.thickness_rate[pos_index][0],
                             output.thickness_rate[neg_index][0],
                             output.direct_area_rate[pos_index][0],
                             output.direct_area_rate[neg_index][0],
                             output.active_fraction_rate[pos_index][0],
                             output.active_fraction_rate[neg_index][0] };
    for (int quantity = 0; quantity < 6; ++quantity) {
      const auto q = static_cast<std::size_t>(quantity);
      const double scale = std::max(std::abs(expected[q]), 1e-30);
      CAPTURE(model_id, quantity, expected[q], actual[q]);
      REQUIRE(std::abs(actual[q] - expected[q]) / scale <= 1e-12);
    }
  }
}

TEST_CASE("LAM RHS composes active-fraction and direct-area loss", "[core][ageing][LAM]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  LamReferenceCell cell;
  auto params = cell.params();
  params.model_mask = core::lam_model_bit(4);
  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(1);
  core::StateArena ydot{ arena.n_rows(), 1 };
  core::RhsViews views{ core::BatchShape::from(arena) };
  views.rebind(arena.raw(), ydot.raw());
  views.zero_derivative();
  core::LamScratch scratch{ 1 };
  auto output = scratch.view();
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    output.thickness_rate[d][0] = -1e-12;
    output.active_fraction_rate[d][0] = -2e-10;
    output.direct_area_rate[d][0] = -3e-6;
  }
  core::addLamRhs(params, views.ydot, layout, output);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    REQUIRE(views.ydot.at(layout.electrode_thickness[d], 0, 0) == -1e-12);
    REQUIRE(views.ydot.at(layout.active_fraction[d], 0, 0) == -2e-10);
    REQUIRE(views.ydot.at(layout.specific_surface_area[d], 0, 0)
            == -3e-6 + 3.0 / params.particle_radius[d] * -2e-10);
  }
}
