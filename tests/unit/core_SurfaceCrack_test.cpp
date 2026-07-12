/**
 * @file core_SurfaceCrack_test.cpp
 * @brief Direct legacy parity for surface-crack mechanisms 1--5.
 *
 * REGISTERED BEFORE FIRST RUN: SEI multiplier, crack rate, and optional negative-D rate
 * match Cell_SPM::CS within 1e-12 relative for every mechanism.
 */

#include "../../src/slide.hpp"
#include "../../src/core/SurfaceCrack.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

using namespace slide;

namespace {

class CrackReferenceCell : public Cell_SPM
{
public:
  core::SurfaceCrackParams params(unsigned model, bool reduce_diffusivity) const
  {
    core::SurfaceCrackParams p;
    p.model_mask = core::surface_crack_model_bit(model);
    p.reduce_negative_diffusivity = reduce_diffusivity;
    p.reference_temperature = T_ref;
    p.electrode_area = geo.elec_surf;
    p.negative_cs_max = electrode[neg].Cmax;
    p.sei_resistivity_area = rsei;
    p.sei_equilibrium_potential = OCVsei;
    p.model1_alpha = csparam.CS1alpha;
    p.model2_alpha = csparam.CS2alpha;
    p.model3_alpha = csparam.CS3alpha;
    p.model4_alpha = csparam.CS4alpha;
    p.model4_max_surface = csparam.CS4Amax;
    p.model5_k = csparam.CS5k;
    p.model5_k_activation = csparam.CS5k_T;
    p.diffusion_exponent = csparam.CS_diffusion;
    return p;
  }

  std::array<double, 3> evaluate(unsigned model, bool reduce_diffusivity,
                                 double ocv_negative, double eta_negative,
                                 double current_stress, double previous_stress,
                                 double interval)
  {
    deg_id.CS_id.add_model(static_cast<DEG_ID::data_t>(model));
    deg_id.CS_diffusion = static_cast<DEG_ID::data_t>(reduce_diffusivity);
    sparam.s_dt = interval;
    sparam.s_dai_update = true;
    sparam.s_lares_update = true;
    sparam.s_dai_n = current_stress;
    sparam.s_dai_n_prev = previous_stress;
    sparam.s_lares_n = current_stress;
    sparam.s_lares_n_prev = previous_stress;
    double multiplier{}, crack_rate{}, diffusion_rate{};
    CS(ocv_negative, eta_negative, &multiplier, &crack_rate, &diffusion_rate);
    return { multiplier, crack_rate, diffusion_rate };
  }
};

constexpr std::size_t core_index(slide::Domain domain)
{
  return core::domain_index(domain == pos ? core::Domain::pos : core::Domain::neg);
}

} // namespace

TEST_CASE("Surface-crack mechanisms match legacy Cell_SPM", "[core][ageing][crack]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr double ocv_negative = 0.17;
  constexpr double eta_negative = -0.04;
  constexpr double current_stress = 12.0;
  constexpr double previous_stress = 3.5;
  constexpr double interval = 60.0;

  for (unsigned model_id = 1; model_id <= 5; ++model_id) {
    for (const bool reduce_diffusivity : { false, true }) {
      CrackReferenceCell cell;
      cell.setT(310.0);
      if (model_id == 5)
        cell.setC({ 0.5, 0.1 });
      REQUIRE(cell.setCurrent(-12.0, false, false) == Status::Success);
      const auto params = cell.params(model_id, reduce_diffusivity);
      REQUIRE(core::validateSurfaceCrackParams(params) == Status::Success);
      const auto expected = cell.evaluate(model_id, reduce_diffusivity, ocv_negative, eta_negative, current_stress, previous_stress, interval);
      if (model_id == 5)
        REQUIRE(expected[1] != 0.0);

      core::BatchBuilder builder;
      const auto layout = core::declareSpmState<NCH>(builder);
      const auto history_layout = core::declareStressHistory(builder);
      auto arena = builder.build(1);
      auto &legacy = cell.getStateObj();
      const auto neg_index = core::domain_index(core::Domain::neg);
      arena.at(layout.temperature, 0, 0) = legacy.T();
      arena.at(layout.sei_thickness, 0, 0) = legacy.delta();
      arena.at(layout.crack_surface, 0, 0) = legacy.CS();
      arena.at(layout.electrode_thickness[neg_index], 0, 0) = legacy.thickn();
      arena.at(layout.specific_surface_area[neg_index], 0, 0) = legacy.an();
      arena.at(layout.diffusion_coefficient[neg_index], 0, 0) = legacy.Dn();
      arena.at(history_layout.previous_dai[neg_index], 0, 0) = previous_stress;
      arena.at(history_layout.previous_laresgoiti_negative, 0, 0) = previous_stress;
      arena.at(history_layout.interval, 0, 0) = interval;

      core::SpmObservableScratch<NCH> observable_scratch{ 1 };
      auto observables = observable_scratch.view();
      for (const auto legacy_domain : { pos, neg }) {
        const auto d = core_index(legacy_domain);
        const auto concentration = cell.getC(legacy_domain);
        for (int node = 0; node < NCH + 2; ++node)
          observables.concentration[d][static_cast<std::size_t>(node)] = concentration[static_cast<std::size_t>(node)];
      }
      DPair surface{};
      REQUIRE(cell.getCSurf(surface, false));
      observables.surface_stoichiometry[core_index(pos)][0] = surface[pos] / 51385.0;
      observables.surface_stoichiometry[neg_index][0] = surface[neg] / 30555.0;
      observables.electrode_ocv[neg_index][0] = ocv_negative;
      observables.negative_entropic_coefficient[0] = 0.0;
      observables.overpotential[neg_index][0] = eta_negative;

      core::SpmStressScratch stress_scratch{ 1 };
      auto stress = stress_scratch.view();
      stress.dai_maximum_hydrostatic[neg_index][0] = current_stress;
      stress.laresgoiti_negative[0] = current_stress;
      const std::array current_density{ cell.I() / params.electrode_area };
      const core::StepCtx ctx{ .i_app = current_density };
      core::SurfaceCrackScratch scratch{ 1 };
      const auto output = scratch.view();
      const core::ConstBatchView state{ core::BatchShape::from(arena), arena.raw() };
      REQUIRE(core::computeSurfaceCrack(params, state, layout, history_layout, ctx, observables, stress, output)
              == Status::Success);

      const std::array actual{ output.sei_multiplier[0], output.crack_surface_rate[0], output.negative_diffusivity_rate[0] };
      for (int quantity = 0; quantity < 3; ++quantity) {
        const auto q = static_cast<std::size_t>(quantity);
        const double scale = std::max(std::abs(expected[q]), 1e-30);
        CAPTURE(model_id, reduce_diffusivity, quantity, expected[q], actual[q]);
        REQUIRE(std::abs(actual[q] - expected[q]) / scale <= 1e-12);
      }

      if (model_id == 1 && !reduce_diffusivity) {
        arena.at(layout.specific_surface_area[neg_index], 0, 0) = 0.0;
        REQUIRE(core::computeSurfaceCrack(params, state, layout, history_layout, ctx, observables, stress, output)
                == Status::Invalid_states);
        arena.at(layout.specific_surface_area[neg_index], 0, 0) = legacy.an();
        arena.at(history_layout.interval, 0, 0) = 0.0;
        REQUIRE(core::computeSurfaceCrack(params, state, layout, history_layout, ctx, observables, stress, output)
                == Status::Invalid_states);
      } else if (model_id == 5 && !reduce_diffusivity) {
        auto explosive = params;
        explosive.model5_k = std::numeric_limits<double>::max();
        observables.surface_stoichiometry[neg_index][0] = 0.8;
        REQUIRE(core::computeSurfaceCrack(explosive, state, layout, history_layout, ctx, observables, stress, output)
                == Status::Numerical_failure);
      }
    }
  }
}

TEST_CASE("Surface-crack parameters classify masks and numeric bounds",
          "[core][ageing][crack][validation]")
{
  CrackReferenceCell cell;
  const auto valid = cell.params(1, false);
  REQUIRE(core::validateSurfaceCrackParams(valid) == Status::Success);

  auto invalid = valid;
  invalid.model_mask = 0;
  REQUIRE(core::validateSurfaceCrackParams(invalid)
          == Status::Invalid_parameters);

  const double quiet_nan = std::bit_cast<double>(UINT64_C(0x7ff8000000000000));
  invalid = valid;
  invalid.F = quiet_nan;
  REQUIRE(core::validateSurfaceCrackParams(invalid)
          == Status::Invalid_parameters);

  invalid = valid;
  invalid.F = 0.0;
  REQUIRE(core::validateSurfaceCrackParams(invalid)
          == Status::Invalid_parameters);
}

TEST_CASE("Surface-crack RHS adds crack SEI and state rates", "[core][ageing][crack]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  CrackReferenceCell cell;
  core::SurfaceCrackRhsParams<NCH> params;
  params.mechanism = cell.params(4, true);
  auto *model = Model_SPM<>::makeModel();
  for (int mode = 0; mode < NCH; ++mode)
    params.negative_input_map[static_cast<std::size_t>(mode)] = model->B[neg](mode);
  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(1);
  core::StateArena ydot{ arena.n_rows(), 1 };
  const auto neg_index = core::domain_index(core::Domain::neg);
  arena.at(layout.electrode_thickness[neg_index], 0, 0) = 2.0;
  arena.at(layout.specific_surface_area[neg_index], 0, 0) = 3.0;
  core::RhsViews views{ core::BatchShape::from(arena) };
  views.rebind(arena.raw(), ydot.raw());
  views.zero_derivative();
  std::array<double, 1> side_current{ 2e-8 }, fraction_rate{};
  std::array<double, 1> multiplier{ 0.25 }, crack_rate{ 4e-10 }, diffusion_rate{ -3e-18 };
  const core::BasicSeiOutput sei{ std::span<double>{ side_current },
                                  std::span<double>{ fraction_rate } };
  const core::BasicSurfaceCrackOutput crack{ std::span<double>{ multiplier },
                                             std::span<double>{ crack_rate },
                                             std::span<double>{ diffusion_rate } };
  core::addSurfaceCrackRhs(params, views.y, views.ydot, layout, sei, crack);
  REQUIRE(views.ydot.at(layout.crack_surface, 0, 0) == crack_rate[0]);
  REQUIRE(views.ydot.at(layout.diffusion_coefficient[neg_index], 0, 0)
          == diffusion_rate[0]);
  const double expected_lost_lithium = side_current[0] * multiplier[0] * params.mechanism.electrode_area * 2.0 * 3.0;
  REQUIRE(std::abs(views.ydot.at(layout.lost_lithium, 0, 0)
                   - expected_lost_lithium)
          <= 1e-13 * std::abs(expected_lost_lithium));
}
