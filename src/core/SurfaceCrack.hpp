/**
 * @file SurfaceCrack.hpp
 * @brief Surface-crack ageing mechanisms 1--5 and additive SPM RHS mapping.
 * @surface internal
 */

#pragma once

#include "AgeingKernel.hpp"
#include "Sei.hpp"
#include "SpmScalarKernels.hpp"
#include "SpmStress.hpp"
#include "SurfaceCrackParams.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>

namespace slide::core {

template <class Real>
struct BasicSurfaceCrackOutput
{
  std::span<Real> sei_multiplier{};
  std::span<Real> crack_surface_rate{};
  std::span<Real> negative_diffusivity_rate{};
};

template <class Real = real_t>
class SurfaceCrackScratch
{
public:
  explicit SurfaceCrackScratch(int n_lanes)
    : storage_{ n_lanes }
  {}

  BasicSurfaceCrackOutput<Real> view()
  {
    return { storage_.field(0), storage_.field(1), storage_.field(2) };
  }

private:
  detail::AgeingScratchStorage<Real, 3> storage_;
};

template <class Real>
[[nodiscard]] slide::Status computeSurfaceCrack(
  const SurfaceCrackParams &p,
  const BasicBatchView<const Real> &state,
  const SpmStateLayout &layout,
  const StressHistoryLayout &history_layout,
  const BasicStepCtx<Real> &ctx,
  const BasicSpmObservables<Real> &observables,
  const BasicSpmStress<Real> &stress,
  BasicSurfaceCrackOutput<Real>
    output)
{
  const int lanes = state.n_lanes();
  const auto L = static_cast<std::size_t>(lanes);
  assert(output.sei_multiplier.size() == L && output.crack_surface_rate.size() == L
         && output.negative_diffusivity_rate.size() == L);
  ctx.assert_valid_for(lanes);
  detail::clear_ageing_fields<Real, 3>(
    lanes,
    { output.sei_multiplier,
      output.crack_surface_rate,
      output.negative_diffusivity_rate });
  const auto neg = domain_index(Domain::neg);
  using std::abs;
  using std::exp;
  using std::pow;
  using std::sqrt;
  const auto exceeds_crack_surface_ceiling = [&](Real crack_surface) {
    return primal_value(crack_surface) > p.model4_max_surface;
  };

  const auto model_status = detail::for_each_enabled_ageing_model_lane<5>(
    p.model_mask, lanes, [&](unsigned model, int lane) {
      const auto i = static_cast<std::size_t>(lane);
      const Real area = state.at(layout.specific_surface_area[neg], 0, lane);
      const Real thickness = state.at(layout.electrode_thickness[neg], 0, lane);
      const Real crack_surface = state.at(layout.crack_surface, 0, lane);
      if (!(primal_value(area) > 0.0 && primal_value(thickness) > 0.0
            && primal_value(crack_surface) >= 0.0))
        return slide::Status::Invalid_states;
      const Real active_surface = area * p.electrode_area * thickness;
      output.sei_multiplier[i] += crack_surface / active_surface;

      if (model == 1 || model == 2) {
        const Real interval = state.at(history_layout.interval, 0, lane);
        if (!(primal_value(interval) > 0.0))
          return slide::Status::Invalid_states;
        const Real current_stress = model == 1
                                      ? stress.laresgoiti_negative[i]
                                      : stress.dai_maximum_hydrostatic[neg][i];
        const Real previous_stress = model == 1
                                       ? state.at(history_layout.previous_laresgoiti_negative, 0, lane)
                                       : state.at(history_layout.previous_dai[neg], 0, lane);
        const Real alpha = model == 1 ? p.model1_alpha : p.model2_alpha;
        output.crack_surface_rate[i] += alpha * sqrt(abs(current_stress - previous_stress) / interval);
      } else if (model == 3) {
        const Real surface = observables.concentration[neg][i];
        const Real centre = observables.concentration[neg][static_cast<std::size_t>(layout.z[neg].rows + 1) * L + i];
        const Real gradient = (surface - centre) / p.negative_cs_max;
        output.crack_surface_rate[i] += p.model3_alpha * gradient * gradient;
      } else if (model == 4) {
        const Real remaining_surface = exceeds_crack_surface_ceiling(crack_surface)
                                         ? Real{}
                                         : Real{ p.model4_max_surface } - crack_surface;
        const Real current = ctx.i_app[i] * p.electrode_area;
        output.crack_surface_rate[i] += p.model4_alpha * remaining_surface * abs(current);
      } else {
        const Real T = state.at(layout.temperature, 0, lane);
        const Real current = ctx.i_app[i] * p.electrode_area;
        const Real ocv_negative = observables.electrode_ocv[neg][i]
                                  + (T - p.reference_temperature)
                                      * observables.negative_entropic_coefficient[i];
        const Real eta_sei = ocv_negative + observables.overpotential[neg][i]
                             - p.sei_equilibrium_potential
                             + p.sei_resistivity_area
                                 * state.at(layout.sei_thickness, 0, lane) * current;
        const Real z_negative = observables.surface_stoichiometry[neg][i];
        Real reaction_rate{};
        if (primal_value(current) <= 0.0) {
          const Real arrhenius = spm_scalar::arrheniusFactor(
            static_cast<Real>(p.reference_temperature),
            T,
            static_cast<Real>(p.Rg));
          const Real temperature_factor =
            exp(p.model5_k_activation * arrhenius);
          if (primal_value(z_negative) < 0.3)
            reaction_rate = Real{ 2 } * p.model5_k * temperature_factor;
          else if (primal_value(z_negative) >= 0.7)
            reaction_rate = p.model5_k * temperature_factor;
        }
        output.crack_surface_rate[i] += p.n_sei * p.F * reaction_rate
                                        * exp(-p.alpha_sei * p.n_sei * p.F / (p.Rg * T) * eta_sei);
      }
      return slide::Status::Success;
    });
  if (model_status != slide::Status::Success)
    return model_status;

  const auto lane_status = detail::for_each_ageing_lane_while_success(
    lanes, [&](int lane) {
      const auto i = static_cast<std::size_t>(lane);
      if (p.reduce_negative_diffusivity) {
        const Real crack_surface = state.at(layout.crack_surface, 0, lane);
        const bool saturated = exceeds_crack_surface_ceiling(crack_surface);
        const Real maximum = saturated ? crack_surface
                                       : Real{ p.model4_max_surface };
        const Real intact_fraction = saturated
                                       ? Real{}
                                       : Real{ 1 } - crack_surface / maximum;
        Real rate_fraction = p.diffusion_exponent
                             * pow(intact_fraction,
                                   p.diffusion_exponent - real_t{ 1 })
                             / maximum * output.crack_surface_rate[i];
        if (primal_value(rate_fraction) > 2e-7)
          rate_fraction = Real{ 2e-7 };
        output.negative_diffusivity_rate[i] = -rate_fraction * state.at(layout.diffusion_coefficient[neg], 0, lane);
      }
      if (!(is_finite_primal(output.sei_multiplier[i])
            && is_finite_primal(output.crack_surface_rate[i])
            && is_finite_primal(output.negative_diffusivity_rate[i])))
        return slide::Status::Numerical_failure;
      return slide::Status::Success;
    });
  return lane_status;
}

template <int NCH>
struct SurfaceCrackRhsParams
{
  SurfaceCrackParams mechanism{};
  std::array<real_t, NCH> negative_input_map{};
};

template <int NCH, class Real>
void addSurfaceCrackRhs(const SurfaceCrackRhsParams<NCH> &p,
                        const BasicBatchView<const Real> &state,
                        BasicBatchView<Real>
                          derivative,
                        const SpmStateLayout &layout,
                        const BasicSeiOutput<Real> &sei,
                        const BasicSurfaceCrackOutput<Real> &crack)
{
  const int lanes = state.n_lanes();
  const auto neg = domain_index(Domain::neg);
  detail::for_each_ageing_lane(lanes, [&](int lane) {
    const auto i = static_cast<std::size_t>(lane);
    const Real extra_side_current = sei.side_reaction_current[i] * crack.sei_multiplier[i];
    for (int mode = 0; mode < NCH; ++mode)
      derivative.at(layout.z[neg], mode, lane) += p.negative_input_map[static_cast<std::size_t>(mode)] * extra_side_current
                                                  / (p.mechanism.n_sei * p.mechanism.F);
    derivative.at(layout.lost_lithium, 0, lane) += extra_side_current * p.mechanism.electrode_area
                                                   * state.at(layout.electrode_thickness[neg], 0, lane)
                                                   * state.at(layout.specific_surface_area[neg], 0, lane);
    derivative.at(layout.crack_surface, 0, lane) += crack.crack_surface_rate[i];
    derivative.at(layout.diffusion_coefficient[neg], 0, lane) += crack.negative_diffusivity_rate[i];
  });
}

} // namespace slide::core
