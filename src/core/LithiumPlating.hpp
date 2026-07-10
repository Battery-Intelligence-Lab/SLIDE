/**
 * @file LithiumPlating.hpp
 * @brief Yang lithium-plating side reaction and additive SPM RHS mapping.
 */

#pragma once

#include "SpmObservables.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <span>

namespace slide::core {

struct LithiumPlatingParams
{
  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n{ 1.0 };
  real_t n_plating{ 1.0 };
  real_t alpha_plating{ 1.0 };
  real_t reference_temperature{ 298.15 };
  real_t electrode_area{};
  real_t sei_resistivity_area{};
  real_t equilibrium_potential{};
  real_t plated_lithium_molar_density{ 10e6 };
  real_t reaction_rate_ref{ 4.5e-10 };
  real_t reaction_rate_activation{ -2.014008e5 };
};

[[nodiscard]] inline slide::Status validateLithiumPlatingParams(
  const LithiumPlatingParams &p)
{
  const std::array values{ p.F,
                           p.Rg,
                           p.n,
                           p.n_plating,
                           p.alpha_plating,
                           p.reference_temperature,
                           p.electrode_area,
                           p.sei_resistivity_area,
                           p.equilibrium_potential,
                           p.plated_lithium_molar_density,
                           p.reaction_rate_ref,
                           p.reaction_rate_activation };
  for (const real_t value : values)
    if (!is_finite(value))
      return slide::Status::Invalid_parameters;
  return (p.F > 0.0 && p.Rg > 0.0 && p.n > 0.0 && p.n_plating > 0.0
          && p.reference_temperature > 0.0 && p.electrode_area > 0.0
          && p.plated_lithium_molar_density > 0.0 && p.reaction_rate_ref >= 0.0)
           ? slide::Status::Success
           : slide::Status::Invalid_parameters;
}

template <class Real>
[[nodiscard]] slide::Status computeLithiumPlating(
  const LithiumPlatingParams &p,
  const BasicBatchView<const Real> &state,
  const SpmStateLayout &layout,
  const BasicStepCtx<Real> &ctx,
  const BasicSpmObservables<Real> &observables,
  std::span<Real>
    side_reaction_current)
{
  const int lanes = state.n_lanes();
  assert(static_cast<int>(side_reaction_current.size()) == lanes);
  ctx.assert_valid_for(lanes);
  const auto neg = domain_index(Domain::neg);
  using std::exp;
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    const Real T = state.at(layout.temperature, 0, lane);
    const Real current = ctx.i_app[i] * p.electrode_area;
    const Real arrhenius = (Real{ 1 } / p.reference_temperature - Real{ 1 } / T) / p.Rg;
    const Real reaction_rate = p.reaction_rate_ref * exp(p.reaction_rate_activation * arrhenius);
    const Real ocv_negative_temperature = observables.electrode_ocv[neg][i]
                                          + (T - p.reference_temperature) * observables.negative_entropic_coefficient[i];
    const Real plating_overpotential = ocv_negative_temperature + observables.overpotential[neg][i]
                                       - p.equilibrium_potential
                                       + p.sei_resistivity_area * state.at(layout.sei_thickness, 0, lane) * current;
    side_reaction_current[i] = p.n_plating * p.F * reaction_rate
                               * exp(-p.n * p.F / (p.Rg * T) * p.alpha_plating * plating_overpotential);
    if (!is_finite_primal(side_reaction_current[i]))
      return slide::Status::Numerical_failure;
  }
  return slide::Status::Success;
}

template <int NCH>
struct LithiumPlatingRhsParams
{
  LithiumPlatingParams mechanism{};
  std::array<real_t, NCH> negative_input_map{};
};

template <int NCH, class Real>
void addLithiumPlatingRhs(const LithiumPlatingRhsParams<NCH> &p,
                          const BasicBatchView<const Real> &state,
                          BasicBatchView<Real>
                            derivative,
                          const SpmStateLayout &layout,
                          std::span<const Real>
                            side_reaction_current)
{
  const int lanes = state.n_lanes();
  assert(static_cast<int>(side_reaction_current.size()) == lanes);
  const auto neg = domain_index(Domain::neg);
  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    const Real plating_current = side_reaction_current[i];
    for (int mode = 0; mode < NCH; ++mode)
      derivative.at(layout.z[neg], mode, lane) += p.negative_input_map[static_cast<std::size_t>(mode)] * plating_current
                                                  / (p.mechanism.n_plating * p.mechanism.F);
    derivative.at(layout.lost_lithium, 0, lane) += plating_current * p.mechanism.electrode_area
                                                   * state.at(layout.electrode_thickness[neg], 0, lane)
                                                   * state.at(layout.specific_surface_area[neg], 0, lane);
    derivative.at(layout.plated_lithium_thickness, 0, lane) += plating_current
                                                               / (p.mechanism.n_plating * p.mechanism.F
                                                                  * p.mechanism.plated_lithium_molar_density);
  }
}

} // namespace slide::core
