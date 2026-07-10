/**
 * @file SurfaceCrack.hpp
 * @brief Surface-crack ageing mechanisms 1--5 and additive SPM RHS mapping.
 */

#pragma once

#include "Sei.hpp"
#include "SpmStress.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace slide::core {

constexpr std::uint8_t surface_crack_model_bit(unsigned model)
{
  return model >= 1 && model <= 5
           ? static_cast<std::uint8_t>(std::uint8_t{ 1 } << (model - 1))
           : std::uint8_t{};
}

struct SurfaceCrackParams
{
  std::uint8_t model_mask{};
  bool reduce_negative_diffusivity{};
  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n_sei{ 1.0 };
  real_t alpha_sei{ 1.0 };
  real_t reference_temperature{ 298.15 };
  real_t electrode_area{};
  real_t negative_cs_max{ 30555.0 };
  real_t sei_resistivity_area{};
  real_t sei_equilibrium_potential{ 0.4 };
  real_t model1_alpha{ 4.25e-5 };
  real_t model2_alpha{ 6.3e-7 };
  real_t model3_alpha{ 2.31e-16 };
  real_t model4_alpha{ 4.3306e-8 };
  real_t model4_max_surface{};
  real_t model5_k{ 1e-18 };
  real_t model5_k_activation{ -127040.0 };
  real_t diffusion_exponent{ 2.0 };
};

[[nodiscard]] inline slide::Status validateSurfaceCrackParams(const SurfaceCrackParams &p)
{
  if (p.model_mask == 0 || (p.model_mask & std::uint8_t{ 0xe0 }) != 0)
    return slide::Status::Invalid_parameters;
  const std::array values{ p.F,
                           p.Rg,
                           p.n_sei,
                           p.alpha_sei,
                           p.reference_temperature,
                           p.electrode_area,
                           p.negative_cs_max,
                           p.sei_resistivity_area,
                           p.sei_equilibrium_potential,
                           p.model1_alpha,
                           p.model2_alpha,
                           p.model3_alpha,
                           p.model4_alpha,
                           p.model4_max_surface,
                           p.model5_k,
                           p.model5_k_activation,
                           p.diffusion_exponent };
  for (const real_t value : values)
    if (!is_finite(value))
      return slide::Status::Invalid_parameters;
  if (!(p.F > 0.0 && p.Rg > 0.0 && p.n_sei > 0.0
        && p.reference_temperature > 0.0 && p.electrode_area > 0.0
        && p.negative_cs_max > 0.0 && p.model4_max_surface > 0.0
        && p.model5_k >= 0.0 && p.diffusion_exponent > 0.0))
    return slide::Status::Invalid_parameters;
  return slide::Status::Success;
}

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
    : storage_(static_cast<std::size_t>(3 * n_lanes)), n_lanes_(n_lanes)
  {
    assert(n_lanes > 0);
  }

  BasicSurfaceCrackOutput<Real> view()
  {
    const auto L = static_cast<std::size_t>(n_lanes_);
    auto storage = std::span<Real>{ storage_ };
    return { storage.first(L), storage.subspan(L, L), storage.subspan(2 * L, L) };
  }

private:
  std::vector<Real> storage_{};
  int n_lanes_{};
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
  std::fill(output.sei_multiplier.begin(), output.sei_multiplier.end(), Real{});
  std::fill(output.crack_surface_rate.begin(), output.crack_surface_rate.end(), Real{});
  std::fill(output.negative_diffusivity_rate.begin(), output.negative_diffusivity_rate.end(), Real{});
  const auto neg = domain_index(Domain::neg);
  using std::abs;
  using std::exp;
  using std::pow;
  using std::sqrt;

  for (unsigned model = 1; model <= 5; ++model) {
    if ((p.model_mask & surface_crack_model_bit(model)) == 0)
      continue;
    for (int lane = 0; lane < lanes; ++lane) {
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
        const Real maximum = std::max(p.model4_max_surface, primal_value(crack_surface));
        const Real current = ctx.i_app[i] * p.electrode_area;
        output.crack_surface_rate[i] += p.model4_alpha * (maximum - crack_surface) * abs(current);
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
          const Real temperature_factor = exp(p.model5_k_activation / p.Rg
                                              * (Real{ 1 } / p.reference_temperature - Real{ 1 } / T));
          if (primal_value(z_negative) < 0.3)
            reaction_rate = Real{ 2 } * p.model5_k * temperature_factor;
          else if (primal_value(z_negative) >= 0.7)
            reaction_rate = p.model5_k * temperature_factor;
        }
        output.crack_surface_rate[i] += p.n_sei * p.F * reaction_rate
                                        * exp(-p.alpha_sei * p.n_sei * p.F / (p.Rg * T) * eta_sei);
      }
    }
  }

  for (int lane = 0; lane < lanes; ++lane) {
    const auto i = static_cast<std::size_t>(lane);
    if (p.reduce_negative_diffusivity) {
      const Real crack_surface = state.at(layout.crack_surface, 0, lane);
      const Real maximum = std::max(p.model4_max_surface, primal_value(crack_surface));
      Real rate_fraction = p.diffusion_exponent
                           * pow(Real{ 1 } - crack_surface / maximum,
                                 p.diffusion_exponent - Real{ 1 })
                           / maximum * output.crack_surface_rate[i];
      if (primal_value(rate_fraction) > 2e-7)
        rate_fraction = Real{ 2e-7 };
      output.negative_diffusivity_rate[i] = -rate_fraction * state.at(layout.diffusion_coefficient[neg], 0, lane);
    }
    if (!(is_finite_primal(output.sei_multiplier[i])
          && is_finite_primal(output.crack_surface_rate[i])
          && is_finite_primal(output.negative_diffusivity_rate[i])))
      return slide::Status::Numerical_failure;
  }
  return slide::Status::Success;
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
  for (int lane = 0; lane < lanes; ++lane) {
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
  }
}

} // namespace slide::core
