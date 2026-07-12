/**
 * @file LithiumPlating.hpp
 * @brief Yang lithium-plating side reaction and additive SPM RHS mapping.
 */

#pragma once

#include "AgeingKernel.hpp"
#include "SpmObservables.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <span>
#include <type_traits>

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

struct LithiumPlatingScales
{
  real_t faradaic_scale{};        //!< n_plating * F
  real_t kinetic_charge{};        //!< n * F
  real_t thickness_denominator{}; //!< n_plating * F * plated molar density
};

template <class Real>
struct BasicLithiumPlatingOutput
{
  std::span<Real> side_reaction_current{};
};

template <class Real = real_t>
class LithiumPlatingScratch
{
public:
  explicit LithiumPlatingScratch(int n_lanes)
    : storage_{ n_lanes }
  {}

  BasicLithiumPlatingOutput<Real> view()
  {
    return { storage_.field(0) };
  }

private:
  detail::AgeingScratchStorage<Real, 1> storage_;
};

[[nodiscard]] inline bool tryLithiumPlatingScales(
  const LithiumPlatingParams &p,
  LithiumPlatingScales &scales) noexcept
{
  LithiumPlatingScales candidate;
  if (!try_multiply_nonnegative(
        p.n_plating, p.F, candidate.faradaic_scale)
      || !try_multiply_nonnegative(
        p.n, p.F, candidate.kinetic_charge)
      || !try_multiply_nonnegative(
        candidate.faradaic_scale,
        p.plated_lithium_molar_density,
        candidate.thickness_denominator))
    return false;
  scales = candidate;
  return true;
}

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
  LithiumPlatingScales scales;
  if (!tryLithiumPlatingScales(p, scales))
    return slide::Status::Invalid_parameters;
  return (p.F > 0.0 && p.Rg > 0.0 && p.n > 0.0 && p.n_plating > 0.0
          && p.reference_temperature > 0.0 && p.electrode_area > 0.0
          && p.plated_lithium_molar_density > 0.0 && p.reaction_rate_ref >= 0.0)
           ? slide::Status::Success
           : slide::Status::Invalid_parameters;
}

/**
 * Compute lane-wise plating current.  A non-success status invalidates the
 * entire output span; lanes completed before the failing lane may be present.
 */
template <class Real>
[[nodiscard]] slide::Status computeLithiumPlating(
  const LithiumPlatingParams &p,
  const BasicBatchView<const Real> &state,
  const SpmStateLayout &layout,
  const BasicStepCtx<Real> &ctx,
  const BasicSpmObservables<Real> &observables,
  BasicLithiumPlatingOutput<Real>
    output)
{
  const int lanes = state.n_lanes();
  assert(static_cast<int>(output.side_reaction_current.size()) == lanes);
  ctx.assert_valid_for(lanes);
  const auto neg = domain_index(Domain::neg);
  using std::exp;
  LithiumPlatingScales scales;
  if (!tryLithiumPlatingScales(p, scales))
    return slide::Status::Numerical_failure;
  const auto lane_status = detail::for_each_ageing_lane_while_success(
    lanes, [&](int lane) {
      const auto i = static_cast<std::size_t>(lane);
      const Real T = state.at(layout.temperature, 0, lane);
      const Real current = ctx.i_app[i] * p.electrode_area;
      const Real arrhenius = (Real{ 1 } / p.reference_temperature - Real{ 1 } / T) / p.Rg;
      const Real reaction_rate = p.reaction_rate_ref
                                 * exp(p.reaction_rate_activation * arrhenius);
      const Real ocv_negative_temperature = observables.electrode_ocv[neg][i]
                                            + (T - p.reference_temperature) * observables.negative_entropic_coefficient[i];
      const Real plating_overpotential = ocv_negative_temperature + observables.overpotential[neg][i]
                                         - p.equilibrium_potential
                                         + p.sei_resistivity_area * state.at(layout.sei_thickness, 0, lane) * current;
      const Real kinetic_factor = exp(
        -scales.kinetic_charge / (p.Rg * T)
        * p.alpha_plating * plating_overpotential);
      if (!(is_finite_primal(reaction_rate)
            && is_finite_primal(kinetic_factor)))
        return slide::Status::Numerical_failure;
      const real_t reaction_primal = primal_value(reaction_rate);
      const real_t kinetic_primal = primal_value(kinetic_factor);
      real_t faradaic_reaction{};
      real_t current_primal{};
      if (!try_multiply_nonnegative(
            scales.faradaic_scale, reaction_primal, faradaic_reaction)
          || !try_multiply_nonnegative(
            faradaic_reaction, kinetic_primal, current_primal))
        return slide::Status::Numerical_failure;
      if constexpr (std::is_same_v<std::remove_cvref_t<Real>, real_t>) {
        output.side_reaction_current[i] = current_primal;
      } else {
        const Real faradaic_reaction_real = scales.faradaic_scale * reaction_rate;
        const Real candidate = faradaic_reaction_real * kinetic_factor;
        // The two checked primal products above use this same order.  Preserve the
        // non-scalar candidate here so derivative components are not discarded.
        output.side_reaction_current[i] = candidate;
      }
      return slide::Status::Success;
    });
  return lane_status;
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
                          const BasicLithiumPlatingOutput<const Real> &output)
{
  const int lanes = state.n_lanes();
  assert(static_cast<int>(output.side_reaction_current.size()) == lanes);
  LithiumPlatingScales scales;
  const bool valid_scales = tryLithiumPlatingScales(p.mechanism, scales);
  assert(valid_scales);
  if (!valid_scales)
    return;
  const auto neg = domain_index(Domain::neg);
  detail::for_each_ageing_lane(lanes, [&](int lane) {
    const auto i = static_cast<std::size_t>(lane);
    const Real plating_current = output.side_reaction_current[i];
    for (int mode = 0; mode < NCH; ++mode)
      derivative.at(layout.z[neg], mode, lane) += p.negative_input_map[static_cast<std::size_t>(mode)] * plating_current
                                                  / scales.faradaic_scale;
    derivative.at(layout.lost_lithium, 0, lane) += plating_current * p.mechanism.electrode_area
                                                   * state.at(layout.electrode_thickness[neg], 0, lane)
                                                   * state.at(layout.specific_surface_area[neg], 0, lane);
    derivative.at(layout.plated_lithium_thickness, 0, lane) += plating_current
                                                               / scales.thickness_denominator;
  });
}

} // namespace slide::core
