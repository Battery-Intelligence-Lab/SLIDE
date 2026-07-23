/**
 * @file Sei.hpp
 * @brief Scalar-generic SEI mechanisms and additive RHS mapping (PLAN.md Phase 1).
 * @surface internal
 */

#pragma once

#include "AgeingKernel.hpp"
#include "Numeric.hpp"
#include "SeiParams.hpp"
#include "SpmObservables.hpp"
#include "SpmScalarKernels.hpp"
#include "SpmState.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>

namespace slide::core {

template <class Real>
struct BasicSeiOutput
{
  std::span<Real> side_reaction_current{}; //!< A/m2
  std::span<Real> active_fraction_rate{};  //!< 1/s
};

template <class Real = real_t>
class SeiScratch
{
public:
  explicit SeiScratch(int n_lanes)
    : storage_{ n_lanes }
  {}

  BasicSeiOutput<Real> view()
  {
    return { storage_.field(0), storage_.field(1) };
  }

private:
  detail::AgeingScratchStorage<Real, 2> storage_;
};

template <class Real>
[[nodiscard]] slide::Status computeSei(const SeiParams &p,
                                       const BasicBatchView<const Real> &state,
                                       const SpmStateLayout &layout,
                                       const BasicStepCtx<Real> &ctx,
                                       const BasicSpmObservables<Real> &observables,
                                       BasicSeiOutput<Real>
                                         output)
{
  const int lanes = state.n_lanes();
  [[maybe_unused]] const auto count = static_cast<std::size_t>(lanes);
  assert(output.side_reaction_current.size() == count
         && output.active_fraction_rate.size() == count);
  assert(observables.negative_entropic_coefficient.size() == count);
  ctx.assert_valid_for(lanes);
  detail::clear_ageing_fields<Real, 2>(
    lanes,
    { output.side_reaction_current, output.active_fraction_rate });

  const auto neg = domain_index(Domain::neg);
  using std::exp;
  const auto model_status = detail::for_each_enabled_ageing_model_lane<4>(
    p.model_mask, lanes, [&](unsigned model, int lane) {
      const auto i = static_cast<std::size_t>(lane);
      const Real T = state.at(layout.temperature, 0, lane);
      const Real delta = state.at(layout.sei_thickness, 0, lane);
      const Real current = ctx.i_app[i] * p.electrode_area;
      const Real arrhenius = spm_scalar::arrheniusFactor(
        static_cast<Real>(p.reference_temperature), T, static_cast<Real>(p.Rg));
      const Real ocv_neg_temperature = observables.electrode_ocv[neg][i]
                                       + (T - p.reference_temperature)
                                           * observables.negative_entropic_coefficient[i];
      const Real eta_neg = observables.overpotential[neg][i];
      const Real film_drop = p.sei_resistivity_area * delta * current;

      Real contribution{};
      if (model == 1) {
        const Real kt = p.model1_k * exp(p.model1_k_activation * arrhenius);
        contribution = p.n_sei * p.F * kt
                       * exp(-p.n_sei * p.F / (p.Rg * T) * p.alpha_sei
                             * (ocv_neg_temperature + eta_neg
                                - p.sei_equilibrium_potential + film_drop));
      } else if (model == 2) {
        const Real kt = p.model2_k * exp(p.model2_k_activation * arrhenius);
        const Real Dt = p.model2_D * exp(p.model2_D_activation * arrhenius);
        const Real kinetics = p.n_sei * p.F * kt
                              * exp(-p.n_sei * p.F / (p.Rg * T) * p.alpha_sei
                                    * (ocv_neg_temperature + eta_neg
                                       - p.sei_equilibrium_potential + film_drop));
        const Real diffusion = delta / (p.n_sei * p.F * Dt);
        contribution = p.electrolyte_reactant_concentration
                       / (Real{ 1 } / kinetics + diffusion);
      } else {
        const bool model3 = model == 3;
        const Real k_ref = model3 ? p.model3_k : p.model4_k;
        const Real k_activation = model3 ? p.model3_k_activation
                                         : p.model4_k_activation;
        const Real D_ref = model3 ? p.model3_D : p.model4_D;
        const Real D_activation = model3 ? p.model3_D_activation
                                         : p.model4_D_activation;
        const Real kt = k_ref * exp(k_activation * arrhenius);
        const Real Dt = D_ref * exp(D_activation * arrhenius);
        constexpr real_t a_L_K = 0.134461;
        const Real first_exponent = model3 ? eta_neg + film_drop : eta_neg;
        const Real first = a_L_K
                           * exp(-p.n_sei * p.F * first_exponent / (p.Rg * T));
        const Real second = p.n_sei * p.F * kt
                            * exp(-p.n_sei * p.F / (p.Rg * T) * p.alpha_sei
                                  * (ocv_neg_temperature
                                     - p.sei_equilibrium_potential));
        const Real third = delta / (p.n_sei * p.F * Dt);
        contribution = first / (Real{ 1 } / second + third);
      }
      output.side_reaction_current[i] += contribution;
      return slide::Status::Success;
    });
  if (model_status != slide::Status::Success)
    return model_status;

  const auto lane_status = detail::for_each_ageing_lane_while_success(
    lanes, [&](int lane) {
      const auto i = static_cast<std::size_t>(lane);
      const Real side_current = output.side_reaction_current[i];
      if (!is_finite_primal(side_current))
        return slide::Status::Numerical_failure;
      if (p.reduce_active_fraction) {
        const Real area = state.at(layout.specific_surface_area[neg], 0, lane);
        const Real thickness = state.at(layout.electrode_thickness[neg], 0, lane);
        if (!(primal_value(area) > 0.0 && primal_value(thickness) > 0.0))
          return slide::Status::Invalid_states;
        const Real molar_flux = ctx.i_app[i] / (area * p.n * p.F * thickness);
        output.active_fraction_rate[i] = -p.porosity_coefficient
                                         * (molar_flux * p.main_molar_volume
                                            + side_current * p.side_molar_volume);
        if (!is_finite_primal(output.active_fraction_rate[i]))
          return slide::Status::Numerical_failure;
      }
      return slide::Status::Success;
    });
  return lane_status;
}

template <int NCH>
struct SeiRhsParams
{
  SeiParams mechanism{};
  std::array<real_t, NCH> negative_input_map{}; //!< Model_SPM B[neg]
};

template <int NCH, class Real>
void addSeiRhs(const SeiRhsParams<NCH> &p,
               const BasicBatchView<const Real> &state,
               BasicBatchView<Real>
                 derivative,
               const SpmStateLayout &layout,
               const BasicSeiOutput<Real> &output)
{
  const int lanes = state.n_lanes();
  const auto neg = domain_index(Domain::neg);
  assert(derivative.n_lanes() == lanes && layout.z[neg].rows == NCH);
  detail::for_each_ageing_lane(lanes, [&](int lane) {
    const auto i = static_cast<std::size_t>(lane);
    const Real side_current = output.side_reaction_current[i];
    const Real active_fraction_rate = output.active_fraction_rate[i];
    for (int mode = 0; mode < NCH; ++mode)
      derivative.at(layout.z[neg], mode, lane) += p.negative_input_map[static_cast<std::size_t>(mode)] * side_current
                                                  / (p.mechanism.n_sei * p.mechanism.F);
    derivative.at(layout.sei_thickness, 0, lane) += side_current
                                                    / (p.mechanism.n_sei * p.mechanism.F * p.mechanism.sei_molar_volume);
    derivative.at(layout.lost_lithium, 0, lane) += side_current * p.mechanism.electrode_area
                                                   * state.at(layout.electrode_thickness[neg], 0, lane)
                                                   * state.at(layout.specific_surface_area[neg], 0, lane);
    derivative.at(layout.active_fraction[neg], 0, lane) += active_fraction_rate;
    derivative.at(layout.specific_surface_area[neg], 0, lane) += Real{ 3 } / p.mechanism.negative_particle_radius * active_fraction_rate;
  });
}

} // namespace slide::core
