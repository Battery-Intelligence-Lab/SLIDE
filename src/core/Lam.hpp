/**
 * @file Lam.hpp
 * @brief Loss-of-active-material mechanisms 1--4 for the v4 SPM core.
 * @surface internal
 */

#pragma once

#include "AgeingKernel.hpp"
#include "LamParams.hpp"
#include "SpmStress.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>

namespace slide::core {

template <class Real>
struct BasicLamOutput
{
  PerDomain<std::span<Real>> thickness_rate{};
  PerDomain<std::span<Real>> active_fraction_rate{};
  PerDomain<std::span<Real>> direct_area_rate{};
};

template <class Real = real_t>
class LamScratch
{
public:
  explicit LamScratch(int n_lanes)
    : storage_{ n_lanes }
  {}

  BasicLamOutput<Real> view()
  {
    BasicLamOutput<Real> output;
    std::size_t cursor = 0;
    for (auto *field : { &output.thickness_rate,
                         &output.active_fraction_rate,
                         &output.direct_area_rate })
      for (const Domain domain : domains) {
        (*field)[domain_index(domain)] = storage_.field(cursor);
        ++cursor;
      }
    return output;
  }

private:
  detail::AgeingScratchStorage<Real, 6> storage_;
};

template <class Real>
[[nodiscard]] slide::Status computeLam(const LamParams &p,
                                       const BasicBatchView<const Real> &state,
                                       const SpmStateLayout &layout,
                                       const StressHistoryLayout &history_layout,
                                       const BasicStepCtx<Real> &ctx,
                                       const BasicSpmObservables<Real> &observables,
                                       const BasicSpmStress<Real> &stress,
                                       BasicLamOutput<Real>
                                         output)
{
  const int lanes = state.n_lanes();
  ctx.assert_valid_for(lanes);
  detail::clear_ageing_fields<Real, 6>(
    lanes,
    { output.thickness_rate[domain_index(Domain::neg)],
      output.thickness_rate[domain_index(Domain::pos)],
      output.active_fraction_rate[domain_index(Domain::neg)],
      output.active_fraction_rate[domain_index(Domain::pos)],
      output.direct_area_rate[domain_index(Domain::neg)],
      output.direct_area_rate[domain_index(Domain::pos)] });
  using std::abs;
  using std::exp;
  using std::sqrt;

  const auto model_status = detail::for_each_enabled_ageing_model_lane<4>(
    p.model_mask, lanes, [&](unsigned model, int lane) {
      const auto i = static_cast<std::size_t>(lane);
      const Real T = state.at(layout.temperature, 0, lane);
      const Real arrhenius = (Real{ 1 } / p.reference_temperature - Real{ 1 } / T) / p.Rg;
      if (model == 1) {
        const Real interval = state.at(history_layout.interval, 0, lane);
        if (!(primal_value(interval) > 0.0))
          return slide::Status::Invalid_states;
        for (const Domain domain : domains) {
          const auto d = domain_index(domain);
          const Real previous = state.at(history_layout.previous_dai[d], 0, lane);
          output.thickness_rate[d][i] += -p.model1_stress_coefficient[d]
                                         * abs(stress.dai_maximum_hydrostatic[d][i] - previous) / interval;
        }
      } else if (model == 2) {
        for (const Domain domain : domains) {
          const auto d = domain_index(domain);
          const Real area = state.at(layout.specific_surface_area[d], 0, lane);
          const Real thickness = state.at(layout.electrode_thickness[d], 0, lane);
          if (!(primal_value(area) > 0.0 && primal_value(thickness) > 0.0))
            return slide::Status::Invalid_states;
          const Real flux = molar_flux_sign(domain) * ctx.i_app[i]
                            / (area * p.n * p.F * thickness);
          const Real absolute_flux = abs(flux);
          const Real temperature_factor = exp(p.model2_activation * arrhenius);
          output.active_fraction_rate[d][i] += p.model2_linear_flux[d] * temperature_factor * absolute_flux
                                               + p.model2_sqrt_flux[d] * temperature_factor * sqrt(absolute_flux);
        }
      } else if (model == 3) {
        const auto pos = domain_index(Domain::pos);
        const Real electrode_potential = p.positive_ocv.eval(observables.surface_stoichiometry[pos][i]);
        const Real dissolution_overpotential = electrode_potential + observables.overpotential[pos][i]
                                               - p.model3_equilibrium_potential;
        const Real kt = p.model3_k * exp(p.model3_k_activation * arrhenius);
        Real dissolution_current = -kt * exp(p.n * p.F / (p.Rg * T) * dissolution_overpotential)
                                   / (p.n * p.F);
        if (primal_value(dissolution_current) < -5e-6)
          dissolution_current = Real{ -5e-6 };
        output.active_fraction_rate[pos][i] += dissolution_current;
      } else {
        for (const Domain domain : domains) {
          const auto d = domain_index(domain);
          output.direct_area_rate[d][i] += -p.model4_area_coefficient[d]
                                           * state.at(layout.specific_surface_area[d], 0, lane);
        }
      }
      return slide::Status::Success;
    });
  if (model_status != slide::Status::Success)
    return model_status;

  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    const auto status = detail::for_each_ageing_lane_while_success(
      lanes, [&](int lane) {
        const auto i = static_cast<std::size_t>(lane);
        if (!(is_finite_primal(output.thickness_rate[d][i])
              && is_finite_primal(output.active_fraction_rate[d][i])
              && is_finite_primal(output.direct_area_rate[d][i])))
          return slide::Status::Numerical_failure;
        return slide::Status::Success;
      });
    if (status != slide::Status::Success)
      return status;
  }
  return slide::Status::Success;
}

template <class Real>
void addLamRhs(const LamParams &p,
               BasicBatchView<Real>
                 derivative,
               const SpmStateLayout &layout,
               const BasicLamOutput<Real> &output)
{
  const int lanes = derivative.n_lanes();
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    detail::for_each_ageing_lane(lanes, [&](int lane) {
      const auto i = static_cast<std::size_t>(lane);
      derivative.at(layout.electrode_thickness[d], 0, lane) += output.thickness_rate[d][i];
      derivative.at(layout.active_fraction[d], 0, lane) += output.active_fraction_rate[d][i];
      derivative.at(layout.specific_surface_area[d], 0, lane) += output.direct_area_rate[d][i]
                                                                 + Real{ 3 } / p.particle_radius[d] * output.active_fraction_rate[d][i];
    });
  }
}

} // namespace slide::core
