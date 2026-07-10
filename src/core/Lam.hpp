/**
 * @file Lam.hpp
 * @brief Loss-of-active-material mechanisms 1--4 for the v4 SPM core.
 */

#pragma once

#include "SpmStress.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <span>
#include <vector>

namespace slide::core {

constexpr std::uint8_t lam_model_bit(unsigned model)
{
  return model >= 1 && model <= 4
           ? static_cast<std::uint8_t>(std::uint8_t{ 1 } << (model - 1))
           : std::uint8_t{};
}

struct LamParams
{
  std::uint8_t model_mask{};
  real_t F{ 96487.0 };
  real_t Rg{ 8.314 };
  real_t n{ 1.0 };
  real_t reference_temperature{ 298.15 };
  PerDomain<real_t> particle_radius{};

  PerDomain<real_t> model1_stress_coefficient{};
  PerDomain<real_t> model2_linear_flux{};
  PerDomain<real_t> model2_sqrt_flux{};
  real_t model2_activation{};
  real_t model3_k{};
  real_t model3_k_activation{};
  real_t model3_equilibrium_potential{ 4.1 };
  PerDomain<real_t> model4_area_coefficient{};
  IndexedPiecewiseLinear positive_ocv{};
};

[[nodiscard]] inline slide::Status validateLamParams(const LamParams &p)
{
  if (p.model_mask == 0 || (p.model_mask & std::uint8_t{ 0xf0 }) != 0
      || !p.positive_ocv.valid())
    return slide::Status::Invalid_parameters;
  const std::array scalars{ p.F,
                            p.Rg,
                            p.n,
                            p.reference_temperature,
                            p.model2_activation,
                            p.model3_k,
                            p.model3_k_activation,
                            p.model3_equilibrium_potential };
  for (const real_t value : scalars)
    if (!is_finite(value))
      return slide::Status::Invalid_parameters;
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    const std::array values{ p.particle_radius[d],
                             p.model1_stress_coefficient[d],
                             p.model2_linear_flux[d],
                             p.model2_sqrt_flux[d],
                             p.model4_area_coefficient[d] };
    for (const real_t value : values)
      if (!is_finite(value))
        return slide::Status::Invalid_parameters;
    if (!(p.particle_radius[d] > 0.0))
      return slide::Status::Invalid_parameters;
  }
  return (p.F > 0.0 && p.Rg > 0.0 && p.n > 0.0
          && p.reference_temperature > 0.0 && p.model3_k >= 0.0)
           ? slide::Status::Success
           : slide::Status::Invalid_parameters;
}

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
    : storage_(static_cast<std::size_t>(6 * n_lanes)), n_lanes_(n_lanes)
  {
    assert(n_lanes > 0);
  }

  BasicLamOutput<Real> view()
  {
    const auto L = static_cast<std::size_t>(n_lanes_);
    auto storage = std::span<Real>{ storage_ };
    BasicLamOutput<Real> output;
    std::size_t cursor = 0;
    for (auto *field : { &output.thickness_rate,
                         &output.active_fraction_rate,
                         &output.direct_area_rate })
      for (const Domain domain : domains) {
        (*field)[domain_index(domain)] = storage.subspan(cursor, L);
        cursor += L;
      }
    return output;
  }

private:
  std::vector<Real> storage_{};
  int n_lanes_{};
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
  const auto L = static_cast<std::size_t>(lanes);
  ctx.assert_valid_for(lanes);
  for (auto *field : { &output.thickness_rate,
                       &output.active_fraction_rate,
                       &output.direct_area_rate })
    for (const Domain domain : domains) {
      assert((*field)[domain_index(domain)].size() == L);
      std::fill((*field)[domain_index(domain)].begin(),
                (*field)[domain_index(domain)].end(),
                Real{});
    }
  using std::abs;
  using std::exp;
  using std::sqrt;

  for (unsigned model = 1; model <= 4; ++model) {
    if ((p.model_mask & lam_model_bit(model)) == 0)
      continue;
    for (int lane = 0; lane < lanes; ++lane) {
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
    }
  }

  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    for (int lane = 0; lane < lanes; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      if (!(is_finite_primal(output.thickness_rate[d][i])
            && is_finite_primal(output.active_fraction_rate[d][i])
            && is_finite_primal(output.direct_area_rate[d][i])))
        return slide::Status::Numerical_failure;
    }
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
    for (int lane = 0; lane < lanes; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      derivative.at(layout.electrode_thickness[d], 0, lane) += output.thickness_rate[d][i];
      derivative.at(layout.active_fraction[d], 0, lane) += output.active_fraction_rate[d][i];
      derivative.at(layout.specific_surface_area[d], 0, lane) += output.direct_area_rate[d][i]
                                                                 + Real{ 3 } / p.particle_radius[d] * output.active_fraction_rate[d][i];
    }
  }
}

} // namespace slide::core
