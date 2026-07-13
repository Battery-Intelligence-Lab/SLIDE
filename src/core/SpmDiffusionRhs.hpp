/**
 * @file SpmDiffusionRhs.hpp
 * @brief The solid-diffusion modal RHS mapping.
 *
 * Owns: `SpmDiffusionRhsParams` and `addSpmDiffusionRhs`. Implements PLAN.md §3.5 (exponential
 * modal propagator) through the shared `spm_scalar::diffusionRate`. Hot: runs once per RHS
 * evaluation, over both electrodes and every mode.
 *
 * Every other mechanism's RHS mapping already had its own header (`Sei.hpp`, `Lam.hpp`,
 * `SurfaceCrack.hpp`, `LithiumPlating.hpp`, `ThermalLumped.hpp`); diffusion's was the one still
 * living inside the pipeline that composes them (MC-1, MC-2).
 * @surface internal
 */

#pragma once

#include "SpmObservables.hpp"
#include "SpmScalarKernels.hpp"

#include <array>
#include <cassert>

namespace slide::core {

template <int NCH>
struct SpmDiffusionRhsParams
{
  PerDomain<std::array<real_t, NCH>> A{};
  PerDomain<std::array<real_t, NCH>> B{};
};

/** Add the two solid-diffusion modal systems using the shared observable reconstruction. */
template <int NCH, class Real>
void addSpmDiffusionRhs(const SpmDiffusionRhsParams<NCH> &p,
                        const BasicBatchView<const Real> &state,
                        BasicBatchView<Real>
                          derivative,
                        const SpmStateLayout &layout,
                        const BasicSpmObservables<Real> &observables)
{
  const int lanes = state.n_lanes();
  assert(derivative.n_lanes() == lanes && derivative.n_rows() == state.n_rows());
  for (const Domain domain : domains) {
    const auto d = domain_index(domain);
    assert(layout.z[d].rows == NCH
           && static_cast<int>(observables.effective_diffusivity[d].size()) == lanes
           && static_cast<int>(observables.molar_flux[d].size()) == lanes);
    for (int mode = 0; mode < NCH; ++mode) {
      const Real A = p.A[d][static_cast<std::size_t>(mode)];
      const Real B = p.B[d][static_cast<std::size_t>(mode)];
      for (int lane = 0; lane < lanes; ++lane) {
        const auto i = static_cast<std::size_t>(lane);
        derivative.at(layout.z[d], mode, lane) += spm_scalar::diffusionRate(
          state.at(layout.z[d], mode, lane),
          observables.effective_diffusivity[d][i],
          A,
          B,
          observables.molar_flux[d][i]);
      }
    }
  }
}

} // namespace slide::core
