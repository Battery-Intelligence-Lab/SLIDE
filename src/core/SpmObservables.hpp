/**
 * @file SpmObservables.hpp
 * @brief v4 core OBSERVABLE-RECONSTRUCTION layer for the SPM (PLAN.md §3.7/§3.12, D-10, D-23).
 *
 * "Observables derived from state, one code path shared by the RHS internals and the Recorder."
 * These are the quantities the legacy Cell_SPM recomputes on demand from its 29 states — surface
 * concentration, then (built on it) overpotential / OCV / terminal voltage. Only time/Ah/Wh are
 * path-dependent; everything here is a pure function of the arena state + operating point, so it
 * is reconstructed lazily instead of stored (D-10).
 *
 * This file lands the KEYSTONE observable: SURFACE CONCENTRATION
 *     c_surf(dom,c) = sum_j C[dom](0,j) * z(j,dom,c)  +  D[dom](0) * molarFlux(dom,c) / Dt(dom,c)
 * (legacy Cell_SPM::calcSurfaceConcentration, Cell_SPM_diffusion.cpp:22-33; the C/D output maps
 * are the same Model_SPM matrices whose surface row carried the nch!=5 bug, §2.2). Overpotential,
 * OCV and V compose on top of it (next increment).
 *
 * SHAPE (mirrors SpectralDiffusion.hpp, deliberately): batch-shared cold-built params in a struct
 * of designated-initialisable named fields (D-15, no loose 15-arg call); a free-function kernel
 * over a rebindable BatchView and StepCtx (D-23); and a lane sweep over whole SoA rows (PC-3).
 * Output scratch is caller-owned and allocated once by the composed batch (PC-1).
 *
 * @date 2026-07-10
 */

#pragma once

#include "BatchView.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <span>

namespace slide::core {

/**
 * Batch-shared parameters for the surface-concentration observable of one SPM archetype.
 * Indexed [pos=0, neg=1] to match slide::Domain. Cold-built once (Phase 1: by the test/factory);
 * immutable in the hot loop. C0/Dout0 are the SURFACE ROW (row 0) of the Model_SPM output maps.
 */
template <int NCH>
struct SurfaceConcentrationParams {
  //!< Physical constants — mirror src/settings/constants.hpp PhyConst (F=96487, Rg=8.314, n=1).
  real_t F{ 96487.0 };  //!< Faraday constant            [C mol⁻¹]
  real_t Rg{ 8.314 };   //!< ideal gas constant          [J mol⁻¹ K⁻¹]
  real_t n{ 1.0 };      //!< electrons in main reaction    [-]
  real_t T_ref{};       //!< Arrhenius reference temperature [K]
  real_t elec_surf{};   //!< electrode surface area (geo.elec_surf) [m²]

  std::array<std::array<real_t, NCH>, 2> C0{}; //!< C[dom](0, :) — surface row of the output map
  std::array<real_t, 2> Dout0{};               //!< D[dom](0)   — surface feedthrough of the output map
  std::array<real_t, 2> D0{};                  //!< st.D(dom): diffusion constant at T_ref [m²/s]
  std::array<real_t, 2> D_T{};                 //!< electrode Arrhenius activation for D
  std::array<real_t, 2> a{};                   //!< st.a(dom): effective surface area
  std::array<real_t, 2> thick{};               //!< st.thick(dom): electrode thickness [m]
  std::array<int, 2> sgn{};                    //!< molar-flux sign: pos=-1, neg=+1
};

/**
 * Reconstruct surface Li concentration for every lane and both electrodes.
 *
 * This is a free-function kernel so the same code is called by the RHS observable stage and
 * by lazy recording. The view is rebound to an integrator's current trial vector before entry.
 * The operation order inside each lane matches legacy Cell_SPM::calcSurfaceConcentration.
 */
template <int NCH, class Real>
void computeSurfaceConcentration(const SurfaceConcentrationParams<NCH> &p,
                                 const BasicBatchView<const Real> &state,
                                 StateSlice zp, StateSlice zn, StateSlice temperature,
                                 const BasicStepCtx<Real> &ctx,
                                 std::span<Real> csp, std::span<Real> csn)
{
  assert(zp.rows == NCH && zn.rows == NCH && temperature.rows == 1);
  const int L = state.n_lanes();
  ctx.assert_valid_for(L);
  assert(static_cast<int>(csp.size()) == L && static_cast<int>(csn.size()) == L);

  const std::array<StateSlice, 2> slice{ zp, zn };
  const std::array<std::span<Real>, 2> out{ csp, csn };
  const std::span<const Real> T = state.row(temperature.row_begin);
  using std::exp;

  for (int d = 0; d < 2; ++d) {
    const Real D0d = p.D0[d], D_Td = p.D_T[d], Dout0d = p.Dout0[d];
    const Real flux_den = p.a[d] * p.n * p.F * p.thick[d];
    const Real sgnd = static_cast<Real>(p.sgn[d]);
    std::span<Real> cs = out[d];

    for (int c = 0; c < L; ++c) {
      const Real Arr = (Real{ 1 } / p.T_ref - Real{ 1 } / T[c]) / p.Rg;
      const Real Dt = D0d * exp(D_Td * Arr);
      const Real molarFlux = sgnd * ctx.i_app[c] / flux_den;

      Real acc{};
      for (int j = 0; j < NCH; ++j)
        acc += p.C0[d][j] * state.at(slice[d], j, c);
      cs[c] = acc + Dout0d * molarFlux / Dt;
    }
  }
}

} // namespace slide::core
