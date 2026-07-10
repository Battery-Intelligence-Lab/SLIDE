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
 * This file lands the KEYSTONE observable: the complete particle concentration profile at
 * {surface, NCH interior nodes, centre}. The surface uses the Model_SPM C/D output map and the
 * centre uses its independent Cc/cc_coeff derivative identity (the historical nch!=5 bug path,
 * §2.2). Overpotential, OCV and V compose on the surface row (next increment).
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
 * Batch-shared parameters for the concentration observable of one SPM archetype.
 * Indexed [pos=0, neg=1] to match slide::Domain. Cold-built once (Phase 1: by the test/factory);
 * immutable in the hot loop. C/Dout contain the surface and interior Model_SPM output rows;
 * Cc/cc_coeff reconstruct the centre node from those outputs.
 */
template <int NCH>
struct SpmConcentrationParams
{
  //!< Physical constants — mirror src/settings/constants.hpp PhyConst (F=96487, Rg=8.314, n=1).
  real_t F{ 96487.0 }; //!< Faraday constant            [C mol⁻¹]
  real_t Rg{ 8.314 };  //!< ideal gas constant          [J mol⁻¹ K⁻¹]
  real_t n{ 1.0 };     //!< electrons in main reaction    [-]
  real_t T_ref{};      //!< Arrhenius reference temperature [K]
  real_t elec_surf{};  //!< electrode surface area (geo.elec_surf) [m²]

  std::array<std::array<std::array<real_t, NCH>, NCH + 1>, 2> C{};
  std::array<std::array<real_t, NCH + 1>, 2> Dout{};
  std::array<real_t, NCH + 1> Cc{}; //!< centre-node derivative map over surface+interior values
  real_t cc_coeff{};
  std::array<real_t, 2> R{};     //!< particle radius [m]
  std::array<real_t, 2> D0{};    //!< st.D(dom): diffusion constant at T_ref [m²/s]
  std::array<real_t, 2> D_T{};   //!< electrode Arrhenius activation for D
  std::array<real_t, 2> a{};     //!< st.a(dom): effective surface area
  std::array<real_t, 2> thick{}; //!< st.thick(dom): electrode thickness [m]
  std::array<int, 2> sgn{};      //!< molar-flux sign: pos=-1, neg=+1
};

/**
 * Reconstruct Li concentration at surface, interior nodes and centre for every lane/electrode.
 *
 * This is a free-function kernel so the same code is called by the RHS observable stage and
 * by lazy recording. The view is rebound to an integrator's current trial vector before entry.
 * The operation order inside each lane matches legacy Cell_SPM::calcSurfaceConcentration.
 */
template <int NCH, class Real>
void computeSpmConcentrations(const SpmConcentrationParams<NCH> &p,
                              const BasicBatchView<const Real> &state,
                              StateSlice zp, StateSlice zn, StateSlice temperature,
                              const BasicStepCtx<Real> &ctx,
                              std::array<std::span<Real>, 2> concentration)
{
  assert(zp.rows == NCH && zn.rows == NCH && temperature.rows == 1);
  const int L = state.n_lanes();
  ctx.assert_valid_for(L);
  constexpr int output_rows = NCH + 2; // surface + NCH interior + centre
  assert(static_cast<int>(concentration[0].size()) == output_rows * L
         && static_cast<int>(concentration[1].size()) == output_rows * L);

  const std::array<StateSlice, 2> slice{ zp, zn };
  const std::span<const Real> T = state.row(temperature.row_begin);
  using std::exp;

  for (int d = 0; d < 2; ++d) {
    const Real D0d = p.D0[d], D_Td = p.D_T[d];
    const Real flux_den = p.a[d] * p.n * p.F * p.thick[d];
    const Real sgnd = static_cast<Real>(p.sgn[d]);

    for (int c = 0; c < L; ++c) {
      const Real Arr = (Real{ 1 } / p.T_ref - Real{ 1 } / T[c]) / p.Rg;
      const Real Dt = D0d * exp(D_Td * Arr);
      const Real molarFlux = sgnd * ctx.i_app[c] / flux_den;

      for (int node = 0; node < NCH + 1; ++node) {
        Real acc{};
        for (int j = 0; j < NCH; ++j)
          acc += p.C[d][node][j] * state.at(slice[d], j, c);
        concentration[d][static_cast<std::size_t>(node) * L + c] = acc + p.Dout[d][node] * molarFlux / Dt;
      }

      Real centre_acc{};
      for (int node = 0; node < NCH + 1; ++node)
        centre_acc += p.Cc[node] * concentration[d][static_cast<std::size_t>(node) * L + c];
      concentration[d][static_cast<std::size_t>(NCH + 1) * L + c] = p.cc_coeff * (centre_acc + molarFlux * p.R[d] / Dt);
    }
  }
}

} // namespace slide::core
