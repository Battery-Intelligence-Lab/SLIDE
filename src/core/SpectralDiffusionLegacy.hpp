/**
 * @file SpectralDiffusionLegacy.hpp
 * @brief Q8 / P1-G0 PARITY kernel: bit-faithful replica of the legacy Cell_SPM
 *        forward-Euler solid-diffusion update, operating on v4 StateArena rows.
 *
 * THIS IS NOT THE PRODUCTION KERNEL. Its single purpose is the Phase-1 parity-drift
 * pilot (PLAN.md §6 Phase 1, §7 Q8): quantify floating-point drift between the legacy
 * scalar update in src/cells/Cell_SPM/Cell_SPM_dstate.cpp and a v4 core kernel that
 * replicates it operation-for-operation. It therefore deliberately mirrors the legacy
 * expression *shape and evaluation order* rather than a vectorised production form.
 *
 * Legacy sources replicated (verified 2026-07-07):
 *  - Derivative: Cell_SPM_dstate.cpp:47-53 (dState_diffusion):
 *      for dom in {pos,neg}:
 *        D         = electrode[dom].Dt(ArrheniusCoeff, st)              // Electrode_SPM.hpp:61-64
 *                  = st.D(dom) * exp(D_T * ArrheniusCoeff)
 *        molarFlux = electrode[dom].molarFlux(i_app, st)               // Electrode_SPM.hpp:55-59
 *                  = sign(dom) * i_app / (st.a(dom) * n * F * st.thick(dom))
 *        for k<nch: d_st.z(k,dom) = (D * A[dom](k) * z(k,dom) + B[dom](k) * molarFlux)
 *      with ArrheniusCoeff = (1/T_ref - 1/st.T()) / Rg                  // Cell_SPM.hpp:119
 *           i_app          = I() / geo.elec_surf                        // Cell_SPM_dstate.cpp:45
 *           sign(pos)=-1, sign(neg)=+1                                  // enum_definitions.hpp:65-68
 *  - Euler: Cell_SPM_dstate.cpp:257-258:  st.z(i) += dt * d_st.z(i)  for i < 2*nch.
 *
 * Ordering invariant (mirrors legacy exactly): the derivatives for ALL modes of BOTH
 * domains are computed into a buffer FIRST (legacy fills the full d_st), and only then
 * are all 2*nch z-states advanced by forward Euler. Because A[dom] is diagonal, mode k's
 * derivative reads only z(k) — so buffered-vs-interleaved is bit-identical — but we
 * buffer regardless to match the legacy control flow literally.
 *
 * Physical constants (F, Rg, n) are embedded to keep slide::core free of legacy headers;
 * their values mirror src/settings/constants.hpp (PhyConst): F=96487, Rg=8.314, n=1
 * (constants.hpp:20,21,23).
 *
 * @date 2026-07-07
 */

#pragma once

#include "StateArena.hpp"

#include <array>
#include <cmath>

namespace slide::core {

/**
 * Bit-faithful legacy diffusion parity kernel for a single batch archetype.
 *
 * @tparam NCH number of positive inner Chebyshev nodes per electrode (legacy settings::nch).
 *
 * All physical inputs are captured ONCE from a legacy Cell_SPM whose degradation and
 * thermal ODEs are blocked (Cell_SPM::setBlockDegAndTherm(true)); under that mode
 * timeStep_CC advances only the diffusion Euler step, so T, D0, a, thick and i_app are
 * constant across steps and D/molarFlux are recomputed each step to identical values.
 */
template <int NCH>
struct SpectralDiffusionLegacyKernel
{
  //!< Physical constants — mirror src/settings/constants.hpp PhyConst.
  real_t F{ 96487.0 }; //!< Faraday's constant [C mol-1]           (constants.hpp:20)
  real_t Rg{ 8.314 };  //!< ideal gas constant [J mol-1 K-1]       (constants.hpp:21)
  real_t n{ 1.0 };     //!< electrons in the main reaction [-]     (constants.hpp:23)

  //!< Scalar operating point (constant under blockDegAndTherm).
  real_t T_ref{}; //!< reference temperature [K]        (Cell_SPM::getTemperatures)
  real_t T{};     //!< cell temperature [K]             (Cell_SPM::getStateObj().T())
  real_t i_app{}; //!< current density I()/elec_surf [A m-2]

  //!< Per-domain constants, indexed [pos=0, neg=1] to match slide::Domain.
  std::array<real_t, 2> D0{};    //!< st.D(dom): diffusion const at reference T
  std::array<real_t, 2> D_T{};   //!< electrode[dom].D_T: Arrhenius activation for D
  std::array<real_t, 2> a{};     //!< st.a(dom): effective surface area
  std::array<real_t, 2> thick{}; //!< st.thick(dom): electrode thickness
  std::array<int, 2> sgn{};      //!< sign(dom): pos=-1, neg=+1

  //!< Diagonal state-space matrices from Model_SPM::makeModel(): A[dom](k), B[dom](k).
  std::array<std::array<real_t, NCH>, 2> A{};
  std::array<std::array<real_t, NCH>, 2> B{};

  /**
   * Advance the z-states of one lane by a single forward-Euler step of dt seconds,
   * replicating Cell_SPM_dstate.cpp:47-53 + :257-258 operation-for-operation.
   *
   * @param arena the state arena holding the z rows.
   * @param zp    slice for the positive-electrode z-modes (NCH rows).
   * @param zn    slice for the negative-electrode z-modes (NCH rows).
   * @param lane  lane index to advance.
   * @param dt    time step [s].
   */
  void step(StateArena &arena, StateSlice zp, StateSlice zn, int lane, real_t dt) const
  {
    const std::array<StateSlice, 2> slice{ zp, zn };

    //!< ArrheniusCoeff — computed once, exactly as legacy calcArrheniusCoeff() (Cell_SPM.hpp:119).
    const real_t ArrheniusCoeff = (1.0 / T_ref - 1.0 / T) / Rg;

    //!< Phase 1: fill the full derivative buffer for BOTH domains (mirrors d_st).
    std::array<real_t, 2 * NCH> dz{};
    for (int d = 0; d < 2; ++d) {
      const real_t D = D0[d] * std::exp(D_T[d] * ArrheniusCoeff);       //!< Electrode_SPM::Dt
      const real_t molarFlux = sgn[d] * i_app / (a[d] * n * F * thick[d]); //!< Electrode_SPM::molarFlux
      for (int k = 0; k < NCH; ++k) {
        const real_t z_k = arena.at(slice[d], k, lane);
        dz[d * NCH + k] = (D * A[d][k] * z_k + B[d][k] * molarFlux); //!< dz/dt = D*A*z + B*j
      }
    }

    //!< Phase 2: forward Euler on all 2*nch z-states (Cell_SPM_dstate.cpp:257-258).
    for (int d = 0; d < 2; ++d)
      for (int k = 0; k < NCH; ++k)
        arena.at(slice[d], k, lane) += dt * dz[d * NCH + k];
  }
};

} // namespace slide::core
