/**
 * @file SpectralDiffusion.hpp
 * @brief v4 core PRODUCTION solid-diffusion kernel (PLAN.md §3.2, §3.5), Phase 1.
 *
 * Vectorised-across-lanes replacement for the legacy per-cell forward-Euler diffusion
 * update (src/cells/Cell_SPM/Cell_SPM_dstate.cpp:47-53,257-258). One CellBatch archetype
 * holds N lanes of identical model composition in a StateArena (SoA, variable-major); this
 * kernel sweeps whole state rows (one variable across all lanes) so the inner loops
 * vectorise over lanes (PC-3). Exactly one kernel call per batch per step (PC-2).
 *
 * Physics (per electrode domain d ∈ {pos, neg}, per Chebyshev mode k, per lane c):
 *   D_eff(d,c)  = D0[d] · exp(D_T[d] · (1/T_ref − 1/T[c]) / Rg)      [m²/s]   (Arrhenius)
 *   flux(d,c)   = sgn[d] · i_app[c] / (a[d] · n · F · thick[d])      [mol m⁻² s⁻¹]
 *   dz/dt(d,k,c)= D_eff(d,c) · A[d][k] · z(d,k,c) + B[d][k] · flux(d,c)
 * The A[d], B[d] state-space diagonals come from Model_SPM::makeModel(); D_T, D0, a, thick,
 * sgn are batch-shared material/geometry constants. T[c] and i_app[c] are the per-lane
 * operating point (each cell has its own temperature and current density).
 *
 * STEPPING MODE (Phase 1): forward Euler — z(d,k,c) += dt · dz/dt(d,k,c). This is the
 * "legacy-Euler mode" the parity gates require (§5.2). The exact exponential modal
 * propagator (D-07, unconditionally stable) is added in Phase 3 as an alternative step
 * routine on the SAME state layout and params; it is intentionally NOT here yet.
 *
 * PER-LANE ARITHMETIC IDENTITY: for a single lane, stepEuler executes the identical
 * operations, in the identical order, as SpectralDiffusionLegacyKernel::step — so a batch
 * of identical lanes reproduces that (already legacy-validated, §7 Q8) kernel bit-for-bit
 * within one TU. That is this file's Phase-1 correctness gate (core_SpectralDiffusion_test).
 *
 * DEFERRED (needs the §3.11 BatchView/StepCtx interface + architecture review, not Phase-1
 * blocking): T[c]/i_app[c] arrive here as plain spans; in the final design a BatchView
 * bundles arena rows + resolved parameter views and a StepCtx carries dt. The math and
 * memory layout below are the durable part; the argument bundling is the thin part that
 * changes when BatchView lands.
 *
 * @date 2026-07-08
 */

#pragma once

#include "StateArena.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <span>
#include <vector>

namespace slide::core {

//!< Batch-shared, cold-built, immutable-in-hot-loop material/geometry + state-space model
//!< for one diffusion archetype. Indexed [pos=0, neg=1] to match slide::Domain.
//!< Populated once at build() from a CellDesign + Model_SPM (Phase 1: by the test/factory).
template <int NCH>
struct DiffusionParams {
  //!< Physical constants — mirror src/settings/constants.hpp PhyConst (F=96487, Rg=8.314, n=1).
  real_t F{ 96487.0 };  //!< Faraday constant           [C mol⁻¹]
  real_t Rg{ 8.314 };   //!< ideal gas constant         [J mol⁻¹ K⁻¹]
  real_t n{ 1.0 };      //!< electrons in main reaction  [-]
  real_t T_ref{};       //!< Arrhenius reference temperature [K]

  std::array<std::array<real_t, NCH>, 2> A{}; //!< state-space diagonal A[dom](k) (Model_SPM)
  std::array<std::array<real_t, NCH>, 2> B{}; //!< state-space input   B[dom](k) (Model_SPM)
  std::array<real_t, 2> D0{};    //!< diffusion constant at T_ref, per domain [m²/s]
  std::array<real_t, 2> D_T{};   //!< Arrhenius activation for D, per domain
  std::array<real_t, 2> a{};     //!< effective surface area, per domain
  std::array<real_t, 2> thick{}; //!< electrode thickness, per domain [m]
  std::array<int, 2> sgn{};      //!< molar-flux sign: pos=-1, neg=+1 (enum_definitions.hpp:65)
};

/**
 * Production diffusion kernel for one CellBatch archetype.
 *
 * @tparam NCH inner Chebyshev nodes per electrode (compile-time; per-batch, not a global — §3.2).
 *
 * Owns a small per-lane scratch (D_eff, flux), allocated exactly once at construction so the
 * hot step performs ZERO heap allocation (PC-1). The scratch is ephemeral (recomputed each
 * step), never part of a StateArena snapshot (PC-7).
 */
template <int NCH>
class SpectralDiffusion
{
public:
  SpectralDiffusion(DiffusionParams<NCH> params, int n_lanes)
    : p_{ params }, n_lanes_{ n_lanes },
      D_eff_(static_cast<std::size_t>(2) * n_lanes),
      flux_(static_cast<std::size_t>(2) * n_lanes)
  {
    assert(n_lanes > 0);
  }

  //!< Rows this component contributes to the batch layout: zp (NCH) then zn (NCH).
  //!< The owner declares them on the BatchBuilder in this order and hands back the slices.
  static constexpr int rows() { return 2 * NCH; }

  /**
   * Advance every lane's z-states by one forward-Euler step of dt seconds.
   *
   * @param arena state arena holding the batch.
   * @param zp    slice for positive-electrode modes (NCH rows).
   * @param zn    slice for negative-electrode modes (NCH rows).
   * @param T     per-lane temperature [K]      (size == n_lanes).
   * @param i_app per-lane current density [A m⁻²] (size == n_lanes).
   * @param dt    time step [s].
   */
  void stepEuler(StateArena &arena, StateSlice zp, StateSlice zn,
                 std::span<const real_t> T, std::span<const real_t> i_app, real_t dt)
  {
    assert(zp.rows == NCH && zn.rows == NCH);
    assert(static_cast<int>(T.size()) == n_lanes_ && static_cast<int>(i_app.size()) == n_lanes_);
    const std::array<StateSlice, 2> slice{ zp, zn };
    const int L = n_lanes_;

    //!< Phase A — per-lane effective diffusion coefficient and molar flux, per domain.
    //!< These are mode-independent, so hoisted out of the k-loop (one exp per lane, not per mode).
    for (int d = 0; d < 2; ++d) {
      real_t *De = D_eff_.data() + static_cast<std::size_t>(d) * L;
      real_t *fl = flux_.data() + static_cast<std::size_t>(d) * L;
      const real_t D0d = p_.D0[d], D_Td = p_.D_T[d];
      const real_t flux_den = p_.a[d] * p_.n * p_.F * p_.thick[d];
      const real_t sgnd = static_cast<real_t>(p_.sgn[d]);
      for (int c = 0; c < L; ++c) {
        const real_t ArrheniusCoeff = (1.0 / p_.T_ref - 1.0 / T[c]) / p_.Rg; //!< == calcArrheniusCoeff()
        De[c] = D0d * std::exp(D_Td * ArrheniusCoeff);                       //!< Electrode_SPM::Dt
        fl[c] = sgnd * i_app[c] / flux_den;                                  //!< Electrode_SPM::molarFlux
      }
    }

    //!< Phase B — modal forward Euler, one state row (all lanes) at a time (SIMD sweep over c).
    for (int d = 0; d < 2; ++d) {
      const real_t *De = D_eff_.data() + static_cast<std::size_t>(d) * L;
      const real_t *fl = flux_.data() + static_cast<std::size_t>(d) * L;
      for (int k = 0; k < NCH; ++k) {
        const real_t Ak = p_.A[d][k], Bk = p_.B[d][k];
        std::span<real_t> z = arena.row(slice[d].row_begin + k); //!< this mode across all lanes
        for (int c = 0; c < L; ++c) {
          const real_t dz = De[c] * Ak * z[c] + Bk * fl[c]; //!< dz/dt = D·A·z + B·j
          z[c] += dt * dz;                                  //!< forward-Euler advance
        }
        //!< Two statements (compute dz, then advance) — not a fused z += dt·(…) — to match the
        //!< legacy-shaped kernel's evaluation structure operation-for-operation, so a batch of
        //!< identical lanes reproduces it bit-for-bit (independent of -O3 FMA contraction).
      }
    }
  }

  const DiffusionParams<NCH> &params() const { return p_; }
  int n_lanes() const { return n_lanes_; }

private:
  DiffusionParams<NCH> p_;
  int n_lanes_{ 0 };
  std::vector<real_t> D_eff_{}; //!< [2·n_lanes] per-domain per-lane diffusion coeff (scratch)
  std::vector<real_t> flux_{};  //!< [2·n_lanes] per-domain per-lane molar flux    (scratch)
};

} // namespace slide::core
