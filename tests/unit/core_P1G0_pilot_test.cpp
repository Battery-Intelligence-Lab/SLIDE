/**
 * @file core_P1G0_pilot_test.cpp
 * @brief Phase-1 P1-G0 parity-drift pilot (PLAN.md §6 Phase 1, §7 Q8).
 *
 * Measures floating-point drift, over 1200 lockstep 1 s forward-Euler steps on ONE cell,
 * between:
 *   (legacy) slide::Cell_SPM::timeStep_CC with degradation+thermal blocked — the reference
 *            solid-diffusion update in src/cells/Cell_SPM/Cell_SPM_dstate.cpp:47-53,257-258;
 *   (core)   slide::core::SpectralDiffusionLegacyKernel on a v4 StateArena — an
 *            operation-for-operation replica of that update.
 *
 * ===========================================================================================
 * REGISTERED BANDS — written BEFORE the first run (CLAUDE.md §3; register the prediction
 * before the run). The kernel is a faithful op-order replica compiled by the same toolchain,
 * so the two updates should be bit-identical.
 *
 *   H0 (primary hypothesis):  max |z_legacy - z_core| == 0.0 EXACTLY, over all 1200 steps and
 *                             all 2*nch = 10 z-modes. Reported via CHECK (non-fatal): a tiny
 *                             nonzero drift (e.g. from differing FMA contraction across TUs)
 *                             still surfaces rather than aborting the run.
 *
 *   Q8 keep-band (decisive):  max relative drift  drift_rel <= 1e-12  -> REQUIRE (fatal).
 *                             This is the criterion that gates keeping the parity approach.
 *
 *   Sanity:                   every compared z (legacy and core) is std::isfinite.
 *
 * A FALSIFIED H0 with drift_rel still <= 1e-12 is a valid, reportable outcome — the bands are
 * NOT to be tuned to fit the result.
 * ===========================================================================================
 *
 * @date 2026-07-07
 */

#include "../../src/slide.hpp"                    // legacy Cell_SPM, Model_SPM, settings, sign, Domain
#include "../../src/core/BatchBuilder.hpp"        // slide::core StateArena/BatchBuilder
#include "../../src/core/SpectralDiffusionLegacy.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

using namespace slide;

TEST_CASE("P1-G0 parity-drift pilot: legacy vs v4 core diffusion kernel", "[core][P1-G0]")
{
  constexpr int NCH = static_cast<int>(settings::nch); // 5
  constexpr int NZ = 2 * NCH;                           // 10 z-modes (zp[0..4], zn[0..4])
  constexpr double dt = 1.0;
  constexpr int nsteps = 1200;

  // ---------------------------------------------------------------------------------------
  // Legacy side: default Cell_SPM, block degradation + thermal, apply 1C (16 A) discharge.
  // ---------------------------------------------------------------------------------------
  Cell_SPM cell; // default constructor: 16 Ah Kokam NMC, SOC 0.5, T = settings::T_ENV
  cell.setBlockDegAndTherm(true);
  const Status st_current = cell.setCurrent(16.0, false, false); // 1C discharge, no V-check, no print
  REQUIRE(st_current == Status::Success);

  // ---------------------------------------------------------------------------------------
  // Capture the operating point from the SAME sources the legacy cell uses. Under
  // blockDegAndTherm, T / D0 / a / thick / I are constant across all steps.
  // ---------------------------------------------------------------------------------------
  auto &s0 = cell.getStateObj();

  double Tenv{}, Tref{};
  cell.getTemperatures(&Tenv, &Tref); // Cell_SPM.hpp:169 — public getter for T_ref

  // elec_surf: geo.elec_surf = Acell*31 = (0.1*0.2)*31 [Geometry_SPM.hpp:15-17], protected
  // (no getter). Literal reproduced with identical operation order for bit-exact i_app.
  const double elec_surf = 0.1 * 0.2 * 31;
  const double i_app = cell.I() / elec_surf; // == 16.0 / elec_surf

  auto *M = Model_SPM<>::makeModel(); // same singleton the legacy cell holds via `M`

  core::SpectralDiffusionLegacyKernel<NCH> kernel;
  kernel.T_ref = Tref;
  kernel.T = s0.T();
  kernel.i_app = i_app;

  // D_T (electrode Arrhenius activation for D) is a protected Electrode_SPM member with no
  // getter — literals read from the Cell_SPM constructor: pos=29000 (Cell_SPM.hpp:46),
  // neg=35000/5=7000 (Cell_SPM.hpp:55).
  kernel.D_T[pos] = 29000.0;
  kernel.D_T[neg] = 35000.0 / 5.0;

  for (auto dom : { pos, neg }) {
    kernel.D0[dom] = s0.D(dom);       // st.D(dom): 8e-14 (pos) / 7e-14 (neg)
    kernel.a[dom] = s0.a(dom);        // st.a(dom)
    kernel.thick[dom] = s0.thick(dom);// st.thick(dom)
    kernel.sgn[dom] = sign(dom);      // enum_definitions.hpp:65 — pos=-1, neg=+1
    for (int k = 0; k < NCH; ++k) {
      kernel.A[dom][k] = M->A[dom](k);
      kernel.B[dom][k] = M->B[dom](k);
    }
  }

  // ---------------------------------------------------------------------------------------
  // Core side: BatchBuilder with zp/zn rows (NCH each), 1 lane; init from legacy INITIAL z.
  // ---------------------------------------------------------------------------------------
  core::BatchBuilder builder;
  const core::StateSlice zp = builder.declare({ "zp", NCH, core::Unit::none });
  const core::StateSlice zn = builder.declare({ "zn", NCH, core::Unit::none });
  core::StateArena arena = builder.build(1);

  for (int k = 0; k < NCH; ++k) {
    arena.at(zp, k, 0) = s0.zp(k); // legacy initial transformed concentration, pos particle
    arena.at(zn, k, 0) = s0.zn(k); // legacy initial transformed concentration, neg particle
  }

  // ---------------------------------------------------------------------------------------
  // Lockstep integration + drift tracking.
  // ---------------------------------------------------------------------------------------
  double max_abs = 0.0;
  double max_rel = 0.0;
  bool all_finite = true;

  auto compare_all = [&](slide::State_SPM &leg) {
    for (int i = 0; i < NZ; ++i) {
      const double zc = (i < NCH) ? arena.at(zp, i, 0) : arena.at(zn, i - NCH, 0);
      const double zl = leg.z(static_cast<size_t>(i)); // z(i) = [i_zp + i]: zp then zn
      if (!std::isfinite(zc) || !std::isfinite(zl)) all_finite = false;
      const double diff = std::abs(zl - zc);
      max_abs = std::max(max_abs, diff);
      const double denom = std::abs(zl);
      const double rel = (denom > 1e-30) ? (diff / denom) : diff;
      max_rel = std::max(max_rel, rel);
    }
  };

  for (int t = 0; t < nsteps; ++t) {
    cell.timeStep_CC(dt);                    // legacy Euler diffusion step
    kernel.step(arena, zp, zn, 0, dt);       // core replica step
    compare_all(cell.getStateObj());
  }

  // ---------------------------------------------------------------------------------------
  // Verdict against the REGISTERED bands.
  // ---------------------------------------------------------------------------------------
  std::printf("P1-G0: max_abs=%.17g  max_rel=%.17g\n", max_abs, max_rel);
  std::printf("P1-G0: drift_rel=%.17g -> Q8 keep-band %s\n",
              max_rel, (max_rel <= 1e-12 ? "PASS" : "FAIL"));

  REQUIRE(all_finite);                 // sanity: all compared z finite
  CHECK(max_abs == 0.0);               // H0 (primary): bit-identical (non-fatal)
  REQUIRE(max_rel <= 1e-12);           // Q8 keep-band (decisive)
}
