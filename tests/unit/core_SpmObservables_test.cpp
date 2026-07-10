/**
 * @file core_SpmObservables_test.cpp
 * @brief Phase-1 correctness gate for the SURFACE-CONCENTRATION observable
 *        slide::core::SurfaceConcentration<NCH> (PLAN.md §3.7/§3.12, D-10).
 *
 * The observable reconstructs surface Li-concentration from arena z-rows via the Model_SPM output
 * maps (the C/D surface row — the nch!=5 bug path, §2.2). Validated TWO independent ways:
 *
 * ===========================================================================================
 * REGISTERED BANDS (written BEFORE first run, CLAUDE.md §3):
 *
 *   PARITY (vs legacy Cell_SPM::getCSurf): a 1-lane batch built from a real Cell_SPM's exact
 *       state reproduces cell.getCSurf bit-for-bit — the compute() operation sequence is the
 *       legacy calcSurfaceConcentration sequence. REQUIRE rel <= 1e-12 (both configs; Debug/-O0
 *       expected exactly 0). OUTCOME 2026-07-10: rel = 0 for both electrodes (bit-identical).
 *
 *   ROUND-TRIP (analytic, heterogeneous batch, independent of legacy): with ZERO current and z
 *       set to a UNIFORM concentration per lane (distinct conc per lane), the surface value must
 *       equal that uniform concentration -> c_surf(c) == conc(c). REQUIRE rel <= 1e-8 (matches
 *       the Chebyshev round-trip tolerance). Non-degeneracy: lane concentrations are distinct.
 *       OUTCOME 2026-07-10: max rel = 1.276e-15; positive-electrode lane spread = 7193.9.
 * ===========================================================================================
 *
 * @date 2026-07-10
 */

#include "../../src/slide.hpp"                 // legacy Cell_SPM, Model_SPM, settings, sign, Domain, pos/neg
#include "../../src/core/BatchBuilder.hpp"
#include "../../src/core/SpmObservables.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace slide;

namespace {

//!< Capture the batch-shared surface-concentration params from the same sources as the P1-G0
//!< pilot / SpectralDiffusion test (protected legacy members reproduced by literal, documented).
template <int NCH>
core::SurfaceConcentrationParams<NCH> captureParams(Cell_SPM &cell, Model_SPM<> &M, double Tref)
{
  auto &s0 = cell.getStateObj();
  core::SurfaceConcentrationParams<NCH> p;
  p.T_ref = Tref;
  p.elec_surf = 0.1 * 0.2 * 31; // geo.elec_surf = Acell*31 (Geometry_SPM.hpp:15-17), protected
  p.D_T[pos] = 29000.0;         // Cell_SPM.hpp:46
  p.D_T[neg] = 35000.0 / 5.0;   // Cell_SPM.hpp:55
  for (auto dom : { pos, neg }) {
    p.D0[dom] = s0.D(dom);
    p.a[dom] = s0.a(dom);
    p.thick[dom] = s0.thick(dom);
    p.sgn[dom] = sign(dom);
    p.Dout0[dom] = M.D[dom](0);
    for (int j = 0; j < NCH; ++j)
      p.C0[dom][j] = M.C[dom](0, j);
  }
  return p;
}

} // namespace

TEST_CASE("SurfaceConcentration<NCH> reproduces legacy Cell_SPM::getCSurf (parity)", "[core][observables]")
{
  constexpr int NCH = static_cast<int>(settings::nch);

  Cell_SPM cell; // 16 Ah Kokam NMC, SOC 0.5
  cell.setBlockDegAndTherm(true);
  REQUIRE(cell.setCurrent(20.0, false, false) == Status::Success); // nonzero -> flux term exercised

  double Tenv{}, Tref{};
  cell.getTemperatures(&Tenv, &Tref);
  auto *M = Model_SPM<>::makeModel();
  auto p = captureParams<NCH>(cell, *M, Tref);

  // Legacy surface concentration (the oracle).
  DPair cs_legacy;
  cell.getCSurf(cs_legacy, false);

  // 1-lane batch built from the cell's exact state.
  auto &s0 = cell.getStateObj();
  core::BatchBuilder builder;
  const core::StateSlice zp = builder.declare({ "zp", NCH, core::Unit::none });
  const core::StateSlice zn = builder.declare({ "zn", NCH, core::Unit::none });
  const core::StateSlice Ts = builder.declare({ "T", 1, core::Unit::K });
  core::StateArena arena = builder.build(1);
  for (int k = 0; k < NCH; ++k) {
    arena.at(zp, k, 0) = s0.zp(k);
    arena.at(zn, k, 0) = s0.zn(k);
  }
  arena.at(Ts, 0, 0) = s0.T();

  const std::array<double, 1> iapp{ cell.I() / p.elec_surf };
  std::array<double, 1> csp{}, csn{};
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  const core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };

  core::computeSurfaceConcentration(p, state, zp, zn, Ts, ctx,
                                    std::span<core::real_t>{ csp },
                                    std::span<core::real_t>{ csn });

  const double relp = std::abs(csp[0] - cs_legacy[pos]) / std::abs(cs_legacy[pos]);
  const double reln = std::abs(csn[0] - cs_legacy[neg]) / std::abs(cs_legacy[neg]);
  std::printf("SpmObservables parity: csp=%.12g (legacy %.12g, rel=%.3e)  csn=%.12g (legacy %.12g, rel=%.3e)\n",
              csp[0], cs_legacy[pos], relp, csn[0], cs_legacy[neg], reln);

  CAPTURE(csp[0], cs_legacy[pos], csn[0], cs_legacy[neg]);
  REQUIRE(std::isfinite(csp[0]));
  REQUIRE(std::isfinite(csn[0]));
  REQUIRE(relp <= 1e-12);
  REQUIRE(reln <= 1e-12);
}

TEST_CASE("SurfaceConcentration<NCH> uniform round-trip on a heterogeneous batch", "[core][observables]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr int L = 8;

  Cell_SPM cell;
  cell.setBlockDegAndTherm(true);
  REQUIRE(cell.setCurrent(0.0, false, false) == Status::Success); // zero current -> flux term = 0
  double Tenv{}, Tref{};
  cell.getTemperatures(&Tenv, &Tref);
  auto *M = Model_SPM<>::makeModel();
  auto p = captureParams<NCH>(cell, *M, Tref);

  const std::array<double, 2> R{ M->Rp, M->Rn };
  const std::array<double, 2> Cmax{ 51385.0, 30555.0 }; // pos/neg (Cell_SPM.hpp:47,56)

  core::BatchBuilder builder;
  const core::StateSlice zp = builder.declare({ "zp", NCH, core::Unit::none });
  const core::StateSlice zn = builder.declare({ "zn", NCH, core::Unit::none });
  const core::StateSlice Ts = builder.declare({ "T", 1, core::Unit::K });
  core::StateArena arena = builder.build(L);
  const std::array<core::StateSlice, 2> slice{ zp, zn };

  // Per-lane distinct uniform concentration (fraction 0.30 .. 0.44 of Cmax), set exactly as the
  // legacy setC: z[zero] = sum_i V[dom](zero,i) * (R*conc*xch(i)); other modes 0.
  std::array<std::array<double, L>, 2> conc{};
  for (int d = 0; d < 2; ++d)
    for (int c = 0; c < L; ++c) {
      const double frac = 0.30 + 0.02 * c;
      conc[d][c] = frac * Cmax[d];
      double zu = 0.0;
      for (int i = 0; i < NCH; ++i)
        zu += M->V[d](M->zero, i) * (R[d] * conc[d][c] * M->xch(i));
      for (int k = 0; k < NCH; ++k)
        arena.at(slice[d], k, c) = (k == M->zero) ? zu : 0.0;
    }

  std::array<double, L> iapp{};
  for (int c = 0; c < L; ++c)
    arena.at(Ts, 0, c) = cell.T() + 1.5 * c;

  std::array<double, L> csp{}, csn{};
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  const core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };
  core::computeSurfaceConcentration(p, state, zp, zn, Ts, ctx,
                                    std::span<core::real_t>{ csp },
                                    std::span<core::real_t>{ csn });

  double max_rel = 0.0, spread = 0.0;
  for (int c = 0; c < L; ++c) {
    const double rp = std::abs(csp[c] - conc[pos][c]) / conc[pos][c];
    const double rn = std::abs(csn[c] - conc[neg][c]) / conc[neg][c];
    max_rel = std::max({ max_rel, rp, rn });
    spread = std::max(spread, std::abs(csp[c] - csp[0]));
  }
  std::printf("SpmObservables round-trip: max_rel=%.3e  lane_spread(csp)=%.6g\n", max_rel, spread);

  CAPTURE(max_rel, spread);
  REQUIRE(spread > 0.0);      // non-degeneracy: lanes genuinely differ
  REQUIRE(max_rel <= 1e-8);   // uniform concentration recovered at the surface
}
