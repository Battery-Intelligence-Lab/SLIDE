/**
 * @file core_SpectralDiffusion_test.cpp
 * @brief Phase-1 correctness gate for the PRODUCTION diffusion kernel
 *        slide::core::SpectralDiffusion<NCH> (PLAN.md §3.2, §3.5).
 *
 * The production kernel sweeps whole SoA state rows across all lanes of one batch. Its
 * Phase-1 oracle is slide::core::SpectralDiffusionLegacyKernel — the op-order replica already
 * validated bit-identical to legacy Cell_SPM in Debug and within rel ≤ 1e-12 of it in Release
 * (P1-G0, §7 Q8). The production kernel executes the identical per-lane operation SEQUENCE
 * (compute dz, then Euler-advance) as the legacy-shaped kernel; both are compiled in THIS TU.
 *
 * ===========================================================================================
 * REGISTERED BANDS — written BEFORE the first run (CLAUDE.md §3); OUTCOME recorded post-run.
 *
 *   H_math (Debug arbiter):  in Debug/-O0, max |z_prod(c) − z_legacyShaped(c)| == 0.0 EXACTLY,
 *                       over all lanes / 2·nch modes / steps -> REQUIRE (fatal). -O0 disables
 *                       FMA contraction and auto-vectorisation, so exact zero here proves the
 *                       vectorised sweep is the SAME operation sequence as the legacy math — a
 *                       reassociation-by-parenthesisation bug would surface at -O0 too.
 *                       OUTCOME 2026-07-08: max_abs = 0 -> CONFIRMED (math identity established).
 *
 *   Q8 band (decisive, both configs):  max relative drift <= 1e-12 -> REQUIRE (fatal). This is
 *                       the band Q8 already fixed for PRODUCTION VECTORISED kernels (they are
 *                       validated by a rel band, NOT the digit-diff — §7 Q8 standing condition
 *                       1). Under -O3/-Ofast the vectorised FMA sweep reassociates vs the scalar
 *                       oracle at ~0.25 ulp/step; the Debug arbiter above rules that a rounding
 *                       artefact (not a bug), and this band bounds it.
 *                       OUTCOME 2026-07-08: Release max_abs = 5.9e-17, max_rel = 3.8e-15 ->
 *                       HOLDS with ~3 orders of margin. (Exact bit-identity is NOT expected at
 *                       -O3 for a vectorised kernel — that was why Q8 chose a rel band here.)
 *
 *   Non-degeneracy:     lanes carry distinct T[c], i_app[c] and initial z; the lane-to-lane
 *                       spread of the final zp[0] must be > 0 (the sweep genuinely runs over
 *                       differing operating points, not N copies of one lane).
 *
 *   Sanity:             every compared z is std::isfinite.
 * ===========================================================================================
 *
 * @date 2026-07-08
 */

#include "../../src/slide.hpp" // legacy Cell_SPM, Model_SPM, settings, sign, Domain
#include "../../src/core/BatchBuilder.hpp"
#include "../../src/core/SpectralDiffusion.hpp"       // production kernel under test
#include "../../src/core/SpectralDiffusionLegacy.hpp" // per-lane oracle (legacy-validated)
#include "../support/RecordedBits.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using namespace slide;

TEST_CASE("SpectralDiffusion rejects invalid lane counts before allocation",
          "[core][diffusion][validation]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  const core::DiffusionParams<NCH> parameters;
  REQUIRE_THROWS_AS(core::SpectralDiffusion<NCH>(parameters, 0),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(core::SpectralDiffusion<NCH>(parameters, -1),
                    std::invalid_argument);
}

TEST_CASE("SpectralDiffusion<NCH> production kernel == legacy-shaped kernel per lane", "[core][diffusion]")
{
  constexpr int NCH = static_cast<int>(settings::nch); // 5
  constexpr int L = 8;                                 // lanes in the batch
  constexpr double dt = 1.0;
  constexpr int nsteps = 600;

  // ---------------------------------------------------------------------------------------
  // Shared model params, captured from the same sources as the P1-G0 pilot.
  // ---------------------------------------------------------------------------------------
  Cell_SPM cell; // 16 Ah Kokam NMC, SOC 0.5, T = settings::T_ENV
  cell.setBlockDegAndTherm(true);
  REQUIRE(cell.setCurrent(16.0, false, false) == Status::Success);

  auto &s0 = cell.getStateObj();
  double Tenv{}, Tref{};
  cell.getTemperatures(&Tenv, &Tref);
  const double elec_surf = 0.1 * 0.2 * 31; // geo.elec_surf, reproduced (protected in legacy)
  auto *M = Model_SPM<>::makeModel();

  core::DiffusionParams<NCH> p;
  p.T_ref = Tref;
  const auto core_index = [](slide::Domain domain) {
    return core::domain_index(domain == pos ? core::Domain::pos : core::Domain::neg);
  };
  p.D_T[core_index(pos)] = 29000.0;       // Cell_SPM.hpp:46
  p.D_T[core_index(neg)] = 35000.0 / 5.0; // Cell_SPM.hpp:55
  for (auto dom : { pos, neg }) {
    const auto d = core_index(dom);
    p.D0[d] = s0.D(dom);
    p.a[d] = s0.a(dom);
    p.thick[d] = s0.thick(dom);
    p.sgn[d] = sign(dom);
    for (int k = 0; k < NCH; ++k) {
      p.A[d][k] = M->A[dom](k);
      p.B[d][k] = M->B[dom](k);
    }
  }

  // ---------------------------------------------------------------------------------------
  // Per-lane heterogeneous operating point + distinct initial state.
  // ---------------------------------------------------------------------------------------
  std::array<double, L> Tl{}, iappl{};
  for (int c = 0; c < L; ++c) {
    Tl[c] = s0.T() + 2.0 * c;                // distinct temperatures (298.15 … +14 K)
    iappl[c] = (16.0 + 1.0 * c) / elec_surf; // distinct current densities
  }

  core::BatchBuilder builder;
  const core::StateSlice zp = builder.declare({ "zp", NCH, core::Unit::none });
  const core::StateSlice zn = builder.declare({ "zn", NCH, core::Unit::none });
  core::StateArena arena = builder.build(L);
  for (int c = 0; c < L; ++c)
    for (int k = 0; k < NCH; ++k) {
      arena.at(zp, k, c) = s0.zp(k) + 1e-4 * c; // distinct initial z per lane
      arena.at(zn, k, c) = s0.zn(k) - 1e-4 * c;
    }

  core::SpectralDiffusion<NCH> prod(p, L);

  // ---------------------------------------------------------------------------------------
  // Per-lane oracle: one legacy-shaped kernel + one 1-lane arena per lane, same initial z.
  // ---------------------------------------------------------------------------------------
  std::vector<core::StateArena> refArena;
  refArena.reserve(L);
  std::array<core::SpectralDiffusionLegacyKernel<NCH>, L> refKern{};
  for (int c = 0; c < L; ++c) {
    core::BatchBuilder rb;
    const core::StateSlice rzp = rb.declare({ "zp", NCH, core::Unit::none });
    const core::StateSlice rzn = rb.declare({ "zn", NCH, core::Unit::none });
    (void)rzp;
    (void)rzn; // same layout as the production batch (zp then zn)
    refArena.push_back(rb.build(1));
    for (int k = 0; k < NCH; ++k) {
      refArena[c].at(zp, k, 0) = arena.at(zp, k, c);
      refArena[c].at(zn, k, 0) = arena.at(zn, k, c);
    }
    auto &rk = refKern[c];
    rk.T_ref = Tref;
    rk.T = Tl[c];
    rk.i_app = iappl[c];
    rk.D_T[pos] = p.D_T[core_index(pos)];
    rk.D_T[neg] = p.D_T[core_index(neg)];
    for (auto dom : { pos, neg }) {
      const auto d = core_index(dom);
      rk.D0[dom] = p.D0[d];
      rk.a[dom] = p.a[d];
      rk.thick[dom] = p.thick[d];
      rk.sgn[dom] = p.sgn[d];
      for (int k = 0; k < NCH; ++k) {
        rk.A[dom][k] = p.A[d][k];
        rk.B[dom][k] = p.B[d][k];
      }
    }
  }

  // ---------------------------------------------------------------------------------------
  // Lockstep integration + exact drift tracking.
  // ---------------------------------------------------------------------------------------
  double max_abs = 0.0;
  double max_rel = 0.0;
  bool all_finite = true;
  const std::span<const core::real_t> Tspan{ Tl.data(), L };
  const std::span<const core::real_t> ispan{ iappl.data(), L };
  test_support::RecordedBits recorded;
  recorded.append(arena.raw());

  auto track = [&](double zc, double zr) {
    if (!std::isfinite(zc) || !std::isfinite(zr)) all_finite = false;
    const double diff = std::abs(zc - zr);
    max_abs = std::max(max_abs, diff);
    const double denom = std::abs(zr);
    max_rel = std::max(max_rel, (denom > 1e-30) ? (diff / denom) : diff);
  };

  for (int t = 0; t < nsteps; ++t) {
    prod.stepEuler(arena, zp, zn, Tspan, ispan, dt);
    recorded.append(arena.raw());
    for (int c = 0; c < L; ++c)
      refKern[c].step(refArena[c], zp, zn, 0, dt);

    for (int c = 0; c < L; ++c)
      for (int k = 0; k < NCH; ++k) {
        track(arena.at(zp, k, c), refArena[c].at(zp, k, 0));
        track(arena.at(zn, k, c), refArena[c].at(zn, k, 0));
      }
  }

  // Non-degeneracy: lanes must have diverged (heterogeneous operating point genuinely exercised).
  double zmin = arena.at(zp, 0, 0), zmax = zmin;
  for (int c = 0; c < L; ++c) {
    zmin = std::min(zmin, arena.at(zp, 0, c));
    zmax = std::max(zmax, arena.at(zp, 0, c));
  }
  const double lane_spread = zmax - zmin;

  std::printf("SpectralDiffusion: max_abs(prod-legacyShaped)=%.17g  max_rel=%.17g  lane_spread=%.17g\n",
              max_abs,
              max_rel,
              lane_spread);

  REQUIRE(all_finite);
  REQUIRE(lane_spread > 0.0); // non-degeneracy: the batch is genuinely heterogeneous
#ifndef NDEBUG
  REQUIRE(max_abs == 0.0); // H_math (Debug arbiter): vectorised sweep == legacy math exactly
#else
  std::printf(
    "SpectralDiffusion: Release/-O3 -> exact bit-identity not expected for a vectorised "
    "kernel; decisive rel band is the gate (see header, §7 Q8)\n");
#endif
  REQUIRE(max_rel <= 1e-12); // Q8 decisive band for production vectorised kernels (both configs)

  CAPTURE(recorded.values, recorded.fnv1a, recorded.mixed);
  REQUIRE(recorded.values == 48080);
#if defined(SLIDE_TEST_HAS_RECORDED_SCALAR_BITS)
#if defined(SLIDE_TEST_IPO) && defined(__FAST_MATH__)
  constexpr auto expected_fnv = UINT64_C(0x785b850afc76b768);
  constexpr auto expected_mixed = UINT64_C(0x3a0476924499700f);
#elif defined(__FAST_MATH__)
  constexpr auto expected_fnv = UINT64_C(0x8f7f609c10bb3d1a);
  constexpr auto expected_mixed = UINT64_C(0x8332970bbbba9a67);
#else
  constexpr auto expected_fnv = UINT64_C(0xd2601e8e7d13249b);
  constexpr auto expected_mixed = UINT64_C(0x5eabfd035f26594c);
#endif
  CHECK(recorded.fnv1a == expected_fnv);
  CHECK(recorded.mixed == expected_mixed);
#endif
}
