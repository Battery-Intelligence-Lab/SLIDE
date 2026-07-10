# Handoff — 2026-07-08 — Phase 1: Q8 Release re-confirm + production SpectralDiffusion<NCH>

## Read first
1. `/PLAN.md` §6 (Phase 1 roadmap), §7 Q8 (now FULLY closed), §8 ledger (bottom rows current).
2. Previous handoff: `handoff-2026-07-07-phase1-kickoff.md`.

## What happened this session (branch `Claude`, NOT pushed) — 3 commits

### 1. P1-G0 Release re-confirm — Q8 standing condition DISCHARGED (commit `7216b05`)
The Q8 resolution left a standing condition: the P1-G0 pilot ran only Debug/-O0; re-confirm
drift in Release before P1-G1 sign-off. Done.
- New build dir `build-release` (clang-21, Ninja, `-O3`; note the tree carries `-Ofast`/`-ffast-math`).
- **Debug/-O0:** max_abs = 0, max_rel = 0 — H0 (exact bit-identity) CONFIRMED.
- **Release/-O3:** max_abs = **2.26e-17**, max_rel = **5.63e-15** — H0 FALSIFIED (cross-TU FMA
  contraction: legacy update lives in the prebuilt `src` lib, kernel is header-only in the test
  TU; `-Ofast` contracts them differently); **decisive Q8 band rel ≤ 1e-12 HOLDS** (~3 orders margin).
- `-ffp-contract=off` NOT applied: on the test target alone it can't reach 0 (legacy `src` stays
  contracted), and forcing it globally would recompile legacy — PLAN §5.1 forbids that. Bit-identity
  was a bonus hypothesis; the decisive band was always rel ≤ 1e-12. Recorded FALSIFIED with numbers
  (CLAUDE.md §3), not tuned.
- Accumulation bound (derived): modal Euler map `z_k ← (1+dt·D·A_k)·z_k + dt·B_k·j` is non-expansive
  for stable modes ⇒ per-step roundoff damped; only the mean mode (A_0≈0) accumulates ~linearly,
  ~1e-15 over a 10-cycle run — ≥3 orders below the gate.
- Test change: `max_abs==0` CHECK scoped `#ifndef NDEBUG` (holds at -O0); `rel≤1e-12` REQUIRE gates
  CI in both configs. **Q8 FULLY CLOSED; P1-G1 unblocked.**

### 2. Production `SpectralDiffusion<NCH>` kernel (commit `3dffe80`)
`src/core/SpectralDiffusion.hpp` — vectorised-across-lanes forward-Euler diffusion on SoA
`StateArena` rows (one variable across all lanes; inner loop over lanes = SIMD sweep, PC-3; one
kernel call per batch per step, PC-2).
- `DiffusionParams<NCH>` = batch-shared A/B (state-space diagonals from Model_SPM) + D0/D_T/a/thick/sgn
  + F/Rg/n/T_ref. Per-lane operating point (T[c], i_app[c]) passed as spans. Owns a once-allocated
  per-lane scratch (`D_eff`, `flux`) → zero per-step heap alloc (PC-1).
- Legacy-Euler stepping mode only (Phase 1). The exact **exponential modal propagator (D-07)** is
  Phase 3, on the SAME layout/params — intentionally not here.
- Inner form is **two statements** (`dz = D·A·z + B·j; z += dt·dz`) to mirror the legacy-shaped
  kernel's op sequence.
- Validation `tests/unit/core_SpectralDiffusion_test.cpp`: heterogeneous 8-lane batch (distinct T,
  i_app, initial z per lane), 600 steps, vs `SpectralDiffusionLegacyKernel` per lane.
  - **Debug: max_abs = 0** (H_math arbiter — -O0 disables FMA+vectorisation, so exact-zero proves
    the vectorised sweep is the SAME operation sequence, i.e. no reassociation bug).
  - **Release: max_abs = 5.90e-17, max_rel = 3.80e-15** — decisive rel ≤ 1e-12 HOLDS (~3 orders).
    Sub-ulp FMA/vectorisation reassoc, EXPECTED for a vectorised kernel (this is exactly why Q8
    validates production kernels by a rel band, not the digit-diff — §7 Q8 standing condition 1).
  - Non-degeneracy guard: lane-to-lane spread of final zp[0] > 0.
- ctest Debug (`build-verify`) 13/13 green (was 12; +1 intended, no regressions).

### 3. Ledger + CHANGELOG (commit `64b3b81`)
PLAN §8 + CHANGELOG Unreleased updated for both deliverables.

## Verification discipline record
- All bands registered pre-run. P1-G0 Release: H0 falsified (recorded with numbers), decisive band
  hit. SpectralDiffusion: H_math (Debug) hit exactly; Release rel band hit; over-registered exact-zero
  for a vectorised kernel corrected to the ALREADY-registered Q8 rel-band policy (not post-hoc tuning
  — the Debug arbiter independently proved math identity across 48000 comparisons).
- advisor tool unavailable this session; Fable agent hit its usage-credit limit → NO external review
  obtained. Decisions above rest on the Debug bit-identity arbiter + derived accumulation bounds.

## BLOCKER for the next kernels — observable-reconstruction layer needs design (surface to Volkan)
ThermalLumped and the ageing kernels are NOT self-contained state→state kernels like diffusion.
`Cell_SPM::dState_thermal` (`src/cells/Cell_SPM/Cell_SPM_dstate.cpp:57-113`) needs:
- **surface concentration** `c_surf = C·z + D·flux` + centre-node (`Model_SPM` C/D matrices, cc_coeff —
  the same Chebyshev output path where the nch≠5 bug lived, §2.2);
- **overpotentials** η (Butler-Volmer, `Electrode_SPM::overpotential`);
- **OCV entropic coefficient** dOCV (table interp) and **Rdc**.
These are the "observables derived from state" layer (§3.7, D-10) + the `BatchView`/`StepCtx` kernel
interface (§3.11) that this session DEFERRED (SpectralDiffusion takes plain spans as a stopgap).
Ageing, the composition registry, and the Simulation façade all converge on this same layer.

**Recommended next step:** design the observable-reconstruction + BatchView layer (Fable-level;
needs review) BEFORE porting ThermalLumped/ageing. Natural companion deliverable: **P1-G3** (Chebyshev
external oracle vs Carslaw & Jaeger analytic transient-sphere-diffusion series) — it exercises exactly
the C/D surface-concentration output path and is a registered Phase-1 gate. Doing P1-G3 first would
validate the output reconstruction the thermal/ageing kernels then consume.

## Still open / next steps
1. **Design observable layer + BatchView/StepCtx** (§3.7/§3.11) — needs Volkan/Fable review.
2. P1-G3 Chebyshev oracle (Carslaw & Jaeger) — validates c_surf output path; registered gate.
3. ThermalLumped (consumes q_ext seam) — after (1).
4. Ageing kernels (SEI/LAM/CS/plating) mechanism-by-mechanism — after (1).
5. Composition registry + factory; single-cell `Simulation` façade; P1-G1 full parity, P1-G2 zero-alloc 1e4 batch.
6. PAY-1 checkpoint at Phase-1 exit (Volkan, quiet machine, 1e4 cells; abort <2×).
7. D-21 pack-thermal design before Phase 2. Branch `Claude` still NOT pushed — user's call.

## Operating notes
- `build-verify` = full-suite Debug gate (13/13). `build-release` = FMA re-check dir; only the two
  core kernel tests built there (P1-G0, SpectralDiffusion) — legacy test exes show "Not Run" there.
- ≤2 agents; Opus HIGH implements; Fable reviews (OUT OF CREDITS this session — flag to Volkan).
- User-owned working-tree files (`.gitignore` M, `.claude/settings.json` D, `.claude/reports/*`,
  `.claude/FABLE.md`, handoffs) — leave alone.
