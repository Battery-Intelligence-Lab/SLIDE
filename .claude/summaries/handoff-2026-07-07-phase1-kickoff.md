# Handoff — 2026-07-07 — Phase 1 kickoff (Q1–Q9 decided, keystone + P1-G0 done)

## Read first
1. `/PLAN.md` — §7 all questions now DECIDED/RESOLVED except none; §8 ledger current.
2. Previous handoff: `handoff-2026-07-07-phase0c-complete.md`.

## What happened this session (branch `Claude`, NOT pushed)

### 1. Fable review + Volkan sign-off (commit `c361279`)
- advisor() endpoint kept 429ing → fallback: Agent tool with `model: "fable"` WORKS (memory note saved:
  `advisor-fable-fallback.md`).
- Fable reviewed §3/§4/§7: all six ASSUMED defaults survive; Volkan accepted all (Q1 f64-only, Q2 keep
  ECM<NRC>, Q3 façade→v5, Q4 KLU optional, Q5 CUDA-first, Q7 PyBaMM latest-at-P7).
- Two NEW gaps found outside the Q-ledger, logged:
  - **Q8** parity band vs FP reassociation → pilot-first rule (Volkan's pick), gate P1-G0 added.
  - **Q9** pack thermal coupling missing from compile() → seam reserved in Phase 1, **D-21 design due
    BEFORE Phase 2 opens** (P2-G1 unpassable until then).
  - R3 archetype-fragmentation assumption written into §3.1 + compile() diagnostic requirement.
- Phase 0C spot-checks passed: `Cell_ECM::SOC()` == `st.SOC()` (Cell_ECM.hpp:77); CHANGELOG P0-C2 complete.

### 2. Phase 1 keystone (commit `471f60a`)
- `src/core/StateArena.hpp` + `BatchBuilder.hpp` (namespace `slide::core`, header-only INTERFACE lib
  `slide_core`; legacy `src` target does NOT link it — strangler §5.1).
- `real_t = double` alias (Q1 insurance). `reserve_thermal_flux()` = the Q9 seam (one `q_ext` [W] row).
- `tests/unit/core_StateArena_test.cpp`: registered bands (exact offsets, stride/alignment, bit-identical
  snapshot round-trip, ZERO-allocation hot path at 1e4 lanes via replaced global new — P1-G2 groundwork).
- Baseline recorded pre-phase: 10/10. After: 11/11, delta = +core test only.

### 3. P1-G0 parity-drift pilot — Q8 RESOLVED (commit `b0c5847`)
- `src/core/SpectralDiffusionLegacy.hpp`: op-order-faithful replica of legacy diffusion Euler
  (Cell_SPM_dstate.cpp:47-53 + :257-258). PARITY kernel, not production.
- Pilot: 1200×1 s lockstep 1C, per-step compare of all 10 z-modes: **max_abs = 0, max_rel = 0 —
  bit-identical** (H0 registered pre-run, CONFIRMED; verified directly, not just agent-reported).
- §5.2 band KEPT at 1e-12. Cheap because modal update is diagonal (no dot products to reassociate).
- **Standing condition**: pilot ran Debug/-O0 → re-confirm drift==0 in Release config BEFORE P1-G1
  sign-off (FMA contraction may differ; if nonzero → `-ffp-contract=off` on parity targets only).

## Verification discipline record
- All registered pre-run: baseline 10/10 (hit), core test 11/11 (hit), pilot H0 drift==0 (hit).
- Agent results verified against real artifacts (test file read, test re-run directly, constants
  cross-checked constants.hpp:20-21).

## Still open / next steps (Phase 1 remainder, PLAN.md §6)
1. Production `SpectralDiffusion<NCH>` (vectorised across lanes) + `ThermalLumped`.
2. `Domain`/`ElectrodeParams`; ageing kernels ported mechanism-by-mechanism.
3. Composition registry + factory; legacy-Euler stepping mode; single-cell `Simulation` façade.
4. Gates: P1-G1 (single-cell parity, band 1e-12 — Release drift re-check first), P1-G2 (zero-alloc 1e4
   batch — groundwork done), P1-G3 (Chebyshev oracle vs Carslaw & Jaeger).
5. PAY-1 checkpoint at Phase-1 exit: Volkan runs 1e4-cell wall-clock on quiet machine; abort <2×.
6. D-21 pack-thermal design before Phase 2.
7. Branch `Claude` still NOT pushed — user's call.

## Operating notes
- ≤2 agents; Opus HIGH implements; Fable (via Agent tool fallback) reviews/architects.
- build-verify/ is the good Debug build dir. User-owned working-tree files (`.gitignore` M,
  `.claude/settings.json` D, `.claude/reports/*`, `.claude/FABLE.md`, handoffs) — leave alone.
