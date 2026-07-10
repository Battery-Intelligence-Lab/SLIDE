# Handoff — 2026-07-07 — Phase 0C COMPLETE (full ctest green)

## Read first
1. `/PLAN.md` — living contract; §8 ledger updated with the Phase 0C DONE row (root causes + evidence).
2. Previous handoff: `handoff-2026-07-07-v4-refactor-kickoff.md` (context for Phase 0 and rules).

## What happened this session (branch `Claude`, NOT pushed)
The "in flight" Phase 0C agent from the previous session had **died mid-configure**
(build-phase0c had no CMakeCache; only the fmt edit in `cmake/Dependencies.cmake`
survived, uncommitted). This session verified and completed everything:

Commits `8f9efed`…`207a3c2` (6):
- `8f9efed` P0-C1: fmt 11.0.2→11.2.0 (clang≥21 consteval blocker). Verified: full 1144-target build green.
- `0259d1c` P0-C2: default `ocv_coefs` (−55 782 V at every SOC) → empty default, `getOCV()` falls back to
  OCV-table interp; polynomial opt-in via `set_ocv_coefs`. KEY EVIDENCE: the full 8-term fit (found in
  b0a1c82's integration test, same commit that introduced the truncated default) spans only [2.53, 3.46] V —
  no polynomial default can span the 2.7–4.2 V window, so table interp is the only green-suite option.
  Test expectations 3.15→3.45 V (3.15 was for the DELETED standalone Cell_Bucket, 2.0–4.3 V ramp).
- `44fcfb5` P0-C3: VERDICT test-stale, not param regression — thickp/thickn recalibrated intentionally in
  147ee4c (Volkan, Aug 2024, changelog'd). CSurf rederived: csurf = x_init·Cmax → 35562.14 / 14694.92,
  closed form matched actuals to all printed digits.
- `5d4b12f` P0-C5+C6: dead `v_now` removed (private method, 1 internal caller); CV energy trapezoid
  (consistent with B3) **+ th.time() was NEVER accumulated in CV** (CV/CCCV time throughput read 0 s) — found
  en route, fixed, covered by registered test (pre-fix: time=0 FAILED, Wh=3.49986 missed 3.5±1e-6 band).
- `9cf7dc9` Module_p test: SPM SOC band 1e-15→5e-5 (SPM SOC = Li-fraction estimate since v3, not coulomb
  count; assertion previously masked by P0-C2 baseline failures). ECM keeps 1e-15.
- `207a3c2` docs: CHANGELOG Unreleased + PLAN.md §8.

## Verification discipline record
- Baseline registered BEFORE run: predicted 4/10 failing {Cell_Bucket, Cell_ECM, Cell_SPM, Module_p} — HIT exactly.
- P0-C6 test written + run pre-fix (failed as registered), then post-fix (all bands hit).
- Final gate: **ctest 10/10 GREEN** in `build-verify` (clang 21.1.8, Ninja, Debug). Delta: 4 failing → 0.
- Two NEW stale assertions surfaced only after P0-C2/C3 unmasked deeper test code (CSurf, SPM SOC) — both
  root-caused to the same 147ee4c/v3 recalibration lineage, not to this session's changes.

## Still open / next steps
1. **Volkan reviews PLAN.md §3/§4/§7 (Q1–Q7)** — Phase 1 remains GATED on this.
2. Phase 1 (PLAN.md §6): StateArena/StateSpec/BatchBuilder, Domain/ElectrodeDesign packs,
   SpectralDiffusion<NCH>, kernels, registry, legacy-Euler parity. Gates P1-G1..G3.
3. P0-C4 (static-sizing in dead Boost path) intentionally left for v4 (dead code).
4. Branch `Claude` still NOT pushed — user's call.
5. Working tree keeps user-owned changes (`.gitignore` modified, `.claude/settings.json` deleted,
   `.claude/reports/*`, `.claude/FABLE.md` untracked) — left alone again.

## Notes
- Advisor (Fable) was rate-limited/unavailable this session; decisions made on primary evidence, logged above.
- `build-verify/` is the current good build dir (ignored via user's `/build-*` gitignore line, uncommitted).
