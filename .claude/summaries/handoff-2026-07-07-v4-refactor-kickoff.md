# Handoff — 2026-07-07 — SLIDE v4 refactor kickoff (Fable session end: quota)

## Read first
1. **`/PLAN.md`** — THE living contract: goals, performance contract §1.1, audit evidence §2 (+§2.5 Phase-0
   outcomes), target architecture §3 (incl. §3.10 coding conventions, §3.11 COMSOL-style parameter system),
   decision log §4 (D-01…D-17, do not re-litigate), migration/verification discipline §5, phased roadmap §6,
   open questions §7, status ledger §8.
2. `.claude/FABLE.md` — Volkan's original requirements (untracked, user's file).
3. Operating rules (also in persistent memory): Fable/architect plans, Opus implements; ≤3 parallel agents;
   PLAN.md updated after every decision; wall-clock benchmarks UNRELIABLE on this machine (parallel jobs).

## State at handoff (branch `Claude`, NOT pushed)
- Commits `7529039`…`91faaf9`: Phase-0 fixes A1–A7 (solver/state) + B1–B5 (data/docs), each with a
  regression test that failed pre-fix. A6 was REFUTED (non-bug, proved by test). `11ac335`: docs (PLAN.md,
  CHANGELOG consolidation, TODO/discussions pointers).
- Working tree also has user-owned changes NOT mine: `.gitignore` modified, `.claude/settings.json` deleted,
  `.claude/reports/*` untracked. Leave them alone.

## IN FLIGHT — Phase 0C agent (Opus, background) may still be running or just finished
- Task: fmt≥11.1 bump in `cmake/Dependencies.cmake` (clang-21 build blocker — tree does NOT compile without
  it or a build-dir patch); root-cause+fix default `ocv_coefs` in `src/cells/Cell_ECM/Cell_ECM.hpp:48`
  (OCV = −55782 V ⇒ Cell_Bucket/Cell_ECM/Module_p failing); `thickp` 70e-6 vs 8.687e-5 test/param mismatch;
  dead `v_now` out-param in `Cycler::setCurrent`; CV/CCCV trapezoid energy consistency.
- Registered target (set before run): ALL 6 test binaries green (baseline: 4 failing pre-existing).
- NEXT SESSION: check `git log` for its commits; verify its report claims by running
  `cmake -B build-verify -DCMAKE_BUILD_TYPE=Debug && cmake --build build-verify && ctest --test-dir build-verify`;
  then update PLAN.md §8 and CHANGELOG if it didn't.

## Next steps (in order)
1. Verify Phase 0C (above). Treat its "COMPLETE" as hypothesis; open the cited evidence.
2. **Volkan reviews PLAN.md §3/§4/§7** (Q1–Q7; esp. Q6: is Nilsu's `setCurrent_analytical_impl` the "Ross
   analytical solution"?). Phase 1 is gated on this review.
3. Phase 1 (PLAN.md §6): core data model — StateArena/StateSpec/BatchBuilder, Domain/ElectrodeDesign packs,
   SpectralDiffusion<NCH>, kernels, registry, legacy-Euler parity mode. Gates P1-G1..G3 (digit-diff parity,
   zero per-step allocations, Chebyshev analytic oracle). Delegate to Opus 4.8 (xhigh via prompt), ≤3 agents,
   disjoint files, own build dirs (build-phaseXY pattern — shared bin/ output dir caused one flake).

## Discipline reminders that mattered this session
- Registered bands BEFORE decisive runs; baseline ctest recorded before fixes; per-fix regression test that
  fails pre-fix; "no-op" claims verified by digit-identical reruns (B5 did revert-arbiter correctly).
- Subagent reports contained real contradictions of the audit (A6 non-bug; B3 worse than audited; B4 swap was
  in the .cpp not DEG_ID.hpp) — verifying against code paid off every time.
