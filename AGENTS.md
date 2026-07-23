# AGENTS.md — SLIDE standing contract for autonomous agents (Codex and any successor)

**`PLAN.md` at the repo root is the single source of truth.** This file is a compression of its
operating rules for agents that load `AGENTS.md` automatically; where the two ever disagree,
PLAN.md wins. Read order at session start: this file → `PLAN.md` §0 (who you are), §0.3 (the
non-stop protocol), §6 (the ladder — your work queue), §8 (current state, last rows). Full
background requirements live in `.claude/FABLE.md`.

## 1. The loop — you never idle, you never ask

You run non-stop until your context or time is exhausted. There is no user to ask mid-run.

1. Take the **FIRST unticked box in PLAN.md §6** whose dependencies are ticked. Finish it, tick it,
   take the next. Milestone boundaries are not pauses.
2. A box is ticked only when: its gates are green (Debug AND Release, plus CUDA where the change
   touches it), `CHANGELOG.md` (Unreleased) is updated for user-visible changes, a §8 row is added,
   and small reviewable commits are made on branch **`Claude`**.
3. Gate stays red after **three genuinely different attempts** (different hypotheses, not retries)?
   Record FALSIFIED/BLOCKED in §8 with the numbers, leave the box unticked, move to the next
   unblocked box. A falsified registered hypothesis is a deliverable, not a failure.
4. Only three things wait for Volkan: overturning a PC/EC/MC/VC contract row, irreversible or
   outward-facing actions (release tags, publishing, deleting data or untracked files), and writing
   outside this repository. Everything else: decide, record, proceed.
5. Design freedom is yours (§0.3(4)): boxes state WHAT and the gate; the HOW is yours. For boxes
   marked "design note first", write `.claude/designs/<box>.md` (problem, options, choice, why,
   invariants touched), add a §4 D-entry marked ASSUMED, and continue without waiting. Improving
   the ladder — splitting boxes, adding gates, redesigning within the invariants — is expected;
   record every change in §8.
6. Session end = exhaustion only. Then: update box states and §8, commit, and leave a ≤10-line
   handoff in `.claude/summaries/handoff-<date>-<topic>.md` pointing at the next box.

## 2. Evidence discipline — the science bar (PLAN.md §5 in full; this is the summary)

- **Register before you run.** Every decisive test gets its quantitative pass/fail band written down
  BEFORE the run. Unregistered runs are exploration, never evidence.
- **Think first, simulate little.** Derivation is free; simulation costs wall-clock. Every test/gate
  simulation stays SHORT (≤ a few hundred steps, seconds of wall-clock, small lane counts unless the
  gate is specifically about scale). Prefer analytic oracles, conservation identities, and
  structural counters over long trajectories.
- **Wall-clock on this machine is NOT evidence** (parallel jobs; D-27). Judge performance by
  structural counters: allocations, iterations, factorisations, bytes per cell.
- **"No-op" / "refactor only" claims require digit-identical outputs** on recorded cases — never
  code inspection. Record the fixture BEFORE refactoring.
- **Every claim traceable**: to a file:line, a command you ran, an artifact you read, or a stated
  assumption. Quote outputs and errors verbatim. Tag load-bearing claims [confirmed] (name the
  evidence) or [inferred] (name what would confirm it). No nameable artifact → you may not write
  the number.
- **Killed ideas stay killed** (§2.2, §4, the bug ledger's refuted rows). Before starting any
  investigation, search PLAN.md §2/§4/§8 and `.claude/reports/` for whether it was already tested.
  Re-opening requires overturning the recorded kill evidence, not forgetting it.
- **A suspected bug becomes a failing SHORT test before it becomes a fix**; the ledger row records
  test, fix, and post-fix evidence. Refuted candidates are recorded refuted.
- **Every equation unit-checked; every approximation names its regime and leading-order error.**
  Derivations of lasting value go in `docs/derivations/` complete enough to re-derive from scratch.
- **Citations are web-verified at the box that registers them** — a DOI is unverified until you
  resolved it to the named paper. Never cite from memory.
- **After substantial work, run your own adversarial pass**: the counterexample, the missing case,
  the hidden assumption, the unit slip. Test gates are code too — try to make your own gate pass
  wrongly (mutation checks) before trusting it.

## 3. Hard rules (violations are never acceptable)

- Never write outside the repository root. `data/` and `results/` are read-only.
- Never delete untracked files or data — propose deletions in a report instead.
- No runtime dependence on repo-relative paths in shipped code.
- Python dependencies via `uv add`/`uv remove` only.
- Optional dependencies (OpenMP, CUDA, KLU, …) stay optional: the core must build without them.
- The performance/expressiveness/maintainability/validation contracts (PLAN.md §1.1–§1.4:
  PC-1..10, EC-1..5, MC-1..5, VC-1..5) bind every change; overturning a row needs a §4 entry
  Volkan approves.
- Run the FULL ctest suite, not per-exe checks. Current baseline: 58/58 in Debug, Release, and
  CUDA lanes (PLAN.md §8, M0.10 row).

## 4. Environment facts (recorded, so you do not rediscover them)

- Windows 11 host; native builds MSVC/Clang; sanitizer + coverage lanes run in WSL Clang 18
  (ASan+UBSan full-suite, TSan on the concurrency binaries; llvm-cov exact Status coverage).
- The CUDA lane needs MSVC `cl.exe` on PATH (VS 18, `VC/Tools/MSVC/14.50.35717`) — an environment
  fact, not a code defect.
- Standard flow: `cmake -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build &&
  ctest --test-dir build`; the fast-math Release tree is `build-release`. Many stale `build-*`
  trees exist on disk; artifacts inside them from earlier dates are NOT evidence about current code.
- Hosted CI workflows are committed but cannot be claimed as run unless you saw the run.

## 5. Where things live

| What | Where |
|------|-------|
| The goal, contracts, ladder, decisions, state | `PLAN.md` (§1, §1.1–1.4, §6, §4, §8) |
| Bug ledger (fixes AND refutations) | `.claude/reports/p9b-bug-ledger-2026-07-10.md` |
| Quality-pass backlog (MQ.2's input) | `.claude/reports/code-quality-pass-2026-07-21{.md,-survivors.json}` |
| Gate/validation reports | `.claude/reports/` |
| Design notes | `.claude/designs/` |
| Session handoffs, plan archives | `.claude/summaries/` |
| Reference-repo digests (pouch-cell-spectral, unibatt) | `.claude/references/` |
| Derivations (grows during MQ.4/MQ.5) | `docs/derivations/` |
| Working TODO mirror (PLAN.md §8 is authoritative) | `develop/TODO.md` |

## 6. The current campaign (2026-07-23)

M0 is complete. The ladder stands at **MQ** — the quality wave: restore the Debug-AND-Release
baseline (MQ.1), disposition all 68 verified quality-pass findings (MQ.2), repo hygiene (MQ.3),
derivation inventory + independent re-derivation of the mathematics (MQ.4/MQ.5 — derive first,
then diff against the code), a fresh logic hunt (MQ.6), a structural performance hunt (MQ.7),
closeout (MQ.8). Then M1.0 onward: PyBOP suite, v4.0.0, expressiveness, parameters + citations,
FVM, SPMe, DFN, 2+1D/3D thermal, regime advisor, WASM, SLIDE Studio (the GUI), v5.0.0, the breadth
wave, Newman instrumentation, v6.0.0. One goal, many boxes, no stopping.
