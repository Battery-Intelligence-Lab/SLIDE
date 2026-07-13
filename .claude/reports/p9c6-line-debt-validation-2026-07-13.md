# M0.10 / 9C-6 dead-code + line-debt sweep — validation (2026-07-13)

Baseline commit `4259caa`. Preregistration (bands, splits, and the prediction that the line count
would *grow*) written before any edit:
`.claude/reports/p9c6-line-debt-preregistration-2026-07-13.md`. Gate **P9-G5**.

## 1. Dead code: there was none. That is the finding.

Every namespace-scope entity in `src/core` is referenced. The 17 that a naive "not named outside
its own header" scan flagged are aggregate members or same-file callees (`ElectrodeData` is a
member of `KernelParams`; `computeSpmSurfaceConcentrations` is called at `SpmObservables.hpp:445`).
`SpectralDiffusionLegacy.hpp` is the deliberate parity oracle — PLAN says KEEP. So no deletions
were manufactured to look productive, and 9C-6 is a line-debt sweep.

## 2. The splits, and the boundary each turned out to have

| Before | After | The boundary that made it real |
|---|---|---|
| `SpmFactory.cpp` 713 | `SpmFactory.cpp` 494 + **`SpmBatch.cpp` 237** | the runtime batch reaches its implementation *only* through type-erased pointers — so `SpmBatch.cpp` does not include `SpmPipeline.hpp` at all, which is what makes the erasure real rather than decorative |
| `Recorder.cpp` 721 | `Recorder.cpp` 293 + **`RecordingFormat.cpp` 451** | recorder + text/columnar sinks vs the CRC-hardened binary format. Writer and reader stayed together: they must agree byte for byte, and separating them would put the header struct and its CRC behind a seam that exists only to satisfy a line count |
| `AsyncRecorder.cpp` 890 | `AsyncRecorder.cpp` 396 + **`AsyncRecordingCodec.cpp` 452** + **`detail/AsyncRecordingFormat.hpp` 116** | transport (slot ring, drain worker) vs codec (shuffle, block layouts, reader); the format they share has one definition |
| `PackSolver.cpp` 910 | `PackSolver.cpp` 610 + **`PackSolverIterative.cpp` 266** | `solveLadder` and `solveRelaxation` never touch `SolverWorkspace::Impl` — they assemble and factorise nothing. That is not a line-count cut; it is the actual difference between the matrix-free strategies and the sparse Newton solve, and it is why they *can* move |
| `SpmPipeline.hpp` 783 | `SpmPipeline.hpp` 745 + **`SpmDiffusionRhs.hpp` 64** | every other mechanism's RHS mapping already had its own header; diffusion's was the one still inside the pipeline that composes them |
| `scripts/status_failure_coverage.py` 1,281 | `status_failure_coverage.py` 135 + **`slide_coverage/{scanner 351, session 278, classify 251, report 347, constants 11}`** | scan → session/profile → classify → report, with the CLI left as the entry point |

**One real duplicate removed:** `checkedAdd`/`checkedMultiply` existed twice, in `Recorder.cpp` and
`AsyncRecorder.cpp`, with identical semantics and different parameter names. They now have one
definition in `detail/CheckedArithmetic.hpp` (49 lines).

## 3. Registered bands — verdicts

**B1 — digit-identical. PASS.** All 53 Catch2 binaries report assertion and test-case counts
identical to the pre-M0.9 baseline, in Debug **and** in a like-for-like fast-math Release build.
CTest **58/58 Debug, 58/58 Release, 58/58 CUDA**; `unit_test_core_CudaSpmBatch` unchanged at
433,671 assertions / 4 cases. The coverage reporter's own `self-test` passes and its `census`
output is byte-identical (`{"conditional": 16, "direct": 354, "total": 370, …}`), and `ruff`
reports no undefined name across the new package — the module split cannot have silently dropped a
reference into the llvm-cov path I cannot run on this machine.

**B2 — line count. PASS against the registered prediction, which contradicted PLAN's default.**
`src/core` went **17,363 → 17,595 (+232, +1.3%)**, inside the registered 1–4% growth band. PLAN
expects a net reduction; splitting a translation unit *adds* lines (a guard, an include block, a
re-declared namespace), and there was no dead code to delete. The MC-1 goal — files under ~700
lines — is met: **only two `src/core` files exceed 700, and both are recorded exceptions** (below).
The −219/−494/−428/−300 reductions in the four oversized TUs are real; the +232 total is the cost
of six new files, and it is stated rather than hidden.

**B3 — the moved code is the code that runs. PASS: 5/5 mutations red.** A mutation in each new
translation unit turns a gate red: halving `dt` in `SpmBatch.cpp` (ExponentialModal RED), dropping
a snapshot in `RecordingFormat.cpp` (Recorder RED), removing `byteShuffle` from
`AsyncRecordingCodec.cpp` (build RED), gutting `solveLadder` in `PackSolverIterative.cpp` (build
RED), and dropping the 64-byte alignment in `detail/CheckedArithmetic.hpp` (Recorder RED).

**First attempt discarded, not reported as evidence.** Three of my initial mutations came back
green — a CRC flipped in `RecordingFormat.cpp` (writer and reader use the same mutated function, so
it is self-consistent), an unexercised overflow in `checkedMultiply`, and a skipped `SpmBatch::rhs`
that `ExponentialModal` never calls. Those prove nothing about the split, so they were replaced
with mutations an existing oracle can actually see, rather than counted as passes.

**B4 — not claimed.** No performance claim. The Linux llvm-cov coverage lane was **not re-run**:
the reporter's structure changed, its behaviour is verified only by its self-test, its census, and
a static undefined-name check. The hosted workflow's path filter now includes
`scripts/slide_coverage/**`, so the lane will trigger on the package — but this session does not
claim it ran.

## 4. MC-1 exceptions, recorded

- **`src/core/SpmPipeline.hpp` (745).** One concept: the compile-time composed RHS pipeline. It is
  a template class whose size comes from composing five optional mechanisms with explicit
  construction proofs; it cannot be split across translation units, and splitting the *class* would
  destroy the zero-overhead composition (D-02). The one genuinely separable piece — the diffusion
  RHS mapping — has been extracted.
- **`src/core/CyclerV2.cpp` (771).** Already a recorded M0.7 cohesion exception (the frozen
  non-IPO runner transaction).

## 5. UNMET — the test-file line debt, for the second milestone running

M0.8 was assigned the oversized test fixtures and grew them instead. `develop/TODO.md` then passed
that debt to M0.10. **M0.10 has not cleared it either**, and I am recording that plainly rather
than letting it disappear:

| Test file | lines |
|---|---|
| `tests/unit/core_ParserAllocation_test.cpp` | 1,229 |
| `tests/unit/core_PackSolver_test.cpp` | 1,198 |
| `tests/unit/core_Experiment_test.cpp` | 1,168 |
| `tests/unit/core_AsyncRecorder_test.cpp` | 1,118 |
| `tests/unit/core_ParameterSet_test.cpp` | 1,007 |

These are untouched. The production-side sweep that P9-G5 gates is done; the test-side sweep is
not, and it now needs an explicit owner rather than a third inherited assignment. It is carried in
PLAN §8 and `develop/TODO.md` as open MC-1 debt.
