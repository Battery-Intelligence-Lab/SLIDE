# M0.10 / 9C-6 dead-code + line-debt sweep — PREREGISTRATION (2026-07-13)

Written **before** any edit. Baseline commit `d79d194`. Gate **P9-G5**: every simplification
digit-identical; core line count recorded before/after in PLAN §8 (net reduction expected; growth
requires written justification).

## 1. Measured baseline (facts)

`src/core` totals **17,363 lines** across 51 files. Files over the MC-1 ~700-line smell:

| File | lines | status |
|---|---|---|
| `src/core/PackSolver.cpp` | 910 | split candidate |
| `src/core/AsyncRecorder.cpp` | 890 | split candidate |
| `src/core/SpmPipeline.hpp` | 783 | split candidate (partial) |
| `src/core/CyclerV2.cpp` | 771 | **already a recorded M0.7 cohesion exception** |
| `src/core/Recorder.cpp` | 721 | split candidate |
| `src/core/SpmFactory.cpp` | 713 | split candidate |
| `scripts/status_failure_coverage.py` | 1,281 | split candidate (PLAN names it) |

Test files over 1,000 lines (M0.8 passed its line debt here): `core_ParserAllocation_test.cpp`
1,229; `core_PackSolver_test.cpp` 1,198; `core_Experiment_test.cpp` 1,168;
`core_AsyncRecorder_test.cpp` 1,118; `core_ParameterSet_test.cpp` 1,007. (`Module_s_test.cpp`
1,138 is legacy v3, out of scope.)

**Dead code: there is none to remove.** A scan of every namespace-scope entity in `src/core`
found no entity that is unreferenced: the 17 flagged by a naive "not named outside its own
header" pass are all aggregate members or same-file callees (e.g. `ElectrodeData` is a member of
`KernelParams`; `computeSpmSurfaceConcentrations` is called at `SpmObservables.hpp:445`).
`SpectralDiffusionLegacy.hpp` is the deliberate parity oracle — PLAN says KEEP. So 9C-6 is a
**line-debt** sweep, not a dead-code sweep, and I register that now rather than manufacturing
deletions to look productive.

## 2. Splits to attempt, and the concept boundary each claims

| File | Split | The boundary |
|---|---|---|
| `SpmFactory.cpp` 713 | → `SpmBatch.cpp` + `SpmFactory.cpp` | the **runtime batch object** (move semantics, `rhs`, `evaluate`, steps, observations) vs the **cold factory/registry** that compiles parameters and builds it |
| `Recorder.cpp` 721 | → `Recorder.cpp` + `BinaryRecording.cpp` | **writing** a recording (recorder + CSV/binary/parquet sinks) vs **reading** one back (the mmap `BinaryRecording` + its layout) |
| `AsyncRecorder.cpp` 890 | → `AsyncRecorder.cpp` + `AsyncRecordingCodec.cpp` | **transport** (slot ring, worker, drain, lifecycle) vs **codec** (byte shuffle/unshuffle, compression, buffer layouts) |
| `PackSolver.cpp` 910 | → `PackSolver.cpp` + `PackSolverStrategies.cpp` | **assembly** (Thevenin linearisation, workspace, configure, dispatch) vs the three **solve strategies** (sparse Newton, ladder, relaxation) |
| `SpmPipeline.hpp` 783 | → extract `SpmDiffusionRhs.hpp` | every other mechanism's RHS mapping already has its own header (`Sei.hpp`, `Lam.hpp`, …); diffusion's is the one that does not. After the extraction the pipeline is still ~720 lines and gets a **cohesion justification**, not a forced second cut. |
| `scripts/status_failure_coverage.py` 1,281 | → package modules | tool discovery/validation vs profile export vs arm mapping vs report emit |

Where a boundary is not real, the file gets a written per-file cohesion justification in the
validation report and PLAN §8 — that is what the box permits, and it is more honest than cutting
a cohesive file in half to satisfy a number.

## 3. Registered bands — judged verbatim after the runs

**B1 — digit-identical. ** Every one of the 53 Catch2 binaries reports assertion and test-case
counts **identical** to the `d79d194` baseline (`baseline-debug.txt` / `baseline-release.txt`, the
same files used for M0.9). CTest 57/57 Debug, 57/57 fast-math Release, 57/57 CUDA.
`unit_test_core_CudaSpmBatch` stays at exactly 433,671 assertions / 4 cases. The coverage-script
split additionally must leave `scripts/status_failure_coverage.py`'s output schema unchanged
(same JSON keys); it is not re-run as a coverage lane here and no coverage claim is made.

**B2 — line count.** Registered **prediction, and it contradicts PLAN's default expectation**:
splitting a translation unit *adds* lines (a new header guard, an include block, a re-declared
anonymous namespace). I predict `src/core` **grows by 1–4%** (17,363 → 17,500–18,050) and I am
registering that now rather than discovering it afterwards. The MC-1 goal is *files under ~700
lines*, not a smaller total; a smaller total would only come from deleting code, and there is no
dead code to delete (§1). Registered target: **no `src/core` file over 700 lines except those
with a written cohesion justification**, and every growth line accounted for.

**B3 — the gates stay real.** The existing structural gates (PC-10, 9C-2, 9C-3, 9C-4, 9C-5) must
all stay green, and the 9C-3 gate's per-file line ceilings must be updated deliberately, not
loosened to accommodate a bad split. At least one mutation per new translation unit must turn a
gate red (a moved function must still be the one that runs).

**B4 — not claimed.** No performance claim. No new coverage claim. TSan not rerun.
