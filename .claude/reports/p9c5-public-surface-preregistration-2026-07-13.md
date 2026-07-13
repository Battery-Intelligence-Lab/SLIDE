# M0.9 / 9C-5 public-surface audit — PREREGISTRATION (2026-07-13)

Written **before** any edit. Baseline commit `df2da52` (tree clean). Every band below is
registered here and judged verbatim in the validation report; post-hoc bands are not allowed
(CLAUDE.md §3).

## 1. Measured baseline (facts, not targets)

A public-consumer translation unit — the headers the Python bindings, the MATLAB MEX, the
benchmarks, and `docs/v4/quickstart-cpp.md` actually include —

```cpp
#include <core/Experiment.hpp>
#include <core/ParameterSet.hpp>
#include <core/ExponentialModal.hpp>
#include <core/ForwardSensitivity.hpp>
#include <core/CudaSpmBatch.hpp>
```

transitively pulls **26 `src/core` headers** (`clang++ -std=c++20 -Isrc -MM`, run 2026-07-13):

```
Experiment EulerLegacy SpmFactory SpmPipeline AgeingKernel Lam SpmStress CompiledCurve
Numeric StateArena SpmScalarKernels SpmObservables BatchView CellDesign SpmState
BatchBuilder LithiumPlating Sei SurfaceCrack ThermalLumped ExponentialModal ParameterSet
ForwardSensitivity CudaSpmBatch AsyncRecorder Recorder
```

Root cause [confirmed, `src/core/SpmFactory.hpp:8`]: the public factory header includes
`SpmPipeline.hpp`, which drags in the entire ageing/observable kernel stack, because
`SpmFactoryInput` stores the four ageing parameter structs by value and `SpmBatch` stores an
`SpmPipelineLayout` by value. The **types** are needed; the **kernels** are not.

Kernel headers presently reachable from the public surface (the eight to eliminate):
`SpmPipeline.hpp`, `SpmObservables.hpp`, `Sei.hpp`, `Lam.hpp`, `SurfaceCrack.hpp`,
`LithiumPlating.hpp`, `ThermalLumped.hpp`, `SpmStress.hpp`.

## 2. Design (chosen under §0.3(4); recorded, not asked)

Three header surfaces, one tag per header (`@surface` in the file's contract comment):

| Surface | Meaning | May be included by |
|---|---|---|
| `api` | supported user API | anything |
| `support` | value/vocabulary types the API's signatures are written in | `api`, `support`, `internal`, tests |
| `internal` | kernels and machinery; implementation detail (MC-5) | `internal` sources, tests, benchmarks — **never** an `api`/`support` header, never a binding |

Enabling split (behaviour-preserving): move the pure-POD types out of the kernel headers.
`SeiParams`, `SurfaceCrackParams`, `LamParams`, `LithiumPlatingParams` → four `*Params.hpp`
`support` headers (MC-2 glob cohesion preserved: `core/Sei*`); `SpmStateLayout`,
`ThermalLumpedLayout`, `StressHistoryLayout`, `SpmPipelineLayout` → one `SpmBatchLayout.hpp`
(`support`), one concept: which rows an SPM batch owns.

## 3. Registered bands — judged verbatim after the runs

**B1 — digit-identical behaviour.** Every pre-existing Catch2 binary reports assertion and
test-case counts **identical** (not ≥) to `baseline-debug.txt`, captured at `df2da52` before
any edit. CTest 57/57 Debug, 57/57 fast-math Release, 57/57 CUDA Release; `CudaSpmBatch`
stays at exactly 433671 assertions / 4 cases. Any new binary is additive and named.
*Failure = the refactor is not a refactor.*

**B2 — MC-5 property, compiler-verified.** After the split, the dependency list of the probe
TU above contains **zero** of the eight kernel headers (registered target: 8 → 0), and no
`api` or `support` header includes an `internal` header. Judged by `clang -MM`, not by eye.
Predicted core-header count of the probe TU after the split: 20 ± 3 (registered; the eight
kernels leave, four `*Params.hpp` + `SpmBatchLayout.hpp` arrive).

**B3 — the gate is real.** Every one of these mutations must turn a gate RED. Registered now:

| # | Mutation | Rule it must trip |
|---|---|---|
| 1 | delete a header's `@surface` tag | R1 (every core header classified) |
| 2 | set `@surface bogus` | R1 (valid value) |
| 3 | `#include "SpmPipeline.hpp"` in `Recorder.hpp` (an `api` header) | R2 (api ⊅ internal) |
| 4 | `#include "SpmObservables.hpp"` in an `api` header | R5 (probe TU dependency scan) |
| 5 | include an `internal` header from `python/bindings.cpp` | R6 (bindings see only `api`) |
| 6 | add a public entity to an `api` header without updating the manifest | R4 (surface manifest) |
| 7 | delete a public entity from an `api` header | R4 (manifest, other direction) |
| 8 | re-add `nLanes()` alongside `n_lanes()` | R7 (one name per concept) |
| 9 | declare a `validate…` function returning `bool` | R7 (validate → Status; valid/is → bool) |
| 10 | add a `getFoo()` accessor | R7 (no `get` prefix in core) |

**B4 — not claimed.** No performance or compile-time claim is made. Wall-clock on this machine
is unreliable (§8 standing note); the MC-5 win is stated as a header-dependency fact, not a
build-time number.

## 4. Naming lexicon to be enforced (one verb per concept)

| Verb | Concept | Current uses |
|---|---|---|
| `build…` | validated input → runnable object | `buildSpmBatch`, `BatchBuilder::build` |
| `compile…` | description → immutable compiled form | `compilePackDescription`, `compileSpectralModel`, `compileThermalLumped` |
| `declare…` | register state rows in the arena | `declareSpmState`, `declareStressHistory`, `declareThermalLumped` |
| `compute…` | pure kernel evaluation into caller storage | `computeSei`, `computeLam`, `computeSpmObservables`, … |
| `add…Rhs` | accumulate into the RHS | `addSeiRhs`, `addLamRhs`, `addSpmDiffusionRhs`, … |
| `validate…` | → `slide::Status` | `validateSeiParams`, `validateElectricalNetlist`, … |
| `valid…`/`is…` | → `bool` | `validSegment`, `is_finite`, `valid_ageing_model_mask` |
| `parse…` | in-memory text → object | `parseLiionpackNetlistCsv`, `parseStrictJson` |
| `load…` | file path → object | `loadLiionpackNetlistCsv` |
| `T::from…` | named constructor | `ParameterSet::fromBpxFile` |
| `step…` | one integration step (stepper API) | `ExponentialModal::step`, `PackStepper::step` |
| `advance…` | internal per-lane primitive | `advanceModal`, `advanceEuler` |

Duplicate-concept names to be killed [confirmed by grep at `df2da52`]:

- lane/row counts: `nLanes()`/`nRows()` (`Recorder.hpp:79,80,139,140`, `AsyncRecorder.hpp:252,253`,
  `CudaSpmBatch.hpp:57`) vs the majority spelling `n_lanes()`/`n_rows()` (`StateArena`, `BatchView`,
  `SpmBatch`, `AgeingKernel`). **Canonical: `n_lanes()`/`n_rows()`.**
- device counters: `cuda::deviceAllocations`/`deviceWideSynchronizations`/`deviceBytes`
  (`CudaSpmData.hpp`) vs `CudaSpmBatch::deviceAllocationCount`/`deviceWideSynchronizationCount`/
  `deviceArenaBytes`. **Canonical: the `…Count` / `deviceArenaBytes` spellings.**
- private lane-check helper, one concept and three spellings: `checked_lane_count`
  (`AgeingKernel.hpp:92`, `SpmObservables.hpp:410`), `checkedLaneCount` (`SpmPipeline.hpp:690`),
  `checkedShape`/`checked_shape` (`SpmObservables.hpp:104`, `StateArena.hpp:196`).
  **Canonical: `checked_lane_count` / `checked_shape`.**

Case convention is NOT swept globally: `snake_case` for vocabulary/predicate helpers and
`camelCase` for service-class methods is the existing split, and a global rename would be a
large diff with no stated mandate. What M0.9 enforces is *one name per concept*, which is what
the box asks for.
