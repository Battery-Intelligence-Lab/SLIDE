# M0.9 / 9C-5 — design note: what the public surface is, and how it stays that way

**Date** 2026-07-13. **Status** ASSUMED (§0.3(4): decided by the implementer, recorded here, not
awaiting Volkan). **Contract touched** MC-5; MC-1, MC-2, MC-3 preserved; PC-10 respected.

## Problem

MC-5 says "public surface minimal and intentional; implementation detail lives in `detail::` or
private headers", and M0.9 asks for API-vs-detail classified per header plus one verb per concept.
Two facts made this more than a labelling exercise [both confirmed at `df2da52`]:

1. `src/core/SpmFactory.hpp:8` included `SpmPipeline.hpp`. Because `SpmFactoryInput` stores the four
   ageing parameter structs **by value** and `SpmBatch` stores the row layout by value, the *types*
   were needed — but the include dragged the whole kernel stack with them. A user translation unit
   including `<core/Experiment.hpp>` compiled `SpmPipeline`, `SpmObservables`, `Sei`, `Lam`,
   `SurfaceCrack`, `LithiumPlating`, `ThermalLumped`, and `SpmStress`: 26 core headers in total.
2. The same concept had several names on the public surface: `n_lanes()`/`nLanes()`,
   `deviceAllocationCount`/`deviceAllocations`, and a private lane check spelled three ways.

## Options considered

**A. Tag headers, change nothing else.** Cheap, and honest about the status quo — but it would
classify a header as `api` while that header still compiles a kernel into every user TU. It labels
the violation instead of removing it. Rejected.

**B. Pimpl the factory input.** Hide `SpmFactoryInput`'s parameter blocks behind a pointer so the
public header needs no parameter types. Removes the kernels *and* the params from the API, but costs
an allocation and an indirection on a cold path, and makes `SpmFactoryInput` no longer an aggregate a
caller can brace-initialise — a real expressiveness loss (EC-2/EC-3 in spirit). Rejected.

**C (chosen). Split the POD types out of the kernel headers.** The parameter structs and the row
layouts are value types; the kernels that consume them are not. Move
`SeiParams`/`SurfaceCrackParams`/`LamParams`/`LithiumPlatingParams` into four `*Params.hpp` headers,
the four layout structs into one `SpmBatchLayout.hpp`, and the ageing model-bit vocabulary into
`AgeingModelMask.hpp`. `SpmFactory.hpp` then includes types only; `SpmFactory.cpp` includes the
pipeline. No allocation, no indirection, aggregates stay aggregates, and MC-2's glob cohesion holds
(`core/Sei*` still finds SEI's params, kernel, and test).

## The classification

Three tiers, one `@surface` tag per header, enforced by `tests/structural/p9c5_public_surface.cmake`:

- **api** — the supported user surface; the only headers a binding or a doc snippet may include.
  16 headers: AsyncRecorder, CellDesign, CudaSpmBatch, EulerLegacy, Experiment, ExponentialModal,
  ForwardSensitivity, NetlistCsv, PackSolver, PackStepper, PackTopology, ParameterSet, Recorder,
  Simulation, SpmFactory, ThreadPool.
- **support** — the vocabulary the API's signatures are written in. Reachable from an api header by
  construction, not intended for direct user inclusion. 12 headers: AgeingModelMask, BatchBuilder,
  BatchView, CompiledCurve, LamParams, LithiumPlatingParams, Numeric, SeiParams, SpmBatchLayout,
  SpmScalarKernels, StateArena, SurfaceCrackParams.
- **internal** — kernels and machinery. Everything else, including all of `src/core/detail/`.

**The load-bearing rule is R3: no api or support header may include an internal header.** Because
every core header is classified (R1) and includes are the only reachability channel, R3 makes the
MC-5 property inductive: no internal header is reachable from the api set. `clang -MM` confirms it
empirically — 21 core headers in a five-api-header TU, zero of them kernels.

## Accepted residuals, recorded not hidden

- `SpmScalarKernels.hpp` is classified **support**, not internal, and therefore does reach a user TU
  (via `CompiledCurve`'s inline evaluator, which calls `linearInterpolate`). Re-implementing
  interpolation inside `CompiledCurve` would duplicate a physics kernel and violate PC-10, which is
  the stronger contract. PC-10 wins; the fact is stated rather than buried.
- `EulerLegacy.hpp` is **api**, not internal, because `CyclerIntegrator::euler_legacy` is a
  user-selectable option and `CyclerV2` holds the stepper by value.
- R4 pins each api header's namespace-scope declaration count. That is a tripwire, not a parser: it
  sees additions and removals, but two simultaneous offsetting edits would net out.

## The naming lexicon (one verb per concept)

`build…` validated input → runnable object · `compile…` description → immutable compiled form ·
`declare…` register arena rows · `compute…` pure kernel evaluation into caller storage ·
`add…Rhs` accumulate into the RHS · `validate…` → `slide::Status` · `valid…`/`is…` → `bool` ·
`parse…` in-memory text → object · `load…` file path → object · `T::from…` named constructor ·
`step…` a stepper's public step · `advance…` internal per-lane primitive.

Killed duplicate spellings: `nLanes()`/`nRows()` → `n_lanes()`/`n_rows()`; `deviceAllocations`,
`deviceWideSynchronizations`, `deviceBytes` → `deviceAllocationCount`,
`deviceWideSynchronizationCount`, `deviceArenaBytes`; `checkedLaneCount`/`checkedShape` →
`checked_lane_count`/`checked_shape`. The gate forbids the dead spellings and requires the canonical
ones, so the decision cannot rot.

**Case style is deliberately NOT swept.** `snake_case` for vocabulary/predicate helpers and
`camelCase` for service-class methods is the existing split; a global rename would be a large,
risky diff with no mandate in PLAN §3.10. M0.9 enforces *one name per concept*, which is what the
box asks for.
