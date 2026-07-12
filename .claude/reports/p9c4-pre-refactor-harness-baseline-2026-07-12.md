# M0.8 / 9C-4 pre-refactor shared-harness baseline (2026-07-12)

Baseline source: `93e75d1`, after M0.7 production and coverage closure code.
This inventory was completed before introducing a shared test helper or moving
any test operation behind it.

## Scope and invariant

PLAN M0.8 factors repeated successful build-batch/run/compare mechanics in
`tests/unit/core_*` into one header-only test harness. It does not move
reference values, physics oracles, tolerances, file I/O, allocation-fault
windows, backend-specific setup, or the operation under test. Every existing
binary must retain at least its baseline assertion and test-case counts;
recorded numeric hashes remain immutable.

The current tree contains **83 physical `buildSpmBatch` calls across 21 of 43
`core_*` test files**. The count is one above the earlier static inventory
because M0.7 added one valid NCH8 construction for the stale-Cycler-scratch
coverage witness. CTest remains 56 tests before the harness work.

## Exact duplicate and semantic groups

- `core_Recorder_test.cpp` and `core_AsyncRecorder_test.cpp` have identical
  whitespace-normalized 279-character `makeBatch(double electrode_area=-1)`
  bodies, serving 12 and 17 static call sites.
- The Recorder/AsyncRecorder allocation binaries repeat six ordinary two-lane
  Kokam construction blocks outside their fault-injection windows.
- Experiment's one-lane Kokam factory tuple is repeated in ParserAllocation,
  but ParserAllocation must remain direct because warm-up and allocation
  ordinals are test inputs.
- P1G4 restart repeats one construction three times; ExponentialModal repeats
  successful Euler/exponential and coarse/fine/adaptive construction tuples.
- ForwardSensitivity `productionTrace` and P7G3 `compareTrace` share one
  semantic pipeline: one-lane build, explicit current-density conversion,
  initial observation, ExponentialModal configuration, positive-grid stepping,
  and sample-major voltage collection.
- Experiment and Recorder repeat explicit current in amperes -> current density
  -> terminal-voltage observation.
- Recorder and AsyncRecorder enqueue/record loops deliberately stay separate:
  cadence, backpressure, snapshot timing, and I/O order are their oracles.

## Baseline assertion/case counts

Debug and ordinary Release match for every M0.8-scoped CPU binary. CUDA is the
validated CUDA Release artifact. Counts include the final M0.7 JSON and
stale-scratch regressions.

```text
AgeingKernel 384/6              AsyncRecorder 463/10
AsyncRecorderAllocation 60/4   CellDesign 10/2
ChebyshevEigenvalues 109/4      ChebyshevTransient 260/3
CompiledCurve 172/7             CudaDisabled 15/1
CudaSpmBatch 433671/4           Experiment 370/14
ExponentialModal 845/6          ForwardSensitivity 3919/4
Lam 57/3                        LithiumPlating 36/4
ModeC 52/6                      NetlistCsv 800/7
P1G0_pilot 4/1                  P1G2_allocation 17/3
P1G4_restart 48/1               P2G1_allocation 14/2
P2G5_ModeB 978/2                P7G3_PyBaMM 6334/2
PackSolver 878/26               PackStepper 214/8
PackTopology 148/9              ParameterSet 2375/8
ParserAllocation 1143/13        Recorder 168/7
RecorderAllocation 33/3         Sei 57/3
Simulation 175/6                SpectralDiffusion 10/2
SpectralModel 423/6             SpmElectrical 39/2
SpmFactory 58/7                 SpmObservables 21/2
SpmPipeline 18/1                SpmScalarKernels 213/4
SpmStress 19/3                  StateArena 626/11
SurfaceCrack 82/3               ThermalLumped 27/2
ThreadPool 605/9
```

CPU Debug totals 22,279 assertions in 227 cases across 42 binaries. Including
the CUDA artifact gives 455,950 assertions in 231 cases across all 43 source
tests. These are per-binary floors, not an aggregate budget: one binary cannot
compensate for a dropped assertion in another.

## Direct-call allowlist

Direct factory calls remain required where construction or special ownership is
the subject:

| Family | Reason to remain direct |
|---|---|
| AgeingKernel | heterogeneous/isolated lanes, recorded bits, failure atomicity |
| CudaSpmBatch | device/reference/backend setup is the oracle |
| ModeC, P2G5 | scale, allocation, and solver-mode construction |
| P1G2, P2G1, ParserAllocation | allocation measurement/fault ordinals |
| PackSolver, PackStepper | heterogeneous topology, aliasing, rollback, workers |
| ParameterSet, SpmFactory | absorption/factory construction is the subject |
| Simulation | build/step/solve, invalid input, and rollback boundaries |

Normal setup outside fault windows may use the harness in RecorderAllocation and
AsyncRecorderAllocation. Experiment may use the helper for ordinary one-lane
setup/voltage observation, but singular controls, callbacks, event roots,
Cycler configuration, and the NCH replacement witness remain direct.

## Preregistered harness contract

One `tests/support/CoreSpmTestHarness.hpp` will own only these mechanics:

1. successful construction of a caller-owned `SpmBatch` from explicit input,
   options, and lane count;
2. terminal-voltage observation from an explicitly tagged span of amperes or
   A/m2, with caller-owned density/voltage scratch;
3. constant-current ExponentialModal traces on a caller-supplied strictly
   increasing time grid, including the initial sample and every partial final
   interval, into caller-owned sample-major output;
4. transparent maximum-absolute and RMS voltage error over exact-size finite
   spans.

It will contain no default model options, chemistry data, C-rate conversion,
reference values, tolerances, file loading, stored batch/stepper/scratch state,
or mutable static object. `KokamSpmFixture.hpp` remains data-only. Expected
traces and pass/fail bands stay in their mirrored `core_*` tests.

## Adversarial preregistration

Before accepting migrations, a source-named harness test must turn red for each
of these controlled mutations:

- hard-code default NCH, lane count, or electrode area;
- flip the positive-discharge sign or omit A-to-A/m2 area conversion;
- drop initial/final samples or one interval, or step at the wrong time;
- skip a nonuniform partial interval;
- compare the actual trace with itself;
- ignore the final metric element (the unique maximum is placed there);
- divide RMS by `N-1` for literal errors `{1, 2, -1}`; expected RMS is
  `sqrt(2)` V and maximum absolute error is 2 V;
- map NaN/Inf to a passing zero;
- perturb one committed PyBaMM sample;
- remove one migrated `REQUIRE`, which must violate the per-binary count floor;
- move helper work inside an allocation-fault window, forbidden by the
  structural allowlist/review gate.

The first functional migration is the exact Recorder/AsyncRecorder batch
factory pair. Trace extraction follows only after the harness's independent
literal forwarding/framing/metric tests are mutation-red.
