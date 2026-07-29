# MQ.2 quality-backlog disposition — preregistration (2026-07-23)

## Scope and source identity

This registration precedes every MQ.2 source/test edit and every decisive MQ.2 run.
The source baseline is clean commit `940716b` (`Close MQ.1 with the restored
three-lane baseline`). MQ.1 established 58/58 in fresh Debug, fast-math Release,
and CUDA lanes at source commit `a49224e`; `940716b` changes evidence and standing
documentation only.

MQ.2 has **71 decisions**:

- the 68 unique IDs in
  `.claude/reports/code-quality-pass-2026-07-21-survivors.json`;
- `surfacecrack-arrhenius-association`, the separately recorded `SurfaceCrack`
  Arrhenius association;
- `clang-ofast-deprecated`, MQ.1's deprecated redundant `-Ofast` observation;
  and
- `benchmark-lp-size-narrowing`, MQ.1's legacy `benchmark_LP_cases`
  `size_t`-to-`int` narrowing diagnostic.

The final disposition table must contain all 68 JSON IDs exactly once and all
three supplemental IDs exactly once. Each row has exactly one state:
`APPLIED`, `REFUTED`, or `DEFERRED — <named ladder owner>`. A row is not
`APPLIED` merely because prose was changed: its registered runtime, structural,
or documentation gate must pass. A deferral names the existing PLAN box that
owns both the work and its gate.

Preliminary triage predicts 60 `APPLIED`, one `REFUTED`, and seven named
deferrals among the 68 JSON rows. The three supplemental rows are intended
`APPLIED`. These are hypotheses, not results; the final counts are derived from
the completed table.

## Anti-balloon decisions registered before implementation

Seven findings are real but belong to an already-defined deeper box:

| Finding | Registered owner | Why it is not a safe MQ.2 cleanup |
|---|---|---|
| `observables-struct-passed-by-value` | MQ.7 structural performance hunt | A 368-byte source type does not prove a copy survives inlining; inspect IR/counters before changing aliasing. |
| `stress-mirror-index-algebra` | MQ.5 stress re-derivation | The legacy node association must be derived before changing index algebra. |
| `mc3-contract-comments-missing` | MQ.3 repo-wide contract-comment hygiene | PLAN MQ.3 expressly owns the full census/content gate; the finding undercounts the debt and its proposed presence-only gate can pass content-free prose. |
| `crack-centre-node-index-arithmetic` | MQ.5 stress/ageing re-derivation | The centre-node association and const-preserving accessor shape must be derived first. |
| `forward-sensitivity-double-observe` | MQ.7 structural performance hunt | Moving the surviving call across an inline boundary can change fast-math hashes; prove the retained evaluation and counters first. |
| `forward-sensitivity-cold-rebuild-per-parameter` | MQ.7 structural performance hunt | The cold compile is seed-independent, but optimizer-context digit identity and factorisation/allocation counters are MQ.7 evidence. |
| `mc3-contract-comments-below-bar` | MQ.3 repo-wide contract-comment hygiene | PLAN MQ.3 expressly owns it as the same repo-wide MC-3 census, not four isolated headers. |

`enqueuesnapshot-success-untested` is registered for attempted refutation:
the CUDA P8-G3 path appears to execute `AsyncRecorder::enqueueSnapshot` through
`CudaSpmBatch.cpp` and bitwise-check the recorded densities/state. It becomes
`REFUTED` only if the current call graph, assertions, and MQ.1 CUDA execution
jointly establish that semantic success coverage. Otherwise it returns to
`APPLIED` with a direct test.

## Explicit decisions for the six named landmines

- `branch-graph-derived-twice`: extract one validated graph derivation now; a
  structural count and topology fixtures must reject either duplicate's return.
- `relaxation-target-triple-role`: give the three quantities independent
  storage now; the registered Mode-C trace must remain bit-identical.
- `crc32-twice`: make the recording subsystem use one owner now; unrelated
  gate-pinned CRC implementations outside recording are not swept into it.
- `fastmath-unsafe-stoichiometry-gate`: add the volatile-bit finite primitive
  before the range comparison now; Release NaN/Inf rejection is the decisive
  behavior gate.
- `substeps-name-contradicts-code`: the intended semantics are **retained** as
  `substeps * dt` with one frozen electrical solve, matching PLAN §3.5. The
  unusually explicit `PackStepper.hpp` contract landed in the 2026-07-21 pass
  and is the API disposition; an executable N-step oracle must pin it. No
  subdivision or rename is silently introduced.
- supplemental `SurfaceCrack` Arrhenius debt: route model 5 through the shared
  scalar helper now. The different association is allowed to change recorded
  low bits, but the legacy comparison must stay within `1e-12` relative and the
  new bits are frozen only after that independent oracle passes.

## Implementation batches

Each batch is a separate reviewable commit after its focused gates are green.
Dependencies run left-to-right within a subsystem; independent subsystem
batches may be prepared in parallel, but builds/tests are serialized.

| Batch | Findings / scope | Registered character |
|---|---|---|
| O1 observable storage/access | `raw-pointer-reimplements-batchview-at`, `scratch-per-lane-magic-21`, `transport-cache-encapsulation-and-lexicon`, `scratch-storage-double-validation`, `lane-index-cast-noise`, `duplicate-status-include` | Valid-input no-op; added bounds/shape assertions and tests. |
| O2 fast-math validity | `fastmath-unsafe-stoichiometry-gate` | Invalid-state behavior fix; valid recorded traces unchanged. |
| O3 observable/stress/thermal tests | `stress-uniform-profile-zero-oracle`, `thermal-untested-validation-branches` | Test-only. |
| A1 ageing scalar ownership | `sei-kinetic-current-duplicated-in-one-function`, supplemental SurfaceCrack Arrhenius | Keep SEI branch/model policy local while moving its repeated scalar expression to the existing scalar-kernel owner; SEI no-op. SurfaceCrack is a separate registered low-bit behavior change landed last. |
| F1 batch factory | `homogeneous-init-lane-loop`, `shared-constants-fanout`, `advance-euler-asserts-only` | Valid-input no-op; direct invalid Euler requests become `Invalid_parameters`. |
| P1 pack-solver algebra | `cell-current-reconstruction-x4`, `branch-drop-and-kcl-current`, `dead-usings-and-misplaced-comment` | Valid-input no-op under exact pack fixtures. |
| P2 Mode-C storage | `relaxation-target-triple-role` | No-op; three cold vectors replace one triple-use buffer without hot allocation. |
| S1 pack stepper cleanup | `packstepper-gather-scatter`, `unique-in-prefix-three-spellings`, `substeps-name-contradicts-code` | No numerical change; `substeps * dt` retained and executable. |
| S2 source-step atomicity | `diagnostics-not-rolled-back`, `untested-rollback-and-advance` | Deliberate failed-attempt diagnostics behavior plus rollback tests. |
| T1 topology derivation | `branch-graph-derived-twice`, `three-parallel-archetype-maps`, `test-gap-ladder-rollback` | Valid-input no-op plus deterministic rollback gate. |
| T2 thermal topology | `assemble-noexcept-and-stale-cold-brief`, `test-gap-isothermal-zero` | Signature/docs plus test-only exact-zero oracle. |
| C1 netlist parsing | `csv-limits-restated-in-prose`, `ascii-digit-scan-helper`, `failsemantic-returns-unused-bool` | Parser no-op on the recorded grammar; rejection offsets/status remain exact. |
| R1 recording common/core | `crc32-twice`, `header-crc-helper`, `csv-parquet-schema-twice`, `state-index-int-arith`, `snapshotview-twice`, `dead-usings-copy-pasted` | Byte-identical encoded files/CSV for recorded fixtures. |
| R2 recording semantics/tests | `crc-zero-rejected`, `shuffle-oracle-gap`, `csv-values-untested`, and disposition of `enqueuesnapshot-success-untested` | Zero-CRC behavior fix; test-only additions; prior shuffle gate retained. |
| Q1 parameter lexicons/helpers | `dup-electrode-name-lexicon`, `dup-segment-locator`, `dup-text-scanner`, `dup-electrode-status-ladder`, `dup-bpx-provenance-string`, `redundant-required-lookups` | Valid-input recorded ParameterSet fingerprints unchanged. |
| Q2 parameter bugs/API | `bug-nan-swallowed-by-sampler`, `bug-constant-probe-aliasing`, `clarity-handrolled-index-loops`, `api-nodiscard-noexcept`, `clarity-untraceable-constants` | Two registered invalid-expression behavior fixes; otherwise exact fingerprints/statuses. |
| E1 experiment parser/control constants | `dup-duration-unit-lexicon`, `clarity-unnamed-tolerances` | Parser/solver no-op; constants named with their units and scale dependence. |
| E2 drive-cycle alignment | `drive-cycle-report-lag` | Deliberate right-endpoint reporting/event behavior fix; applied interval current remains left-endpoint. |
| E3 simulation contract | `test-gap-solve-atomicity` | Test-only mid-run failure contract. |
| G1 CUDA cleanup/tests | `cuda-kernel-loop-invariant-and-open-coded-sign`, `cuda-kernel-unreachable-isfinite-and-duplicated-failure-arm`, `test-gap-untested-rejection-arms` | CUDA valid-output no-op plus test-only validation/cadence/inline-executor gates. |
| B1 build diagnostics | supplemental deprecated `-Ofast`; supplemental `benchmark_LP_cases` narrowing | Build-policy no-op: retain Release `-O3` + `-ffast-math`; make the constant extent type match the legacy `int` API. |

The already-landed rows (`pc10-arrhenius-handrolled`,
`quadratic-ladder-branch-lookup`, `shuffle-oracle-gap`, and
`test-gap-exponent-notation`) are not reimplemented. Their original evidence is
linked in the final table and their current gates still ride the final suites.

## Registered behavior gates

These bands are fixed before the corresponding test or source change.

1. **Fast-math surface stoichiometry.** Start from an otherwise-valid state
   with positive diffusion, then inject opaque bit-built qNaN, +Inf, and -Inf
   into a modal `z` row. Each returns `Invalid_states` in Debug and fast-math
   Release; finite stoichiometry at/below 0 and at/above 1 is also rejected. A
   valid interior control retains the existing recorded observable/ageing bits.
   The old comparison need not fail on this compiler—the decisive gate is the
   new structural requirement for the exact compact predicate
   `if(!is_finite_primal(z_surface)||!(primal_value(z_surface)>0.0&&primal_value(z_surface)<1.0))`
   appearing once in the `z_surface` rejection slice, plus its removal mutation.
2. **SurfaceCrack Arrhenius.** Each single model-bit mask
   `{1u, 2u, 4u, 8u, 16u}` still satisfies the independent
   `core_SurfaceCrack` legacy comparison at relative error `<= 1e-12`; the
   model-5 old/new relative
   delta is `<= 1e-13`; models 1–4 remain bit-identical; the Dual tangent retains
   a new registered temperature-Dual gate: seed `T=310 K` with tangent 1, compare
   model-5 crack-surface-rate value/tangent against the double value and a
   centered `h=1e-4 K` finite difference, and require
   `abs(dual_tangent-fd) <= 1e-8*max(abs(fd),1e-30)`. Valid Status and assertion
   counts do not fall. Relative error uses the existing oracle denominator
   `max(abs(expected), 1e-30)`; if old and new are both exactly zero, their
   relative delta is defined as zero. Both double and Dual value/tangent paths
   are covered. Recorded scalar hashes may change and are accepted only after those
   independent bands pass. Restoring the private associated expression must
   make the PC-10 structural gate red.
3. **Direct internal Euler validation.** `SpmPipeline::advanceEuler` validates
   before constructing `lanePeriod` or mutating a `z` row: the state view has
   exactly the pipeline row count and lane count, `ctx.i_app.size()` and
   `terminal_voltage.size()` each equal `n_lanes`, every current and
   `ctx.time` and `ctx.dt` are finite, and the explicit `dt` is finite and
   strictly positive. Wrong
   state rows/lanes, empty/short/long current or output spans, NaN/Inf current,
   zero/negative/NaN/Inf explicit `dt`, NaN/Inf `ctx.dt`, or non-finite time returns
   `Invalid_parameters` with state and output byte-identical. Valid
   public-wrapper fixtures remain exact.
4. **Mode-C/source-step rollback.** The test forces the direct caller attempt to
   fail, at least one hidden source step to succeed, and a later hidden step to
   fail. Two contracts and mutations are separate: (a) deleting restoration of
   `solution_`/`has_solution_` changes the published cell-current, node-voltage,
   terminal-voltage, or warm-start result and must turn red; (b) deleting the
   new diagnostics restoration changes iterations, residual, constraint,
   numeric-factorisation count, or `source_steps` and must turn red.
   Restored diagnostics describe the caller-requested failed direct attempt,
   with `source_steps == 0`; they do not describe the prior accepted solve.
5. **CRC value zero.** A syntactically valid recording whose computed CRC is
   exactly zero is accepted; corruption is still rejected. The deterministic
   fixture is a 3,584-byte zero vector with the preregistered 64-byte header
   `534c4944455245430100000004030201400000000000000070000000c700000087030000000000004000000000000000000e000000000000000e000000000000`
   at offset 0 and little-endian `uint64_t{3584}` at offset 64. Its CRC
   self-check must equal exactly zero. No unbounded search is permitted.
6. **Parameter expression non-finites.** The quarter-point BPX expression
   `x + 0/(x-0.25)` fails evaluation at exactly `x=0.25`, reaches the sampler
   as qNaN, and is rejected as `Invalid_parameters` with diagnostic
   `BPX function is non-finite on stoichiometry [0,1]` and the caller's
   sentinel `ParameterSet` byte-for-byte unchanged, rather than being swallowed
   by min/max. Separately, a direct `sampleParameterCurve` fixture returns an
   opaque bit-built +Inf at exactly the existing `x=0.75` probe and must return
   an empty/invalid curve; this is the infinity coverage because the BPX
   evaluator converts its own non-finite failures to qNaN before sampling. A fast-math-safe
   boolean records the non-finite probe; no optimizer-visible infinity sentinel
   is introduced. Finite sampled expressions retain fingerprints.
7. **Constant-probe aliasing.** The concrete cubic
   `3.3e-14+x*(x-0.5)*(x-1)` aliases the old `{0,0.5,1}` probe set. The new
   `{0,0.25,0.5,0.75,1}` check compares **all five** samples, rejects the
   state-dependent alias as `Invalid_parameters` with diagnostic
   `state-dependent BPX diffusivity is not supported by the constant-D SPM composition`
   and atomic caller output. The comparison tolerance stays exactly
   `64*epsilon*max(min_normal,abs(samples[0]))`; a true constant retains its
   scalar representation and exact bits. Restoring the old three probes, or
   ignoring either new sample, must turn the fixture red.
8. **Drive-cycle alignment.** For the existing US06 fixture with
   `(t,I)={(0,1),(0.5,-1),(1,0)}` and `dt=0.25`, the first right-endpoint drive
   samples at output indices 9–12 are exactly `{0,-1,-0.5,0} A`. A custom
   termination `1 A - current` on a linear `0 -> 2 A` drive cycle is located at
   local time `0.5 s` within `1e-12`, and reports `1 A` within `1e-12`.
   The implementation must call `driveCurrent` in all three post-advance
   resolution ladders: full `dt`, bisection `middle`, and accepted `high`;
   omitting any arm must turn the output/event test red. A second identical
   batch advanced once for exactly `0.5 s` at constant `0 A` is the exact
   comparator for the event run and must finish with a byte-identical state:
   reporting changes, applied interval current does not.
9. **Simulation mid-run failure.** An extreme-current first step returns a
   non-Success status, publishes a one-sample truncated solution with
   `solution.termination == status`, preserves the time/voltage lane-shape
   invariant, and leaves elapsed time exactly zero because the failed Euler step
   rolls back.
10. **CUDA host rejection.** Wrong-size current, NaN current, zero/negative/NaN/
    infinite `dt`, and non-finite time each return `Invalid_parameters` and leave
    downloaded device state byte-identical.
11. **CUDA recorder cadence.** Configure cadence 3, three ring slots, and
    `AsyncBackpressurePolicy::block`. Step the GPU before each accepted enqueue;
    steps 1..6 write only 3 and 6. Repeating 6 is rejected before `finish`, step
    9 succeeds, then `finish` reports exactly three snapshots and zero
    backpressure thinning.
12. **Single-worker executor.** A callback throwing at index 1 visits all three
    indices exactly once and returns `Unknown_problem`; an
    `Invalid_parameters` at index 1 followed by a throw at index 2 preserves the
    lowest-index `Invalid_parameters`.

## Registered no-op and structural gates

- Existing embedded `RecordedBits` fixtures are the primary digit oracle:
  `core_AgeingKernel` (all-mask trace), `core_Experiment` (parser trace),
  `core_ParameterSet` (Chen/input/BPX traces), and `core_CudaSpmBatch`
  (CPU/device trace). Expected constants stay unchanged except the explicitly
  registered SurfaceCrack and drive-cycle behavior rows.
- Before the homogeneous factory rewrite, an NCH=12, thermal+all-ageing,
  nine-lane fixture records the SHA-256 of the entire padded state arena,
  derivative arena, and one valid post-step state in Debug and Release. All
  three digests must be identical afterward; padding remains exactly zero.
- Pack solver/topology/stepper refactors retain every exact/bitwise assertion in
  their existing unit binaries. Before the first edit, the touched binaries'
  assertion/test-case summaries are captured from the MQ.1 Debug and Release
  trees; after each batch, no pre-existing case or assertion may disappear.

### P0 pack-algebra oracle addendum

A read-only review fixed the following oracle inputs before the first P1/P2
source edit. Traces A–C were specified before their first test run; this
repository addendum makes that earlier registration durable. The initial D0
candidate was registered before its run after adversarial review found that
A–C use only dyadic resistances. Its mutation was subsequently falsified and
the revised Trace D below is registered before its first evidentiary run.

1. **Trace A — exact high-dynamic-range all-mode solve.** Compile
   `parallel(4, affine)` with `H = 0x1p53`,
   `E = {H, 1, H-1, -H}`, `R = {1,1,1,1}`, and applied current
   `H-8`. Fresh sparse, ladder, and relaxation solvers each use tolerance
   `1e-12`, maximum four iterations, and must finish in exactly two.
   Terminal voltage is bitwise `2`, cell currents are bitwise
   `{H-2,-1,H-3,-H-2}`, and node voltages are bitwise `{2,0}`.
   An independent branch walk must produce exact positive-zero KCL residuals
   at every node. Relaxation drift is positive zero and gain is exactly one.
2. **Trace B — resistor arm.** Compile
   `parallel(1, affine, link R=0.5)` with cell `E=2`, cell `R=0.5`,
   and applied current `1`. A cold sparse solve takes two iterations and
   publishes cell current `1`, terminal voltage `1`, and node voltages
   `{1,0,1.5}` bitwise. Reusing that solution as a warm relaxation start with
   gain one takes exactly one iteration, retains those bits, reports positive
   zero drift and constraint bound bitwise `1e-12`, and satisfies independent
   exact-zero KCL.
3. **Trace C — sparse damping cadence.** Compile `parallel(2, affine)` with
   `E={8,0}`, `R={1,1}`, and zero applied current. With tolerance `1e-12`
   and maximum three iterations, sparse solve must take exactly three:
   half damping, root, confirmation. It publishes voltage `4`, currents
   `{4,-4}`, nodes `{4,0}`, and independent exact-zero KCL.
4. **Trace D0 — first non-dyadic candidate, FALSIFIED.** A test-only affine
   adapter recorded every current span passed to `linearizeThevenin` before
   copying its constant `E/R` outputs. Fresh sparse, ladder, and relaxation
   solvers used `parallel(4, affine)`,
   `E={4,4.1,3.9,4.2}`, `R={0.1,0.2,0.15,0.3}`, applied current `2`,
   tolerance `1e-12`, and maximum four iterations; a damped sparse solver used
   `parallel(2, affine)`, `E={4,0}`, `R={0.3,0.7}`, zero applied current, and
   maximum three iterations. Baseline records were 19 doubles for each
   all-mode path and 15 for the damped path. All three configurations produced
   the same four FNV recurrences
   `{6ff139e00d3a135a,4cd7e2c4e437d0a6,8941772cd3334768,a943f04490d81eb6}`
   and mixed recurrences
   `{32a19f0565ed72b7,c7540d7f77909867,536ea92d297f75c9,a002fa5404b9b09b}`.
   However, changing all four reconstruction sites from `x/R` to
   `x*(1/R)` left the Debug gate green at 26/26: every realized
   numerator/resistance pair rounded identically. Production was restored.
   D0 is therefore exploration and falsification evidence, not an identity
   gate, and its hashes are not retained as acceptance criteria.
5. **Trace D — bit-separating callback and publication bits.** Before its
   first evidentiary build/run, the replacement is fixed as follows. Fresh
   sparse, ladder, and relaxation solvers use `parallel(2, affine)`, `E={0,0}`,
   `R={0.1,0.1}`, applied current `0.1`, tolerance `1e-12`, and maximum four
   iterations; each must take exactly two. The analytic terminal voltage is
   `-0.005`. Direct reconstruction `0.005/0.1` has bits
   `3fa9999999999999`, while the forbidden reciprocal multiplication has bits
   `3fa999999999999a`.

   The damped sparse path uses `parallel(2, affine)`, `E={0.5,0}`,
   `R={0.1,0.1}`, zero applied current, the same tolerance, and maximum three
   iterations. Its first undamped voltage is `0.25`; current magnitude `2.5`
   selects damping `0.8`, producing voltage `0.2`. The first reconstructed
   positive current is direct `(0.5-0.2)/0.1` with bits
   `4007ffffffffffff`, versus reciprocal-multiply bits
   `4008000000000000`. The solve must take exactly three iterations.

   For all four paths, two independent `RecordedBits` recurrences hash, in
   order, every callback-current frame, final cell currents, final node
   voltages, and one scalar frame containing terminal voltage, residual norm,
   constraint drift, constraint bound, and relaxation gain. Expected hashes
   are captured separately for Debug, fast-math Release/IPO-off, and
   Release/ThinLTO/CUDA-host builds through deliberate zero placeholders. The
   all-mode records contain 13 doubles each; the damped record contains 15.
   No hash may be reblessed during P1 or P2.

Revised Trace D is the decisive digit-identity gate for non-dyadic division at
all four current-reconstruction sites: sparse undamped, sparse damped, ladder,
and relaxation. A final solution-only hash is insufficient because a later
iteration could wash out a changed reconstructed current; callback frame two
is therefore part of the record.

### P1/P2 pack-algebra implementation addendum

The implementation shape and structural gates below are registered after the
read-only review and before the first P1/P2 production-source edit.

P1 gives the duplicated algebra five owners in `PackSolverInternal.hpp`.
Every scalar operand passed across a helper boundary is a `const real_t &`;
that spelling is load-bearing under finite-math because `Numeric.hpp`
documents that a by-value argument may acquire `nofpclass` and let Clang fold
the bitwise finiteness gate away.

- `branchDrop` owns the ordered positive-node minus negative-node
  subtraction.
- `BranchAffine` carries resistance and source. `branchAffine` must branch on
  kind before indexing the cell arrays; a resistor returns its own resistance
  and positive zero source.
- `branchCurrentNumerator` owns `drop-source`, and `branchCurrentOut` owns
  `numerator/resistance`. They deliberately do not validate: sparse retains
  its runtime numerator rejection while relaxation retains its documented
  assertion.
- `cellCurrentFromDrop` owns the checked discharge-positive
  `(ocv-drop)/resistance` reconstruction.

In compacted production sources the exact token census is:

| Token | owner | `PackSolver.cpp` | `PackSolverIterative.cpp` | total |
|---|---:|---:|---:|---:|
| `cellCurrentFromDrop(` | 1 | 2 | 2 | 5 |
| `branchDrop(` | 1 | 4 | 2 | 7 |
| `branchAffine(` | 1 | 2 | 2 | 5 |
| `branchCurrentNumerator(` | 1 | 1 | 1 | 3 |
| `branchCurrentOut(` | 1 | 1 | 1 | 3 |

The gate also forbids every old raw cell reconstruction, the six raw branch
drops, and the two cell/resistor affine ternaries. It slices sparse undamped
and damped bodies and requires one reconstruction-helper call in each. The
dead `<cstring>` include and three dead namespace `using` declarations are
forbidden. The Kahan/fast-math rationale moves immediately above its pragma
block; the false claim that the compensated helper serves the sparse solve is
deleted.

P2 replaces four multiply-used vectors with one private
`RelaxationScratch` containing exactly seven cold-sized vectors:
`diagonal`, `diagonal_compensation`, `rhs`, `rhs_compensation`, `target`,
`residual`, and `residual_compensation`. This is an increase of three vectors,
all allocated during transactional `configure`; the candidate aggregate is
published by one no-throw move only after all fallible work. Hot traffic stays
at exactly six fills: diagonal and its compensation plus RHS and its
compensation at entry, then residual and its compensation immediately before
KCL. `target` is not filled because every node is unconditionally assigned.
No vector construction is permitted inside `solveRelaxation`.

Within compacted `solveRelaxation`, indexed-use counts are exactly
`diagonal[`: 4, `diagonal_compensation[`: 2, `rhs[`: 6,
`rhs_compensation[`: 4, `target[`: 4, `residual[`: 5, and
`residual_compensation[`: 4. All four old member names are forbidden
repository-wide. The private P9C5 `PackSolver` two-space member-anchor count
changes from 85 to 83 solely because four private member declarations become
one private struct declaration plus one member.

P1 and P2 must retain every Trace D hash in all three configurations without
reblessing. A reciprocal-multiply mutation must make the non-dyadic record
red. Replacing diagonal compensation or residual storage with `target` must
make the indexed-role gate red; removing the residual reset must make the
six-fill gate red and the repeated relaxation record red. All mutations are
reversed with exact source hashes before acceptance.

### S1 PackStepper ownership and substeps addendum

The following implementation shape, focused counts, and mutations are
registered after read-only review and before the first S1 source or test edit.

`substeps-name-contradicts-code` retains its already decided public meaning and
name: `substeps=N` performs `N` full-`dt` advances at
`time + substep * dt` under one electrical solve and one thermal assembly. It
does not divide `dt`, rename the argument, or re-solve inside the loop.

An exact SHORT test uses two independently built one-cell, one-lane,
isothermal NCH=5 batches, binary-exact nonzero current `8.0 A`, `N=4`, and
`dt=0.125 s`. One stepper receives one call with `substeps=4`; the other
receives four calls with `substeps=1` at times `{0, 0.125, 0.25, 0.375}`.
Both published single-cell currents must equal exactly `8.0`, both elapsed-time
rows must equal exactly `0.5`, and the complete arenas must have equal size and
be byte-identical. The case has exactly 16 assertions, so the frozen
PackStepper floor becomes 230 assertions / 9 cases from 214/8. Dividing the
inner `dt` by `substeps` must turn its elapsed-time and arena comparisons red.

`unique-in-prefix-three-spellings` gets one cold owner in
`PackTopologyInternal.hpp`:

```cpp
template <class Range, class Projection = std::identity>
  requires std::ranges::random_access_range<const Range>
        && std::ranges::sized_range<const Range>
[[nodiscard]] constexpr bool firstOccurrence(
  const Range &range, std::size_t index, Projection projection = {})
```

It returns false out of range, otherwise compares the projected current value
only with `[begin, current)`. There is exactly one owner, two PackSolver
consumers (archetype string and projected Thevenin identity), and one
PackStepper pointer consumer. The duplicate-archetype unit fixture is
strengthened so duplicate names are its only defect: both cells use that name
and two distinct valid one-lane views are supplied. The existing aliased-view
fixture independently covers the projected identity consumer. Making the
owner always return true must turn both numerical gates red; removing only the
PackStepper guard must turn the structural gate red because downstream
validation otherwise masks it.

`packstepper-gather-scatter` gets exactly two private owners:
`gatherStates(std::span<real_t>) const` and
`scatterStates(std::span<const real_t>)`. The former owns the one
state-to-contiguous `memcpy` loop and the latter the inverse. Internal
save/restore and public checkpoint/restore retain their current validation,
solver-vector, diagnostic, heat, and invalidation behavior while calling these
owners. Compacted `PackStepper.cpp` therefore contains three `gatherStates(`
tokens, three `scatterStates(` tokens, and exactly two `std::memcpy(` tokens.
Reversing either copy direction or wiring public checkpoint to the internal
buffer must turn the exact helper/call gate red and the existing two-batch
checkpoint oracle red.

The structural gate additionally pins:

- `firstOccurrence(` as owner 1 / PackSolver 2 / PackStepper 1, the three
  negated rejection tokens, and absence of all three former prefix searches;
- one gather and one scatter declaration in the private header slice, exact
  helper bodies and all four caller arguments, and the 3/3/2 token census;
- one `for(int substep...)` loop, two exact
  `time + substep * dt, dt` call suffixes, and exactly one solver call and one
  thermal assembly before the loop with zero of either inside it;
- the raw public-header phrases ``substeps * dt`, NOT by `dt` `` and
  `held frozen`.

The two new private methods change the PackStepper two-space P9C5 member anchor
33→35; its namespace/API anchor remains 4. The focused acceptance band in
Debug, fast-math Release/IPO-off, and the host-C++ CUDA tree is PackStepper
230/9, PackSolver 949/30, and the aggregate structural gate 1/1, with no
recorded hash reblessed. The existing zero-accepted-step allocation binary
must remain 14/2. This is a public-contract clarification without a signature
or behavior change; add an Unreleased CHANGELOG entry because the executable
contract was not recorded when the prose first landed.

#### S1 oracle amendment after the exact-current comparator was falsified

The first pre-production run falsified two overstrong assumptions in the
original 16-assertion shape. A one-cell ladder solve published
`8.00000000000001066 A`, not exact `8.0 A`; after four separate PackStepper
calls the comparator published `8.00000000000002309 A`. Re-solving after each
state advance therefore changed low current bits and the full-arena `memcmp`
returned 1. The run passed 13/16 assertions. No production file had changed.

The public semantics and exact arena criterion remain registered; only the
comparator is corrected before its next run. The replacement compares one
PackStepper `substeps=4` call with four direct `EulerLegacy::step` calls on the
independently built batch, using the PackStepper's one published current
divided by electrode area as the frozen current-density input. This directly
isolates the inner-loop contract without introducing four additional
electrical solves. The published current vector must contain exactly one
value, and KCL conservation must satisfy the already supplied solver band
`abs(I_cell - 8.0 A) <= current_tolerance = 1e-10 A`; no post-hoc tighter band
is introduced. Both elapsed-time rows remain exact `0.5`, and equal arena size
plus byte identity remain decisive. One PackStepper configure assertion is
replaced by one EulerLegacy configure assertion, while the two exact-current
assertions become size and KCL-band assertions, so the registered count stays
16 and the accepted PackStepper floor remains 230/9. Dividing the direct or
inner full `dt` by `substeps` must still turn elapsed time and arena identity
red.

### S1.1 moved-owner and independent-oracle hardening addendum

This addendum is registered after S1's read-only adversarial review and three
independent read-only design audits, before the first S1.1 test or production
edit and before any S1.1 binary is run. It owns one newly discovered
high-severity bug and the three test-strength gaps recorded by S1. None is
retroactively attributed to S1.

#### Explicit move contract and failing-test-first boundary

`PackSolver` and `PackStepper` currently rely on compiler-generated moves.
Their owning vectors and `SolverWorkspace::impl_` move away, but plain scalar
validity flags remain true. A moved-from sparse solver can therefore pass
`solveImpl`'s configured gate and dereference a null `workspace_.impl_`;
`PackStepper::solveElectrical` and `step` can reach the same invalid owner.
`SolverWorkspace` also continues to report its copied validity and
factorisation counters after its implementation pointer has moved. This is
the same invalid-owner class excluded at the public mutation seam by P9-B17,
not a previously dispositioned survivor.

The fix is an explicit, public rule-of-five boundary:

- `PackSolver` and `PackStepper` each declare a default constructor, delete
  copy construction/assignment, and declare no-throw move
  construction/assignment;
- `SolverWorkspace` retains its existing public declarations but replaces
  both defaulted moves with self-guarded memberwise moves;
- every owned member is transferred in declaration order, assignment accepts
  `x = std::move(x)` as an exact no-op, and all validity flags are transferred
  with `std::exchange(..., false)`;
- a moved-from solver reports an invalid, zero-counter workspace, zero
  workers, empty/zero solution and diagnostics, and rejects `solve` and
  `setRelaxationGain` with `Invalid_parameters`;
- a moved-from stepper reports zero checkpoint size/workers, empty
  solution/diagnostics/heat spans, and rejects `checkpoint`, `restore`,
  `solveElectrical`, `step`, and `stepExponential` with
  `Invalid_parameters`;
- the destination retains exact solution, diagnostics, workspace counters,
  worker count, checkpoint bytes, thermal publications, relaxation gain, and
  all state needed for the next solve/step; a moved-from source remains
  destructible and reconfigurable; and
- concurrent move versus `solve`/`step` remains outside the thread-safety
  contract. Only quiescent ownership transfer is registered.

Compile-time tests require `SolverWorkspace`, `PackSolver`, and `PackStepper`
to remain no-throw move-constructible and move-assignable; both outer owners
must remain non-copyable. Member-level no-throw proofs cover the complete
container/owner sets so the outer `noexcept` spelling cannot conceal future
throwing member drift. Five explicit public special-member declarations raise
the P9C5 two-space member anchors `PackSolver` 83 -> 88 and `PackStepper`
35 -> 40; namespace/API anchors remain 16 and 4.

Four separate SHORT cases make constructor and assignment failures
independent:

1. PackSolver move construction, 15 assertions / 1 case: configure and
   sparse-solve one affine cell, snapshot every public result and workspace
   counter, move-construct, prove exact destination continuity, then require
   the source's `setRelaxationGain(0.5)` to reject before probing the
   crash-prone `solve` path. The destination solves again to identical bits.
2. PackSolver move assignment, 20 / 1: move the same source into an already
   configured, differently shaped two-cell destination. In addition to the
   source/destination invariants, the old affine callback count must stay
   fixed while the transferred callback advances exactly once on the next
   warm solve.
3. PackStepper move construction, 38 / 1: advance a heterogeneous two-batch
   thermal pack beside an independent control, snapshot checkpoint bytes,
   solution, diagnostics, heat, workspace counters, and workers, then
   move-construct. The moved-from checkpoint uses its own reported size and
   must reject before any crash-prone call; the destination's next step must
   remain bit-identical to the control.
4. PackStepper move assignment, 44 / 1: move the same fixture into an already
   configured one-cell destination, prove the old external arena remains
   untouched, apply all moved-from checks independently, and advance beside
   the control with exact checkpoint and observable equality.

The fatal pre-fix discriminators are deliberately non-UB:
`PackSolver::setRelaxationGain(0.5)` currently returns `Success` because
`configured_` was copied; zero-sized `PackStepper::checkpoint` currently
returns `Success` for the same reason. Each appears before a moved-from
`solve`/`step` probe. Both constructor and assignment live in separate test
cases so Catch2 continues to the second red boundary after the first fatal
assertion. No signal/SEH crash is treated as evidence.

The move-only deltas are therefore PackSolver +35 assertions / +2 cases and
PackStepper +82 / +2. Starting from the S1 floors, their intermediate
move-oracle floors are 984/32 and 312/11.

#### Independent checkpoint, thermal, and prefix oracles

The checkpoint-layout case is not a round trip. It builds two one-lane
archetypes with different raw arena sizes (isothermal NCH=5 and thermal
NCH=12), fills them with distinct exact-integer sentinels, and requires the
public checkpoint to equal an independently concatenated
`batch0.raw || batch1.raw`. It then restores a separately generated wire
vector and independently compares each batch with its correct slice. The
case is exactly 11 assertions / 1 case and makes paired wrong
gather/scatter permutations red.

The frozen-thermal case is analytic and unit-checked. A two-lane NCH=5
thermal batch begins at `T0={300,310} K` in two series cells joined by
`G=2 W K^-1`, with zero applied current, no boundary, and zero convective
surface term. Both electrode OCV curves are replaced by the same constant
zero curve and both entropic curves are empty. The ladder solution therefore
has bit-exact zero current, while reaction, reversible, ohmic, and
environmental heat are all bit-exact zero; the test also checks that both
generated-heat-energy states remain zero. The initial link heat is

```text
q0 = G (310 - 300) = 20 W
Ccell = rho cp V = 1626 kg m^-3 * 750 J kg^-1 K^-1 * 1e-4 m^3
      = 121.95 J K^-1
```

One `substeps=4`, `dt=0.25 s` PackStepper call must publish
`cellExternalHeat={+20,-20} W`, bit-exact zero cell currents and generated
heat energies, and final temperatures within `1e-10 K` of
`{300 + 4*dt*q0/Ccell, 310 - 4*dt*q0/Ccell}`. This is exactly 9 assertions /
1 case. Reassembling after each Euler substep instead follows a contracting
temperature-difference recurrence and differs by about `2e-3 K`, far outside
the registered band.

Prefix coverage becomes non-adjacent:

- a compile-time `{7,11,7}` probe requires `firstOccurrence` true at indices
  0 and 1, false at duplicate index 2, and false out of range at 3;
- the duplicate-archetype and Thevenin-identity fixtures become valid
  three-entry `[a,b,a]` / `[shared,middle,shared]` inputs without changing
  PackSolver's assertion count; and
- the PackStepper pointer-alias fixture becomes
  `[shared,middle,shared]`, adding only the middle batch's successful build
  assertion.

The independent-oracle delta is therefore PackStepper +21 assertions / +2
cases and zero PackSolver or PackTopology runtime-count change. The complete
S1.1 focused acceptance floors, in Debug, fast-math Release/IPO-off, and the
host-C++ CUDA tree, are:

```text
PackSolver       984 assertions / 32 cases
PackStepper      333 assertions / 13 cases
P2-G1 allocation  14 assertions /  2 cases
structural         1 assertion  /  1 case
```

No existing pack trace or factory/ageing hash may be reblessed. The aggregate
structural gate pins all public declarations, source-reset exchanges,
self-move guards, independent prefix consumers, and the frozen
solve/assembly placement. The old exact S1 gather/scatter owners remain
unchanged.

#### S1.1 controlled mutations

At a clean committed test/implementation boundary, at minimum:

1. copy rather than exchange `PackSolver::configured_`;
2. default-move or copy `SolverWorkspace::valid_`;
3. copy rather than exchange `PackStepper::configured_`;
4. remove one move-assignment self guard;
5. omit one destination-critical owner (`thevenin_`/workspace/executor or a
   PackStepper checkpoint/stepper owner);
6. pair the wrong checkpoint gather and scatter order;
7. reassemble and republish thermal heat inside every substep; and
8. make `firstOccurrence` accept the non-adjacent duplicate or out-of-range
   index

must each turn its registered behavioral, compile-time, or structural gate
red. Every mutation is reversed to exact pre-mutation hashes before the
three-configuration acceptance run.

### S2 source-step rollback exact fixture

The following analytic fixture is fixed before its first test implementation
or run. One scripted affine batch has two parallel one-ohm cells and emits
these OCV rows by callback:

```text
{2,0}, {2,0}, {0,0}, {0,0}, {0,0}, {4,0}, {0,0}
```

With ladder mode and tolerance `0.5`, an initial zero-current solve with at
most two iterations publishes bitwise currents `{1,-1}`, nodes `{1,0}`, and
terminal voltage `1`. The test then calls `invalidate()` so those published
values remain but the internal warm-state flag is false. A caller-requested
`I=8`, one-iteration solve must fail: its direct attempt produces `{4,4}` and
maximum change `4`; hidden source stepping then succeeds at `I=0` with
`{0,0}`, succeeds at `I=1` with `{0.5,0.5}`, and fails at `I=2` with
`{3,-1}` and maximum change `2.5`. Exactly six callbacks have occurred.

Rollback must restore the three published solution fields and the false
warm-state flag. Its diagnostics deliberately describe the failed direct
caller attempt, exactly:

```text
iterations=1, numeric_factorizations=0, symbolic_factorizations=0,
jacobian_refreshes=0, source_steps=0, residual_norm=4,
constraint_drift=0, constraint_bound=0, relaxation_gain=0
```

Without diagnostics rollback the observable residual is `2.5`, describing
the abandoned later hidden solve. A final zero-current, one-iteration solve
uses the seventh `{0,0}` row and must succeed; this is the independent probe
that `has_solution_` returned to false rather than warming from `{1,-1}`.
Every expected scalar is binary-exact.

Three mutations are independently decisive: removing only solution-field
restoration turns exactly those solution checks red; removing only diagnostics
restoration changes residual `4` to `2.5`; removing only the warm-state restore
makes the final one-iteration probe fail. Production changes are limited to a
POD diagnostics snapshot beside the existing rollback snapshot and its restore
before workspace invalidation; no allocation or header/member change is
authorized.

- Recording refactors retain byte-identical encoded headers, payloads, CRCs, and
  CSV text on the existing fixtures. Added CSV value tests parse every data cell
  and require `max_digits10` round-trip equality.
- `substeps` is pinned by comparing one `N`-substep call with `N` explicit
  full-`dt` calls under one frozen electrical current: elapsed time is exactly
  `N*dt` and state is bit-identical.
- Extend `tests/structural/p9c_architecture.cmake` so
  `src/core/PackTopology.cpp` has exactly one `BranchGraph buildBranchGraph(`
  definition, exactly two call sites (three `buildBranchGraph(` tokens total),
  and exactly two adjacency `push_back` statements, both inside that owner.
  Reintroducing either former graph derivation must make the count gate red.
- Extend `tests/structural/p9c_architecture.cmake` so
  `src/core/detail/RecordingFormatCommon.hpp` has exactly one definition each
  of `std::uint32_t crc32(`, `endian_marker`, and
  `allocationFailureStatus(`, while `RecordingFormat.cpp`,
  `Recorder.cpp`, and `detail/AsyncRecordingFormat.hpp` have zero definitions
  of those three owners. Call-site uses are not forbidden.
- Extend `tests/structural/pc10_single_source.cmake` to load compacted
  `SpmScalarKernels.hpp`, `Sei.hpp`, and `SurfaceCrack.hpp`.
  `SpmScalarKernels.hpp` has exactly one `seiKineticCurrent` definition;
  `Sei.hpp` has exactly five `spm_scalar::activatedValue(` calls, exactly three
  `spm_scalar::seiKineticCurrent(` calls, and zero definitions (four
  `seiKineticCurrent` name tokens total across owner and consumer).
  `SurfaceCrack.hpp` has exactly one `spm_scalar::arrheniusFactor(` call and
  zero occurrences of the compact old association
  `p.model5_k_activation/p.Rg*(Real{1}/p.reference_temperature-Real{1}/T)`.
  Reintroducing any eliminated private expression must make the gate red.

### A1 SEI owner amendment after the registered exact gate falsified a function boundary

The function-shaped `seiKineticCurrent` extraction above is retained as the
original registration, but its digit-identity hypothesis is **FALSIFIED** for
the Release/ThinLTO capture mode. Three genuinely different implementations
all produced the same changed all-ageing fingerprint
`01ada7c263317e7f / b2d2672a19e42c12` instead of
`0be580cf849e57e1 / 4c6451971b74f789`, and changed the factory initial
derivative SHA-256 from
`3f5de3d1146e51b5d8d5755130753f329944a7c9cada50e5909aa348e76ed178`
to
`3c712d4b1d5c625c29f8cb2b43468da904f103847ad3c371cd6601e243a4c2b6`:

1. the preregistered `<Real, Scalar>` function with `const Real &` operands;
2. the same function forced `always_inline`; and
3. a fully deduced, by-value, `auto`-returning function.

Keeping the five `activatedValue` calls while restoring only the three kinetic
expressions restored both exact fixtures, isolating the optimizer-sensitive
boundary. Therefore, before trying the replacement, A1 registers the same
expression-macro remedy already required by the modal kernel:

- `SpmScalarKernels.hpp` owns exactly one
  `SLIDE_SPM_SEI_KINETIC_CURRENT` definition;
- `Sei.hpp` contains exactly three calls and no definition (four name tokens
  total across owner and consumer);
- the macro accepts the caller's `exp` token so double and ADL `Dual` retain
  the existing caller expression tree;
- the five `spm_scalar::activatedValue(` counts and all old-private-expression
  prohibitions remain as registered;
- Debug, Release, and Release/ThinLTO all-ageing and factory fingerprints must
  remain unchanged. No changed hash is authorized.

This amendment changes the structural spelling, not the finding's disposition
or behavior contract. It is registered before the first macro build/run.
- The final disposition-table census is executable: JSON ID count = 68, unique
  ID count = 68, and the table matches that set exactly. Its supplemental set
  must equal exactly
  `{surfacecrack-arrhenius-association, clang-ofast-deprecated,
  benchmark-lp-size-narrowing}`—not merely have size three. Every state parses
  as one of the three allowed forms and every deferral names its registered
  PLAN owner.

## Build-policy gate

In a fresh Clang Release tree after B1:

- generated C++ commands contain **zero** `-Ofast`;
- Release commands retain `-O3` (from CMake's Release defaults) and explicit
  `-ffast-math`;
- build output contains zero `-Ofast is deprecated` diagnostics and zero
  `benchmark_LP_cases.cpp` `size_t`-to-`int` narrowing diagnostics; and
- all recorded fast-math hashes remain at their preregistered values unless a
  behavior-change row above explicitly owns the delta.

No claim is made that GCC deprecates `-Ofast`; its spelling is reviewed
separately and is not changed merely because Clang 21 warns.

## Mutation/adversarial gate

At minimum, one controlled mutation must turn red for each new behavior/test
family: fast-math finiteness, SurfaceCrack PC-10 ownership, direct Euler invalid
span, source-step rollback, zero CRC, parameter non-finite sampling, constant
aliasing, drive-cycle right-endpoint event resolution, Simulation truncation,
CUDA host rejection, CUDA cadence/monotonicity, single-worker exception
selection, topology shared derivation, and recording shared CRC. Mutations run
only at a clean committed batch boundary. Hash every touched file before the
mutation; after reversing it, the hashes must match exactly, `git diff
--exit-code HEAD -- <touched-files>` must return zero, and `git status
--porcelain` must be empty. `git diff --check` is retained only as a whitespace
check and is never restoration evidence.

The adversarial review also tries to make the final ledger pass with a duplicate,
missing, or unknown ID and with an unowned deferral; each mutation must fail the
census.

## Final gates

After all MQ.2 batches:

1. fresh tree-local Debug, fast-math Release/IPO-off, and CUDA Release with
   host-C++ ThinLTO each discover exactly 58 CTest tests and pass 58/58;
2. CPU lanes select `CudaDisabled`; CUDA selects `CudaSpmBatch`;
3. no test command resolves outside its own build tree except the three
   intentional current-source structural scripts;
4. every pre-existing Catch2 test case remains, all registered hash deltas are
   explained, and no unrelated recorded hash moves;
5. the CUDA binary has exactly six test cases (baseline four plus host-validation
   and recorder-cadence cases), passes all assertions, and its pre-existing
   433,671 assertions remain a floor rather than a replacement target;
6. a final no-op rebuild in every lane reports no work;
7. `CHANGELOG.md`, PLAN §8, and `develop/TODO.md` record applied/refuted/deferred
   counts and name every deferred owner; and
8. no timing, sanitizer, coverage, hosted-CI, installed-package, Linux, macOS,
   or device-LTO claim is inferred from these gates.
