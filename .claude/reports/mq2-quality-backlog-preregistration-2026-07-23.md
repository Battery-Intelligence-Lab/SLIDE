# MQ.2 quality-backlog disposition — preregistration (2026-07-23)

## Scope and source identity

This registration precedes every MQ.2 source/test edit and every decisive MQ.2 run.
The source baseline is clean commit `940716b` (`Close MQ.1 with the restored
three-lane baseline`). MQ.1 established 58/58 in fresh Debug, fast-math Release,
and CUDA lanes at source commit `a49224e`; `940716b` changes evidence and standing
documentation only.

MQ.2's original baseline has **71 decisions**:

- the 68 unique IDs in
  `.claude/reports/code-quality-pass-2026-07-21-survivors.json`;
- `surfacecrack-arrhenius-association`, the separately recorded `SurfaceCrack`
  Arrhenius association;
- `clang-ofast-deprecated`, MQ.1's deprecated redundant `-Ofast` observation;
  and
- `benchmark-lp-size-narrowing`, MQ.1's legacy `benchmark_LP_cases`
  `size_t`-to-`int` narrowing diagnostic.

The final disposition table must contain all 68 JSON IDs exactly once and all
three original supplemental IDs exactly once. Post-baseline discoveries are
tracked in the append-only stable-ID registry added below; the current combined
set has 84 IDs without rewriting the original 71-ID baseline. Each row has exactly one state:
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

#### Post-baseline supplemental registry

The following stable IDs were assigned on 2026-07-29 after S1.1 acceptance.
This is a census-only amendment: it changes no frozen oracle, expected value,
pass band, implementation, or prior disposition. It prevents findings
discovered after the original 68+3 baseline from disappearing behind a
changing count:

| Stable ID | Finding | Disposition / owner |
|---|---|---|
| `moved-owner-validity-after-defaulted-move` | moved `PackSolver`/`PackStepper` sources retain validity after their owners move away | APPLIED by S1.1 at `574cfaf` |
| `checkpoint-roundtrip-oracle-symmetry-gap` | the prior checkpoint round trip could accept paired wrong layouts | APPLIED test-only by S1.1; independent layout mutation red |
| `thermal-assembly-placement-oracle-gap` | structural placement had no heterogeneous numerical cadence oracle | APPLIED test-only by S1.1; analytic recurrence mutation red |
| `prefix-duplicate-adjacency-oracle-gap` | prefix fixtures covered only adjacent duplicates | APPLIED test-only by S1.1; non-adjacent/OOB mutation red |
| `eigen-strong-inline-redefinition-warning` | a fresh Debug build reports the pre-existing Eigen macro-redefinition diagnostic | DEFERRED with explicit owner MQ.2 B1 build diagnostics |
| `recorder-schema-default-vector-debug-oom-terminates` | R1's initially frozen bare schema-vector constructor allocates a Debug iterator proxy inside a `noexcept` constructor | APPLIED at `5a462b6`; RecorderAllocation 33/3 is green without an escaped exception |
| `recording-owner-contract-comments-stale` | R1 moved CRC ownership but two leading production comments still attribute all CRC ownership/co-location to `RecordingFormat.cpp` | APPLIED at `a7a4236`; the leading comments name the split owners truthfully |
| `r1-self-containment-gate-overclaim` | 9C-5 compiler evidence was described as proving every direct standard include rather than only load-bearing standalone dependencies | APPLIED at `a7a4236`; claim narrowed and lexical/compiler roles separated |
| `r1-helper-internal-linkage-gate-gap` | exact helper bodies were pinned but their required anonymous-namespace linkage was not | APPLIED at `a7a4236`; mutation 14 makes the linkage assertion red |
| `r1-linux-debug-byte-oracle-unregistered` | R1's complete byte oracle treated every non-Release build as the registered Windows Debug fixture, so the fresh Linux/Clang 18 coverage build selected the wrong exact tuple | APPLIED at `7b8b1f8` after exact `ad700ca` old-production reproduction |
| `p9c5-script-policy-unset` | the aggregate `cmake -P` gate uses `IN_LIST` without setting CMP0057, which the fresh retained Linux CMake 3.31.10 invocation rejects before 9C-5 evaluation | APPLIED at `7b8b1f8`; direct and aggregate 9C-5 entry paths pass |
| `status-exception-manifest-stale-after-packsolver-split` | the exact Status exception manifest retains six pre-S1.1/S2 PackSolver/SpectralModel source identities, so the authoritative reporter refuses to evaluate complete fresh profiles | APPLIED at `7b8b1f8`; six identities relocated one-for-one and exact coverage passes |
| `instrumented-workflow-test-count-stale` | both hosted instrumented workflows still require 57 discovered test commands after the committed suite reached 58 | DEFERRED with explicit owner MQ.2 B1 build diagnostics; hosted workflows are not claimed run |

This post-baseline registry is append-only during MQ.2. A new finding must
receive a stable ID here before any disposition claim. At wave closeout, the
disposition table must match both the original baseline sets and this
post-baseline set exactly.

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

#### S2 executable count and acceptance amendment

This count amendment is registered on 2026-07-29 before the first S2 test
edit and before any S2 binary run. It changes none of the analytic values or
production limits above.

One test named `failed source stepping restores caller-attempt diagnostics`
with tags `[core][pack][solver][source-stepping][rollback][MQ.2][S2]` has
exactly ten assertions / one case:

1. the existing `compile(...)` helper requires successful compilation of the
   two-parallel-cell topology;
2. solver configuration succeeds;
3. the initial zero-current, two-iteration ladder solve succeeds;
4. all three initially accepted solution fields bit-match
   `{cell_current={1,-1}, node_voltage={1,0}, terminal_voltage=1}`;
5. after `invalidate()`, the caller's `I=8`, one-iteration solve returns
   `Numerical_failure`;
6. the scripted callback count is exactly six;
7. all three published solution fields still bit-match the accepted snapshot;
8. all nine diagnostics fields bit-match the registered failed direct attempt,
   including `iterations=1`, `source_steps=0`, and `residual_norm=4`;
9. a final zero-current, one-iteration solve succeeds; and
10. the callback count is exactly seven.

The last successful status is the independent warm-flag discriminator: if
rollback leaves `has_solution_` true, the seventh `{0,0}` row begins from
`{1,-1}`, changes by `1 > 0.5`, and cannot converge in one iteration.

The current PackSolver floor is 984 assertions / 32 cases, so the S2 floor is
exactly **994 / 33**. In Debug, fast-math Release/IPO-off, and the retained
host-C++ CUDA tree, the focused acceptance set is PackTopology 148/9,
PackSolver 994/33, PackStepper 333/13, P2-G1 allocation 14/2, and structural
aggregate 1/1. The restored implementation must also pass the unfiltered
58/58 CTest suite in all three trees.

The frozen old-production run is expected to reach 9/10 assertions and fail
only the exact diagnostics comparison with residual `2.5` instead of `4`.
No oracle is reblessed if another assertion fails. After the local
snapshot/restore fix is committed, three independent mutations remove only:

- the three solution-field restores;
- `diagnostics_ = rollback_diagnostics`; and
- `has_solution_ = rollback_has_solution`.

They must make assertions 7, 8, and 9 red respectively, with inverse patches
and exact source-hash restoration before the three-configuration acceptance
run.

### T1 topology derivation exact fixtures and structural ownership

This amendment is registered on 2026-07-29 before the first T1 test/source
edit and before any T1 binary run. It owns exactly
`branch-graph-derived-twice`, `three-parallel-archetype-maps`, and
`test-gap-ladder-rollback`.

The no-op oracle compiles this exact description:

```text
series({
  cell({archetype="zeta", thermal=false}),
  parallel({
    cell({archetype="alpha", thermal=true}),
    cell({archetype="zeta", thermal=false})
  })
})
```

One test named `topology derivation preserves exact graph and batch ordering`
with tags `[core][pack][compile][topology][MQ.2][T1]` has exactly eight
assertions / one case:

1. compilation succeeds;
2. the three complete cell records are exactly
   `{s00,zeta,batch=1,lane=0,false}`,
   `{s01.p00,alpha,batch=0,lane=0,true}`, and
   `{s01.p01,zeta,batch=1,lane=1,false}`;
3. `batch_archetypes` is exactly `{"alpha","zeta"}`;
4. electrical scalars are exactly `node_count=3`, terminals `0/1`,
   `connected=true`, `index1_candidate=true`, and
   `series_parallel_ladder=true`;
5. the complete ordered branch records are exactly
   `{0,2,cell,cell=0,R=0}`, `{2,1,cell,cell=1,R=0}`, and
   `{2,1,cell,cell=2,R=0}`;
6. nodal sparsity is exactly
   `{(0,0),(0,2),(1,1),(1,2),(2,2)}`;
7. ladder metadata is exactly offsets `{0,1,3}`, cells `{0,1,2}`, and
   nodes `{0,2,1}`; and
8. the empty thermal graph has exactly `cell_count=3`, `boundary_count=0`,
   no edges/incidents/edge scratch, offsets `{0,0,0,0}`, and three zero
   endpoint-scratch values.

Every scalar and sequence is equality-tested directly; no implementation
helper or derived expectation is shared with production. This fixture must
pass against the old production source before extraction, and remain
digit-identical afterward.

One test named `import finalization clears every ladder vector after an
orientation mismatch` with tags
`[core][pack][import][ladder][rollback][MQ.2][T1]` has exactly seven
assertions / one case. It compiles a two-cell series ladder, checks the
initial ladder flag, reverses only the second branch, requires imported
finalization to succeed as a connected non-ladder, requires the final ladder
flag false, and independently requires each of `ladder_offsets`,
`ladder_cells`, and `ladder_nodes` empty. The test is always enabled; it is
not hidden behind `assert`, Debug-only configuration, or an exceptional path.

The existing `electrical validator rejects independently corrupted metadata`
case gains exactly one direct validator assertion: setting the first branch's
positive endpoint to `node_count` returns `Invalid_parameters`. This is
deliberately separate from `finalizeImportedPackTopology`, whose own endpoint
precheck otherwise masks whether the shared graph builder is called too soon.

The current PackTopology floor is 148 assertions / 9 cases. These additions
are exactly +16 assertions / +2 cases, so the frozen T1 floor is
**164 / 11**. In Debug, fast-math Release/IPO-off, and the retained host-C++
CUDA tree, the final focused acceptance set is PackTopology 164/11,
PackSolver 994/33, PackStepper 333/13, P2-G1 allocation 14/2, and aggregate
structural 1/1, followed by unfiltered 58/58 in all three trees.

The production extraction is file-local and preserves all public types and
statuses:

- exactly one `BranchGraph` owns adjacency and sorted/unique nodal sparsity;
- exactly one `buildBranchGraph(` definition serves the compiler and direct
  validator, for three name tokens total, and its body owns exactly the two
  undirected adjacency `push_back` statements; the builder documents/asserts
  its already-validated endpoint precondition and retains the existing
  `sizeof(CompiledElectricalBranch) > 3` reserve-overflow proof;
- exactly one `isConnectedFrom(` definition serves those same two callers,
  also for three name tokens total;
- the validator retains its initial node/terminal/branch-count allocation
  guards (including `node_count - 1 <= branches.size()`), then all endpoint,
  kind, resistance, and cell-identity checks in its original branch loop,
  completes that loop and the complete cell census, and only then calls
  `buildBranchGraph`; sparsity comparison still precedes connectivity
  comparison;
- exactly one `BatchSlot` record owns batch index, next lane, and thermal
  composition, and `assignBatchLocations` uses exactly one sorted
  `std::map<std::string, BatchSlot>` with one `try_emplace`, one `at`, and one
  post-increment of `next_lane`;
- the three former parallel maps for `batches`, `thermal`, and `next_lane`
  have zero definitions; and
- mixed thermal composition for one archetype remains rejected before any
  locations are published, while lexicographic batch order and encounter-order
  lanes remain unchanged.

`tests/structural/p9c_architecture.cmake` loads compacted
`src/core/PackTopology.cpp`, counts the owners/callers and one-map vocabulary,
forbids all three old map declarations, slices `buildBranchGraph` to prove
both adjacency insertions live there, and slices
`validateElectricalNetlist` to require this order:

```text
branch loop -> endpoint guard -> complete-cell guard -> buildBranchGraph
            -> sparsity comparison -> isConnectedFrom
```

The validator slice also proves the initial size-amplification guard precedes
the branch loop. This preserves the P9-B18 fix: neither an input-sized graph
allocation nor adjacency indexing occurs before its respective hostile-input
bound has been established. The extraction claims published-value identity,
not allocation-attempt identity; moving the reserve into the shared cold
builder may change allocation scheduling.

The same gate slices `assignBatchLocations`, requires exactly one
`std::map<` token in that slice, and requires this order:

```text
try_emplace -> mixed-thermal rejection -> sorted batch publication
            -> slots.at -> location assignment with next_lane++
```

This structurally proves the validation pass completes before any batch/lane
location is published; the existing mixed-composition behavior case proves
the rejection remains live.

The test-only boundary deliberately leaves this future structural gate red
against old production while the 164/11 behavior gate is green. No structural
expectation is weakened or reblessed to fit the implementation.

At the clean production boundary, at least these independent mutations are
required before final acceptance:

1. reintroduce either former adjacency/sparsity derivation outside
   `buildBranchGraph` (structural owner/count gate red);
2. make `isConnectedFrom` accept the disconnected direct-validator fixture;
3. drop or alter a sparsity contribution (exact fixture or validator red);
4. move the validator's graph call before its endpoint loop (ordering gate
   red before any unsafe binary run);
5. change `next_lane++` to `++next_lane` (exact location fixture red);
6. accept mixed thermal composition for one archetype (existing imported
   topology case red); and
7. omit each of the three orientation-failure ladder clears in turn (the
   corresponding independent rollback assertion red).

Every mutation is inverse-patched to the exact registered source/test hashes,
then `git diff --exit-code HEAD -- <touched-files>` and an empty porcelain
status prove restoration before the three-configuration acceptance run.

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
  benchmark-lp-size-narrowing}`—not merely have size three. A separate
  post-baseline set must equal exactly
  `{moved-owner-validity-after-defaulted-move,
  checkpoint-roundtrip-oracle-symmetry-gap,
  thermal-assembly-placement-oracle-gap,
  prefix-duplicate-adjacency-oracle-gap,
  eigen-strong-inline-redefinition-warning,
  recorder-schema-default-vector-debug-oom-terminates,
  recording-owner-contract-comments-stale,
  r1-self-containment-gate-overclaim,
  r1-helper-internal-linkage-gate-gap,
  r1-linux-debug-byte-oracle-unregistered,
  p9c5-script-policy-unset,
  status-exception-manifest-stale-after-packsolver-split,
  instrumented-workflow-test-count-stale}` at this boundary and must grow by
  explicit registry amendment, never by count alone. Every state parses as
  one of the three allowed forms and every deferral names its registered PLAN
  owner.

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

## T2 thermal-topology amendment (registered 2026-07-29)

This amendment is registered before the first T2 production, unit-test, or
structural-gate edit and before any T2 binary run. The prior-art audit found no
FALSIFIED, REFUTED, killed, or deferred record for either
`assemble-noexcept-and-stale-cold-brief` or `test-gap-isothermal-zero`. This
batch owns exactly those two original survivor IDs.

### Executable isothermal oracle

Add one always-on test named
`isothermal thermal graph assembles exact zero heat` with tags
`[core][pack][thermal][oracle][MQ.2][T2]`. It compiles three thermal cells in
parallel, boundaries `a` and `b`, and exactly these five links:

```text
p00 --1.50 W/K-- p01
p01 --2.50 W/K-- p02
p00 --0.50 W/K-- p02
p00 --0.75 W/K-- a
p02 --3.25 W/K-- b
```

The canonical endpoint degrees are `{3,2,3,1,1}` and the exact incident
offsets are `{0,3,5,8,9,10}`. The case has exactly six executed assertions:

1. `compilePackDescription` returns `Success`;
2. one grouped shape assertion requires exactly five edges and those exact
   offsets;
3. `assemble` returns `Success`;
4. every cell `q_ext` value compares exactly equal to `0.0`;
5. every boundary-heat value compares exactly equal to `0.0`; and
6. every published `edge_flux` value compares exactly equal to `0.0`.

Before assembly, the test poisons all three output cells, both boundary
outputs, every published and trial edge-flux slot, every trial endpoint-heat
slot, and every incidence byte with distinct nonzero sentinels. Thus a
pre-zeroed compile result or omitted calculation/reset/publication cannot
make the oracle pass vacuously. Every cell and boundary temperature is the
same finite `300.0 K`. `tests/unit/core_PackTopology_test.cpp` adds its direct
`<algorithm>` include for `std::ranges::all_of`; it must not depend on a
transitive include.

The exact-equality band is analytic and unit-checked:

```text
Delta T = 300 K - 300 K = exactly 0 K
q_edge = G [W/K] * Delta T [K] = exactly 0 W
```

For finite positive conductance, the subtraction and multiplication produce
zero exactly; either signed zero compares equal to `0.0`, and fixed-order
endpoint accumulation of only signed zeros remains zero. This does not
authorize an exact general non-isothermal conservation sum: separately
rounded endpoint totals can prevent exact cancellation.

The current PackTopology floor is 164 assertions / 11 cases. This addition is
exactly +6 assertions / +1 case, so the frozen T2 floor is **170 / 12**. No
extreme-temperature repeat and no extra Catch assertion is included.

At the oracle-only boundary, old production must pass PackTopology 170/12.
The future structural gate described below must be red, with its first
diagnostic reporting that the exact declaration token occurs 0 times rather
than 1. This deliberately separates the old-production numerical oracle from
the future signature/comment contract. No `static_assert(noexcept(...))` is
added at that boundary because it would prevent the old-production numerical
oracle from compiling; the exact structural signature gate is the executable
owner of that guarantee.

### Signature and leading-header contract

The production change is restricted to the exception specification and the
leading contract comment:

- add `noexcept` to both the `CompiledThermalGraph::assemble` declaration and
  definition, with no body, arithmetic, ordering, storage, or status change;
- replace the stale all-cold brief with a leading MC-3 contract of six lines
  or fewer that names ownership, `PLAN.md` section 3.4 with D-19/D-21, and the
  hot/cold split; and
- state specifically that description authoring and
  `compilePackDescription()` are cold, while
  `CompiledThermalGraph::assemble()` is hot once per PackStepper step attempt.

The full current body has been inspected: it performs span/vector
`operator[]` access, finite/shape guards, arithmetic on `real_t = double`,
`std::fill`, and `std::copy` over pre-sized storage, then returns a
`slide::Status`. It neither allocates nor calls a throwing user operation.
Adding `noexcept` therefore documents and enforces the PC-4 no-throw boundary
at the type level without changing any successful or rejected numerical
path. It does **not** enforce PC-1: `noexcept` cannot observe a successful
allocation. The independent P2-G1 coupled-thermal allocation test remains the
owner of the zero-allocation claim.

`tests/structural/p9c_architecture.cmake` loads compacted
`src/core/PackTopology.hpp` and requires exactly one complete declaration
token ending:

```text
std::span<real_t>boundary_heat)noexcept;
```

It requires exactly one complete definition token in compacted
`src/core/PackTopology.cpp` ending:

```text
std::span<real_t>boundary_heat)noexcept{
```

The gate separately reads the raw header, removes whitespace while retaining
comments, requires exactly one each of the new brief, owner, PLAN-section,
cold, and hot fragments, and requires zero copies of the former
`cold-compiled electrical/thermal topology` brief. This T2 edit fixes only
PackTopology's stale comment; it does not claim to close the repo-wide MC-3
debt owned by MQ.3.

### Registered mutations and acceptance

At the clean implementation boundary, at least these independent mutations
must turn red before final acceptance:

1. remove both exception specifications so production still compiles:
   behavior remains green but the exact declaration/definition structural
   counts are red;
2. remove `noexcept` from only one side: the declaration/definition mismatch
   must fail to compile, while the structural gate also rejects the missing
   token;
3. restore or damage any one owner/PLAN/hot/cold comment fragment: the
   raw-header structural gate is red;
4. bias the edge temperature difference by `+1.0`: at least the exact-zero
   edge-flux assertion is red;
5. omit publication from trial edge flux to `edge_flux`: the poisoned
   published-flux sentinel makes assertion 6 red while cell and boundary
   zeros remain green;
6. omit cell-heat publication: the poisoned `q_ext` values make assertion 4
   red; and
7. omit boundary-heat publication: the poisoned boundary values make
   assertion 5 red.

Every touched header/source/test/gate file is SHA-256 anchored at the clean
implementation commit before mutations. Each mutation is inverse-patched
individually; the exact hashes, `git diff --exit-code HEAD -- <touched-files>`,
and an empty porcelain status are required before the next mutation and
before final acceptance. No oracle value may be reblessed to fit production.

The final focused acceptance set in Debug, fast-math Release/IPO-off, and the
retained host-C++ CUDA tree is PackTopology 170/12, PackSolver 994/33,
PackStepper 333/13, P2-G1 allocation 14/2, and aggregate structural 1/1.
Each lane then runs the unfiltered 58/58 suite, including the real CUDA test
only in the CUDA tree, followed by a no-op rebuild. `CHANGELOG.md`, PLAN
section 8, the validation report, and `develop/TODO.md` receive the final
applied census only after all gates and mutations are green.

## C1 netlist-parsing amendment (registered 2026-07-29)

This amendment is registered at clean HEAD `e5e0b9d` before the first C1
unit-test, structural-gate, source, or header edit and before any C1 binary
run. The prior-art and survivor-ledger audit found no FALSIFIED, REFUTED,
killed, applied, or deferred record for:

- `csv-limits-restated-in-prose`;
- `ascii-digit-scan-helper`; and
- `failsemantic-returns-unused-bool`.

These are original survivor rows 32--34. Each verifier verdict is REAL,
behavior-preserving, low-risk, and contract-compatible. C1 owns exactly these
three IDs. It does not absorb the analogous BPX/Experiment message
duplication, change the public import limits, broaden the liionpack grammar,
or reopen the killed sparse-label, zero-ohm-row, or dense-index-overflow
candidates.

The clean pre-C1 anchors are:

```text
516DFCD1141304D5B858F5A24404E9179F57997F2D63E29452D3CDE870E2BE5A  510  src/core/NetlistCsv.cpp
55D8B05F8BFA9A065C39746B22195D31B4F6283E04F38DCA2B4C76E24A6E2B11   54  src/core/NetlistCsv.hpp
D2F79313FD54AA321B89C57DC5B296F1219682EC5DFD79AD55BD83C8090F5B57  405  tests/unit/core_NetlistCsv_test.cpp
07B30AF084297B59DA1952699CD4F79F0A27452AE87E20654313A89CD7AC0004  978  tests/structural/p9c_architecture.cmake
```

The final T2 three-lane logs independently record the unchanged pre-C1
NetlistCsv floor as 821 assertions / 9 cases and ParserAllocation as
1143/13. No executable was invoked to establish this registration.

### Frozen grammar and exact diagnostic oracle

The accepted number grammar is an invariant, not an implementation choice:

```text
value    := "-"? integer fraction? exponent?
integer  := "0" | [1-9][0-9]*
fraction := "." [0-9]+
exponent := ("e" | "E") ("+" | "-")? [0-9]+
```

The complete spelling must be consumed and the subsequent finite
`std::from_chars(..., std::chars_format::general)` conversion must succeed.
A leading `+`, a leading zero followed by another digit, a missing integer,
an empty fraction, and an empty exponent remain rejected. Existing accepted
scientific-notation values remain exactly `1e-05`, `250.0`, and `1500.0`.

The existing exponent-rejection case adds `+4.2`, `01`, `.5`, and `1.` and
changes its predicate from non-Success to exact `Invalid_parameters`. Each
new spelling executes the two existing unpublished-output checks plus the
outer status check: exactly 12 new assertions and no new case. The already
present `1e`, `1e+`, `1e-`, `1e+x`, and `1.e5e5` probes remain.

Add one always-on test named
`liionpack CSV reports exact first-failure diagnostics`, tagged
`[core][pack][netlist][csv][diagnostic][MQ.2][C1]`. It contains exactly
twelve fixtures. Before every parse, `row`, `offset`, and `message` are
independently poisoned with nonzero/nonmatching values. Every fixture then
executes exactly four assertions: exact `Invalid_parameters`, logical row,
absolute byte offset, and byte-exact message.

| Triggering input/path | Row | Offset | Exact message |
|---|---:|---:|---|
| data row has three fields under a four-field header | 2 | 0 | `liionpack CSV row width differs from header` |
| descriptor `X0` | 2 | 0 | `unsupported liionpack descriptor` |
| a second `V0` descriptor | 3 | 0 | `duplicate liionpack descriptor` |
| node label `+1` | 2 | 0 | `invalid liionpack node label` |
| `V0` has identical endpoints | 2 | 0 | `liionpack element has identical endpoints` |
| value `+4.2` | 2 | 0 | `invalid liionpack element value` |
| row-3 resistor has value `-0.1` | 3 | 0 | `liionpack resistance must be non-negative` |
| a row-3 zero-ohm wire contracts the row-4 current-source endpoints | 4 | 0 | `liionpack ideal wire shorts an element or terminal source` |
| a header has 33 one-byte fields | 1 | 64 | `CSV row exceeds 32 columns` |
| a 65,537-byte unquoted first field | 1 | 65537 | `CSV field exceeds 65536 bytes` |
| a quoted first field contains 65,537 payload bytes | 1 | 65538 | `CSV field exceeds 65536 bytes` |
| `V0,1,0,4"2` under the canonical 23-byte header | 2 | 31 | `quote inside unquoted CSV field` |

The last four offsets are derived from the reader cursor, not copied from
production expectations. The 33rd one-byte field begins after 32
`field,` pairs, at absolute byte 64. An unquoted oversized field has consumed
65,537 bytes when its limit guard fires; the quoted version has additionally
consumed the opening quote. The canonical header occupies bytes `[0,23)`;
the offending quote is row-relative byte 8 and therefore absolute byte 31.
It is rejected before cursor advancement.

The existing oversized-file path poisons its diagnostic immediately before
the 4 MiB + 1 load and gains exact row 0, offset 0, and
`liionpack CSV exceeds 4194304 bytes` assertions: +3. The existing
100,001-data-row path likewise starts poisoned, replaces its message
substring check with exact equality, and adds exact logical row 100002 and
offset 0: +2. Its row follows from one header plus the first disallowed data
row.

The new exact arithmetic is therefore:

```text
821 + (4 spellings * 3) + (12 diagnostics * 4) + 3 + 2 = 886 assertions
9 + 1 = 10 test cases
```

No assertion is added to the already complete accepted-topology and
scientific-value fixtures, and no poison-output assertion is removed.
`parallel_fixture`, the exact scientific values, sentinel topology
atomicity, the independent imported-topology validator, P9-G2 deterministic
fuzz behavior, and P9-B38's `sizeof(array)-1` embedded-NUL extent remain live
owners. At the oracle-only boundary, unchanged production must pass
NetlistCsv exactly **886/10**.

### One message owner, one ASCII grammar owner, one semantic-status owner

The production change is restricted to `src/core/NetlistCsv.cpp`,
`src/core/NetlistCsv.hpp`, and `CHANGELOG.md`.

Four adjacent file-local `constexpr std::string_view` values own the four
limit diagnostics next to the existing numeric limits:

```text
msg_too_large     = "liionpack CSV exceeds 4194304 bytes"
msg_too_wide      = "CSV row exceeds 32 columns"
msg_field_large   = "CSV field exceeds 65536 bytes"
msg_too_many_rows = "liionpack CSV exceeds 100000 rows"
```

One numeric `static_assert` pins 4,194,304 bytes, 100,000 rows, 32 columns,
and 65,536 bytes/field to those literal spellings. Each literal occurs
exactly once in the translation unit. Complete owner-plus-consumer token
counts are respectively 3, 2, 3, and 2. Consumers remain exact:

- direct input assigns `msg_too_large`;
- the bounded-file aggregate explicitly constructs
  `std::string{msg_too_large}`;
- the 33-column path calls `fail(msg_too_wide)`;
- both quoted and unquoted field paths call `fail(msg_field_large)`; and
- the first disallowed data row assigns `msg_too_many_rows`.

The public header's readable `4 MiB, 100,000 rows, 32 columns, and
65,536 bytes/field` prose deliberately remains; it is a public contract, not
a replaceable C++ owner. Analogous strings in BPX and Experiment remain
outside C1.

Two file-local `constexpr` helpers transcribe, without widening, the current
ASCII grammar:

```cpp
[[nodiscard]] constexpr bool isAsciiDigit(char value) noexcept;
constexpr std::size_t scanDigits(
  std::string_view text, std::size_t &cursor) noexcept;
```

`isAsciiDigit` has the exact inclusive `'0'`/`'9'` bounds and a comment
explaining why locale-sensitive `std::isdigit` is not substituted.
`scanDigits` advances while that predicate holds and returns exactly
`cursor - begin`. It is deliberately **not** `[[nodiscard]]`: the
already-proven nonzero integer branch legitimately discards its count.

The compact translation-unit census is exactly four `isAsciiDigit(` tokens
(definition, scanner, parseValue leading-zero guard, descriptor consumer)
and four `scanDigits(` tokens (definition plus three parseValue consumers).
The descriptor call explicitly converts its `unsigned char` loop value with
`static_cast<char>(value)`.

Within a sliced `parseValue`, the structural gate requires exactly three
scanner calls, one direct digit-predicate call, the unchanged `'-'`-only
leading-sign token, the unchanged `'0'` special case, the unchanged
`'1'..'9'` branch, two zero-count checks for fraction/exponent, and the final
full-consumption check. The old raw
`text[cursor] >= '0' && text[cursor] <= '9'` spelling has zero occurrences.
This pins helper extraction and grammar simultaneously; replacing the
`'1'..'9'` branch with a generic digit predicate is not accepted as proof of
the no-leading-zero contract.

`failSemantic` becomes exactly:

```cpp
[[nodiscard]] slide::Status failSemantic(
  NetlistCsvDiagnostic &diagnostic,
  std::size_t row,
  std::string_view message)
{
  diagnostic.row = row;
  diagnostic.message.assign(message);
  return slide::Status::Invalid_parameters;
}
```

All eight consumers directly return the helper: seven pass `row`, and the
post-contraction check passes `element.row`. Thus the exact census is nine
`failSemantic(` tokens, eight `return failSemantic(` tokens, split 7/1.
There are exactly 12 remaining compact
`return slide::Status::Invalid_parameters;` tokens (the helper plus eleven
unrelated direct exits), down from 19; together with the unchanged allocation
failure return, the explicit failure-return census is 13 rather than 20.
`CsvReader::fail` remains a bool because `CsvReader::next` consumes it as a
bool.

The public `offset` comment is corrected without changing outputs: it states
that the value is the source cursor at or immediately after a
reader-detected syntax/limit failure and is zero when no cursor position is
available. The structural gate requires the new compact comment and forbids
the stale unrestricted `source byte at or immediately after the failure`
wording. It also pins the existing exact `CsvReader::fail` body.

Semantic, direct-bound, archetype, and allocation failures retain offset
zero. No synthetic row-start offset is introduced. The unit matrix makes
that zero executable from poisoned state, while the reader cases make
nonzero cursor offsets executable.

### Future structural gate and old-production boundary

Append C1 to `tests/structural/p9c_architecture.cmake`, after loading compact
`NetlistCsv.cpp` and `.hpp`. The gate requires:

- the four exact message-owner declarations, numeric pin, once-only literals,
  3/2/3/2 reference counts, and exact five consumer forms above;
- exact `>` guards for input bytes, data rows, and both field paths, plus the
  pre-field `fields.size() >= max_csv_columns` column guard;
- an ordered parse-entry slice containing `diagnostic = {}` before
  `csv.size() > max_csv_bytes` and the direct `msg_too_large` publication;
- exactly one `diagnostic.offset` token (the allocation-failure reset) and
  one `diagnostic_.offset` token in the exact `CsvReader::fail` body, so a
  new semantic/direct-limit offset write is structural-red;
- exact helper bodies, the 4/4 global helper census, absence of
  `[[nodiscard]]` on `scanDigits`, and the sliced parseValue/descriptor
  grammar constraints above;
- the exact Status-returning `failSemantic` body and 9/8/7/1 ownership
  census;
- the 12 direct Invalid-parameters returns; and
- the new truthful header offset comment with the stale wording absent.

Against unchanged production, this future gate must fail first at the first
message owner: expected one exact `msg_too_large` declaration, found zero.
All preceding included structural suites must remain green. No future
structural expectation may be weakened or reblessed to accommodate an
implementation.

The oracle-only commit may change only
`tests/unit/core_NetlistCsv_test.cpp` and
`tests/structural/p9c_architecture.cmake`. In Debug it must produce:

```text
unit_test_core_NetlistCsv                  886 assertions / 10 cases, PASS
structural_test_core_9C2AgeingKernel       RED only at C1 msg_too_large 0/1
```

The production source/header hashes must still equal the anchors above.
Release, CUDA, and full-suite claims are withheld at the oracle-only
boundary.

### Registered mutations and final acceptance

At a clean committed implementation boundary, independently run at least
these mutation families:

1. drift any one numeric limit or change its registered `>`/`>=` guard;
2. alter one limit-message byte or route a consumer through the wrong owner;
3. remove the quoted or unquoted field-limit consumer independently;
4. change either inclusive ASCII digit bound or reintroduce one raw digit
   loop;
5. break `scanDigits` advancement or its returned count;
6. admit a leading `+`, a leading zero, an empty fraction, or an empty
   exponent (structural-red is required even where `from_chars` also rejects
   the broadened spelling);
7. make `failSemantic` return `Success`;
8. change its row/message behavior, synthesize a nonzero offset, or move the
   entry diagnostic reset after the direct-size guard;
9. reintroduce one discarded two-line
   `failSemantic(...); return Invalid_parameters;` consumer;
10. damage the header offset contract or restore its stale unrestricted
    wording; and
11. change the quote guard or the independently derived byte-31 expectation.

Every touched source/header/test/gate file is SHA-256 anchored at the clean
implementation commit. Each mutation is inverse-patched individually; exact
hashes, `git diff --exit-code HEAD -- <touched-files>`, and an empty porcelain
status are required before the next mutation and before final acceptance. No
oracle count, status, row, offset, message, or accepted value may be reblessed
to fit production.

The final focused set in Debug, fast-math Release/IPO-off, and the retained
host-C++ CUDA tree is NetlistCsv 886/10, ParserAllocation 1143/13,
PackTopology 170/12, PackSolver 994/33, PackStepper 333/13, P2-G1 allocation
14/2, and aggregate structural 1/1. Every lane then runs the unfiltered
58/58 suite, including the real CUDA test only in the CUDA tree, followed by
a no-op rebuild. `clang-format --dry-run --Werror`, `git diff --check`, and an
adversarial gate review precede disposition.

Only after those gates are green do the three IDs become APPLIED. That would
move the original census from 31 APPLIED / 0 REFUTED / 8 named deferrals /
32 pending to 34 / 0 / 8 / 29, and the combined 76-ID census from
35 / 0 / 9 / 32 to 38 / 0 / 9 / 29. `CHANGELOG.md`, PLAN section 8, the
validation report, `AGENTS.md`, and `develop/TODO.md` receive that final
census only at closeout.

## R2 recording-semantics amendment (2026-07-29)

This amendment is registered at clean HEAD `ac2cf03` before the first R2
production, unit-test, or structural-gate edit and before any R2 binary run.
It specializes the already-registered behavior gate 5 without changing its
decision: a correctly computed zero-valued synchronous header CRC is legal.
The prior-art audit searched PLAN sections 2, 4, and 8, this report, the
validation report, the survivor census, and the recording tests. It found no
later FALSIFIED, BLOCKED, REFUTED, or deferred record for
`crc-zero-rejected` or `csv-values-untested`.

R2 owns exactly these remaining decisions:

- `crc-zero-rejected`: APPLIED only after the frozen correct-zero fixture
  changes from `Invalid_parameters` to `Success`;
- `csv-values-untested`: APPLIED only after every field of all three
  nondegenerate CSV rows is parsed independently and reproduced bitwise; and
- `enqueuesnapshot-success-untested`: attempted REFUTATION through the
  existing real-CUDA P8-G3 path and the registered density-flag mutation.

`shuffle-oracle-gap` is already APPLIED at `cf64d074`; R2 does not add,
rewrite, or reclassify its independent 3-by-4 shuffle and inverse oracle.
Runtime CRC brute force, a forged snapshot count, a first-row-only CSV check,
NaN skipping, a round-trip-only shuffle check, and either minor-version
operator change remain killed. No format version, schema spelling, current or
state representation, successful nonzero-CRC behavior, or allocation
contract changes in R2.

The clean pre-edit anchors are:

```text
24E0B5FA53E0F4AC9D61F612DA549C0BDA2C47F52FD1E9CA9E7A7840D9730A03  436  src/core/RecordingFormat.cpp
10CEA3E821121549879670A4BE8D42D231D2BB16A8D753A88F9297D073F600DB  323  src/core/Recorder.cpp
D60BD14655C6E6321B26602305FB406E2344F54D73FD7317DE2D1542A080E1FB  391  src/core/AsyncRecorder.cpp
79BABE5E9E7E0EB8E322BB670105C5E48A42838DC10A4895D24EAE0635846D3D  619  src/core/CudaSpmBatch.cpp
E69ABFBFAA3454585DD435138CCB1FD485DD2300732A1442C149445DC70FC313  551  tests/unit/core_Recorder_test.cpp
3358FAE048F211E23E71FD220F519CC0FA27A0F5565C63FA7918A9D98E97FE79  352  tests/unit/core_CudaSpmBatch_test.cpp
6A4438E3E4BC7C772D5F015DA15AF0EEF8A6F85A2D9E5266AE41896F4EFAEC20 1791  tests/structural/p9c_architecture.cmake
```

### Correct-zero and corrupted-header oracle

The already-registered deterministic fixture remains the oracle; it is not
replaced by a newly searched or post-output value. It is a 3,584-byte zero
vector with this exact 64-byte header at offset zero:

```text
534c4944455245430100000004030201400000000000000070000000c700000087030000000000004000000000000000000e000000000000000e000000000000
```

The fields are little-endian `SLIDEREC`, major 1, minor 0, endian marker
`0x01020304`, header size 64, stored CRC 0, rows 112, lanes 199, stride 903,
snapshots 0, offset-table position 64, data offset 3,584, and file size 3,584.
The sole table entry at offset 64 is little-endian `uint64_t{3584}`. The test
builds those committed bytes directly and independently requires the
test-side CRC-32/ISO-HDLC result to equal zero; it performs no search.

The empty layout is fully legal: rows and lanes are positive, stride is at
least lanes, dimensions fit `int`, the one-entry table ends at byte 72,
3,584 is 64-byte aligned, `table[0] == data_offset`, and with zero snapshots
the first and final offset are both the file size. Successful
`BinaryRecording::open` must therefore publish `valid()==true`, `size()==0`,
`n_rows()==112`, `n_lanes()==199`, and `stride()==903`.

The recomputation witness changes only the stored CRC field from zero to one.
An independent header-CRC helper zeros that field and still obtains exactly
zero, while the stored value is exactly one. All magic, version, endian,
dimension, offset, table, and file-size facts remain legal, so
`Invalid_parameters` can come only from
`headerCrc(header) != header.header_crc32`. Opening this corrupt copy after
the valid fixture must fail atomically while the valid mapping and metadata
remain published. A rows-byte mutation is not used: changing only the stored
CRC gives the recomputation oracle no downstream alternate rejection.

The new zero-CRC case has exactly 14 assertions: two independent CRC facts,
two calls to the existing `writeBytes` helper at two assertions each, the
valid-open requirement, five valid metadata assertions, corrupt-open
rejection, and post-failure validity. Both files are constructed and written
before the valid-open requirement. Thus unchanged production reaches the
seventh new assertion and fails exactly there; fixed production reaches all
14.

The only allowed production delta is deletion of
`|| header.header_crc32 == 0` from `BinaryRecording::open`. The generic
writer assignment, reader recomputation, `minor > format_minor`, and every
later layout check remain exact. The R1 architecture assertion at the former
zero guard is expressly superseded: it is renamed for R2 and changes its
expected count from one to zero. The gate retains exactly one
`headerCrc(header) != header.header_crc32`, exactly one writer assignment,
and exactly one synchronous minor comparison, and forbids a replacement
`headerCrc(header) == 0` or `headerCrc(header) != 0` special case.

### Independent all-row CSV semantics

The existing `P6-G1 CSV and mmap recordings preserve snapshots` case retains
its three snapshots but changes accepted steps from the degenerate loop
indices to the exact sequence `{7, 11, 19}`. Simulation times remain
`{0, 1, 2}` seconds and currents remain `{8, -4}` A. This is test-fixture
hardening, not a production-output reblessing: R1's opaque platform-specific
CSV fingerprints live in `core_RecordingFormatCommon_test.cpp` and are
unchanged. The nondegenerate steps are required because serializing the CSV
loop index instead of `recorded.accepted_step` would otherwise remain green.

Before each `Recorder::record`, the test independently retains one typed
expected row from the live batch:

- the named accepted step and literal simulation time;
- both current densities computed once from current divided by electrode
  area;
- all 58 live state values, read as 29 `StateArena::row(row)[lane]` pairs,
  never through `Recorder::snapshot`, `snapshotRow`, or padded raw indexing;
  and
- both voltages from the shared test harness's production-observation call on
  that live state and density at the `0.0` observation time used by
  `SpmBatch::terminalVoltageAt`.

One grouped assertion freezes the fixture geometry at 29 rows, two lanes,
and stride eight. The arena therefore owns 232 padded values per snapshot,
but each CSV row must expose only the 58 live state values. The two density
lanes differ at every row, and rows after the first have lane-distinct states;
the first row alone is deliberately insufficient.

A file-local parser uses `std::string_view` slicing and
`std::from_chars`, not `std::stod`, streams, Recorder schema helpers, snapshot
accessors, or terminal-voltage accessors. It requires exactly one header and
three newline-terminated data rows with no trailing bytes, exactly 64
nonempty fields per data row, complete base-10 integer consumption for field
zero, and complete finite scientific-double consumption for fields 1 through
63. One assertion requires the complete parse. The three typed rows then
make exactly 192 comparisons:

```text
3 rows * (1 accepted step + 1 time + 2 currents
          + 29*2 live states + 2 terminal voltages) = 192
```

Every parsed real is compared by `std::bit_cast<uint64_t>`, not tolerance.
`Recorder.cpp` emits scientific format with `max_digits10` digits after the
decimal point, i.e. 18 significant digits for `double`, so exact recovery is
the intended contract. Reducing precision to `digits10` or at least
`max_digits10 - 2` must turn this oracle red; `max_digits10 - 1` is not a
registered mutation because 17 significant digits may still round-trip.

The CSV delta is exactly 227 assertions and zero cases: one grouped geometry
assertion, three 11-assertion independent voltage observations, one complete
parser assertion, and the 192 comparisons above. The existing case therefore
passes exactly 257 assertions / 1 case, and Recorder rises from 275/9 to
502/9 before the zero-CRC case. With that 14/1 case included, final Recorder
is exactly **516 assertions / 10 cases**.

### Existing CUDA seam and attempted refutation

The current call graph is sufficient only if the registered mutation proves
it dynamically. `CudaAsyncRecorder::enqueue` records the device current
density and state into pinned storage. Its drain worker obtains those spans
and calls the sole production caller of
`AsyncRecorder::enqueueSnapshot`; that function delegates to
`enqueueValues(..., true)`, whose distinguishing expression divides by one
rather than electrode area. Worker failure propagates through
`CudaAsyncRecorder::finish`.

The existing real-CUDA P8-G3 case supplies two nonzero densities for six
snapshots. Chen2020's electrode area is 0.1027, not one. The case requires
enqueue and finish success, decodes all six blocks, and compares both current
densities plus complete state exactly. R1's retained CUDA evidence executed
this real test at 433,671 assertions / 4 cases on the RTX 4000 Ada; it was
not the optional-off `CudaDisabled` binary.

At a clean R2 implementation boundary, change only the literal `true` passed
by `enqueueSnapshot` to `false`, rebuild the real CUDA target under the full
VS 18 x64 environment, and run exactly that test. The registered result is
six failed dynamic density checks:

```text
test cases:       4 |       3 passed | 1 failed
assertions: 433671 | 433665 passed | 6 failed
```

Accepted steps, states, counts, and statuses stay green. If and only if that
mutation is unexpectedly green after proving the correct CUDA tree and rebuilt
source were used, `enqueuesnapshot-success-untested` becomes APPLIED with the
already-reviewed CPU fallback: paired codec-none recorders compare
`enqueue(current_A)` against
`enqueueSnapshot(current_A/area)` bitwise, plus the finite-density/tiny-area
quotient-overflow discriminator. When the six CUDA assertions turn red, the
finding is REFUTED as a stale test inventory and no duplicate CPU round-trip
test is added.

### Frozen boundary, mutations, and acceptance bands

The oracle-only boundary may change only
`tests/unit/core_Recorder_test.cpp`,
`tests/structural/p9c_architecture.cmake`, and evidence records. All production
hashes above remain exact. Unchanged production must report:

```text
CSV case:       257 assertions / 1 case, all pass
Recorder:       509 reached assertions / 10 cases,
                508 pass / 1 fail; 9 cases pass / 1 fails
zero-CRC case:  7 reached assertions, the valid-open requirement alone fails
AsyncRecorder:  471 / 12, all pass
Recorder alloc: 33 / 3, all pass
Async alloc:    60 / 4, all pass
```

The aggregate structural gate passes every preceding suite and then fails at
the superseded zero-guard policy: expected zero occurrences, found one. This
is the frozen old-production-red behavior boundary; no source edit is allowed
before it is committed.

After the one-line implementation, independently execute these mutations from
a clean committed boundary and restore exact SHA-256 anchors between them:

1. restore `header.header_crc32 == 0`;
2. replace it with a computed-CRC zero special case;
3. delete or bypass `headerCrc(header) != header.header_crc32`;
4. serialize CSV loop index instead of accepted step, and separately force
   all serialized times to zero;
5. duplicate or reverse current lanes;
6. use lanes rather than stride in state selection, and separately reverse
   live state lanes;
7. duplicate or reverse voltage lanes;
8. reduce precision to `digits10` or `max_digits10 - 2`;
9. change `enqueueSnapshot`'s density flag from `true` to `false` in the real
   CUDA lane; and
10. weaken the test-side three-row bound or omit a field family. The exact
    516/10 count and structural test-source assertions must turn red rather
    than silently accepting a smaller oracle.

Changing the fixture's stored CRC away from zero, changing snapshots from
zero, changing `table[0]`, resealing the corrupt fixture, or changing either
minor-version operator is not an alternate implementation. Those mutations
either make the self-check/layout oracle red or violate a killed boundary.

Final focused gates in Debug, fast-math Release/IPO-off, and host-C++ CUDA are:

```text
Recorder               516 / 10
AsyncRecorder           471 / 12
RecorderAllocation       33 / 3
AsyncRecorderAllocation  60 / 4
aggregate structural      1 / 1
```

The retained CUDA tree additionally passes the restored real
`CudaSpmBatch` 433,671/4. Each configuration then passes unfiltered 58/58
and a no-op rebuild. `clang-format --dry-run --Werror`, `git diff --check`,
and an adversarial gate review are mandatory. Removing one boolean clause
adds or removes no Status-returning statement, so the registered WSL Clang 18
coverage census remains 359 lexical / 329 active = 323 measured + 6 exact
exceptions, 30 inactive, and zero uncovered/unmapped; it is rerun rather than
assumed.

If all registered gates pass, R2 closes two pending findings as APPLIED and
one as REFUTED, with no new supplemental ID. The original 71-ID census moves
from 40 APPLIED / 0 REFUTED / 8 deferrals / 23 pending to
**42 / 1 / 8 / 20**. The combined 84-ID census moves from
51 / 0 / 10 / 23 to **53 / 1 / 10 / 20**. `shuffle-oracle-gap` remains in its
existing APPLIED count and is not counted again.

## R1 post-native-acceptance Linux coverage hardening (registered 2026-07-29)

This amendment is registered at clean `a7a4236` after all fourteen R1
mutation families and all three retained native acceptance lanes passed, but
before any response to the first fresh WSL coverage run. It does not change
an R1 production file, format byte, CRC policy, schema, Status arm, or any of
the three previously frozen Windows tuples.

The first WSL Clang/LLVM 18.1.3 Debug coverage build was fresh and complete:
272/272 build actions, 55 measured test binaries plus one production anchor,
and an immediate no-op rebuild. Its unfiltered CTest was deliberately not
accepted at 56/58:

1. Recorder reached 268/275 assertions. File sizes, canonical header,
   endian, independently recomputed file/block/raw/payload CRCs, distinct
   raw-versus-payload CRCs, eager metadata/current/state, reader open, and
   dimensions were green. Exactly the seven Windows-Debug tuple assertions
   were red. The observed candidate Linux tuple is:

   ```text
   CSV    (4072,   580175175469508733,  8482670618437841473)
   SLREC  (3904, 16236475239795828522, 17177906431584795385)
   SLCMP  (3936, 17797356430432497358, 14598319551377944930)
   block CRCs (raw,payload):
     (4273982167, 2835542638), (2456069782, 2077085578)
   ```

2. The aggregate structural test passed 9C-2, 9C-3, and 9C-4, then CMake
   script mode rejected `stem IN_LIST P9C5_API` at
   `p9c5_public_surface.cmake:104`: CMP0057 was unset and the old-policy
   parser reported `Unknown arguments specified`. The separate compiler
   self-containment test passed.
3. All 55 expected nonempty `.profraw` groups exist, but the authoritative
   reporter stopped before export with
   `stale exception site:
   src/core/PackSolver.cpp:351:27:Invalid_parameters`. No coverage result is
   claimed from those profiles.
4. A read-only workflow audit found exact `57`-command checks in
   `.github/workflows/core-status-coverage.yml` and
   `.github/workflows/core-sanitizers.yml`; the fresh discovery artifact has
   58 commands. This is separately deferred to B1 and is not part of the
   local R1 repair.

Before adding a Linux tuple, prove it against the unchanged old R1
production at oracle commit `ad700ca`. Materialize only that commit's four
recording `.cpp` files and `detail/AsyncRecordingFormat.hpp` under the
already ignored `build-r1-coverage` tree, compile them with the exact
coverage-tree Clang 18 Debug command, replace those four members in a copy of
the current `libslide_core.a`, and relink a separate copy of the current
Recorder test objects. `git diff --name-only ad700ca..a7a4236` confirms that
no other fixture-producing production source changed during R1. The
old-production sentinel must reach exactly 268/275 and fail only the same
seven tuple assertions while printing every candidate value above. Any
different byte, CRC, assertion count, or failure site falsifies the proposed
Linux branch; do not bless the current output. In that fork,
`r1-linux-debug-byte-oracle-unregistered` remains OPEN and a second,
production-regression stable ID is appended before any source fix. Restore
current production to all thirteen old-source values first; only then may
the platform branch be added.

Only after that proof, split the compile-time oracle explicitly: Windows
Release/IPO-on, Windows Release/IPO-off, Linux Clang 18 Debug, and Windows
Debug. The three existing Windows tuples remain byte-for-byte unchanged, the
new branch selects exactly the old-production Linux tuple above, and an
unknown platform/toolchain/configuration must stop at compilation rather
than silently reuse another branch. WSL Recorder must return to 275/9
without changing its count or any native result. Restoring the prior generic
all-Debug fallback is a registered mutation and must reproduce 268/275 with
the same seven failures.

The `IN_LIST` spelling entered 9C-5 at `d0e82e6` without a local policy; the
failure therefore precedes the appended MQ.2 R1 block and is not caused by
its tokens. Because 9C-5 is also independently runnable, establish CMP0057
NEW inside `p9c5_public_surface.cmake`, not transitively in the aggregate.
Both direct `cmake -P tests/structural/p9c5_public_surface.cmake` with
`SLIDE_SOURCE_DIR` and the focused WSL aggregate CTest must pass; removing
the policy line must reproduce the CMP0057 red result in both entry paths.
The policy is already implied by the repository's CMake 3.31 minimum and
does not weaken a structural assertion.

Refresh the six exception identities without changing their classes,
reasons, statement hashes, or exception count:

```text
PackSolver.cpp:351 -> PackSolver.cpp:329  Invalid_parameters
PackSolver.cpp:588 -> PackSolver.cpp:563  Numerical_failure
PackSolver.cpp:652 -> PackSolver.cpp:625  Invalid_states
PackSolver.cpp:796 -> PackSolverIterative.cpp:153  Numerical_failure
PackSolver.cpp:887 -> PackSolverIterative.cpp:233  Invalid_states
SpectralModel.hpp:314 -> SpectralModel.hpp:315  Numerical_failure
```

The first three and SpectralModel retain their exact context hashes. The two
cross-file moves take current scanner context hashes
`9f4768e915b52c3effd0107d19fe39c049b73777a1ee98fab992b23dd1c9fdec`
and
`f081f954b131b5760033986c39739261085f98c12700a9d1a0c5c5f06096fe66`.
Restoring any stale identity must make the reporter refuse the manifest,
rather than silently reducing coverage.

The manifest was last validated at `93e75d1`. Its refresh is a one-to-one
identity relocation: the six classes/reasons and six statement hashes stay
exact, four context hashes stay exact, and only the two cross-file context
hashes change as listed. No exception is added, removed, merged, or
reclassified.

After the three local corrections, rebuild the complete coverage target,
require an immediate no-work dry run, rerun the scanner self-test and exact
359/329 census, prepare a new session, run unfiltered CTest at 58/58, and
run the reporter at exactly:

```text
359 lexical = 343 direct + 16 conditional
329 active = 323 measured + 6 structural exceptions
30 inactive; zero uncovered and zero unmapped
Numerical_failure 63; every other Status count unchanged
```

The compiled-test change invalidates the first profile session; no old
profile is reused for final evidence. Native Debug, fast-math Release/IPO
off, and retained host-C++ CUDA focused gates plus 58/58 are rerun after the
test/gate change. The three new local IDs become APPLIED only after all of
those gates are green. With the workflow ID deferred to B1, R1 closeout then
has the unchanged original census 40 APPLIED / 0 REFUTED / 8 named deferrals
/ 23 pending. The combined prerepair census is 84 IDs: 38 APPLIED /
0 REFUTED / 10 named deferrals / 36 pending. Closing the thirteen R1 IDs
makes it 51 / 0 / 10 / 23.

B1 must later change both hosted instrumented-workflow discovery pins from
57 to 58. Each instrumented configuration must freshly discover exactly 58
CTest commands; coverage must still classify 55 measured tests plus one
anchor. A controlled 58-to-57 mutation must make each workflow assertion
red. No hosted result is claimed unless an observed Actions run supplies it.

## R1 recording-common amendment (registered 2026-07-29)

This amendment is registered at clean HEAD `1b40b84` before the first R1
unit-test, CMake-list, structural-gate, production, or CHANGELOG edit and
before any R1 binary run. The prior-art audit found no FALSIFIED, REFUTED, or
deferred record for:

- `crc32-twice`;
- `header-crc-helper`;
- `csv-parquet-schema-twice`;
- `state-index-int-arith`;
- `snapshotview-twice`; or
- `dead-usings-copy-pasted`.

The frozen batch table classifies all six as byte-identical recording
refactors. R1 does not change a file-format version, CRC acceptance policy,
schema byte, payload byte, state/current selection, Status, or allocation
contract. In particular:

- the synchronous reader's `header_crc32 == 0` rejection remains; the
  `crc-zero-rejected` behavior decision belongs to R2;
- synchronous `minor > format_minor` and compressed
  `minor != format_minor` remain distinct and unchanged;
- `format_major` and `format_minor` remain independently owned by the two
  formats;
- the independent test-side CRC implementations remain test oracles;
- `ParameterSet.cpp` and `CyclerV2.cpp` keep their structurally pinned
  allocation-failure owners;
- `BinaryRecording::snapshot`, which indexes mapped interleaved records,
  remains separate from the two eager SoA readers; and
- R2 still owns independent per-field CSV semantic parsing/coverage and the
  disposition of `enqueuesnapshot-success-untested`. R1 treats its two data
  rows only as an opaque whole-file byte-identity oracle needed to execute
  `snapshotRow`, and does not call `enqueueSnapshot`.

The old-source anchors are:

```text
651C7F81FDE8C0274167DBD503001B5182E5372D93B0AE40071F18A4E1757756  293  src/core/Recorder.cpp
F5107C41F0841964DB13D2393D7B643C87B37A9557D1A66471BABEC4537DEBDE  451  src/core/RecordingFormat.cpp
E8A52D4319FD5198C343AFD3783FE8856F0035B39EB51E76F2D7B54E88CA0F8F  396  src/core/AsyncRecorder.cpp
10620233E09B29588859C01ECC295024EB69216192E2F96EEFE17FC1AEDF58BA  452  src/core/AsyncRecordingCodec.cpp
73A0EB170C4ADDB08816062D86FF1005036C99A89566C2B60C36068894F1C2FC  116  src/core/detail/AsyncRecordingFormat.hpp
E69ABFBFAA3454585DD435138CCB1FD485DD2300732A1442C149445DC70FC313  551  tests/unit/core_Recorder_test.cpp
17C04C42F70A77A111C262EFB6E6B668BF94D8B011E497DA3A4A5E3EDF145D04 1188  tests/unit/core_AsyncRecorder_test.cpp
407417B23C675ED6885D2C03ED4DFC50AAF856C7F7EA427EFB7C5E33BA63672B  245  tests/unit/CMakeLists.txt
B0C595C7BB0F9080C54420AD08BB491BC4CD604966949B7210D5A115F8A55143 1265  tests/structural/p9c_architecture.cmake
86532D4C48AC8367D6428766B44670996CB944444DD9D2BA32E2882954E5299E  326  CHANGELOG.md
```

The current lexical Status census is independently reproduced before edits:

```text
total 361 = 345 direct + 16 conditional
active 331, inactive 30
Invalid_parameters 236, Invalid_states 50,
NotImplementedYet 10, Numerical_failure 65
```

The active classification is from WSL Clang 18 against the retained coverage
tree. Consolidating exactly three recording-local
`allocationFailureStatus` definitions into one removes exactly two active,
direct `Numerical_failure` arms. No exception row moves.

### Frozen byte and eager-index oracle

Add `tests/unit/core_RecordingFormatCommon_test.cpp` as a second source of
the existing `unit_test_core_Recorder` target. This keeps the CTest census at
58, mirrors the new recording-common owner, and does not grow the already
1,188-line AsyncRecorder test. The file adds exactly two cases and 53
assertions, taking Recorder **222/7 -> 275/9**. AsyncRecorder remains
**471/12**, RecorderAllocation remains **33/3**, and
AsyncRecorderAllocation remains **60/4**.

The first case adds exactly 25 assertions:

1. configure a two-lane Recorder of capacity two;
2. for two snapshots, retain the factory's physically valid live state,
   poison every padding slot with distinct exact finite integer doubles, set
   distinct exact-quarter elapsed times, and use total currents equal to
   `{1,-2}` then `{3,-4}` times the exact electrode area;
3. directly require each public snapshot's accepted step, time, current
   density, and complete padded state against the independently retained
   inputs; and
4. write the complete two-row CSV, require its bytes to begin with an
   independently constructed canonical header, and freeze its exact byte
   count plus both recurrences in one grouped assertion. This executes both
   `snapshotRow` calls while remaining an opaque byte-identity check; R2
   retains independent parsing and per-field semantic ownership; and
5. write the synchronous binary file, require its exact byte count and two
   independent whole-file 64-bit recurrences, require the literal endian
   marker, and recompute the stored header CRC with the test-side CRC owner.

The second case adds exactly 28 assertions. It writes the same two
deterministic snapshots through `AsyncRecorder::enqueue` with codec `none`
and blocking backpressure, freezes the exact complete `.slcmp` byte count and
the same two whole-file recurrences, independently verifies the file header
and both block-header CRCs plus the endian marker, opens the file, and
requires both eager snapshots' metadata/current/state directly against the
retained inputs. It does not compare against `Recorder::snapshot`, so a
paired mutation of both eager readers cannot remain self-consistent and
green.

The byte fingerprint records the exact byte count plus:

- conventional byte-wise 64-bit FNV-1a, seeded with
  `14695981039346656037`, applying XOR with each unsigned byte and then
  multiplying by `1099511628211`; and
- a second full-width recurrence seeded with `0x6a09e667f3bcc909`,
  applying
  `rotl(mixed ^ (byte + 0x9e3779b97f4a7c15), 17)` and then multiplying by
  `0xbf58476d1ce4e5b9` for each unsigned byte. Unsigned wrap is intentional.

An initial run with deliberately impossible digest sentinels is explicitly
**exploration**, not evidence. It may only print the old-production byte
counts and recurrences. Those nine values (CSV, `.slrec`, and `.slcmp`) are
then written into this amendment and the test, committed test-only, and rerun
against unchanged production.

A pre-freeze adversarial audit identified one configuration-sensitivity risk:
the complete CSV includes production terminal voltages formatted at
`max_digits10`, and the retained Release/CUDA-host configurations use
fast-math. Therefore the sentinel capture sequence is run against unchanged
production in Debug, fast-math Release/IPO-off, and the host-C++ CUDA tree
before any constant is frozen. Each exploratory run must compile, execute
exactly 275 assertions / 9 cases, fail exactly 7 assertions solely at the
one grouped CSV and six binary fingerprint sentinels, and leave both block
CRC/data-CRC checks green. It also prints each block's independently checked
stored raw/payload CRC pair, which must be unequal. If any of the nine
fingerprint values differ between configurations, the test freezes one exact
tuple for each active recorded-fixture configuration branch; it must not
accept an unordered set of known values. The existing
`slide_configure_recorded_scalar_fixture` helper is applied to the Recorder
target so `SLIDE_TEST_RELEASE` and `SLIDE_TEST_IPO` name those branches
without inferring configuration from fast-math. A non-asserting CRC `CAPTURE`
addition and these repeated sentinel runs remain exploration, not evidence.

The three captures met that band exactly: each reached 275 assertions / 9
cases, with 268 assertions green and exactly the seven impossible digest
sentinels red; every other assertion passed. All five production hashes
remained equal to the anchors above before and after every capture. The
frozen values are `(bytes, FNV-1a, mixed)`:

```text
Debug:
  CSV    (4072, 9687440020754251757, 4395482851448939133)
  SLREC  (3904, 5616759404010358714, 2046236715092195141)
  SLCMP  (3936,  674564965144750450, 6807551341330539846)
  block CRCs (raw,payload):
    (4136208346, 118461792), (2589125531, 3586173060)

Release, IPO off:
  CSV    (4072,  9836986554487851297, 10195804735363157979)
  SLREC  (3904,  7644568223896275882,  9507010167800238048)
  SLCMP  (3936, 12852073895150042970,  2881486469091071908)
  block CRCs (raw,payload):
    (1129188246, 4052193948), (798320599, 591966072)

Release, IPO on, host-C++ CUDA tree:
  CSV    (4072,   113020753288656565, 1220644351662519581)
  SLREC  (3904,  4980822913020495434, 5810136181434859582)
  SLCMP  (3936, 15459230203036584454, 6543818133775003568)
  block CRCs (raw,payload):
    (1528093825, 3638640079), (936358080, 170786859)
```

The differing raw CRCs prove that fast-math/IPO changed factory-state bytes;
this is a toolchain fixture distinction, not a format difference. In every
configuration each raw CRC differs from its payload CRC. The two block
assertions freeze the applicable exact pairs while continuing to recompute
raw and payload CRCs independently.

No production or structural-gate edit may begin until the frozen
old-production run passes exactly 275/9 and all five production source hashes
still match the anchors above.

### One format-fact owner

The already frozen filename `src/core/detail/RecordingFormatCommon.hpp`
wins over the survivor prose's earlier `RecordingCommon.hpp` spelling. Its
leading MC-3 contract names ownership of `endian_marker`, `crc32`,
`headerCrc`, and `allocationFailureStatus`, PLAN §3.7, and its cold
file-format/allocation-translation role. It is `@surface internal`,
self-contained through `types/Status.hpp`, and contains exactly:

```cpp
inline constexpr std::uint32_t endian_marker = 0x01020304U;

inline std::uint32_t crc32(std::span<const std::byte> bytes)
{
  std::uint32_t crc = 0xffffffffU;
  for (const auto byte : bytes) {
    crc ^= std::to_integer<std::uint8_t>(byte);
    for (int bit = 0; bit < 8; ++bit)
      crc = (crc >> 1U) ^ (0xedb88320U & (0U - (crc & 1U)));
  }
  return ~crc;
}

template <class Header>
[[nodiscard]] inline std::uint32_t headerCrc(Header header)
{
  static_assert(std::is_trivially_copyable_v<Header>);
  header.header_crc32 = 0;
  return crc32(std::as_bytes(std::span{ &header, 1 }));
}

inline slide::Status allocationFailureStatus() noexcept
{
  return slide::Status::Numerical_failure;
}
```

No `[[nodiscard]]` or `noexcept` is added to the moved `crc32`; its signature
and body remain the existing async spelling. `headerCrc` is the one generic
by-value zero-then-hash owner. The three packed headers are already
trivially copyable and exactly 64 bytes.

`RecordingFormat.cpp`, `Recorder.cpp`, and
`detail/AsyncRecordingFormat.hpp` include the common header and contain zero
definitions of the three moved facts. `RecordingFormat.cpp` directly calls
`headerCrc` in its writer and reader. The reader retains the explicit
`header.header_crc32 == 0` guard and compares
`headerCrc(header) != header.header_crc32`; no `std::exchange` spelling
remains. The async header's two type-specific CRC wrappers are removed, and
the async writer/reader directly call the one `headerCrc` owner for file and
block headers. Across the five production files plus the common header, the
exact compact token censuses are:

```text
crc32(                    6 = one definition, one headerCrc call, four payload calls
endian_marker             8 = one definition, three using declarations, four header consumers
allocationFailureStatus( 13 = one definition, twelve catch consumers
headerCrc(                7 = one definition, six header consumers
fileHeaderCrc( / blockHeaderCrc(  0 / 0
```

The compact per-file matrix is also frozen so the same global counts cannot
pass after moving a call to the wrong format or raw/payload role. Columns are
`crc32(` / `endian_marker` / `allocationFailureStatus(` / `headerCrc(`:

```text
detail/RecordingFormatCommon.hpp  2 / 1 / 1 / 1
RecordingFormat.cpp               0 / 3 / 2 / 2
Recorder.cpp                      0 / 0 / 6 / 0
detail/AsyncRecordingFormat.hpp   0 / 0 / 0 / 0
AsyncRecorder.cpp                 2 / 2 / 2 / 2
AsyncRecordingCodec.cpp           2 / 2 / 2 / 2
detail/SnapshotIndexing.hpp       0 / 0 / 0 / 0
```

`RecordingFormat.cpp`, `Recorder.cpp`, and
`detail/AsyncRecordingFormat.hpp` directly include the common header.
`Recorder.cpp` and `AsyncRecordingCodec.cpp` directly include the snapshot
header. Every other matrix entry is pinned to zero as well as every positive
entry being pinned to its exact consumer. In particular, the async writer
contains exactly one ordered pair
`.raw_crc32 = crc32(raw)` / `.payload_crc32 = crc32(payload)`. The reader
checks `block.payload_crc32 != crc32(payload_view)` before unshuffle and
`block.raw_crc32 != crc32(raw)` after unshuffle. The structural gate pins
that payload-read -> payload-CRC -> unshuffle -> raw-CRC order.

The common header contains zero `format_major` and `format_minor` tokens.
`RecordingFormat.cpp` and `AsyncRecordingFormat.hpp` each retain exactly one
independent `format_major = 1` and `format_minor = 0` definition and their
existing distinct reader comparisons.

### One eager-snapshot owner

`src/core/detail/SnapshotIndexing.hpp` is separate from the format-fact
header: combining unrelated file-format and reader-indexing concepts would
violate MC-1. Its leading contract names `snapshotView`, PLAN §3.7, and cold
allocation-free eager-reader indexing, and carries `@surface internal`. It
includes `Recorder.hpp` in the permitted internal-to-api direction and owns
exactly:

```cpp
[[nodiscard]] inline SnapshotView snapshotView(
  std::size_t index,
  std::size_t lanes,
  std::size_t state_values,
  std::span<const std::uint64_t> steps,
  std::span<const real_t> times,
  std::span<const real_t> currents,
  std::span<const real_t> states) noexcept
{
  return { .accepted_step = steps[index],
           .time = times[index],
           .current_density = currents.subspan(index * lanes, lanes),
           .state = states.subspan(
             index * state_values, state_values) };
}
```

Only `Recorder.cpp` and `AsyncRecordingCodec.cpp` include it. Each eager
member retains its distinct caller-side assertion and directly returns
`snapshotView` with the same seven logical arguments. Thus the exact global
`snapshotView(` census is three: one owner and two consumers.
`BinaryRecording::snapshot` retains its mapped-byte implementation and does
not include or call this helper.

### One schema order and widened row index

The anonymous namespace in `Recorder.cpp` owns one
`recorderColumnNames(int rows, int lanes)` function. It returns the canonical
order:

```text
accepted_step, time_s,
current_density_lane{lane}_A_m2 for every lane,
state_r{row}_lane{lane} in row-major/live-lane order,
terminal_voltage_lane{lane}_V for every lane
```

Each of the eight quoted schema fragments occurs exactly once, inside that
owner. The CSV sink joins the returned names with commas. The Arrow-only
Parquet branch obtains the same vector once and consumes it monotonically
with one `std::size_t column{}` cursor across its five append sites, then
asserts `column == names.size()`. The exact
`recorderColumnNames(` census is three: owner, CSV, Parquet. No Arrow-enabled
runtime is available in the retained configurations, so R1 claims exact CSV
runtime bytes and structurally proves Parquet name/order consumption; it
does not claim a Parquet runtime.

#### Post-freeze Debug iterator-proxy correction

The first production-WIP Debug run exposed one new supplemental finding,
`recorder-schema-default-vector-debug-oom-terminates`, before any R1
implementation commit. It is a real public-boundary exception-safety defect
in the initially frozen schema-owner spelling, not a reason to weaken the
existing allocation oracle:

- `RecorderAllocation`'s first two cases passed independently at 5/1 and
  15/1. The third, already-frozen public-I/O case reached the fail-next
  allocation after `writeBinary` and then opened a hidden Microsoft Visual
  C++ Runtime abort dialog. Command timeouts that left orphaned dialog
  processes are exploration and are not evidence.
- A disposable, allocation-free diagnostic mutation established that the
  failing allocation was the 16-byte MSVC Debug iterator proxy created by
  `std::vector<std::string>`'s default constructor. The diagnostic test file
  was restored byte-for-byte to SHA-256
  `D8FA21C68BA93468F1194C15E5F1AD40167898BF0A9354BF0F724F226639AC01`
  and is not part of the implementation.
- [confirmed] In the installed VS 18 STL, `_DEBUG` selects
  `_ITERATOR_DEBUG_LEVEL == 2` (`yvals.h`), `vector()` is conditionally
  `noexcept` on the nothrow-default-constructible allocator and calls
  `_Alloc_proxy` (`vector`), and that proxy allocates one object when iterator
  debugging is enabled (`xmemory`). `std::allocator`'s default constructor is
  `noexcept`. Therefore a `bad_alloc` from this proxy cannot reach
  `Recorder::writeCsv`'s catch. The pre-R1 code did not have this site: its
  first product vector used a non-`noexcept` count constructor.

The corrected schema owner starts with the explicit canonical prefix:

```cpp
std::vector<std::string> names{ "accepted_step", "time_s" };
```

The initializer-list constructor is not `noexcept`; a storage or Debug-proxy
allocation failure therefore propagates to the unchanged public catch.
`reserve()` after a default construction is forbidden because the unsafe
proxy allocation has already happened. The architecture gate must require
this exact prefix constructor and reject a bare default constructor. This
changes no successful CSV/Parquet byte, name, order, Status, or allocation-free
hot path. The decisive band remains RecorderAllocation **33/3**, including
`matching_failure_triggered == true`, no escaped exception, and
`Status::Numerical_failure` for the CSV failure. Recorder remains **275/9**
with all five old-production fingerprints exact.

The bare-default-constructor mutation is structural-red and is not executed
as a behavioral mutation on Windows because its expected behavior is an
interactive CRT abort, not a bounded Catch failure. The already-frozen
allocation test is the behavioral acceptance oracle for the safe owner.
Only after all R1 gates pass is
`recorder-schema-default-vector-debug-oom-terminates` APPLIED. It is outside
the original 71-ID census. The post-mutation adversarial review registered
three further R1 IDs below, growing the combined registry from 76 to 80 IDs.

The R1 self-containment requirement is also made persistent rather than
satisfied by one-off commands. 9C-5 keeps its 28-header api/support
classification unchanged and compiles
`detail/RecordingFormatCommon.hpp` and `detail/SnapshotIndexing.hpp` in a
separate two-header internal list. Removing either new header's load-bearing
project dependency (`types/Status.hpp` or `Recorder.hpp`) must make that
compiler gate red; the headers are not reclassified as public. Exact direct
standard-include ownership remains a lexical architecture assertion, because
a standalone compiler cannot distinguish a direct include from one supplied
transitively by another header.

#### Post-mutation adversarial hardening registration

An independent read-only audit of committed implementation `5a462b6` was
performed after the first 13 mutation families and before the changes in this
subsection. It found no format, CRC, schema, eager-index, Status, include-cycle,
or initializer-list defect. It did find three bounded quality/gate issues,
assigned stable IDs in the append-only registry above:

1. `Recorder.cpp` still says the binary format's header, CRC, and mmap all
   live in `RecordingFormat.cpp`, and `RecordingFormat.cpp` still says its
   co-location avoids a CRC seam. The comments must instead distinguish its
   packed-header/layout/mapping owner from `RecordingFormatCommon.hpp`'s
   shared CRC/endian/allocation facts.
2. The 9C-5 claim is narrowed to the load-bearing project include, as stated
   above; the already-executed mutation removed `SnapshotIndexing.hpp`'s
   `Recorder.hpp` include and failed standalone compilation with
   `unknown type name 'SnapshotView'`.
3. The architecture gate must pin `recorderColumnNames` and `snapshotRow`
   inside the one anonymous namespace immediately preceding
   `Recorder::configure`. Exact bodies alone do not prevent those helpers
   acquiring external linkage.

These are no-op comment/test corrections. Before final acceptance, remove the
anonymous-namespace open/close around the two helpers as controlled mutation
14; the new linkage assertion must turn red. Restore the exact final gate and
source hashes, then rerun the full focused and three-lane gates. No runtime
oracle, Status count, output byte, or census baseline is reblessed.

The same anonymous namespace owns one:

```cpp
[[nodiscard]] std::span<const real_t> snapshotRow(
  const SnapshotView &snapshot,
  int row,
  int rows,
  int stride,
  int lanes) noexcept;
```

It asserts `0 <= row && row < rows && 0 < lanes && lanes <= stride`, then
returns a live-lane subspan beginning at
`static_cast<std::size_t>(row) * static_cast<std::size_t>(stride)`.
CSV iterates that span; Parquet indexes its live lane after the helper.
The exact `snapshotRow(` census is three. The structural gate pins both
ordered calls exactly as
`snapshotRow(recorded, row, rows_, stride_, lanes_)` and
`snapshotRow(snapshot(i), row, rows_, stride_, lanes_)`, and Recorder.cpp
contains zero instances of a cast applied after `row * stride_` or
`row * stride_ + lane`.

### Dead-using boundary and future structural gate

`AsyncRecorder.cpp` removes only its dead
`detail/CheckedArithmetic.hpp` include and `checkedAdd`,
`checkedMultiply`, and `compressionBound` using-declarations.
`AsyncRecordingCodec.cpp` retains all four because its layout code consumes
them. The common-header migration replaces the two async CRC-wrapper usings
with the live `headerCrc` using; every other live format name remains.

Append R1 after C1 in `tests/structural/p9c_architecture.cmake`. It loads
compacted versions of the five production files and both future headers,
while preserving raw leading comments for the MC-3 checks. Before detailed
checks, it treats a missing future header as an empty string and must fail
old production first at:

```text
9C-2 MQ.2 R1 RecordingFormatCommon owner: expected 1 occurrences ..., found 0
```

The final gate pins:

- the exact common-header bodies and 6/8/13/7/0/0 combined censuses above;
- zero moved definitions in all old files;
- independent version ownership and both unchanged reader comparisons,
  including the zero-CRC guard;
- the exact `snapshotView` body, two includes, two caller assertions, two
  ordered seven-argument calls, and the untouched mapped reader;
- one schema owner, once-only literals, 3 owner/caller tokens, exact CSV join,
  monotonically consumed Parquet cursor, and final size assertion;
- one widened `snapshotRow` owner, the two exact ordered five-argument
  consumers, and zero late-cast old index spellings; and
- absence of the three dead AsyncRecorder usings/include while the
  corresponding AsyncRecordingCodec dependencies remain.

At the oracle-only boundary, only
`tests/unit/core_RecordingFormatCommon_test.cpp`,
`tests/unit/CMakeLists.txt`, and
`tests/structural/p9c_architecture.cmake` may differ. Old production must pass
Recorder exactly **275/9**, AsyncRecorder **471/12**, RecorderAllocation
**33/3**, and AsyncRecorderAllocation **60/4**, while the aggregate
structural test passes every preceding suite and then fails only at the
registered common-owner 0/1 boundary. Production hashes must still equal the
anchors above. No Release, CUDA, full-suite, or production claim is made at
that boundary.

### Registered mutations and final acceptance

At a clean committed implementation boundary, independently run at least
these mutation families:

1. alter the shared CRC polynomial/body or reintroduce a private CRC owner;
2. make `headerCrc` zero a wrong byte/field or bypass it at one synchronous
   or compressed header consumer;
3. change the shared endian marker or merge either format-version owner;
4. make the shared allocation mapper return a different Status;
5. change one schema fragment/order or make either sink bypass the common
   names/cursor;
6. change `snapshotRow` to multiply in `int`, use `lanes` as its stride,
   swap the same-typed `stride_`/`lanes_` arguments at either consumer, or
   reintroduce either old late-cast index;
7. damage the shared `snapshotView` offset/extent body;
8. swap one same-typed eager-reader argument or make either caller bypass
   the helper;
9. reintroduce one dead AsyncRecorder using/include or remove the live codec
   counterpart; and
10. alias the async writer's raw CRC to the shuffled payload, make either
    reader comparison use the wrong field/bytes, or pair a writer/reader
    raw-payload swap. For both deterministic codec-none blocks, the test
    independently hashes the retained unshuffled current-plus-padded-state
    bytes and the exact on-disk shuffled payload; each one grouped assertion
    requires a valid block-header CRC, both stored data CRCs, and
    `raw_crc32 != payload_crc32`, preserving the registered +28 case count.
11. remove or weaken the synchronous zero-CRC guard or alter either minor
    comparison. This mutation is structural-red; R1 must not execute or
    reclassify R2's behavior change.
12. replace the schema owner's initializer-list prefix with a bare default
    `std::vector<std::string>` construction. The architecture gate must turn
    red; do not execute the known interactive-abort path as mutation evidence.
13. remove one load-bearing project dependency include from either new
    internal header. The separate 9C-5 internal standalone-compilation loop
    must turn red.
14. remove the anonymous namespace around `recorderColumnNames` and
    `snapshotRow`. The architecture linkage assertion must turn red even
    though exact bodies and runtime outputs are unchanged.

Every touched source/header/test/gate file is SHA-256 anchored at the clean
implementation commit. Each mutation is inverse-patched individually; exact
hashes, `git diff --exit-code HEAD -- <touched-files>`, and an empty porcelain
status are required before the next mutation and before final acceptance.
No byte fingerprint, schema spelling, CRC, Status, output value, count, or
hash may be reblessed.

The final focused set in Debug, fast-math Release/IPO-off, and the retained
host-C++ CUDA tree is Recorder 275/9, AsyncRecorder 471/12,
RecorderAllocation 33/3, AsyncRecorderAllocation 60/4, and aggregate
structural 1/1. Every lane then runs the unfiltered 58/58 suite, including
the real CUDA test only in the CUDA tree, followed by a no-op rebuild.
`clang-format --dry-run --Werror`, `git diff --check`, header self-containment,
and an adversarial gate review precede disposition.

Because the Status owner moves, the WSL Clang 18 exact Status-coverage lane
is also mandatory after implementation. Its registered census is:

```text
359 lexical = 343 direct + 16 conditional
329 active = 323 measured + 6 unchanged structural exceptions
30 inactive; zero uncovered and zero unmapped
Numerical_failure 63; every other Status count unchanged
```

Those gates are green at `7b8b1f8`, so the six original R1 IDs, the four
earlier R1 supplemental findings, and the three local coverage-hardening IDs
are APPLIED. This moves the original census from
34 APPLIED / 0 REFUTED / 8 named deferrals /
29 pending to 40 / 0 / 8 / 23, and the combined 76-ID census from
38 / 0 / 9 / 29 through the prerepair 84-ID census
38 / 0 / 10 / 36 to the closeout census 51 / 0 / 10 / 23. The fourth new
coverage-audit ID remains explicitly deferred to B1. `CHANGELOG.md`,
PLAN section 8, the
validation report, `AGENTS.md`, and `develop/TODO.md` receive that final
census only at closeout.
