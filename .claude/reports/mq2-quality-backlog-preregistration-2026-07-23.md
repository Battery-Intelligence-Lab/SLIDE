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
repository addendum makes that earlier registration durable. Trace D is
registered here before its first build/run after adversarial review found that
A–C use only dyadic resistances and therefore cannot freeze the
`x/R` versus `x*(1/R)` fast-math seam.

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
4. **Trace D — non-dyadic callback and publication bits.** A test-only affine
   adapter records every current span passed to `linearizeThevenin` before
   copying its constant `E/R` outputs. First, fresh sparse, ladder, and
   relaxation solvers run `parallel(4, affine)` with
   `E={4,4.1,3.9,4.2}`, `R={0.1,0.2,0.15,0.3}`, applied current `2`,
   tolerance `1e-12`, and maximum four iterations; each must take exactly two.
   Second, a fresh damped sparse solver uses `parallel(2, affine)`,
   `E={4,0}`, `R={0.3,0.7}`, zero applied current, the same tolerance, and
   maximum three iterations; it must take exactly three. For each of the four
   paths, two independent `RecordedBits` recurrences hash, in order, every
   callback-current frame, final cell currents, final node voltages, and one
   scalar frame containing terminal voltage, residual norm, constraint drift,
   constraint bound, and relaxation gain. Expected hashes are captured
   separately for Debug, fast-math Release/IPO-off, and
   Release/ThinLTO/CUDA-host builds through deliberate zero placeholders
   before any pack source edit. The all-mode records contain 19 doubles each;
   the damped record contains 15. No hash may be reblessed during P1 or P2.

Trace D is the decisive digit-identity gate for non-dyadic division at all
four current-reconstruction sites: sparse undamped, sparse damped, ladder, and
relaxation. A final solution-only hash is insufficient because a later
iteration could wash out a changed reconstructed current; callback frame two
is therefore part of the record.
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
