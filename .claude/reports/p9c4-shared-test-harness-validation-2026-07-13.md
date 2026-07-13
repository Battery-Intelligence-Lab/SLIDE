# M0.8 / 9C-4 shared test harness — validation (2026-07-13)

Baseline commit: `93e75d1` (pre-harness), inventory in
`.claude/reports/p9c4-pre-refactor-harness-baseline-2026-07-12.md`.
Work under validation: `6cf44c4..31f5675` (harness foundation + four migration
commits). The harness landed in a prior session; **no gate had been run on it**.
This report is that gate.

## 1. What the harness owns

`tests/support/CoreSpmTestHarness.hpp` (346 lines) owns exactly the four
preregistered mechanics: successful `SpmBatch` construction from explicit
input/options/lanes; terminal-voltage observation behind the `CurrentA` /
`CurrentDensityApm2` unit tags with caller-owned scratch; constant-current
`ExponentialModal` traces on a caller-supplied strictly increasing grid
(initial sample plus every partial final interval); and a transparent
maximum-absolute / RMS voltage error. It holds no default options, chemistry,
C-rate conversion, reference values, tolerances, file I/O, stored batch/stepper
state, or mutable static. Reference values and tolerances remain in the mirrored
`core_*` tests. Structural enforcement: `tests/structural/p9c4_shared_test_harness.cmake`,
run inside `structural_test_core_9C2AgeingKernel`.

Migrated: `core_{Recorder,AsyncRecorder,RecorderAllocation,AsyncRecorderAllocation,
Experiment,ExponentialModal,P1G4_restart,ForwardSensitivity,P7G3_PyBaMM}_test.cpp`.
Direct `buildSpmBatch` calls remain only in the preregistered allowlist families
(construction/allocation-ordinal/topology/backend subjects).

## 2. Suites

| Lane | Result |
|---|---|
| Native Debug (Clang 18, Ninja) | 57/57 |
| Native Release fast-math (Clang 18) | 57/57 |
| CUDA Release (`build-cuda4`) | 57/57 |
| WSL Clang 18 ASan+UBSan+LSan | 57/57, no finding (139.14 s) |

Baseline CTest was 56; the suite is 57 because the harness ships its own
source-named self-test (`unit_test_core_CoreSpmTestHarness`). No test was
removed or merged.

## 3. Assertion floors — mechanically compared, not eyeballed

All 43 pre-existing binaries were run and compared against the frozen floors.
**No binary is below its floor; no test-case count dropped.** 39 binaries are
digit-identical; four increased because harness helpers add their own admission
`REQUIRE`s inside the caller's loop:

| Binary | floor → measured |
|---|---|
| Experiment | 370/14 → 436/14 |
| ForwardSensitivity | 3919/4 → 16172/4 (12512 before the §7 grid assertion) |
| P7G3_PyBaMM | 6334/2 → 7274/2 |
| Recorder | 168/7 → 222/7 |
| CudaSpmBatch (CUDA lane) | 433671/4 → 433671/4 |
| CoreSpmTestHarness | new → 184/8 |

Migrated binaries report identical counts in Debug, Release, and the CUDA
configuration.

### Floors re-frozen (M0.10 and later must not drop below these)

The old floors now carry slack, which falsified one preregistered mutation
(§5, F1). The binding floors are hereby re-registered at the measured Debug
counts above; every other binary keeps its baseline value.

## 4. Mutation battery — the harness gate is real

Registered before running: every mutation must turn a gate red. Results
(`RED` = gate caught it):

| # | Mutation | Caught by | Verdict |
|---|---|---|---|
| 1 | hard-code lane count (ignore `n_lanes`) | self-test | RED (5 failures) |
| 2 | hard-code `nch`/options | self-test | RED (6) |
| 3 | flip the positive-discharge sign | self-test | RED (3) |
| 4 | omit the A → A/m² area conversion | self-test | RED (3) |
| 5 | drop the initial trace sample | self-test | RED (1) |
| 6 | drop the final interval | self-test | RED (6) |
| 7 | uniform `dt` (skip nonuniform partial intervals) | self-test | RED (6) |
| 8 | compare the actual trace with itself | self-test | RED (16) |
| 9 | ignore the final metric element | self-test | RED (9) |
| 10 | RMS divided by `N-1` | self-test | RED (2) |
| 11 | naive metric (no IEEE preflights at all) | self-test | RED (15) |
| 12 | step at the wrong absolute time | structural gate | RED |
| 13 | default lane count on the build seam | structural gate | RED |
| 14 | harness build inside an allocation-fault window | allocation binary | RED |
| 15 | perturb one committed PyBaMM sample by 50 mV | P7G3 | RED (1) |

Mutation 10 confirms the registered literal metric: errors `{1, 2, -1}` give
maximum 2 V and RMS `sqrt(2)` V.

Mutation 15 detail: a **1 mV** perturbation is correctly GREEN — it lies inside
the registered 1C band (11.577 mV max / 2.028 mV RMS). Only a perturbation
exceeding the registered band may turn the gate red, and 50 mV does.

### Masked mutations (reported, not excused)

Removing any **single** IEEE guard from `computeVoltageError` stays green:
`canEvaluateAbsoluteDifference`, the `is_finite(error)` check, and the
`canSquare`/`canAddSquare` preflights are three independent layers, each
sufficient on its own to reject NaN/Inf. Individually they are therefore
redundant. The property itself is covered: mutation 11, which replaces the
whole metric with the naive implementation, fails 15 assertions including the
NaN and Inf rejection cases. This is defence in depth, not dead code, but the
per-guard mutations cannot discriminate and are recorded as masked.

## 5. FALSIFIED preregistrations (killed, with numbers)

**F1 — "removing one migrated `REQUIRE` violates the per-binary count floor."**
False. Harness helpers *raise* the counts well above the baseline floors
(Recorder 222 vs floor 168; ForwardSensitivity 12512 vs 3919). Removing one
assertion leaves the binary far above its floor, so the gate cannot see it.
Remedy applied: floors re-frozen at the measured counts (§3). Residual
limitation, recorded rather than hidden: a per-binary aggregate floor still
cannot tell a lost *oracle* assertion from a gained *admission* assertion.

**F2 — "moving helper work inside an allocation-fault window is forbidden by the
structural gate."** False. `p9c4_shared_test_harness.cmake` contains no
allocation-window rule; with a harness build moved inside the fault window of
`core_RecorderAllocation_test.cpp` the structural gate stays GREEN. The
protection exists, but comes from a different mechanism: the allocation binary's
own ordinals turn RED. The preregistered mechanism is wrong; the invariant holds.

**Note on mutation 12.** The wrong-absolute-time mutation cannot be detected by
any state or voltage oracle, because the exact modal SPM update is autonomous:
no kernel reads `StepCtx.time` (it is only constructed and forwarded —
`SpmPipeline.hpp:347,385,511,566`), and `elapsed_time` accumulates `dt`, not `t`.
The structural gate pins the `start_time_s` token, which is what turns red.

**F3 — the MC-1 fixture-reduction sub-goal is UNMET.** `develop/TODO.md` assigned
M0.8 to "factor repeated fixtures out of the >700-line PackSolver, AsyncRecorder,
and Experiment tests". Measured effect of the migration:

| File | before → after |
|---|---|
| core_AsyncRecorder_test.cpp | 1075 → 1118 (+43) |
| core_Experiment_test.cpp | 1091 → 1168 (+77) |
| core_Recorder_test.cpp | 520 → 551 (+31) |
| core_PackSolver_test.cpp | 1128 → 1128 (untouched; direct-call allowlist) |
| core_P7G3_PyBaMM_test.cpp | 114 → 104 (−10) |

Migrated tests net **+144 lines**; adding the harness, its self-test, and the
structural gate, M0.8 adds ~1,190 test lines overall. The shared *mechanics* are
real (one idiom, MC-4), but caller-owned scratch spans cost more lines at each
call site than the local lambdas they replaced, so the oversized files grew.
What PLAN §6 M0.8 actually demands — one helper, no assertion-count drop — is
met; the line debt is not, and is carried to M0.10 (its registered owner) rather
than quietly dropped.

## 6. Independent checks

- `requireSpmBatch` returns `core::SpmBatch` **by value**. This is sound:
  `SpmBatch`'s move constructor `std::exchange`s every pointer and function
  pointer (`SpmFactory.cpp`), `StateArena` owns its storage through a
  `unique_ptr` so the buffer address is stable across a move, no implementation
  holds a back-pointer to the owning batch, and `ExponentialModal::configure`
  caches only scalars and sizes (`ExponentialModal.hpp:30-37`) — the harness
  configures the stepper only after the batch has reached its final address.
- The self-test verifies the harness against an **independent manual pipeline**
  built in the test itself (`core_CoreSpmTestHarness_test.cpp:194-231`), not
  against the harness.

## 7. Adversarial review and what it changed

An independent reviewer was tasked to break the migration. Its four highest-value
attacks were refuted with evidence (by-value seam, allocation windows, oracle
weakening, circularity — §6 and the checks it cites). Three findings survived and
were verified against the sources before acting:

**R1 (fixed) — the structural blacklist was narrower than its own comment.** It
forbade the single token `core::ExponentialModalstepper`, so a differently named
stepper member, a `EulerLegacy` member, or a Catch matcher tolerance would have
passed a gate that claims the harness owns "no steppers, references, or
tolerances". Added forbidden tokens `EulerLegacy`, `stepper_`, `SpmBatchbatch_`,
`Catch::Matchers`, `WithinAbs`, `WithinRel`, and verified each turns the gate RED
by injection. Accepted residual limitation, recorded not hidden: a bare numeric
threshold literal inside the harness would still pass — the clause is enforced
against a token list, not against the property.

**R2 (fixed) — the integration grid became an implicit consequence of its
inputs.** `core_ForwardSensitivity_test.cpp` previously handed `dt` straight to
`stepper.step`; after migration it accumulates `time += dt` into a grid and the
harness re-derives `dt_s = t[s] − t[s−1]`. In IEEE-754, `fl(fl(t+dt) − t) ≠ dt`
in general, so a future non-dyadic `sample_step` would silently integrate a
different problem, and the `2e-12` primal gate is far too loose to notice an
ulp-level shift. A bitwise `REQUIRE(sample_time[s] - sample_time[s-1] == dt)` now
registers the invariant; it passes, so the present grids are exact
(ForwardSensitivity 12512 → 16172 assertions).

**R3 (documented) — the configured-batch ownership invariant was unstated.**
`Recorder` (`Recorder.hpp:105`) and `CyclerV2` (`Experiment.hpp:158`) store a raw
`SpmBatch *` at configure time, so a configured batch must never be moved or
reseated; `requireSpmBatch` makes batch prvalues easy to produce. No live bug (a
prvalue cannot bind to `SpmBatch &`, and every call site names a local first), so
the seam's doc comment now states the invariant.

Not acted on: `31f5675` also reflowed an unrelated initializer in
`core_Experiment_test.cpp` (formatting churn in a migration commit — noted, not
reverted), and `canAddSquare(max, min_subnormal)` is deliberately conservative
(it rejects a sum that would in fact be finite); that conservatism is pinned by
the self-test and is contractual.

## 8. Not claimed

TSan was not rerun for M0.8 and no new TSan claim is made. Hosted CI workflows
are committed but not claimed as run. Wall-clock timings on this machine remain
unreliable and are not used as evidence.
