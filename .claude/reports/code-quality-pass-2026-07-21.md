# Implementation & code-quality pass — 2026-07-21

Baseline commit `353c905`. Method: a nine-way fan-out review over **disjoint** subsystems of `src/core`
(no two reviewers read the same file), each finding then handed to an independent adversarial verifier
whose default was to refute. **117 findings raised, 68 survived, 49 refuted or contract-blocked.**
Raw survivors: `.claude/reports/code-quality-pass-2026-07-21-survivors.json`.

## 0. Baseline, recorded before any edit

| Quantity | Value |
|---|---|
| CTest (Debug) | **58/58** |
| Catch2 binaries | **53**, 0 failing |
| Total assertions | **36,568** |
| `src/core` warnings under the project's own policy | **1** (`-Wshadow`) |

Baseline artifact: per-binary assertion/test-case counts captured for all 53 binaries and used as the
digit-identity oracle for every change below.

## 1. A measurement error found while taking the baseline

Running the 53 binaries directly gave **25 failures**, while CTest reported 58/58. Cause: `add_test`
sets `WORKING_DIRECTORY "bin"` (`cmake/Coverage.cmake:160`), and `src/settings/slide_paths.hpp` resolves
`data/` as `"../.." / "data"` **against the current working directory**. CTest happens to launch from
`<build>/bin`, which is exactly two levels below the repo, so the tests passed *by coincidence of layout*.
A developer running `./unit_test_core_Sei.exe` got `File ../..\data\Kokam_OCV_C.csv could not be opened`.

The header has always had a `#ifdef SLIDE_ROOT_DIR` escape hatch — **nothing in the build system ever
defined that macro**, so the `#ifdef` branch was dead code and every build silently took the fallback.

**Registered before the fix:** defining the macro should take the repo-root run from 25 failures to 0
while leaving CTest at 58/58 and every assertion count unchanged; anything else reverts.
**Result: 25 → 0, CTest 58/58, all 53 counts digit-identical.** Confirmed.

## 2. Changes landed, each with its evidence

| # | Change | Evidence |
|---|---|---|
| 1 | `SLIDE_ROOT_DIR` defined at directory scope (excluded from `SLIDE_CORE_ONLY`, so wheel/WASM/consumer builds stay relocatable) | §1 registered prediction confirmed |
| 2 | `PathVar::results`/`data` `static` → `inline const` | nothing assigns to them; no namespace-scope object is constructed from them, so there is no static-init-order exposure (checked, not assumed) |
| 3 | `project_warnings` linked `PRIVATE` to `slide_core` | it previously had **no library consumer at all** — `slide_core` and all 53 `tests/unit` binaries compiled with zero `-W` flags. 22/22 core TUs now rebuild **0 warnings, 0 errors** |
| 4 | `CyclerV2.cpp` inner `derivative` → `current_derivative` | the only `-Wshadow` in core: a scalar d(current)/dt shadowing the arena derivative-row span declared ~120 lines above |
| 5 | Boost include marked `SYSTEM` (both branches) | Eigen and range-v3 already were; Boost accounted for 61 of 187 unique warning sites in the test tree |
| 6 | **PC-10 fix:** `Sei`/`Lam`/`LithiumPlating` call `spm_scalar::arrheniusFactor` | nine sites already called it, including `CudaSpmRuntime.cu`; all 53 binaries digit-identical after |
| 7 | Byte-shuffle hand-computed + property oracles | see §3 |
| 8 | `parseValue` scientific-notation and truncated-exponent tests | see §3 |
| 9 | **Ladder detection O(cells²) → O(cells)**: `cell → first branch` index replaces a `find_if` over every branch, per cell | complexity argument verified at `PackTopology.cpp:256`; all 53 binaries digit-identical after |

**Digit-identity result:** 51 of 53 binaries byte-for-byte identical to baseline; the two that changed are
exactly the two I added tests to (`AsyncRecorder` 463→471 assertions / 10→12 cases, `NetlistCsv`
800→821 / 7→9). CTest 58/58 throughout.

### Why the PC-10 fix is digit-identical for `Dual`, provably

The concern was that `arrheniusFactor(static_cast<Real>(T_ref), T, static_cast<Real>(Rg))` might change
the derivative component versus the hand-rolled `(Real{1}/T_ref − Real{1}/T)/Rg`, since `Dual/Dual` and
`Dual/double` need not agree bit-for-bit. **`Dual` declares only `operator/(Dual, Dual)`** (`Dual.hpp:58`)
— there is no `Dual/double` overload — so the hand-rolled form was *already* performing the implicit
`double → Dual` conversion that the cast now writes explicitly. Same operations, same order.

`SurfaceCrack.hpp` is deliberately **excluded**: it writes the factor as `k_act/Rg · (1/T_ref − 1/T)`,
a different association, so routing it through the helper would change the last bits. Recorded as
remaining PC-10 debt rather than silently reassociated.

## 3. Two tests that could not fail, and the proof they now can

- **`byteShuffle`** was covered only by a round-trip test — accurately titled "is a bitwise involution" —
  which passes whenever shuffle and unshuffle share *the same wrong permutation*. Added a hand-computed
  3×4 byte-plane expectation and a property oracle (after shuffling smoothly-varying doubles the plane
  holding the shared exponent byte must be a constant run, the low plane must not be, and the output must
  be a permutation of the input multiset). **Demonstrated:** reversing the plane order in *both*
  directions — still a perfect mutual inverse — leaves the old test green (10 passed) and turns both new
  tests red (4 assertions).
- **`parseValue`'s exponent branch** was never exercised, though pandas writes connection resistances as
  `1e-05`, so a real liionpack `netlist.to_csv()` is precisely the untested input. The branch turned out
  to be **correct** — this is a test gap, not a bug, and is reported as such. Disabling the branch turns
  the new test red.

## 4. A measurement I got wrong twice, and the correction

I recorded one CTest run at **57/58** and went looking for the cause. Two things then went wrong, and
both are worth recording because they are ordinary mistakes that produce confident-looking numbers.

**First, confounding.** The 12-iteration hunt ran while I was building and sweeping 53 binaries
concurrently. Concurrent `ninja` + `ctest` on Windows is a known failure mode already recorded in the
M0.9 handoff, so nothing from that loop is evidence.

**Second — and worse — the detector itself was wrong.** Both loops classified a run as failing with

```sh
if echo "$out" | grep -q "tests failed"; then ...
```

but CTest's **success** line is `100% tests passed, 0 tests failed out of 58`, which *contains*
`tests failed`. Every run matched. The reported "6/12" and the clean re-run's "8/8 had a failure"
summary are both artifacts of that predicate, not observations. The per-iteration output the same
script recorded shows the truth plainly:

```
--- iter 1 ---   100% tests passed, 0 tests failed out of 58
...
--- iter 8 ---   100% tests passed, 0 tests failed out of 58
```

**What is actually true:** 8/8 clean runs are **58/58 green**, as were the earlier consecutive runs.
The single 57/58 observation was real (it printed `98% tests passed, 1 tests failed out of 58`) but has
never reproduced across roughly a dozen subsequent runs, and the failing test's name was never captured.
It is recorded as **one unexplained, unreproduced failure** — not as a characterised flake, and not as
a clean bill of health either.

## 5. Found and NOT fixed — carried as debt

59 verified findings remain unapplied. The ones worth an owner, roughly by value:

- **The branch→{adjacency, nodal sparsity, BFS connectivity} derivation is written twice, line for line**
  (`PackTopology.cpp:180-207` and `374-433`).
- **`PackStepper::substeps` does not subdivide `dt`** — it repeats `dt` N times, so `step()` advances
  `substeps·dt`. **Partially addressed**: the declaration now documents it (verified: the electrical
  solve and thermal assembly do sit outside the loop, so the current really is frozen across substeps).
  The name still contradicts the code; renaming it is a breaking API change and needs its own decision.
- **`relaxation_target_` is reused for three unrelated quantities** in one function, including as the
  Kahan compensation array for a differently-named buffer.
- **`crc32`, `allocationFailureStatus`, `endian_marker` are each defined twice** inside the recording
  subsystem — the exact smell `detail/CheckedArithmetic.hpp` was created to kill.
- **The cell constitutive relation `I = (E − ΔV)/R` with its finiteness gauntlet appears four times**
  across three pack-solver kernels.
- **The PyBaMM electrode name lexicon is written three times** (Chen2020 table, `toSpmInput` suffixes,
  BPX reader).
- **`ForwardSensitivity` rebuilds the seed-independent cold model once per parameter** and recomputes
  the identical primal trajectory P times, discarding it; `solveOne` also evaluates the full observable
  twice per sample.
- **A fast-math-unsafe floating-point validity gate** in `SpmObservables.hpp:463-468`, in a TU compiled
  `-ffast-math` — the construct `Numeric.hpp` itself documents as unreliable there.
- Untested arms: CUDA `step`/`enqueue` rejection paths, `BatchExecutor`'s inline path, pack-solver
  source-stepping rollback, CSV sink *values* (the sink is checked only by header prefix and newline
  count), `enqueueSnapshot`'s success branch.

## 6. What is and is not claimed

**Claimed (Debug, this machine, verified by artifact):** CTest 58/58; 53/53 binaries pass; 51 digit-identical
to baseline and the two deltas accounted for; 22/22 core TUs compile warning-clean under the project policy;
every mutation named above observed red and the source restored (`git diff --stat` clean after each).

**Not claimed:** Release and CUDA lanes were **not** re-run in this pass — the digit-identity evidence is
Debug-only, and the project's usual Debug-AND-Release standard is therefore not met here. TSan not rerun.
The llvm-cov lane not rerun. No performance claim is made about any change; the O(cells²) finding above is
a complexity argument from the code, not a measurement.
