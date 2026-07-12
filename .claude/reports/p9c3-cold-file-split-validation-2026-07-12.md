# M0.7 / 9C-3 cold-file split validation (2026-07-12)

## Result

**PASS.** The 1,496-line fixture-baseline `ParameterSet.cpp` is now
separated into value absorption/SPM compilation, BPX mapping/file I/O,
expression parsing/evaluation, strict JSON, and three narrow detail seams.
The 1,063-line fixture-baseline `Experiment.cpp` is now a parser TU beside one
cohesive Cycler runner TU and a shared semantic seam. Production behavior is
digit-identical in all three recorded-bit configurations, CTest remains exactly
56 tests, and all sources are owned by the one centralized core source list.

The functional split is commit `10b7b36`; the structural ownership gate is
commit `57a442b`; final coverage-boundary closure is `93e75d1`. The frozen
fixtures were captured at `e360218`, before any
production definition moved, and the immediate extraction baseline was
`2bd6ec5`, after the allocation-atomicity fixes that the split preserves.

## Architecture delivered

- `ParameterSet.cpp` owns canonical names, value-semantic set/update/query,
  built-in Chen2020 absorption, and compilation to `SpmFactoryInput`.
- `BpxParameterReader.cpp` alone owns BPX JSON paths, BPX-to-ParameterSet
  mapping, bounded files, public BPX entry points, and atomic publication.
- `BpxExpression.cpp` owns the stack-local AST, scanner, parser, and evaluator.
  Its only two external detail seams evaluate caller-provided samples or
  produce one adaptively sampled `OCVCurve`; no AST lifetime crosses the TU.
- `StrictJson.cpp` owns the resource-bounded strict grammar and UTF-8 scanner.
  Its one value-semantic DOM seam parses into a local candidate and publishes
  by statically no-throw move only after complete success.
- `detail/ParameterCurve.hpp` is the one implementation of curve validation and
  D-16 adaptive priority refinement shared by built-in and BPX absorption.
- `Experiment.cpp` owns `Experiment::parse`; `CyclerV2.cpp` owns all eleven
  runner methods and their transaction scratch. `detail/ExperimentSemantics.hpp`
  is the single inline definition of normalization and semantic admission.
- `cmake/SlideCoreTarget.cmake` lists all six moved production TUs exactly once,
  so ordinary, fuzz, coverage, core-only, and nested-consumer builds cannot
  acquire different implementations.
- `p9c3_cold_file_split.cmake` proves ownership, the internal dependency DAG,
  public-header isolation, exact fixture wiring, source-list uniqueness,
  transactional JSON publication, and physical line ceilings. It is included
  by the existing 9C architecture CTest, preserving the 56-test suite count.

## Mechanical-equivalence proof

Whitespace-normalized comparisons against `2bd6ec5` were exact:

| Extracted region | Compared tokens/chars | Result |
|---|---:|---|
| ParameterSet core methods | 12,217 / 12,217 | exact |
| curve validation | 327 / 327 | exact |
| built-in OCP functions | 340 / 340 | exact |
| adaptive sampler | 2,011 / 2,011 | exact |
| BPX expression AST/parser/evaluator | 5,745 / 5,745 | exact |
| strict JSON parser | 7,629 / 7,629 | exact before transactional wrapper |
| BPX mapping/file methods | 10,860 / 10,860 | exact after seam normalization |
| Experiment parser definition | 6,795 / 6,795 | exact after `detail::` qualification |
| eleven Cycler definitions | 29,529 / 29,529 | exact after `detail::` qualification |

Standalone Clang C++20 syntax checks passed for all four detail headers. The
aggregate structural gate passes directly and through CTest.

## Frozen behavior

### Parameter absorption and compilation

| Trace | Values | Debug/ThinLTO FNV / mix | Release FNV / mix | Release/ThinLTO FNV / mix |
|---|---:|---|---|---|
| Chen2020 values | 11,537 | `a5d211a63dd77522` / `8ccfefe72e03094e` | `e5c3dc16e932ad03` / `b4e4cfcce5167740` | same as Release |
| `SpmFactoryInput` | 11,635 | `134b5b4650709d1a` / `1e63d406b3d98297` | `705555174f1b2e23` / `d9357b2ee6af55e4` | same as Release |
| functional BPX values | 4,135 | `04bdc30e3e6ca305` / `e75964cae8725037` | `d3f393daa406ad80` / `f7d10dea65f3c8f5` | same as Release |

Chen2020 retains 53 entries, 159 metadata strings, and 3,491 metadata bytes;
its metadata hashes are `8ea717ef0459a7b7` / `745379e447872b1f`.
The functional BPX case retains 35 entries, 105 strings, and 1,906 bytes;
its metadata hashes are `f4a1c008bf74be9a` / `fa654b7fe7b9ab89`.

### Experiment parser and runner

The eight-segment parser fixture retains 96 numeric/integer values and the
configuration-independent hashes `1e404a3be7365df5` / `94a28984eed0159b`.
Sources and drive-cycle names are also compared in order.

The independently hand-built runner fixture retains 503 values:

| Configuration | Runner FNV / mix |
|---|---|
| Debug/ThinLTO | `c779e41caac8339c` / `32dec619a41e0324` |
| Release | `e9616b6ce3131bb6` / `b72cecb529e7c3c0` |
| Release/ThinLTO | `9d97787b3976978b` / `4c7bcd0565164aa2` |

The runner fixture bypasses the parser and covers power, rest, current, a
sign-changing drive table at interior and breakpoint samples, sample periods,
and a dyadic event root. Parser/runner common-mode errors therefore cannot
explain an exact pass.

## Adversarial findings and fixes

1. **Missing runner ownership.** The first independent integration audit found
   `CyclerV2.cpp` absent from the centralized source list. A relink would have
   left all runner methods unresolved. It was added exactly once before any
   integrated model test was accepted.
2. **Non-atomic strict-JSON reuse (P9-B54).** The extracted seam initially
   parsed into caller storage. A successful reuse retained stale object members,
   and malformed input partially mutated the prior DOM. The SHORT red witness
   failed at assertion 4/4 (`2 == 1`). Candidate parsing plus no-throw move
   publication passes 17/17 and is structurally pinned. BPX had always supplied
   a fresh local, so no public BPX fixture changed.
3. **Incomplete seam proof.** The first structural draft named the two BPX
   seams but did not prove that the reader consumed them. Exact reader call
   counts, definition/declaration counts, private-parser exclusions, and the
   JSON transaction shape are now mandatory.
4. **Allocation-count hypothesis falsified.** The pre-split ParserAllocation
   total was 1,138. A transactional `StrictJsonValue` candidate creates one
   extra cold map-node allocation, so the exhaustive BPX sweep gains one
   occurrence times five preservation assertions: 902 -> 907 and full
   ParserAllocation 1,138 -> 1,143. The seven 9C-3 allocation cases remain
   exactly 108 assertions. No accepted-step allocation path changed.
5. **Harness working-directory error rejected.** Directly launching the
   Experiment binary from the repository root made twelve cases abort before
   model construction because legacy data paths are configured relative to the
   CTest working directory. That run is recorded as non-evidence; CTest and the
   correct build/bin working directory pass the same binary at 352/352.
6. **Initial exact-coverage hypothesis falsified.** The 369-site M0.6 census
   was stale after the pre-split allocation hardening: the first fresh session
   found 380 lexical / 350 active sites and failed with 337 measured, seven
   exceptions, and six uncovered. Five were duplicate literal returns for
   `length_error` catches whose `bad_alloc` twins were covered. Equivalent
   handlers now delegate to one mutation-tested mapper per TU, while the gate
   structurally proves every catch and both callback rethrows. A valid NCH8
   batch replacing a configured NCH5 batch directly covers the remaining stale
   scratch guard. The honest final census is 369 lexical / 339 active, expected
   as 333 measured plus six unchanged structural exceptions.
7. **ThinLTO proof option repaired.** An older IPO build tree regenerated with
   recorded bits disabled. Its portable assertions passed but totals were lower
   by exactly the guarded six ParameterSet and four Experiment hash checks. The
   tree was explicitly reconfigured with recorded bits enabled; all exact checks
   then passed. No hash was inferred or replaced.

Earlier controlled mutations, all reverted before extraction, independently
turned red for uncaught drive registration, missing-drive late discovery,
post-step OOM rollback, reentrant registration, and state-only event rollback.
The split preserves those transaction bodies exactly.

The final coverage-boundary changes are independently mutation-red: returning
`Success` from the stale-scratch guard passes only 17/18 assertions; mapping the
Parameter allocation boundary to `Invalid_parameters` fails five of 31 focused
assertions; the equivalent Cycler mutation fails eight of 64. All three mutants
were rebuilt, observed red at their intended Status check, reverted explicitly,
and followed by clean exact-fixture runs in all three configurations.

## Test evidence

The new JSON regression raises ParameterSet from 2,358 assertions in seven
cases to 2,375 in eight. Experiment is 370 in fourteen cases after the
18-assertion valid NCH5-to-NCH8 stale-scratch regression. The
intentional cold JSON candidate raises ParserAllocation to 1,143 in thirteen
cases; its 9C-3 transaction slice remains 108 in seven.

| Lane | Result | Notes |
|---|---:|---|
| Windows Clang Debug/ThinLTO | 56/56 | parallel CTest 17.20 s; ModeC 16.57 s |
| Windows Clang Release/fast-math | 56/56 | parallel CTest 6.82 s; ModeC 6.16 s |
| Windows Clang CUDA/fast-math/ThinLTO | 56/56 | parallel CTest 6.31 s; CUDA batch 1.83 s; ModeC 5.75 s |
| WSL Clang 18 ASan+UBSan | 56/56 | 129.66 s; ModeC 128.91 s; zero ASan, LSan, or UBSan finding |
| WSL Clang/LLVM 18 exact Status coverage | 56/56 | 256.90 s; 369 lexical; 339 active = 333 measured + 6 exceptions; 30 inactive; zero uncovered/unmapped |
| Windows core-only optional-off | 1/1 | dependency rebuild, 0.09 s, immediate no-op |
| Windows nested external consumer optional-off | 1/1 | dependency rebuild, 0.09 s, immediate no-op |

Native Debug, Release, and CUDA immediate rebuilds each report
`ninja: no work to do`. TSan was not repeated: M0.7 changes no synchronization,
worker lifetime, or shared concurrent ownership, so no new TSan claim is made.
Hosted CI is not claimed before a push.

The final WSL sanitizer rebuild performed 51 changed-state actions in 1,214.5 s
on the concurrent mounted-filesystem host; this is not a performance claim. It
proved ASan+UBSan flags on all six split TUs and the changed Experiment test,
runtime references in both core and legacy archives, 54 build-local executable
tests plus two structural scripts, and then an immediate no-work rebuild.

## Archive and ownership evidence

LLVM 21.1.8 found 18 members in both Debug and Release `slide_core` archives,
with the six split objects exactly once. Defined public ownership is exactly ten
ParameterSet core methods plus two BPX APIs, one `Experiment::parse`, and eleven
Cycler methods. Each archive has exactly one strong definition of
`evaluateBpxExpressionSamples`, `sampleBpxExpressionCurve`, and
`parseStrictJson`, and zero global AST, `JsonParser`, scanner, parse-node, or
evaluator helpers. Retained private helpers are local `t` symbols.

Some emitted Release inline helpers are printed as `T` by `llvm-nm` on COFF;
`llvm-readobj --sections` proves their sections carry `IMAGE_SCN_LNK_COMDAT`,
so they are ordinary inline COMDATs rather than duplicate strong ownership.
The rebuilt core-only and nested-consumer archives independently show the same
three seams, zero private-parser globals, 18 members, and one copy of every
split object.

Final archive SHA-256 identities at `93e75d1` are
`10b72209a7fc0a18abafe5cd05163006a2a815ef7dd6f11592b1f653b5cee481`
(Debug), `03c94ccce8a3d85c902022fd25c2c7cacb5490882597127fdfb61a040d434e41`
(Release), `4007298bc92431f21b84a8360101a9949326cb2fabcc0440f84690cbe2b86d8d`
(core-only), and
`0198e681e5655b2516aefbbeb13654fb961486a484149176f0436af097c632b2`
(nested consumer).

## Maintainability accounting

The immediate `2bd6ec5` pre-extraction ParameterSet cluster was 1,507 physical
/ 1,436 nonblank lines. Its seven cohesive implementation/detail files are
1,749 / 1,641 (+242 / +205). The Experiment cluster was 1,211 / 1,166 and is
now 1,270 / 1,215 (+59 / +49). The growth is explicit contract comments,
declaration seams, namespace/include boundaries, and transactional JSON reuse;
no physics or runner algorithm was duplicated.

Every Parameter-side file is below 500 physical lines and every detail header
below 200. `CyclerV2.cpp` is 771 physical lines, a reviewed MC-1 exception:
it is one cohesive runner transaction, and an additional TU boundary would
change same-TU non-IPO inlining and could perturb the frozen Release trace.
The gate uses a documented ceiling of 800 and M0.10 will revisit it only after
a stable runner-transaction interface exists.

Artifacts:

- `.claude/reports/p9c3-pre-refactor-cold-fixtures-2026-07-12.md`
- `.claude/reports/p9c3-status-coverage-refresh-2026-07-12.md`
- `.claude/reports/p9c3-status-coverage-refresh-2026-07-12.json`
- `.claude/reports/p9c3-cold-file-split-validation-2026-07-12.md`
