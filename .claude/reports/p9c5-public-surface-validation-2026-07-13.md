# M0.9 / 9C-5 public-surface audit — validation (2026-07-13)

Baseline commit `df2da52`. Preregistration (bands, mutations, lexicon) written before any edit:
`.claude/reports/p9c5-public-surface-preregistration-2026-07-13.md`. Design note and the
rejected alternatives: `.claude/designs/m0-9-public-surface.md`. Every band below is the one
registered there, judged verbatim.

## 1. What changed

`SpmFactory.hpp:8` included `SpmPipeline.hpp`, so the ageing kernel stack compiled into every
user translation unit. Split the value types out of the kernels — four `*Params.hpp` headers,
one `SpmBatchLayout.hpp`, one `AgeingModelMask.hpp` — and the factory header now names types
only; `SpmFactory.cpp` includes the pipeline. `SpmPipelineLayout` → `SpmBatchLayout` (it is what
`SpmBatch::layout()` returns; "Pipeline" is an internal word). Every core header carries an
`@surface api|support|internal` tag; `tests/structural/p9c5_public_surface.cmake` enforces the
classification, the api set, each api header's declaration count, what bindings may include, and
the naming lexicon.

## 2. Registered bands — verdicts

**B1 — digit-identical behaviour. PASS.**

| Lane | Result |
|---|---|
| Debug, 53 Catch2 binaries | assertion + test-case counts **identical** to the `df2da52` baseline |
| Release fast-math, 53 binaries | **identical** to a `df2da52` Release baseline built with the same options |
| CTest Debug | 57/57 |
| CTest Release (fast-math, `-ffast-math` via `cmake/StandardProjectSettings.cmake:59`) | 57/57 |
| CTest CUDA (`build-cuda4`, sm_89) | 57/57 |
| `unit_test_core_CudaSpmBatch` | 433671 assertions / 4 cases — the frozen M0.8 number, unchanged |

**B2 — MC-5 property, compiler-verified. PASS at the registered target, with one claim corrected.**
The probe TU (`Experiment` + `ParameterSet` + `ExponentialModal` + `ForwardSensitivity` +
`CudaSpmBatch`) went from **26 core headers to 21**, and all eight registered kernel headers left:
`SpmPipeline`, `SpmObservables`, `Sei`, `Lam`, `SurfaceCrack`, `LithiumPlating`, `ThermalLumped`,
`SpmStress` (registered target 8 → 0; header count 20 ± 3).

**Corrected overclaim.** An earlier draft of this report, the CHANGELOG, and the commit message
said the TU contains "zero kernels". That is **false**. `SpmScalarKernels.hpp` — 364 lines that
call themselves "the sole owner of SPM scalar expression trees" — is still in every user TU, via
`SpmFactory.hpp` → `LamParams.hpp` → `CompiledCurve.hpp:11`, because `IndexedPiecewiseLinear`'s
inline evaluator calls `spm_scalar::linearInterpolate` (`CompiledCurve.hpp:151`). This was a known,
recorded design decision (PC-10 outranks MC-5: the alternative is duplicating interpolation), but
the *claim* did not match it. The rule R3 enforces is the honest one — **no header classified
`internal` is reachable from the api set** — and the gate's own comment now says exactly that
instead of "must not compile a kernel".

R3 makes this an invariant rather than a measurement: every core header is classified (R1), and
no api/support header may include an internal one, so no internal header is reachable from the
api set by induction. `clang -MM` is the empirical confirmation of the proof, not its substitute.

**B3 — the gate is real. PASS: 13/13 mutations RED** (10 registered + 3 added during the run,
noted as such). Each mutation was applied to the committed tree, the gate run, the file restored
from an in-memory snapshot.

| # | Mutation | Rule | Verdict |
|---|---|---|---|
| 1 | delete a header's `@surface` tag | R1 | RED |
| 2 | `@surface public` (invalid tier) | R1 | RED |
| 3 | `#include "SpmPipeline.hpp"` in `Recorder.hpp` (api) | R3 | RED |
| 4 | `#include "SpmObservables.hpp"` in `Experiment.hpp` (api) | R3 | RED |
| 5 | `SeiParams.hpp` (support) includes `AgeingKernel.hpp` (internal) | R3 | RED |
| 6 | `python/bindings.cpp` includes `core/SpmPipeline.hpp` | R5 | RED |
| 7 | add a public struct to `Simulation.hpp` | R4 | RED |
| 8 | rename a public type (`NetlistCsvDiagnostic` → `…Renamed`) | R4 | RED **after a fix** |
| 9 | re-add `nLanes()` beside `n_lanes()` | R6 | RED |
| 10 | a `validate…` returning `bool` | R6 | RED |
| 11 | a `getFoo()` free function | R4 | RED |
| 12 | retag an api header `internal` | R2 | RED |
| 13 | a `getFoo()` **method** (indented; invisible to R4) | R6 | RED **after a fix** |

**B4 — no performance or compile-time claim is made.** None is.

## 3. Gate holes found by the battery, and fixed

**H1 — R4 could not see a renamed public type.** The type check used `string(FIND)`, and
`NetlistCsvDiagnostic` is a substring of `NetlistCsvDiagnosticRenamed`, so mutation 8 passed a
green gate while the anchor count stayed put. Fixed: R4 now matches the *declaration*
(`(class|struct|enumclass|using)<Type>[{:=;<]`), and mutation 8 turns red.

**H2 — the `get`-prefix rule never fired.** It was checked against whitespace-stripped text,
where a leading `get` is glued to its return type (`inlineintgetWorkerCount`) and the
`[^A-Za-z_]get[A-Z]` anchor cannot match. Mutation 11 was being caught by R4's count, not by the
naming rule — and an indented `getFoo()` *method* would have been caught by nothing. Fixed: the
rule now runs against comment-stripped text that keeps its whitespace, and mutation 13 (which
only that rule can catch) turns red.

Both holes were found because the mutation was registered before the run. Neither would have
been visible from a green gate.

## 3b. Adversarial review — seven more evasions, and one real defect

An independent reviewer was tasked to break this work. It did, and everything below is fixed and
verified red by injection. Mutations 14–20 are its attacks, added to the battery (**20/20 RED**).

**A1 (fixed, HIGH) — R4 pinned almost nothing.** The anchor count matched column-0 declarations
only, so every class member, enumerator, and struct field was unpinned; and **no support header
was pinned at all**, even though the `*Params` structs *are* the public parameter surface — a new
field in `SeiParams` silently changes the layout of `SpmFactoryInput`, an api type. The reviewer
added a public method to `SpmBatch`, an enumerator to `SpmComposition`, and a field to `SeiParams`,
and the gate stayed green all three times. R4 now pins **both** anchor levels (namespace-scope and
member-level, at two-space indent) across **api and support** headers — mutations 14, 15, 16 RED.

**A2 (fixed, HIGH) — three working routes back into the kernels.** An api header could include
`"../core/SpmPipeline.hpp"` (the `../` escape hatch added for `types/Status.hpp` whitelisted a path
that climbs back into `src/core`), or `<core/SpmPipeline.hpp>` (R3 matched quoted includes only);
and a binding could include `"core/detail/StrictJson.hpp"` (R5's regex could not match a nested
path). All three passed a green gate. Fixed and now RED — mutations 17, 18, 19.

**A3 (fixed, MEDIUM, a real code defect) — two of my new headers were not self-contained.**
`SeiParams.hpp` and `SurfaceCrackParams.hpp` use `slide::Status` but included neither
`../types/Status.hpp` nor anything that provides it:

```
SeiParams.hpp:65:29: error: no type named 'Status' in namespace 'slide'
SurfaceCrackParams.hpp:49:29: error: no type named 'Status' in namespace 'slide'
```

They compiled only because `SpmFactory.hpp` includes `LamParams.hpp` (→ `CompiledCurve.hpp` →
`../types/Status.hpp`) alphabetically first. I hit this exact error while writing
`LithiumPlatingParams.hpp`, fixed it there, and did not check its three siblings. Any include
reorder, IWYU pass, or a user writing `#include <core/SeiParams.hpp>` first would have broken the
build. Fixed in all three; and a new **compiler-verified** gate,
`structural_test_core_9C5HeaderSelfContained`, now compiles each of the 28 api/support headers
alone (`-fsyntax-only`). A token gate could never have caught this; the compiler catches it every
time. CTest is now 58.

**A4 (fixed, LOW) — R6 forbade one syntactic shape, not the dead name.** `intnLanes(` missed
`auto nLanes()`, `a->nLanes()`, and `std::size_t nLanes()`. The bare tokens `nLanes(` / `nRows(`
are now forbidden outright (no legitimate use exists), `SpmPipelineLayout` is forbidden too, and
`n_rows()` is required where it was previously only implied — mutation 20 RED.

**A5 (fixed, LOW-MED) — `p9c2` was weakened by include-presence.** After the params moved, the gate
checked only that each mechanism header *includes* its params header. A mechanism could have
stopped consuming the shared mask and stayed green. It now re-pins mask **use**:
`for_each_enabled_ageing_model_lane<N>(p.model_mask,` in each of SEI, surface-crack, and LAM. The
earlier claim that the 9C-2 invariant was "unchanged" was slightly overclaimed; this makes it true.

**Refuted / accepted, not acted on:** `SpmFactory.hpp:195`'s `friend struct SpmBatchFactoryAccess;`
lets any user define that struct and reach `SpmBatch`'s private constructor. Real, but **pre-existing**
(not introduced by M0.9) and out of scope for this box — recorded here so it is not lost.

## 4. FALSIFIED / corrected during the run — recorded, not hidden

**F1 — "the `*_fast` / `*_ipo` binaries in `build-release/bin` are gate evidence." False.** They
are **stale artifacts from 2026-07-12**, not in the current Ninja build graph and never rebuilt
(`unit_test_core_Experiment_fast.exe` is dated 2026-07-12 07:39; today's binaries are 19:07).
One of them fails an assertion. It says nothing about current code: the live Release lane already
compiles with `-ffast-math`, and the live `unit_test_core_Experiment` passes the very test case
the stale binary fails (`108 assertions in 1 test case`). **OPEN for a later box:** these stale
binaries sit in a build directory and will keep polluting any glob-based measurement; a `*_fast`
probe that once failed a non-finite guard under some flag set deserves its own investigation, but
it is not M0.9's and no claim is made about it here.

**F2 — "Release assertion counts must match the Debug baseline." False, and the first reading was
wrong.** Three binaries differed (LithiumPlating 36→34, P1G0_pilot 4→3, SpectralDiffusion 10→9).
Two are pre-existing Debug/Release differences — a `df2da52` Release build reproduces 34 and 3
exactly. The third was **not** a Debug/Release difference and **not** IPO: `SpectralDiffusion`
gave 7 in a fresh baseline build and 9 in mine, and the cause was the build **option**
`SLIDE_ENABLE_RECORDED_SCALAR_BITS`, ON in `build-release` and OFF in a default configure. With
the option matched, the `df2da52` baseline gives **9**, identical to mine — and the two extra
assertions are the fast-math bit-hash checks at `tests/unit/core_SpectralDiffusion_test.cpp:235`,
which **pass**, so the split is bit-exact under fast-math. The band was only closed by building
the baseline like-for-like; comparing Release against a Debug baseline would have manufactured a
regression that does not exist.

**F3 — the CUDA lane needed no code fix but did need an environment one.** `nvcc` failed with
`Cannot find compiler 'cl.exe' in PATH`; MSVC (VS 18 Community,
`VC/Tools/MSVC/14.50.35717/bin/Hostx64/x64`) must be on PATH for `build-cuda4`. Recorded so the
next session does not misread it as a code break.

## 5. What the 9C-2 gate lost, and what replaced it

`p9c2_ageing_kernel.cmake` previously required `ageing_model_bit<N>(model)` and
`valid_ageing_model_mask<N>(p.model_mask)` *inside* `Sei.hpp` / `SurfaceCrack.hpp` / `Lam.hpp`.
Those tokens moved with the parameter blocks, so the gate now checks them in the `*Params.hpp`
headers and additionally requires each mechanism header to include its own params header. The
invariant it defends — one shared mask vocabulary, not four hand-rolled ones — is unchanged and
still enforced. What is *not* enforced any more: that the mask check sits in the same file as the
kernel that consumes it. That coupling was incidental, and nothing in PLAN requires it.

## 6. Not claimed

No performance, compile-time, or binary-size claim. TSan not rerun. The hosted CI workflows are
committed but not claimed as run. `SpmScalarKernels.hpp` remains reachable from a user TU by
design (PC-10 beats MC-5 there — see the design note); the api TU therefore still compiles the
scalar physics library, and that is stated rather than buried. R4 is a tripwire, not a parser: two
offsetting edits within one header could net out its declaration count.
