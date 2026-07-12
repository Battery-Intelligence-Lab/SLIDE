# M0.6 / 9C-2 ageing-kernel validation (2026-07-12)

## Result

**PASS.** `src/core/AgeingKernel.hpp` now owns construction-allocated checked
scratch plus force-inlined, accepted-step-allocation-free lane traversal and
pipeline-stage transaction, with mask/model traversal and clearing for
mechanisms that accumulate numbered alternatives. SEI, surface cracking, LAM, and lithium
plating retain four named scalar-generic equation bodies and their independent
legacy parity oracle. The validated source/test/CMake state is commit
`6fe960d25c180cead2a1e371e70fcdf4363151eb`.

## Architecture delivered

- `Sei.hpp`, `SurfaceCrack.hpp`, `Lam.hpp`, and `LithiumPlating.hpp` use one
  scratch and lane-sweep vocabulary while keeping their physics local.
- Model masks are checked by the same `ageing_model_bit`, full-mask, and
  optional-mask definitions in both mechanisms and `SpmFactory`.
- Numbered mechanisms traverse models in ascending order and each model
  traverses lanes in ascending order. This is a numerical contract, not an
  implementation detail.
- `SpmPipeline` invokes SEI, crack, LAM, and plating through one transaction:
  disabled is a no-op, compute failure returns before publication, and apply
  runs only after success. The established diffusion -> thermal -> SEI -> crack
  -> LAM -> plating order is structurally pinned.
- A disabled SEI stage explicitly clears the SEI scratch consumed by surface
  cracking, so no prior evaluation can leak through that dependency.
- Lithium plating now has the same named output/scratch boundary as the other
  mechanisms. The real `Cell_SPM_degradation.cpp` legacy equations remain raw
  and structurally forbidden from including the production scaffold.
- Shared traversals are forcibly inlined. Besides removing hot-path calls,
  this preserves opaque-reference finite checks under Clang fast-math; a
  merely-inline mutation made an existing Release overflow case return
  `Success`. LLVM 21.1.8 `llvm-nm --defined-only --demangle` found zero defined
  `clear_ageing_fields`, lane/model traversal, or `evaluate_ageing_stage`
  helpers in `build-release/bin/core/slide_core.lib`; the core-only and nested
  consumer archives give the same result.

Construction-time hardening found adjacent instances of the same extent bug.
Ageing, stress, observable, transport, and pipeline scratch plus `StateArena`
now reject invalid extents before conversion, signed padding/multiplication, or
allocation and guard element and byte representability. These are cold-path
exceptions; the accepted-step path remains allocation-free and Status-based.

## Frozen trace and provenance

The fixture was committed at `70f0cc6` against production state `716a1a1`,
before the scaffold existed. It enables all SEI (`0x0f`), surface-crack
(`0x1f`), and LAM (`0x0f`) models with plating, SEI porosity loss, and crack
diffusivity coupling across seven heterogeneous lanes. Current patterns A, B,
A fingerprint the initial arena, every derivative row, and terminal voltage;
the repeated A result must be bit-identical, independently checking scratch
reset. The trace contains 1,077 doubles.

`SLIDE_ENABLE_RECORDED_SCALAR_BITS` is OFF by default and CMake restricts it to
x64 Windows Clang 21.1.8. Reproduction also needs a compatible CRT and ISA
because Release uses `-march=native`. The hashes are provenance checks, not
portable mathematical truth. CUDA evidence used CUDA 13.0, architecture sm_89,
on the development-host RTX 4000 Ada.

| Configuration | Values | FNV-1a | Independent mix | Final result |
|---|---:|---:|---:|---|
| Debug/ThinLTO (`-O0`) | 1,077 | `87119b1b6fa83b81` | `55296fa07292792f` | exact |
| Release/fast-math | 1,077 | `dc8f92a59dd8d67f` | `59085638e06725d3` | exact |
| Release/fast-math/ThinLTO | 1,077 | `0be580cf849e57e1` | `4c6451971b74f789` | exact |

Portable toolchains still enforce value count, finite outputs, A/B/A reset,
independent legacy parity (including an explicit nonzero crack-model-5
witness), Dual finite differences, lane isolation, stage atomicity, allocation,
and structural gates without adopting a new hash.

## Adversarial findings and fixes

1. **Saturated crack sensitivity.** Both maximum-surface sites selected
   `primal_value(crack_surface)`, which retained the right value but discarded
   the selected Dual tangent above the ceiling. The old path produced a
   `-8e-9` model-4 tangent while a centered finite difference was zero. The
   piecewise implementation now retains the scalar-generic selected operand
   and spells the saturated zero fraction explicitly, avoiding a fast-math
   fractional power of a rounded negative zero-neighbour.
2. **Dual zero-exponent boundary.** The generic derivative formula evaluated
   `0 * pow(0, -1) * 0` when the crack diffusivity exponent was one, producing a
   NaN tangent although the primal was finite. `pow(Dual, 0)` now returns the
   exact constant `(1, 0)` and exponent one returns the original Dual.
3. **Unchecked scratch extents.** Several constructors converted or multiplied
   signed lane counts in vector initializers and asserted only afterward.
   `SpmTransportCache` additionally evaluated `2 * lanes` as signed `int`, and
   `SpmPipeline` padded with signed `n_lanes + 7`. Zero, negative, and oversized
   regressions now fail before allocation; `StateArena` also checks padded
   stride, elements, and bytes.
4. **Degenerate model-5 parity.** The previous direct legacy comparison could
   pass with a zero crack rate. Its state is now chosen so the expected and
   production rates are finite and nonzero.
5. **Missing call-site atomicity proof.** Direct mechanism failures and the
   generic transaction test did not prove the real pipeline wiring. A full
   `SpmBatch::evaluate` regression now reaches an SEI `Invalid_states` failure
   after upstream observables, proves the Status is returned, and proves no SEI
   or cumulative derivative is published.
6. **Ambiguous failure-site regions.** Directly returning an inline traversal
   whose lambda contained literal failure returns made exact LLVM source-site
   attribution ambiguous. SEI, crack, and plating now bind the traversal result
   to `lane_status` before returning it. The change is arithmetic-neutral, all
   three recorded hashes remain exact, and the coverage lexer now resolves 369
   failure arms without overlapping return shapes.

Controlled mutations distinguished the intended contracts:

- replacing the shared SEI traversal with a raw loop failed the structural gate;
- reversing ascending model order changed the Release FNV hash to
  `68364d42ddb458f9` and the independent mix to `6f217e9f9acf6596`;
- aliasing a LAM scratch field failed the field-mapping structural gate;
- reversing the Dual power derivative sign failed the finite-difference gate;
- retaining the primal-only saturated maximum failed at `8e-9` versus zero;
- retaining the generic zero-exponent formula produced a NaN tangent;
- accepting a zero transport-cache extent failed the constructor regression;
- swapping the expected crack/LAM stage order made the new ordered structural
  check fail at the first out-of-order token;
- removing forced inlining made the existing Release surface-crack overflow
  regression return `Success` rather than `Numerical_failure`.

## Test evidence

`unit_test_core_AgeingKernel` now contains 384 assertions in six cases on the
recorded-bits capture lanes and 382 assertions when recorded bits are disabled.
Besides the frozen trace it checks at least one Dual finite-difference point in
each of the four mechanism bodies, both sides of the crack ceiling including
exponent one, seven-lane versus seven isolated
evaluations, masks/scratch/mappings/traversal, invalid extents, transaction
semantics, and the integrated failed-stage boundary. The warmed 257-lane full
ageing RHS remains at exactly zero allocations.

| Lane | Result | Notes |
|---|---:|---|
| Windows Clang Debug/ThinLTO | 56/56 | serial CTest, 26.24 s; ModeC 18.06 s |
| Windows Clang Release/fast-math | 56/56 | serial CTest, 13.15 s; ModeC 6.16 s |
| Windows Clang CUDA/fast-math/ThinLTO | 56/56 | serial CTest, 14.29 s; CUDA batch 1.54 s; ModeC 6.47 s |
| WSL Clang 18 ASan+UBSan | 56/56 | serial CTest, 160.29 s; ModeC 146.79 s; no ASan, LSan, or UBSan finding |
| WSL Clang/LLVM 18 exact Status coverage | 56/56 | 287.16 s; ModeC 274.49 s; 54 fresh test profiles + one whole-archive anchor; 369 lexical, 339 active = 332 measured + 7 structural exceptions, 30 inactive, zero uncovered/unmapped |
| Windows core-only optional-off | 1/1 | dependency rebuild + immediate no-op verification; CTest 0.22 s |
| Windows nested external consumer optional-off | 1/1 | dependency rebuild + immediate no-op verification; CTest 0.08 s |

TSan was not repeated. M0.2 passed its ThreadPool/AsyncRecorder/PackStepper
lane 3/3 on the earlier M0.2 source state. M0.6 changes no synchronization,
worker lifecycle, or shared ownership, so no new TSan claim is made.
Hosted workflows are updated for 56 discovered tests and 54 executable coverage
profiles plus one whole-archive anchor, but no hosted run is claimed before
push.

The authoritative sanitizer build performed 82 final-state actions in
1,408.52 s and a subsequent no-work check in 2.39 s. Its only diagnostics were
12 pre-existing `units.hpp` double-promotion warnings in PAY1/PAY2 benchmark
compiles; no M0.6 source emitted a warning.

## Maintainability accounting

The four mechanism headers together move from 882 to 895 nonblank lines; the
13-line net growth buys the named plating output, saturated-branch fix, and
unambiguous exact-coverage site boundaries.
The shared scaffold is 163 nonblank lines. `SpmPipeline.hpp` is now 746
nonblank lines because checked construction adds explicit lane/row/element/byte
proofs beside the compile-time composition. That remains one orchestration
concept but crosses the MC-1 review threshold; this is recorded debt for the
M0.10 line-reduction/dead-code gate rather than hidden as a successful size
reduction in M0.6.

Artifacts:

- `.claude/reports/p9c2-pre-refactor-ageing-fixture-2026-07-12.md`
- `.claude/reports/p9c2-status-coverage-refresh-2026-07-12.md`
- `.claude/reports/p9c2-status-coverage-refresh-2026-07-12.json`
