# MQ.2 quality-backlog validation (2026-07-23)

## Status and provenance

This is the accumulating post-run evidence for the preregistration in
`.claude/reports/mq2-quality-backlog-preregistration-2026-07-23.md`.
MQ.2 remains open: this report records completed batch boundaries but does not
yet claim the final 71-row census or the final three-lane gate.

The registered source floor is commit `5481496` (`Freeze MQ.2 pre-edit
fixtures`). The completed boundaries recorded so far are:

| Commit | Batch | Character |
|---|---|---|
| `ec7ac09` | O2 fast-math validity | invalid-state behavior hardening |
| `025108f` | O3 stress/thermal oracles | test-only |
| `ac15536` | F0 factory byte oracle | oracle-only, before factory source edits |
| `b9c4db1` | F1 shared factory constants | exact no-op; homogeneous rewrite deferred |
| `43c4953` | A1.0 surface-crack oracles | oracle-only, before ageing source edits |
| `a699da5` | A1.1 SEI scalar ownership | exact no-op; function-boundary hypothesis falsified |
| `e34f548` | A1.2 SurfaceCrack Arrhenius | registered low-bit association correction |
| `e034d23` | F1 direct fused-Euler validation | invalid-input behavior hardening |
| `83ab7be` | O1 observable storage/access | exact no-op plus checked cold boundaries |
| `6c7ab74` | P0 pack-algebra oracles | oracle-only, before P1/P2 source edits |
| `d4e1a2d` | A1 structural-anchor correction | test-gate baseline repair, no source change |
| `8630207` | P1 pack-solver algebra | exact no-op plus executable ownership gate |
| `bf185dd` | P2 Mode-C relaxation storage | exact no-op plus seven-role scratch ownership |
| `62f1587` | S1 PackStepper ownership | exact no-op plus executable full-`dt`/frozen-solve contract |

The retained build trees are `build-mq1-debug` (Clang 21.1.8 Debug/ThinLTO),
`build-mq1-release` (Clang 21.1.8 fast-math Release, IPO off), and
`build-mq1-cuda-vsenv` (Clang host-C++ Release/ThinLTO plus the configured CUDA
lane). The last tree requires the complete VS 18 x64 `vcvars64.bat`
environment. The deprecated Clang `-Ofast` warning remains the registered B1
input; no timing evidence is inferred from these runs.

## O2 — fast-math-safe surface stoichiometry

Files:

- `src/core/SpmObservables.hpp`
- `tests/unit/core_SpmElectrical_test.cpp`
- `tests/structural/p9c_architecture.cmake`

The structural gate and the opaque qNaN/+Inf/-Inf test were added before the
production predicate. On the old production source, the structural command

```text
cmake -DSLIDE_SOURCE_DIR=C:/D/git/SLIDE -P tests/structural/p9c_architecture.cmake
```

failed with:

```text
required token is absent: #include"Numeric.hpp"
```

Both modified-test/old-source binaries already rejected the three injected
values on this compiler, so no runtime RED is claimed. The registered decisive
RED was the structural requirement for the fast-math-safe predicate.

The production rejection slice now contains exactly one compact
`!is_finite_primal(z_surface)` guard before the open-interval comparison.
Focused Debug and fast-math Release builds and direct runs each passed
`45 assertions in 2 test cases`, up from the frozen `39 / 2`. The valid control
has positive diffusion in both domains and returns `Success`; opaque bit-built
qNaN, +Inf, and -Inf modal states each return `Invalid_states`.

Removal mutation: deleting only `!is_finite_primal(z_surface)||` made the
structural gate fail:

```text
expected 1 occurrences ... found 0
```

After the exact inverse patch, the three touched-file SHA-256 values were:

```text
0574ae338710d13be65230e6b1cb7c44beb4781299e02adfd6f551100eaf82c6  src/core/SpmObservables.hpp
ef0c81322c2e17ad9dec27836f8d174383f40563e9a083ec126e5a543bff5726  tests/unit/core_SpmElectrical_test.cpp
82a15436949eae47ad495e103738a4db24cf8a95570dd4f58247990c58bb755b  tests/structural/p9c_architecture.cmake
```

The restored structural gate and `git diff --check` passed. The old runtime
behavior is deliberately not generalized beyond this compiler/configuration.

## O3 — independent stress and thermal oracles

Files:

- `tests/unit/core_SpmStress_test.cpp`
- `tests/unit/core_ThermalLumped_test.cpp`

No production file changed. In both Debug and fast-math Release:

| Binary | Frozen floor | Result |
|---|---:|---:|
| `unit_test_core_SpmStress` | 19 / 3 | 40 assertions / 4 cases |
| `unit_test_core_ThermalLumped` | 27 / 2 | 44 assertions / 2 cases |

The four-lane stress oracle uses distinct lane concentrations so node-major
stride mixing is observable. Its registered normalized `1e-10` bands passed
with these maxima:

| Mode | uniform profile | constant-offset difference |
|---|---:|---:|
| Debug | `1.4e-16` | `5e-17` |
| fast-math Release | `1.3e-16` | `6e-17` |

The test states its blind spots: a zero/offset invariant cannot validate the
material factor's magnitude or sign, nor the deliberate legacy node
association. Those remain covered by legacy replay pending MQ.5 derivation.

The thermal additions exercise conductance overflow with atomic parameter
clearing, a finite control followed by separate internal/external opaque-qNaN
rejections, and an adiabatic zero-heat limit. Adversarial review found that the
first adiabatic draft reused `T == T_env`; before commit it was strengthened to
`T = 310 K` with `T_env = 300 K`. Both lanes still pass exact
`dT/dt == 0` and `d(t_thermal)/dt == 1`.

## F0 — padded NCH=12 factory byte fixture

Files:

- `tests/unit/core_SpmFactory_test.cpp`
- `tests/unit/CMakeLists.txt`

This oracle landed before any factory production edit. It constructs nine
lanes with thermal plus every ageing family enabled, including nonzero LAM
coefficients and deliberately conflicting mechanism-local copies of
factory-owned shared constants. Geometry is pinned to 50 rows, 9 live lanes,
stride 16, and 800 doubles per arena.

The test hashes three separate 6,400-byte domains with Boost.Hash2 SHA-256:
the full padded state immediately after construction, the full padded
derivative immediately after construction, and the full padded state after
one successful `1e-3 s` exponential step with nine distinct currents. It also
checks every padding value in state and derivative, before and after the step,
by `uint64_t` bits against positive zero.

Each configuration first ran with three placeholder digests. The intentional
capture RED was:

```text
1414 assertions: 1411 passed, 3 failed
```

Only the three digest comparisons failed. The captured constants, ordered
initial state / initial derivative / post-step state, are:

```text
Debug
f9c0df074580cad37eab525a428eb75641626e1398d3f3fb71eeffe7c10278fa
d7746f27ecd3de2efccc65458bc429421b80d45f3984bb2fe4f28d3d2c01d476
d5de6e0fb2fc2b550728922066c349a3678032421cc3a7962ecca2d412896110

fast-math Release, IPO off
fec245fb6f6751b570f02edbcb00c2f84394b7c28e5652e5d151c464f9461f18
ad6f936b1087ff2a8647c6f1f69b5ef2e03ac374f7ce9c540f03f1b4095f7a08
119011a365f6bb438e63d13aafa3faafcaa3f8bcf65cae0df153a161182b94f4

fast-math Release, host ThinLTO / CUDA tree
a1c372842363e6109d609a217067b51e6af24a0a77480236f58f99c8d76d7955
3f5de3d1146e51b5d8d5755130753f329944a7c9cada50e5909aa348e76ed178
5f72e3a5fbc7ff60d7c63d21f7ee7b11a9299b1aaba7930c3c0ac64b1af53705
```

These are per-configuration pre/post identities, not a cross-configuration
equality claim. The final fixture selects them through the configured
`SLIDE_TEST_RELEASE` and `SLIDE_TEST_IPO` definitions, not compiler incidental
macros.

Final direct-binary results after the permanent constants were installed:

```text
Debug:                           1472 assertions in 8 test cases
fast-math Release, IPO off:      1472 assertions in 8 test cases
Release host ThinLTO/CUDA tree:  1472 assertions in 8 test cases
focused [MQ.2][recorded] case:   1414 assertions in 1 test case
```

`clang-format --dry-run --Werror` and `git diff --check` passed. No factory
source had changed at this evidence boundary.

## F1 — shared factory constants and falsified homogeneous rewrite

The repeated factory propagation of `F`, `Rg`, reference temperature, and the
mechanism-dependent shared fields now has one anonymous-namespace owner.
`sei_resistivity_area` is read from the already-compiled electrical parameter
block, and the negative spectral input map is bound once by reference. Each
helper call remains after its whole-struct mechanism copy, so the factory-owned
values still win over the deliberately conflicting F0 sentinels.

The final shared-constant-only source passed both exact fixture families in all
three configurations:

| Binary | Debug | Release | Release/ThinLTO CUDA tree |
|---|---:|---:|---:|
| `unit_test_core_SpmFactory` | 1472 / 8 | 1472 / 8 | 1472 / 8 |
| `unit_test_core_AgeingKernel` | 384 / 6 | 384 / 6 | 384 / 6 |

All nine F0 SHA-256 constants and all three mode-specific 1,077-value
AgeingKernel hashes remained unchanged.

The registered hypothesis that the homogeneous initialization loop could be
removed as a digit-identical cleanup was **FALSIFIED** after three different
implementations:

1. **Compute each lane-independent value once, then fill live rows.** Debug
   stayed exact, but fast-math Release changed the F0 hashes to
   `b36c14f510de590f7702d9f9a58b543df536aea7d783e93ed05aa651d9a54907`,
   `b37a3718ce1bcb5a1f794ea57ab2954219203902fdadf5a0be7dfdd83266bcc4`,
   and
   `163eac6fd30c298db82fb77f0f6c9a4859ae37597d09b4046e291731f70d616e`.
2. **Cache the first lane's scalar results inside the runtime lane loop.**
   Fast-math Release changed the negative zero mode by one ULP
   (`bfdf03f288bd0ed6` to `bfdf03f288bd0ed7`), with F0 hashes
   `2fa3461c9e33075ea387a0efccbdcfd2feb1ca0a9071350f1d6cb09c7666780e`,
   `f438aecdadd12fdf19cc81f16840acc5a0f1893cdadb673818359d5696d6c677`,
   and
   `e4e641ca37f611fff5b612dc263090f0907c682a33286399f0c5e67e2f9ce1fc`.
3. **Keep the lane-zero arithmetic block verbatim and copy all later arena
   rows from lane zero.** The NCH=12 F0 hashes passed, but the independent NCH=5
   fast-math all-ageing fixture moved from
   `dc8f92a59dd8d67f / 59085638e06725d3` to
   `8ce61773119858d1 / 7aaf4294740b7321` (the 1,077-value count stayed fixed).

No changed hash was accepted. The initialization source was restored exactly,
and the final Debug/Release/ThinLTO runs above are green. The finding
`homogeneous-init-lane-loop` therefore changes from preliminary `APPLIED` to
`DEFERRED — MQ.7 structural performance hunt`: that box must first establish
whether this cold-path operation count warrants an optimizer-sensitive
numerical change. This changes the predicted original-68 distribution to
59 applied, one refuted, and eight named deferrals; it does not alter the
71-row census.

## F1 — direct fused-Euler validation

Commit `e034d23` changes the direct base-archetype fused-Euler boundary from a
Debug-only output-size assertion plus implicit assumptions to a Release-active
transactional guard. Before constructing a state view for lane-period
detection, it now requires the pipeline's exact row/lane counts, exact
current/output span sizes, finite `ctx.time`, `ctx.dt`, explicit `dt`, and every
current, plus strictly positive explicit `dt`.

The preregistered old-source SHORT witness was red: a zero explicit `dt`
returned `Success` and changed terminal output, so 6/8 assertions passed and
two failed. The permanent table has 21 invalid rows: undersized and oversized
state row/lane counts; empty, short, and long current/output spans; opaque
qNaN/+Inf/-Inf currents; zero, negative, qNaN, and +Inf explicit steps; and
qNaN/+Inf context step/time. Every row requires `Invalid_parameters`, a
byte-identical complete candidate arena, and all three terminal sentinels
unchanged. The valid control additionally requires a changed arena and both
published voltages to replace their finite sentinels, so returning `Success`
without doing work cannot satisfy the gate.

Final clean direct-binary results:

| Binary | Debug | fast-math Release | host-ThinLTO/CUDA tree |
|---|---:|---:|---:|
| `unit_test_core_SpmPipeline` | 114 / 2 | 114 / 2 | 114 / 2 |
| `unit_test_core_SpmFactory` | 1472 / 8 | 1472 / 8 | 1472 / 8 |
| `unit_test_core_AgeingKernel` | 384 / 6 | 384 / 6 | 384 / 6 |

The factory and all-ageing exact fixtures retained their previously recorded
hashes. Three adversarial mutations were independently red:

1. changing exact row/lane comparisons from `!=` to `<` accepted both
   oversized shapes and failed 6/96 focused assertions, including state and
   output atomicity;
2. returning `Success` immediately after validation failed 3/96 assertions
   because neither arena nor terminal output changed;
3. deleting the strict-positive explicit-step predicate failed 4/96
   assertions: zero `dt` returned `Success` and published output, while
   negative `dt` reached later physics and mutated state before rejection.

Moving the guard below `lanePeriod` made the deliberately null empty-current
span reach indexed period detection and hang the mutated process; it was
terminated, the mutation was reversed explicitly, and the final source hashes
were restored to:

```text
SpmPipeline.hpp             6786685528480BA65129136896D8823ADCB457A9ABCCDFEBD844D8F633C770FF
core_SpmPipeline_test.cpp   63159E07680EE783377CF7D952543C7A18329C210142F38FEBB6F4E59DB5BD95
```

`clang-format --dry-run --Werror` passes the changed test file and
`git diff --check` passes. The pre-existing formatting debt later in
`SpmPipeline.hpp` remains outside this batch; no whole-header formatting claim
is made.

## O1 — observable storage, cache ownership, and lane lexicon

Commit `83ab7be` applies all six O1 findings. `AgeingScratchStorage` and
`SpmObservableScratch` now share the construction-only
`detail::checked_lane_extent` owner while retaining their distinct,
byte-identical invalid-count and overflow messages. Observable scratch sizing
names the six per-domain and nine shared fields, guards every subspan before
formation with a non-underflowing bound, and asserts exact final consumption.
The concentration kernel obtains its first modal value through
`BatchView::at`; the SEI model body binds its unsigned lane index once.

`SpmTransportCache` now exposes only `n_lanes`, `try_load`, and `store`; its
seven buffers and indexing are private. The cache operations own only exact
key comparison and byte transport. Both Arrhenius/activation/denominator/flux
arithmetic branches remain in `computeSpmTransportLane`, pinned structurally.
The public `TheveninBatchView::lanes()` spelling becomes `n_lanes()` and is
called out in Unreleased. Duplicate includes are rejected only across
classified api/support headers; the three known non-core duplicates are not
silently claimed fixed.

The first attempt to include the new detail owner exposed a stale 9C-3 rule
that described every top-level core header as public. A proposed broad
internal-header exemption was rejected by adversarial review because it would
let unrelated cold detail seams leak into `SpmPipeline.hpp`. The final gate
preserves the original global prohibition and removes exactly one allowlisted
`CheckedLaneExtent.hpp` edge from each of `AgeingKernel.hpp` and
`SpmObservables.hpp` before scanning.

Focused clean results:

| Binary | Debug | fast-math Release | host-ThinLTO/CUDA tree |
|---|---:|---:|---:|
| `unit_test_core_SpmObservables` | 69 / 3 | 69 / 3 | 69 / 3 |
| `unit_test_core_SpmElectrical` | 48 / 3 | 48 / 3 | 48 / 3 |
| `unit_test_core_Sei` | 57 / 3 | 57 / 3 | 57 / 3 |
| `unit_test_core_AgeingKernel` | 384 / 6 | 384 / 6 | 384 / 6 |
| `unit_test_core_PackSolver` | 878 / 26 | 878 / 26 | 878 / 26 |
| `unit_test_core_SpmFactory` | 1472 / 8 | 1472 / 8 | 1472 / 8 |

The pre-existing SpmObservables numerical bands remain unchanged (Debug
legacy maximum relative error exactly zero; Release/CUDA
`2.665e-15`; heterogeneous round-trip `9.027e-15` Debug and
`8.533e-15` Release/CUDA). The 1,077-value ageing hashes and all nine padded
factory SHA-256 values remain at their preregistered constants.

Two new Release-visible oracles close the actual storage seams. A three-lane
scratch spans exactly
`L * (2 * (NCH + 2) + 2 * 6 + 9)` values from the first negative
concentration value through the end of `total_heat`. A two-lane cache fixture
requires cold miss, unchanged-key hit, and each of temperature, reference
diffusivity, specific area, and thickness invalidation to reproduce the
cache-free diffusivity and molar-flux bits in both domains.

Adversarial mutations:

1. Reintroducing the duplicate Status include fails R3 with the exact repeated
   header named.
2. Restoring either `lanes()` or whitespace-obfuscated `lanes ()` fails R6.
3. Restoring `state.raw()` addressing fails the exact zero-count gate.
4. Taking `L-1` values for `total_heat` fails the Release scratch case
   (one of two reached assertions fails).
5. Adding an unreachable scalar-kernel call inside the cache fails the
   cache-slice ownership gate.
6. Swapping the cached diffusivity and denominator initially demonstrated
   that the existing ageing 15/1 and factory 1414/1 recorded cases were
   **non-discriminating**: both stayed green. After registering the direct
   cache oracle, that same mutation fails 30/48 assertions. It was then
   explicitly reversed and all three clean configurations passed 69/3.

The final key source hashes are:

```text
CheckedLaneExtent.hpp       ED7A33A1B4638D2EB306F818247B88B20BEF360221284C49A438E366FEA6A552
AgeingKernel.hpp            FE89EF1AD5470996608A94EC7234341017FF7BE687128B5B05E40692212C95AB
SpmObservables.hpp          65BA7D54F8C9FEBDBF68F34EB19EDF7F74510CE5D0B3E28F4B345C15FBE6FD4E
Sei.hpp                     0B706DE2BFDB5C4B708812B240A24EEC139B8B04CE3F5F1E8E2B9A047E5B0C20
PackSolver.hpp              5A05552D245E0DAD96B70BEF313BD7984AD39BBE3C1BDE4FD53B42E7DABB23F3
PackSolver.cpp              A5155224C2BAC572BF4CD785859988A3C8E23BCC5AB6CC097E4C214B32C088E6
```

Both structural binaries pass, changed C++ files pass
`clang-format --dry-run --Werror`, and `git diff --check` is clean.

## A1 — ageing scalar ownership and SurfaceCrack association

### A1.0 — association-sensitive oracles before source edits

Commit `43c4953` extends the existing direct legacy test before either ageing
source changes. The production outputs for models 1–4 are recorded as 24
doubles (three outputs, two diffusivity settings, four models) with two
independent bit recurrences:

| Configuration | FNV-1a | mixed |
|---|---|---|
| Debug | `26095ff417f83ca0` | `af0313f19f95e3dc` |
| fast-math Release, IPO off | `4f691333ce8072a6` | `823f0fffde2c3d18` |
| Release/ThinLTO CUDA tree | `dc30ce9466aa1ae3` | `699574956434b620` |

Each pair was first captured through two deliberate zero-placeholder
failures: Debug passed 93/95 assertions and each Release mode passed 92/94.
After installing only the observed constants, the three modes passed 95/4,
94/4, and 94/4. This freezes current production, not legacy: an exploratory
bitwise comparison found that Debug model 3 already differs from legacy by one
ULP in crack rate and two ULPs in the reduced-diffusivity rate, so that wrong
identity target was not retained.

The model-5 fixture independently assembles both associations at `T = 310 K`,
uses a nonzero low-stoichiometry charging branch, and evaluates the complete
crack-side reaction rather than only its exponential. In Debug, the old and
shared results are respectively
`0x3e05bf8574544a7a` and `0x3e05bf8574544a77`, a relative delta of
`4.9e-16` against the registered `1e-13` band. The corrected production result
equals the shared result in Debug. The rate is approximately
`6.3295389e-10`, so a zero branch cannot satisfy the test. In fast-math Release
the two scalar associations round identically; this is why the strict
bit-inequality witness is limited to the explicitly recorded compiler mode,
while the numerical oracle is portable.

Seeding temperature as `Dual{310, 1}` gives a derivative of approximately
`-1.2401496e-10`, equal at the printed precision to the centered
`h = 1e-4 K` finite difference and within the registered `1e-8` relative
band. The Release Dual value differs from its double value by only
`1.63e-15` relative, within the registered `1e-13` band.

### A1.1 — exact SEI single ownership

The five activation expressions now call `spm_scalar::activatedValue`. The
three repeated kinetic-current expressions use the one
`SLIDE_SPM_SEI_KINETIC_CURRENT` owner in `SpmScalarKernels.hpp`; branch policy
and the two distinct driving potentials remain local to `computeSei`.

The originally preregistered ordinary function extraction was
**FALSIFIED**, not silently reblessed. Debug and IPO-off Release stayed exact,
but Release/ThinLTO changed:

```text
AgeingKernel:
  0be580cf849e57e1 / 4c6451971b74f789
  -> 01ada7c263317e7f / b2d2672a19e42c12

Factory initial derivative:
  3f5de3d1146e51b5d8d5755130753f329944a7c9cada50e5909aa348e76ed178
  -> 3c712d4b1d5c625c29f8cb2b43468da904f103847ad3c371cd6601e243a4c2b6
```

Three distinct implementations produced those same changed values: the
preregistered `<Real, Scalar>` function with reference operands, a forced
`always_inline` form, and a fully deduced by-value `auto` form. Keeping the
five activation-owner calls while restoring only the kinetic expressions
restored both fingerprints, isolating the inline boundary. The macro remedy
was then registered before its first build/run; it mirrors the modal kernel's
already-established optimizer-sensitive one-source pattern and passes the
caller `exp` token so ADL still selects `Dual`.

No hash was accepted. With the caller-expanded owner, `core_Sei` remains 57/3,
`core_AgeingKernel` 384/6, and `core_SpmFactory` 1472/8 in Debug, fast-math
Release, and Release/ThinLTO. Every pre-A1 ageing/factory fingerprint remains
unchanged.

The kinetic exponent is dimensionless: `F/(Rg*T)` has units `1/V` and its
driving potential has units `V`. The Arrhenius exponent is likewise
dimensionless:
`(J/mol) * ((1/K) / (J/(mol K))) = 1`.

### A1.2 — registered SurfaceCrack association correction

Model 5 now calls `spm_scalar::arrheniusFactor` directly and evaluates
`Ea * ((1/Tref - 1/T)/Rg)`. The `2*k` versus `k` branch multiplication and the
final crack-side Butler–Volmer ordering are deliberately untouched. The
independent oracle above passes after the source change, every model mask
remains within `1e-12` of legacy, and the frozen production bits for models
1–4 remain exact.

Focused final results:

| Binary | Debug | fast-math Release | Release/ThinLTO CUDA tree |
|---|---:|---:|---:|
| `unit_test_core_SurfaceCrack` | 95 / 4 | 94 / 4 | 94 / 4 |
| `unit_test_core_Sei` | 57 / 3 | 57 / 3 | 57 / 3 |
| `unit_test_core_AgeingKernel` | 384 / 6 | 384 / 6 | 384 / 6 |
| `unit_test_core_SpmFactory` | 1472 / 8 | 1472 / 8 | 1472 / 8 |

Although a low-bit change is observable in the dedicated Debug fixture, none
of the existing all-ageing or padded-factory fingerprints moves for their
inputs, so no recorded constant was reblessed.

Four adversarial mutations turned red:

1. restoring one private activation expression made the structural count
   report four rather than five calls;
2. restoring one private kinetic expression made the structural count report
   two rather than three macro calls;
3. giving model 1 the film-free driving potential failed direct legacy parity
   with relative error `0.74232267435407617`;
4. restoring the old SurfaceCrack association made the required shared-owner
   token disappear from the structural gate.

Each mutation was explicitly reversed. Final source SHA-256 values are:

```text
SpmScalarKernels.hpp          4B54927A59CA88D9863AA4325AA1F59089961A03FA1F66986775D12E704117F8
Sei.hpp                      9B5BA40BF4BB8D96E03F0CA4D0342AB316AA3E16FAA30340FA8014CB45388214
SurfaceCrack.hpp             A759F6E01E28271BB0342CD5C5EC0F0508A26F8D9490CC24361CFAAB391ADBAB
core_SurfaceCrack_test.cpp   CE5BB09FB355F7B34193826D675F3C377CB06290FFA2246367B85FF698DDF9B3
pc10_single_source.cmake     641AE57005A8BFFA3309F965E85B3182630F6B09574B0AFB706B3FE933711C43
```

The final structural gate, changed-file formatting checks, and
`git diff --check` pass.

## P0 — pack-algebra oracle reconciliation

### Registration-only commits and dirty-tree disposition

The four commits called out by PLAN MQ.2 were audited against their exact
diffs. They changed only this batch's preregistration report:

| Commit | Durable registration | Disposition |
|---|---|---|
| `c2940a2` | exact pack-algebra Traces A–D | retained as P0/P1/P2 preregistration |
| `87021a5` | P1/P2 helper, token, storage, fill, and mutation gates | retained as future P1/P2 preregistration |
| `39e6a00` | scripted source-step rollback fixture | retained as future S2 preregistration |
| `b1625fb` | D0 falsification and revised bit-separating Trace D | retained; D0 hashes remain non-acceptance evidence |

They are therefore not implementation commits and apply no survivor
disposition. The unexplained working tree was a coherent P0 oracle boundary:
`tests/unit/CMakeLists.txt` enabled the established recorded-scalar fixture for
PackSolver, and `tests/unit/core_PackSolver_test.cpp` added Traces A–D plus an
independently coded KCL branch walk. It landed without production changes at
`6c7ab74`.

The tracked PackSolver test grew from 1,198 to 1,510 lines. This is recorded
debt, not a new silent inheritance: PLAN M1.0 already owns the split and its
source count is updated to the new floor. P0 is oracle-only; none of
`cell-current-reconstruction-x4`, `branch-drop-and-kcl-current`,
`dead-usings-and-misplaced-comment`, or `relaxation-target-triple-role` is
called APPLIED here.

### Three-configuration exact gate

Before the first P1/P2 production edit, each retained build reported no work
for the current source and the direct PackSolver binary passed:

| Configuration | P0 `[MQ.2]` | Complete PackSolver |
|---|---:|---:|
| Debug/ThinLTO | 71 assertions / 4 cases | 949 / 30 |
| fast-math Release, IPO off | 71 / 4 | 949 / 30 |
| Release/ThinLTO CUDA tree, host C++ | 71 / 4 | 949 / 30 |

The committed pre-P0 floor was 878/26, so the delta is exactly the four
registered cases and 71 assertions. Trace D retained its registered value
counts `{13,13,13,15}` and its two frozen recurrences in all three
configurations. These bit records are regression provenance, not a portable
mathematical oracle; Traces A–C and the independent KCL walk supply the
analytic checks.

### Adversarial sensitivity

At clean commit `6c7ab74`, each of the four raw cell-current reconstruction
sites was changed separately from direct division to reciprocal
multiplication. The Debug Trace D case went red every time:

| Mutated path | Trace index | Mutated FNV-1a / mixed | Result |
|---|---:|---|---|
| sparse undamped | 0 | `67674e77f46ba2b8` / `c61c787dcc8a9396` | 24/26 passed; both hashes failed |
| sparse damped | 3 | `0664dfa152e6be22` / `5db2a6502d3dc5d2` | 24/26; both hashes failed |
| ladder | 1 | `ea9ba4b556eddad4` / `80602aa74e05476d` | 24/26; both hashes failed |
| relaxation | 2 | `2f3c710b8f14bdba` / `50827d75714062b9` | 24/26; both hashes failed |

This overturns the D0 gate weakness without reusing or reblessing D0: revised
Trace D observes the intermediate callback frame that the old dyadic fixture
missed.

Two oracle mutations were also red. Flipping the independent cell-branch KCL
sign failed 3 of 71 assertions (68/71 passed), while omitting callback frame
two changed the first record from 13 to 11 values and failed at 14/15 reached
assertions. After every mutation, exact restoration was verified against:

```text
PackSolver.cpp             A5155224C2BAC572BF4CD785859988A3C8E23BCC5AB6CC097E4C214B32C088E6
PackSolverIterative.cpp    382057BED8F21A95DD08F970C5210181341EA9CA43914D17425DA8E0BF8F8856
core_PackSolver_test.cpp   19ED0074E61849B27CF5CB2650D1E027A2CC7A6D46EF205DF9C869A41883703F
```

The restored Debug binary passes 949/30. `clang-format --dry-run --Werror`,
the shared-harness structural gate, and `git diff --check` pass. No full-suite,
timing, sanitizer, coverage, hosted-CI, installed-package, device-code, Linux,
or macOS claim is made at this oracle boundary.

## P1 — single-owner pack affine algebra

### Pre-source structural RED and stale-gate repair

The registered P1 gate was added before the production refactor. Its first
aggregate run stopped on an older contradiction rather than the intended P1
absence:

```text
9C-5 R4: SpmScalarKernels.hpp declares 123 member-level entities,
the pinned surface has 122.
```

A1's registered `SLIDE_SPM_SEI_KINETIC_CURRENT` owner had added the one anchor,
but `p9c5_public_surface.cmake` was never updated. Commit `d4e1a2d` changes
only that exact 122→123 baseline. With the stale gate repaired, unchanged
pre-P1 production then failed at the intended boundary:

```text
MQ.2 P1 cell-current owner: expected 1 occurrences of
cellCurrentFromDrop(, found 0
```

This corrects the earlier broad “structural gate passed” wording for A1: its
PC-10 gate passed, but the aggregate public-surface anchor was stale until
`d4e1a2d`.

### Implementation and structural ownership

Commit `8630207` applies `cell-current-reconstruction-x4`,
`branch-drop-and-kcl-current`, and `dead-usings-and-misplaced-comment`.
`PackSolverInternal.hpp` now has five algebra owners, with floating inputs
crossing the inline boundary by `const real_t &`:

| Owner token | Header | Sparse TU | Iterative TU | Total |
|---|---:|---:|---:|---:|
| `cellCurrentFromDrop(` | 1 | 2 | 2 | 5 |
| `branchDrop(` | 1 | 4 | 2 | 7 |
| `branchAffine(` | 1 | 2 | 2 | 5 |
| `branchCurrentNumerator(` | 1 | 1 | 1 | 3 |
| `branchCurrentOut(` | 1 | 1 | 1 | 3 |

`BranchAffine` branches on kind before touching the cell spans, uses designated
field returns, and supplies a positive-zero resistor source. The gate pins each
complete owner body, argument order, sparse runtime numerator rejection,
relaxation's assertion policy, both sparse reconstruction slices, and broad
absence of raw branch-current assignments/divisions. The dead `<cstring>` and
three dead `using` declarations are absent. The compensated-sum explanation is
checked in a raw-comment view immediately above its strict-FP pragmas; it no
longer claims a nonexistent sparse consumer. The header's MC-3 contract now
names its PLAN §3.4 affine-algebra ownership.

The first fast-math build exposed one additional warning: a local `finite`
lambda was used only inside `assert` and became unused in Release. The lambda
now lives directly inside the assertion, preserving the Debug check without a
Release declaration. Subsequent builds emit no such warning. The redundant
`-Ofast` deprecation remains visible and deliberately owned by B1.

### Exact three-configuration result

No recorded constant moved:

| Configuration | Complete PackSolver | Structural aggregate |
|---|---:|---:|
| Debug/ThinLTO | 949 assertions / 30 cases | 1/1 |
| fast-math Release, IPO off | 949 / 30 | 1/1 |
| Release/ThinLTO CUDA tree, host C++ | 949 / 30 | 1/1 |

Trace B therefore covers the resistor's unified `drop - +0.0` path, and Trace D
retains all four exact callback/publication records without reblessing. A final
targeted build in each tree reported `ninja: no work to do`.

### Adversarial gate checks

At clean `8630207`, changing only the centralized cell-current division to
reciprocal multiplication made the structural owner gate red and made every
Trace D path fail: 18/26 assertions passed, with the eight hash comparisons
failing at indices 0–3. The mutated hashes were the already derived
bit-separating values:

```text
index 0  67674e77f46ba2b8 / c61c787dcc8a9396
index 1  ea9ba4b556eddad4 / 80602aa74e05476d
index 2  2f3c710b8f14bdba / 50827d75714062b9
index 3  0664dfa152e6be22 / 5db2a6502d3dc5d2
```

Four structural-only mutations also turned red: making a scalar helper input
by value, reversing the branch-drop subtraction, swapping the same-typed OCV
and resistance spans at one call, and restoring a false sparse-consumer Kahan
comment. During adversarial review, two initially vacuous gate designs were
fixed before commit: comment-free compact input could not police the comment,
and semicolon-rich bodies passed through variadic CMake lists did not prove
contiguity. The final gate reads comments separately and uses fixed-arity exact
block counts.

Every mutation was reversed to exact hashes:

```text
PackSolverInternal.hpp      BE8E44AFE23F81A1B174B9113FEE77FB6EF3638AD7660FCACA87BEACE0322187
PackSolver.cpp              8449CA29E593C585E2917A6F1F90E7B1A136EB8121AA25D74F30217CD1F44B3C
PackSolverIterative.cpp     2F9D343527B4A8989FB54D4E406E392A6775FDE20C42DC39D5D987A45F431662
p9c_architecture.cmake      20CEC02299F8437C5C98629D6596CCA97210A52BFC0F56C96827AF4F92F9880D
```

The restored Debug/Release/host-C++ trees each re-passed 949/30.
`clang-format --dry-run --Werror` and `git diff --check` pass. No full-suite,
timing, sanitizer, coverage, hosted-CI, installed-package, CUDA-device, Linux,
or macOS claim is made for P1.

## P2 — seven-role Mode-C relaxation storage

### Pre-source structural RED

The P2 gate was added and run before any production storage changed. Against
P1 production it failed at the intended first owner:

```text
MQ.2 P2 exact relaxation scratch owner: expected 1 occurrences of
structRelaxationScratch{...};RelaxationScratchrelaxation_{};
found 0
```

The gate bounds the owner inside `PackSolver`'s private slice, bounds
`solveRelaxation` through its unique return tail, scans every production
`.hpp`/`.cpp`/`.cu` under `src/core` for the four retired names, and checks
allocation/proof/publication order numerically before accepting the contiguous
no-throw member-publication block.

### Implementation and structural ownership

Commit `bf185dd` applies `relaxation-target-triple-role`. One private
`RelaxationScratch` now owns:

| Role | Indexed uses in `solveRelaxation` |
|---|---:|
| `diagonal` | 4 |
| `diagonal_compensation` | 2 |
| `rhs` | 6 |
| `rhs_compensation` | 4 |
| `target` | 4 |
| `residual` | 5 |
| `residual_compensation` | 4 |

`configure` constructs all seven node-sized vectors in one local aggregate
before member publication, proves the aggregate nothrow-move-assignable, and
publishes it once. This adds three cold vector allocations and objects. The
hot solve still performs exactly six fills: four coefficient/RHS buffers at
entry and two KCL-residual buffers immediately before assembly. `target` is
assigned on every node and is never filled; no vector is constructed or
resized in `solveRelaxation`.

The structural gate pins every same-node sum/compensation pair, both terminal
signs, the complete target loop, the two exact fill blocks, and the final
residual read. It rejects alternate fill/copy/allocation APIs rather than
counting only `std::fill`. The private two-space PackSolver anchor changes
85→83 exactly because four members became one nested-struct declaration plus
one aggregate member; its namespace anchor remains 16.

### Exact three-configuration result

No recorded constant moved:

| Configuration | Complete PackSolver | Structural aggregate |
|---|---:|---:|
| Debug/ThinLTO | 949 assertions / 30 cases | 1/1 |
| fast-math Release, IPO off | 949 / 30 | 1/1 |
| Release/ThinLTO CUDA tree, host C++ | 949 / 30 | 1/1 |

The existing Debug hot/cold boundaries also remain green:

```text
P4-G3 100k-cell Mode C is sub-GB and allocation-free: 7 assertions / 1 case
PackSolver late allocation failure preserves prior configuration: 33 / 1
```

### Adversarial gate checks

Six distinct mutations turned the registered boundary red:

1. Reusing `target` as diagonal compensation made the indexed-role gate find
   0/2 diagonal-compensation uses.
2. Removing only the residual reset made the fill gate find 5/6 and made the
   repeated relaxation record fail 2/71 assertions. Only Trace-D index 2
   changed, from `3168649734a681a2 / 5e3b17711acbdbf4` to
   `f44ba3fce50c5a12 / e5500e469212d7b6`.
3. Swapping positive/negative diagonal-compensation indices preserved all
   aggregate counts but failed the exact positive-node association.
4. Moving the complete scratch candidate and nothrow proof below the intact
   publication block failed: `all scratch allocation and its no-throw proof
   must precede member publication`.
5. Adding `std::ranges::fill(relaxation_.target, 0.0)` preserved the six
   `std::fill` calls but failed the alternate-mutation prohibition.
6. Restoring `relaxation_diagonal_` only in an unrelated production comment
   in `PackSolverValidation.cpp` failed the production-wide legacy-name scan.

Every mutation was reversed to exact hashes:

```text
PackSolver.hpp              D41F0B7E2E823B9729E7517069AA37925724CBD43D1A790652B595093E317586
PackSolver.cpp              B532DFDFAAF90DEE495AA482C95A06D140BCC679FDB3A8C9A6C2CA124B77E8B6
PackSolverIterative.cpp     1A8E3B40EC056D6D62B50333E858711AF4A05954D07B29F78AED015201446CDD
p9c_architecture.cmake      5DBC885D9D76BB379630D51D28F5E582D1F5D527A28EFF19BFA0D9E54C273E7D
p9c5_public_surface.cmake   EAFE6A20E500066FE93D2C71FE9A8820E050B1CA68D4AB98601126280BA60DBD
```

The restored three configurations each re-passed 949/30 and 1/1.
`clang-format --dry-run --Werror` and `git diff --check` pass. No full-suite,
timing, sanitizer, coverage, hosted-CI, installed-package, CUDA-device, Linux,
or macOS claim is made for P2.

## S1 — PackStepper ownership and executable substeps contract

### Reconciliation and implementation

The seven-file WIP named by PLAN was first audited without editing either
frozen unit oracle. The refactor is behavior-preserving by direct
correspondence:

- the three old prefix searches become one cold
  `detail::firstOccurrence` owner with the same prefix, equality, and
  short-circuit semantics;
- the four state-copy loops become one state-to-contiguous
  `gatherStates` owner and one contiguous-to-state `scatterStates` owner,
  retaining their original offsets, copy order, validation, invalidation,
  solver snapshot, diagnostics, and heat transactions; and
- `stepImpl` retains one electrical solve and one thermal assembly before
  `substeps` full-`dt` advances. The public spelling remains `substeps`
  because its documented meaning is the established multirate window, not
  `dt` subdivision.

Independent adversarial review found no S1 blocker. Commit `62f1587` applies
`packstepper-gather-scatter`, `unique-in-prefix-three-spellings`, and
`substeps-name-contradicts-code`. The compact source census is exactly the
registered one: `firstOccurrence(` owner/PackSolver/PackStepper = 1/2/1,
`gatherStates(`/`scatterStates(`/`std::memcpy(` = 3/3/2, and the
PackStepper private-member anchor is 35 while its namespace/API anchor stays
4. No numerical oracle or recorded hash was changed or reblessed.

### Positive three-configuration and full-suite gates

Before the implementation commit and again after all mutations were
reversed, the focused gate passed without a count delta:

| Configuration | PackStepper | PackSolver | P2-G1 allocation | Structural aggregate |
|---|---:|---:|---:|---:|
| Debug/ThinLTO | 230 assertions / 9 cases | 949 / 30 | 14 / 2 | 1 / 1 |
| fast-math Release, IPO off | 230 / 9 | 949 / 30 | 14 / 2 | 1 / 1 |
| Release/ThinLTO CUDA tree, host C++ | 230 / 9 | 949 / 30 | 14 / 2 | 1 / 1 |

The amended 16-assertion S1 case therefore proves the registered contract:
one `substeps=4`, `dt=0.125 s` call publishes one current within
`1e-10 A` of `8 A`; four independent direct `EulerLegacy::step` calls use
that same frozen current density; both elapsed-time rows are exactly
`0.5 s`; and the complete arenas are byte-identical.

As an additive gate beyond the focused registration, restored commit
`62f1587` passed the unfiltered CTest suite **58/58** in all three retained
trees. The CUDA tree was built under the complete VS 18 x64 developer
environment and its suite included the real `unit_test_core_CudaSpmBatch`,
not the disabled shim. This makes no timing claim: reported durations are
discarded under D-27.

### Adversarial mutation sensitivity

Nine controlled Debug mutations ran only after the implementation was a
clean committed boundary. Every one turned its registered gate red:

| Mutation | Numerical result | Structural result |
|---|---|---|
| divide both production inner steps by `substeps` | S1 stopped at 12/13 reached assertions; elapsed time was `0.125`, expected `0.5` | full-`dt` token count 0/2 |
| divide only the independent direct-Euler comparator by `substeps` | S1 stopped at 13/14; comparator elapsed time was `0.125`, expected `0.5` | not applicable; this is an oracle mutation |
| make `firstOccurrence` always true | duplicate archetype 17/18, aliased Thevenin 1/2, aliased PackStepper 2/3 | exact owner absent |
| remove only the PackStepper prefix guard | downstream PackSolver validation may mask it, as registered | PackStepper consumer count 0/1 |
| reverse gather direction | two-batch P2-G1 stopped at 6/7 | exact gather owner absent |
| reverse scatter direction | P2-G1 stopped at 14/15; repeated current changed from `{8.93033154179052069, 11.06966845820940115}` to `{8.93475168745069226, 11.06524831254933616}` | exact scatter owner absent |
| make public `checkpoint` gather into the internal buffer | P2-G1 stopped at 13/14 | public gather call absent |
| add an unreachable electrical solve inside the substep loop | not run; structural placement is the registered arbiter | solve count 2/1 |
| add an unreachable thermal assembly inside the substep loop | not run; structural placement is the registered arbiter | assembly count 2/1 |

The production-`dt`, gather, scatter, and public-buffer mutations also make
the exact helper/body gate red, so neither numerical nor structural
sensitivity is being inferred from the other.

After every inverse patch, the touched file matched its pre-mutation hash.
The final restoration anchors were:

```text
473E690A73D3575B683FB424223236A2B51FDE05C3C615E28C98FB8C16677A07  src/core/PackTopologyInternal.hpp
A870C40AA4043E82AA30A0A3CA1969C38DD68AB0E470E9FEFC2A2B536D31DD22  src/core/PackSolver.cpp
DBAB230FE124ECE13967A21B56861CF1996EEB3D2BA7E3FCBC45D2EB6EC9B553  src/core/PackStepper.hpp
A9B72A22A7E723C3A17AAC7FDC8B386313EFF396C8472F44CE7292E6F8F916B3  src/core/PackStepper.cpp
2D2A678801A4DA84368116DC2FDDF19F2A389308F5D49B537F34A8FAD4C5159A  tests/unit/core_PackStepper_test.cpp
A2B8971BC4D803C6307BDAD75A4CEC0BA37F6D49D3234B4328C6C69911B79C0E  tests/structural/p9c_architecture.cmake
```

At the clean restoration observation, `git diff --exit-code HEAD --` over
all touched files returned zero, `git status --porcelain` was empty,
`clang-format --dry-run --Werror` passed on the four changed C++ files, and
`git diff --check` passed.

### Adversarial findings carried forward

The review also found a separate, pre-existing high-severity ownership bug:
implicit moves leave `PackSolver::configured_`/`has_solution_` and
`PackStepper::configured_` true in objects whose vectors and workspace were
moved away. A moved-from `PackStepper::solveElectrical` can consequently
reach a null moved-from solver workspace. No prior PLAN, ledger, survivor,
or refutation row owns it. It is not misreported as an S1 regression or
silently fixed without a failing test; the next supplemental boundary will
preregister explicit no-throw move invariants, moved-from rejection, and
destination-continuity oracles before editing production code.

Three test-strength observations are also retained rather than hidden by the
green gate: add an independent checkpoint-layout oracle (the existing
round-trip is symmetric), a numerical heterogeneous-thermal frozen-assembly
oracle (the present S1 case is isothermal), and non-adjacent `[a,b,a]`
prefix-duplicate cases (the exact structural owner currently supplies the
stronger proof). A fresh full Debug rebuild additionally exposed the
pre-existing `EIGEN_STRONG_INLINE` redefinition warning; B1's later
build-diagnostics batch must disposition it alongside the already registered
`-Ofast` and legacy narrowing warnings.

No sanitizer, coverage, hosted-CI, installed-package, Linux, macOS, or
performance claim is made for S1. MQ.2 remains open; the three S1 survivor
IDs are APPLIED, taking the running census from 21 to 24 APPLIED and leaving
39 of the original 71 decisions pending before newly discovered supplemental
work.

## S1.1 — moved-owner and independent-oracle hardening

### Frozen pre-production RED

The complete S1.1 contract was registered at `b622272` before a test or
production edit. Commit `a747cc1` then froze four independent move cases, the
independent checkpoint and analytic frozen-thermal cases, non-adjacent prefix
fixtures, compile-time move traits, and the future public-surface gate.
The first build found only a test-syntax defect: five grouped Catch2
conjunctions lacked the extra parentheses required by Catch's decomposer.
Test-only commit `b468cd7` corrected that grouping without changing an
assertion, expected value, or production file.

Production remained byte-identical to the restored S1 boundary:

```text
D41F0B7E2E823B9729E7517069AA37925724CBD43D1A790652B595093E317586  src/core/PackSolver.hpp
A870C40AA4043E82AA30A0A3CA1969C38DD68AB0E470E9FEFC2A2B536D31DD22  src/core/PackSolver.cpp
DBAB230FE124ECE13967A21B56861CF1996EEB3D2BA7E3FCBC45D2EB6EC9B553  src/core/PackStepper.hpp
A9B72A22A7E723C3A17AAC7FDC8B386313EFF396C8472F44CE7292E6F8F916B3  src/core/PackStepper.cpp
```

`git diff --exit-code 80fba19 -- src/core` returned zero. The corrected
Debug test targets then built, and the tag-filtered old-production run was
red without executing a null dereference:

```text
PackSolver [S1.1]
test cases:  2 |  0 passed | 2 failed
assertions: 19 | 17 passed | 2 failed

move-constructed ... line 566:
REQUIRE(source.setRelaxationGain(0.5) == Status::Invalid_parameters)
with expansion: 0 == 'x'

move-assigned ... line 627:
REQUIRE(source.setRelaxationGain(0.5) == Status::Invalid_parameters)
with expansion: 0 == 'x'
```

Here `0` is `Status::Success`; both constructor and assignment therefore
reproduce the copied-`configured_` bug independently. Their fatal assertions
precede the moved-from sparse `solve` call that could reach the null
workspace implementation.

```text
PackStepper [S1.1]
test cases:  4 |  2 passed | 2 failed
assertions: 80 | 70 passed | 10 failed
```

The independent checkpoint-layout and frozen-thermal cases are the two
passing cases. Both move forms fail only safe observations before their
fatal checkpoint discriminator: solution and diagnostics remain stale,
workspace age/symbolic-factorisation counters remain `1`, and
`checkpoint(empty)` returns `Status::Success` rather than
`Invalid_parameters`. The fatal assertion again prevents the later
`solveElectrical`/`step` calls from reaching the null moved-from solver.

The direct prefix compile/runtime boundary independently remains green:

```text
unit_test_core_PackTopology
All tests passed (148 assertions in 9 test cases)
```

The aggregate structural gate is independently red before source edits:

```text
9C-5 R4: PackSolver.hpp declares 83 member-level entities, the pinned
surface has 88.
```

This is the intended future rule-of-five boundary, not a reblessed surface.
No Release/CUDA/full-suite, mutation, sanitizer, coverage, timing, or
performance claim is made by this pre-fix run.

### S1.1 implementation and restored acceptance

Commit `574cfaf` implements the frozen contract without changing any S1.1
oracle. `SolverWorkspace` now transfers its implementation and factorisation
storage while exchanging its validity, age, and factorisation counters to
their default values. `PackSolver` and `PackStepper` now expose explicit
non-copyable, no-throw move boundaries, transfer every owner in declaration
order, guard self-assignment, and leave the source observably unconfigured.
Member-level compile-time proofs cover both move operations and the default
construction of the `PackSolution`/`PackSolveDiagnostics` reset temporaries.

The `PackSolver` ownership implementation was split into the new
`src/core/PackSolverOwnership.cpp` rather than growing the sparse solver past
MC-1's 700-line boundary. Physical line counts at the accepted commit are
642 for `PackSolver.cpp`, 100 for `PackSolverOwnership.cpp`, and 401 for
`PackStepper.cpp`. `cmake/SlideCoreTarget.cmake` owns the new source exactly
once, and the aggregate structural gate pins that source wiring, all public
special-member declarations, all complete move bodies, the source-reset
exchanges, and the no-throw reset-temporary proof. The frozen P9C5 public
member anchors pass at PackSolver 88 and PackStepper 40; no surface count was
reblessed after implementation.

The restored focused gate passed at the exact registered floors:

| Configuration | PackTopology | PackSolver | PackStepper | P2-G1 allocation | Structural aggregate |
|---|---:|---:|---:|---:|---:|
| Debug/ThinLTO | 148 assertions / 9 cases | 984 / 32 | 333 / 13 | 14 / 2 | 1 / 1 |
| fast-math Release, IPO off | 148 / 9 | 984 / 32 | 333 / 13 | 14 / 2 | 1 / 1 |
| Release/ThinLTO CUDA tree, host C++ | 148 / 9 | 984 / 32 | 333 / 13 | 14 / 2 | 1 / 1 |

The exact commands, all from the repository root, were:

```powershell
# Debug/ThinLTO
cmake --build build-mq1-debug --parallel 2
ctest --test-dir build-mq1-debug --verbose -R "PackSolver|PackStepper|PackTopology|P2G1_allocation|structural_test_core_9C2AgeingKernel"
ctest --test-dir build-mq1-debug --output-on-failure -j1

# fast-math Release, IPO off
cmake --build build-mq1-release --parallel 2
ctest --test-dir build-mq1-release -V --no-tests=error -j1 -R '^(unit_test_core_PackSolver|unit_test_core_PackStepper|unit_test_core_PackTopology|unit_test_core_P2G1_allocation|structural_test_core_9C2AgeingKernel)$'
ctest --test-dir build-mq1-release --output-on-failure -j1

# Release/ThinLTO CUDA tree, host C++; each command used this full VS wrapper
cmd.exe /d /s /c 'call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul && cmake --build build-mq1-cuda-vsenv --config Release'
cmd.exe /d /s /c 'call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul && ctest --test-dir build-mq1-cuda-vsenv -R "PackSolver|PackStepper|PackTopology|P2G1_allocation|structural_test_core_9C2AgeingKernel" -V -j1'
cmd.exe /d /s /c 'call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul && ctest --test-dir build-mq1-cuda-vsenv --output-on-failure -j1'
```

The unfiltered CTest suite then passed **58/58** in every tree. The CUDA lane
ran under the complete VS 18 Community `vcvars64.bat` environment, retained
Clang as the host C++ compiler, had `SLIDE_WITH_CUDA=ON`, and passed the real
`unit_test_core_CudaSpmBatch`. The Debug full build completed 61/61 steps
and that incremental build output contained no warning or error. This does
not close B1's separately registered pre-existing
`EIGEN_STRONG_INLINE` warning. No wall-clock or performance claim is made.

Two command-wrapper incidents are recorded so they cannot be mistaken for
product evidence. In Release and CUDA, an initial five-second tool wrapper
expired while its child build continued. A concurrent Release retry reported
`ninja: error: failed recompaction: Permission denied`; process inspection
showed the first build still owned that tree. After each child exited, the
identical authoritative build command returned exit zero with
`ninja: no work to do.`, the focused and full gates above passed, and both
tracked and staged diffs remained empty. The transient Release regeneration
also reported optional JNI, Java, and SWIG packages absent. None of these
preliminary wrapper observations is counted as a compiler or test result.

### S1.1 adversarial mutation sensitivity

Eight independent Debug mutations ran only after `574cfaf` was a clean
committed boundary:

| Registered mutation | Red evidence |
|---|---|
| copy rather than exchange move-constructed `PackSolver::configured_` | PackSolver `[S1.1]`: 27/28 reached assertions; moved-source `setRelaxationGain(0.5)` returned `Success` |
| copy rather than exchange move-constructed `SolverWorkspace::valid_` | PackSolver `[S1.1]`: 34/35; `isDefaultWorkspace(source.workspace())` was false |
| copy rather than exchange move-constructed `PackStepper::configured_` | PackStepper `[S1.1]`: 92/93; moved-source empty `checkpoint` returned `Success` |
| remove the `PackStepper` assignment self guard | exact structural move-body owner found 0/1 instead of 1/1 |
| omit the move-constructed destination `SolverWorkspace` owner | exact structural PackSolver move-body owner found 0/1 instead of 1/1 |
| pair a reverse batch-segment order in both checkpoint gather and scatter | independent layout case: 8/11; serialized bytes and both independent restore slices failed |
| reassemble and republish thermal coupling inside every Euler substep | analytic thermal case: 6/9; heat changed and each final temperature missed by `0.0020062352608079 K`, versus `1e-10 K` |
| accept the last non-adjacent duplicate and out-of-range prefix index | PackTopology compilation stopped at the constexpr `{7,11,7}` `static_assert` |

The self-guard mutation is adjudicated by its registered structural gate. An
over-broad exploratory runtime invocation was also allowed to reach the
deliberately corrupted self-move and exited with Windows access-violation
status `-1073741819` (`0xC0000005`) after 73/74 reached assertions; Catch2
rendered this as `SIGSEGV - Segmentation violation signal`. That termination
is not used as evidence and was not repeated. The omitted-workspace mutation's
behavioral command produced no verdict before its 64-second wrapper timeout
and left its test child alive. Its executable path was verified as this
workspace's mutation binary and the process was stopped; only the safe,
immediate structural failure is counted.

Every mutation was reversed with an inverse patch before the next one. The
final restoration hashes are:

```text
470B50C0E6BFBCDD0970C615236AA5B95E0C4A99957E10B3B1185BE3DAFA1F7E  src/core/PackSolver.cpp
9A44E2BBE3BE1AE413212761BBA0511C7541E13E9F1B0146971A6BB812EF63F5  src/core/PackSolverOwnership.cpp
2E228259545DC875E73D3562C088A456EDD04251621C2B8567DEA9659088D01E  src/core/PackStepper.cpp
473E690A73D3575B683FB424223236A2B51FDE05C3C615E28C98FB8C16677A07  src/core/PackTopologyInternal.hpp
3E970FBB935BD6BD715F3F84E6ADB6999A641B522BBB4154ECD2B21B77B1F724  tests/structural/p9c_architecture.cmake
```

Immediately before these evidence-only edits,
`git diff --exit-code HEAD` and the final per-lane clean checks all returned
zero. The original 71-decision census is unchanged by this supplemental
boundary: 24 APPLIED, 0 REFUTED, 8 named deferrals, and 39 pending.

The exact post-baseline registry separately records
`moved-owner-validity-after-defaulted-move`,
`checkpoint-roundtrip-oracle-symmetry-gap`,
`thermal-assembly-placement-oracle-gap`, and
`prefix-duplicate-adjacency-oracle-gap` as APPLIED by S1.1.
`eigen-strong-inline-redefinition-warning` is DEFERRED with owner B1. The
combined working census is therefore 76 stable IDs: 28 APPLIED, 0 REFUTED,
9 named deferrals, and 39 pending. The closeout gate compares the original
and post-baseline ID sets independently. MQ.2 remains open; S2 source-step
rollback is next.

## S2 — source-step rollback

### Frozen old-production RED

The exact 10-assertion / 1-case count amendment was committed at `181a41a`
before the first S2 test edit or run. Commit `4a9d29a` then froze only the
scripted rollback oracle. Its first run reached the registered 9/10 boundary;
test-only commit `d94a0ba` added a zero-assertion Catch2 `CAPTURE` so the
observed residual would be printed rather than supported only by derivation.
No expected value or assertion changed.

Production remained exactly at the accepted S1.1 hash:

```text
470B50C0E6BFBCDD0970C615236AA5B95E0C4A99957E10B3B1185BE3DAFA1F7E  src/core/PackSolver.cpp
5B33D2A0A7CBB43D1559E45E808947F4E306CE6D6496D42D6BE19F4B4E70B713  tests/unit/core_PackSolver_test.cpp
```

The exact Debug command was:

```powershell
cmake --build build-mq1-debug --parallel 2 --target unit_test_core_PackSolver
.\build-mq1-debug\bin\Debug\unit_test_core_PackSolver.exe 'failed source stepping restores caller-attempt diagnostics'
```

It failed only the nonfatal diagnostics comparison:

```text
solver.diagnostics().residual_norm := 2.5
test cases:  1 | 1 failed
assertions: 10 | 9 passed | 1 failed
```

The initial exact solution, six-callback boundary, restored solution, final
one-iteration cold solve, and seventh callback all passed. Thus unchanged
production already restores the solution fields and false warm-state flag;
only `diagnostics_` leaks the abandoned hidden step's residual `2.5` instead
of describing the caller-requested direct attempt at residual `4`.

No Release, CUDA, full-suite, mutation, sanitizer, coverage, timing, or
performance claim is made by this failing-test-first boundary.
