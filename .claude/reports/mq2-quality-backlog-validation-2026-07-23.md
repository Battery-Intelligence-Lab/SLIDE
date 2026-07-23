# MQ.2 quality-backlog validation (2026-07-23)

## Status and provenance

This is the accumulating post-run evidence for the preregistration in
`.claude/reports/mq2-quality-backlog-preregistration-2026-07-23.md`.
MQ.2 remains open: this report records completed batch boundaries but does not
yet claim the final 71-row census or the final three-lane gate.

The registered source floor is commit `5481496` (`Freeze MQ.2 pre-edit
fixtures`). The first three completed boundaries are:

| Commit | Batch | Character |
|---|---|---|
| `ec7ac09` | O2 fast-math validity | invalid-state behavior hardening |
| `025108f` | O3 stress/thermal oracles | test-only |
| `ac15536` | F0 factory byte oracle | oracle-only, before factory source edits |
| `b9c4db1` | F1 shared factory constants | exact no-op; homogeneous rewrite deferred |

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
