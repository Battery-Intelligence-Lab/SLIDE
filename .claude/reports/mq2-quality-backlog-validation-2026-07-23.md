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
