# M0.7 / 9C-3 pre-refactor cold-path fixtures (2026-07-12)

Production source baseline: `e360218` (`ParameterSet.cpp` 1,496 lines,
`Experiment.cpp` 1,063 lines).  These fixtures were captured before moving any
production definition.  Recorded bits are enabled only on the pinned x64 Windows
Clang 21.1.8 capture host; portable semantic gates remain active everywhere.

## ParameterSet fixture

The fixture independently frames ordered names, value kinds, provenance, scalar
values, curve spans, and every numeric/curve field of the resulting
`SpmFactoryInput`.  A functional BPX document exercises addition, subtraction,
multiplication, division, unary signs, parentheses, `exp`, `tanh`, `cosh`,
right-associative `x**2**3`, and the precedence of `-x**2`.  Exact diagnostics
and atomic publication are checked for malformed JSON, an unsupported expression,
and state-dependent diffusivity.

| Trace | Values | Debug/ThinLTO FNV / mixed | Release FNV / mixed | Release/ThinLTO FNV / mixed |
|---|---:|---|---|---|
| Chen2020 values | 11,537 | `a5d211a63dd77522` / `8ccfefe72e03094e` | `e5c3dc16e932ad03` / `b4e4cfcce5167740` | same as Release |
| `SpmFactoryInput` | 11,635 | `134b5b4650709d1a` / `1e63d406b3d98297` | `705555174f1b2e23` / `d9357b2ee6af55e4` | same as Release |
| functional BPX values | 4,135 | `04bdc30e3e6ca305` / `e75964cae8725037` | `d3f393daa406ad80` / `f7d10dea65f3c8f5` | same as Release |

Chen2020 has 53 ordered entries and 159 metadata strings (3,491 bytes),
fingerprinted as `8ea717ef0459a7b7` / `745379e447872b1f`.  The functional BPX
case has 35 entries and 105 strings (1,906 bytes), fingerprinted as
`f4a1c008bf74be9a` / `fa654b7fe7b9ab89`.

All three binaries passed 2,358 assertions in 7 cases.  Launching the three
configuration binaries concurrently first reproduced a pre-existing shared-temp
file race: one run failed 3 assertions after another process truncated or removed
the same path.  The test now reserves a process-safe temporary directory
atomically; the identical concurrent command passes 2,358/2,358 in all three
configurations.

## Experiment fixture

The parser fixture frames 12 numeric/integer fields for each of the eight expanded
segments (96 values), while source strings and drive-cycle names are compared
separately in order.  Its hash is configuration-independent:
`1e404a3be7365df5` / `94a28984eed0159b`.

The runner fixture is manually constructed, so parser and runner cannot share a
false positive.  It covers power, rest, a sign-changing drive cycle sampled at
interior points and exact breakpoints, current control, per-segment sample periods,
and a dyadic event root inside a larger step.  It records initial/final arenas,
time, voltage, current, segment ownership, and solution metadata (503 values).

| Configuration | Runner FNV / mixed |
|---|---|
| Debug/ThinLTO | `c779e41caac8339c` / `32dec619a41e0324` |
| Release | `e9616b6ce3131bb6` / `b72cecb529e7c3c0` |
| Release/ThinLTO | `9d97787b3976978b` / `4c7bcd0565164aa2` |

Each targeted parser run passed 81 assertions and the runner passed 13.  Each full
binary passed 351 assertions in 13 cases.  The pre-run arithmetic prediction of
373 assertions was falsified: it incorrectly treated the parser case's full 81
assertions as newly added; only 5 parser assertions were new, so the correct
derivation is the 333 baseline + 5 parser + 13 runner = 351.

## Suite preservation

Debug/ThinLTO, Release, and Release/ThinLTO each enumerate exactly 56 CTests.
No production source was changed while capturing these fixtures.
