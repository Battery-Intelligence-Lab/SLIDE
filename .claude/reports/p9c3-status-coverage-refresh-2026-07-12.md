# P9-G3 Status-failure coverage

Generated: 2026-07-12T19:17:30+00:00

Result: **PASS**

## Prior falsification retained

> The first M0.7 run at `57a442ba6ab4f34ed286054bfc200151ed514fd0` was a genuine **FAIL** and is retained here rather than overwritten: **380 lexical, 350 active, 337 covered, 7 excepted, 6 uncovered, 30 inactive**. All 56 tests passed in 266.06 s and produced 54 profile groups, so this was a Status-arm coverage failure, not a test-session failure.

The failed exact sites were:

- `src/core/CyclerV2.cpp:107:27:Numerical_failure`
- `src/core/CyclerV2.cpp:130:25:Numerical_failure`
- `src/core/CyclerV2.cpp:354:27:Numerical_failure`
- `src/core/CyclerV2.cpp:763:25:Numerical_failure`
- `src/core/ParameterSet.cpp:223:25:Numerical_failure`
- `src/core/ParameterSet.cpp:374:25:Numerical_failure`

The initial 369/339/333+6 preregistration was therefore falsified before commit `93e75d1`: the M0.7 file split added no Status semantics, while earlier allocation-transaction work had added eleven literal sites. The preserved failed artifacts had SHA-256 `400258ec5088c980ec509f2fa4105004a016f7a4dd490d3ae5a540d60f087ff9` (JSON) and `ace349fe8645f1f3666ff3fc59ef21de023661982ac25b4c1c5b62b48c7a3200` (Markdown).

## Fresh validation session

- Commit: `93e75d123b93e8fe737379e45b3de6bfc1888745`
- Tests: 56 / 56 passed in 256.90 s
- Coverage evidence: 54 fresh profile files in 54 unique test groups; 54 test binaries plus one anchor
- Compile-command SHA-256: `22fb1680a9994eaf0cf2c2f11da524fa0964c2f6a6608a724f6ee01cc7b35ea4`
- Target-manifest SHA-256: `6fd083b573f2010a433c49a1adae094c5f84290f659152f5c49efc2007a3c2ae`
- Source identities: CyclerV2 `e7f98a0d92db28a5cd920c5dae6544672f9a34e17e6a12feee213880573927a0`; ParameterSet `de10d0f2d40eb9a78a49a1eda9e28b19dcaaf1f8795785de6507d52fe1fd0167`; Experiment test `7960ace96117d6e09cb747e4df284da4042e99e6c2fdee782be7766173528b13`

- Lexical failure arms: 369
- Active optional-off arms: 339
- Measured covered arms: 333
- Structural exceptions: 6 / 10
- Inactive optional branches: 30
- Uncovered/unmapped active arms: 0

| File | Active | Covered | Excepted | Inactive |
|---|---:|---:|---:|---:|
| `src/core/AsyncRecorder.cpp` | 33 | 33 | 0 | 3 |
| `src/core/BpxParameterReader.cpp` | 14 | 14 | 0 | 0 |
| `src/core/CompiledCurve.hpp` | 15 | 15 | 0 | 0 |
| `src/core/CudaSpmBatch.cpp` | 9 | 9 | 0 | 19 |
| `src/core/CyclerV2.cpp` | 25 | 25 | 0 | 0 |
| `src/core/EulerLegacy.hpp` | 3 | 3 | 0 | 0 |
| `src/core/Experiment.cpp` | 10 | 10 | 0 | 0 |
| `src/core/ExponentialModal.hpp` | 5 | 5 | 0 | 0 |
| `src/core/ForwardSensitivity.cpp` | 9 | 9 | 0 | 0 |
| `src/core/Lam.hpp` | 8 | 8 | 0 | 0 |
| `src/core/LithiumPlating.hpp` | 6 | 6 | 0 | 0 |
| `src/core/NetlistCsv.cpp` | 20 | 20 | 0 | 0 |
| `src/core/PackSolver.cpp` | 48 | 43 | 5 | 0 |
| `src/core/PackSolverValidation.cpp` | 2 | 2 | 0 | 0 |
| `src/core/PackStepper.cpp` | 7 | 7 | 0 | 0 |
| `src/core/PackTopology.cpp` | 34 | 34 | 0 | 0 |
| `src/core/ParameterSet.cpp` | 8 | 8 | 0 | 0 |
| `src/core/Recorder.cpp` | 24 | 24 | 0 | 8 |
| `src/core/Sei.hpp` | 6 | 6 | 0 | 0 |
| `src/core/Simulation.cpp` | 4 | 4 | 0 | 0 |
| `src/core/SpectralModel.hpp` | 9 | 8 | 1 | 0 |
| `src/core/SpmFactory.cpp` | 13 | 13 | 0 | 0 |
| `src/core/SpmObservables.hpp` | 1 | 1 | 0 | 0 |
| `src/core/SpmPipeline.hpp` | 7 | 7 | 0 | 0 |
| `src/core/SpmStress.hpp` | 3 | 3 | 0 | 0 |
| `src/core/SurfaceCrack.hpp` | 6 | 6 | 0 | 0 |
| `src/core/ThermalLumped.hpp` | 3 | 3 | 0 | 0 |
| `src/core/ThreadPool.cpp` | 4 | 4 | 0 | 0 |
| `src/core/ThreadPool.hpp` | 3 | 3 | 0 | 0 |

## Structural exceptions

- `src/core/PackSolver.cpp:351:27:Invalid_parameters` — **defensive-only**: All dimensions come from already-materialized validated vectors and node_count is bounded by INT_MAX. Deterministic allocation exhaustion reaches bad_alloc; this catch is retained for standard-library or allocator implementations that report an exhausted representability limit with length_error.
- `src/core/PackSolver.cpp:588:27:Numerical_failure` — **defensive-only**: After successful SparseLU factorization and a finite right-hand side, a solve-info failure requires an Eigen backend failure or internal corruption. Public rounded-singular inputs exercise the preceding factorization-failure boundary instead.
- `src/core/PackSolver.cpp:652:29:Invalid_states` — **defensive-only**: The damped state is a convex interpolation between finite pre-correction and full-correction affine states, so its node and cell values remain finite. The combined runtime guard is retained to preserve atomicity against staged floating-point rounding anomalies.
- `src/core/PackSolver.cpp:796:29:Numerical_failure` — **platform-specific**: Every reciprocal conductance is first required to be finite and positive, and connectedness supplies an incident branch for every non-reference node. This guard remains for FP environments that accept a positive subnormal comparison but flush it during compensated diagonal accumulation.
- `src/core/PackSolver.cpp:887:27:Invalid_states` — **defensive-only**: All relaxed branch currents and earlier compensated KCL additions are finite before the final terminal injection. Failure here requires a representability-limit cancellation pattern not constructible by the bounded validated pack witnesses, so the atomicity guard is retained.
- `src/core/SpectralModel.hpp:314:29:Numerical_failure` — **platform-specific**: The registered fixed-order matrices have a validated distinct real spectrum, which implies an independent exact eigenbasis. This guard preserves atomic output if a finite-precision Eigen backend nevertheless returns a numerically rank-deficient approximate basis.

## Gate failures

None.
