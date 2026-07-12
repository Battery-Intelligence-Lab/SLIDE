# P9-G3 Status-failure coverage

Generated: 2026-07-12T13:51:18+00:00

Result: **PASS**

- Lexical failure arms: 369
- Active optional-off arms: 339
- Measured covered arms: 332
- Structural exceptions: 7 / 10
- Inactive optional branches: 30
- Uncovered/unmapped active arms: 0

| File | Active | Covered | Excepted | Inactive |
|---|---:|---:|---:|---:|
| `src/core/AsyncRecorder.cpp` | 33 | 33 | 0 | 3 |
| `src/core/CompiledCurve.hpp` | 15 | 15 | 0 | 0 |
| `src/core/CudaSpmBatch.cpp` | 9 | 9 | 0 | 19 |
| `src/core/EulerLegacy.hpp` | 3 | 3 | 0 | 0 |
| `src/core/Experiment.cpp` | 33 | 33 | 0 | 0 |
| `src/core/ExponentialModal.hpp` | 5 | 5 | 0 | 0 |
| `src/core/ForwardSensitivity.cpp` | 9 | 9 | 0 | 0 |
| `src/core/Lam.hpp` | 8 | 8 | 0 | 0 |
| `src/core/LithiumPlating.hpp` | 6 | 6 | 0 | 0 |
| `src/core/NetlistCsv.cpp` | 20 | 20 | 0 | 0 |
| `src/core/PackSolver.cpp` | 48 | 43 | 5 | 0 |
| `src/core/PackSolverValidation.cpp` | 2 | 2 | 0 | 0 |
| `src/core/PackStepper.cpp` | 7 | 7 | 0 | 0 |
| `src/core/PackTopology.cpp` | 34 | 34 | 0 | 0 |
| `src/core/ParameterSet.cpp` | 24 | 23 | 1 | 0 |
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
- `src/core/ParameterSet.cpp:185:25:Numerical_failure` — **defensive-only**: ParameterSet::set receives already-constructed strings and value containers and requests no derived length beyond their representable sizes. The catch preserves Status atomicity if a conforming string or map implementation nevertheless reports its resource limit with length_error.
- `src/core/SpectralModel.hpp:314:29:Numerical_failure` — **platform-specific**: The registered fixed-order matrices have a validated distinct real spectrum, which implies an independent exact eigenbasis. This guard preserves atomic output if a finite-precision Eigen backend nevertheless returns a numerically rank-deficient approximate basis.

## Gate failures

None.
