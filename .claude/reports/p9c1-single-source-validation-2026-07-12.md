# M0.5 / 9C-1 single-source physics validation (2026-07-12)

## Result

**PASS.** The modal update, SPM observable algebra, curve interpolation, and
spectral-diffusion scalar expressions now have one production definition in
`src/core/SpmScalarKernels.hpp`. CPU `double`, CUDA `double`, and forward-mode
`Dual` paths instantiate those definitions without changing their established
operation order. The independent `SpectralDiffusionLegacy.hpp` oracle remains a
deliberately raw formula copy and is structurally forbidden from including the
shared header.

## Architecture delivered

- `SpmPipeline.hpp` uses the shared modal update, Arrhenius, activated
  diffusivity, molar-flux, and diffusion-rate definitions.
- `SpmObservables.hpp` and `CompiledCurve.hpp` use the shared scalar observable
  and interpolation algebra.
- `CudaSpmRuntime.cu` uses the same host/device-capable scalar expressions.
- `ForwardSensitivity.cpp` instantiates the same expressions for `Dual`; its
  terminal-voltage observation now also differentiates the entropic OCV term.
- `SpectralDiffusion.hpp`, the remaining production consumer found by the
  adversarial audit, uses the same Arrhenius, diffusivity, flux, and RHS
  definitions.
- `tests/structural/pc10_single_source.cmake` scans every production consumer,
  rejects raw duplicate expressions, checks the expected call topology, and
  protects the independent legacy oracle from common-mode reuse.

The modal and diffusion hot loops use statement/expression forms supplied by
the shared header. This is deliberate: moving the same scalar expression behind
an ordinary inline function changed Clang fast-math loop IR and recorded bits;
`always_inline` did not restore them. Inspection of the optimized
`SpectralDiffusion` binary after the final form still found packed AVX
`vmulpd`/`vaddpd` instructions.

## Recorded-bit gate and provenance

`SLIDE_ENABLE_RECORDED_SCALAR_BITS` is OFF by default. CMake restricts opt-in to
x64 Windows with Clang 21.1.8. Successful reproduction additionally depends on
a compatible CRT/ISA because Release uses `-march=native`; CUDA also depends on
the recorded CUDA 13.0 / sm_89 RTX 4000 Ada environment, which CMake does not
hardware-guard. These hashes are provenance checks, not portable mathematical
oracles.

| Trace | Values | Configurations checked | Result |
|---|---:|---|---|
| CPU Chen2020, four heterogeneous lanes, NCH=5/8/12 | 944 / 1,136 / 1,392 | Debug/ThinLTO, Release fast-math, Release fast-math/ThinLTO | all six fingerprints per configuration match the frozen pre-refactor rows |
| Dual P7-G2, 61 times/voltages and 61x10 tangents | 732 | Debug/ThinLTO, Release build with strict-FP producer, Release/ThinLTO build with strict-FP producer | both fingerprints match in all three configurations |
| CUDA Chen2020, four heterogeneous NCH=12 lanes | 1,388 | CUDA 13.0, sm_89 RTX 4000 Ada | both fingerprints match |

The complete CPU/Dual/CUDA hash table and capture history are in
`p9c1-pre-refactor-scalar-fixtures-2026-07-12.md`. The missed
`SpectralDiffusion` production consumer was frozen immediately before its own
migration as a whole 600-step heterogeneous trace:

| Configuration | Values | FNV-1a | Independent mix | Result |
|---|---:|---:|---:|---|
| Debug/ThinLTO (`-O0`) | 48,080 | `d2601e8e7d13249b` | `5eabfd035f26594c` | exact |
| Release/fast-math | 48,080 | `8f7f609c10bb3d1a` | `8332970bbbba9a67` | exact |
| Release/fast-math/ThinLTO | 48,080 | `785b850afc76b768` | `3a0476924499700f` | exact |

On WSL Ubuntu 24.04 with Clang 18.1.3, the recorded hashes were disabled and
the portable `SpectralDiffusion`, `ExponentialModal`, and `ForwardSensitivity`
gates passed. A detached historical comparison showed exact Dual and CPU NCH=8
hashes but different CPU NCH=5/12 hashes across the Linux code-generation
boundary; those post-refactor Linux values were correctly not adopted as new
baselines.

## Adversarial checks and bugs found

1. **Dual entropic OCV omission.** A temperature difference-of-differences
   regression separated the production tangent from a centred finite-difference
   arbiter (approximately 0.0048961 versus 0.008). Compiling the total entropic
   curve and sharing the cell-OCV expression fixed the intentional off-reference
   output transition while preserving the frozen reference-temperature P7-G2
   trace.
2. **ThinLTO plating overflow misclassification.** Optimisation converted a
   derived-scale overflow path to `Success`. Bit-safe primal classification and
   overflow preflight now validate `n_plating*F`, `n*F`, and
   `(n_plating*F)*plated_density` without overflowing the guard itself. Runtime
   and parameter failures retain their distinct Status taxonomy; exact coverage
   reaches all six active `LithiumPlating.hpp` failure arms.
3. **Negative spectral lane count.** The signed count previously converted to a
   huge allocation before an assertion. Zero and negative counts now throw
   `invalid_argument` before allocation.

Controlled mutations proved the gates distinguish the intended defects:

- a private raw diffusion copy failed PC-10 structurally;
- reversing a shared diffusion input sign failed the independent legacy oracle
  with maximum relative error near 2.0;
- swapping Arrhenius reference/current temperature wiring failed that oracle
  with maximum relative error near 0.5564;
- reversing the terminal-resistance sign failed CPU, Dual, and CUDA recorded
  fixtures;
- restoring a raw modal `1/6` coefficient failed the structural gate.

## Validation matrix

| Lane | Result | Notes |
|---|---:|---|
| Windows Clang Debug/ThinLTO | 54/54 | serial CTest, 21.22 s; ModeC 13.82 s |
| Windows Clang Release/fast-math | 54/54 | serial CTest, 9.49 s; ModeC 5.04 s |
| Windows Clang CUDA/fast-math/ThinLTO | 54/54 | serial CTest, 13.68 s; CUDA batch 2.07 s; ModeC 5.88 s |
| WSL Clang 18 ASan+UBSan | 54/54 | full suite, 163.56 s; ModeC 150.74 s; no sanitizer finding |
| WSL Clang/LLVM 18 exact Status coverage | 54/54 | 53 per-test profiles + one whole-archive anchor; 369 lexical, 339 active, 332 measured + 7 structural exceptions, 30 inactive, zero uncovered/unmapped |
| Windows core-only optional-off | 1/1 | ordinary-C++ consumer precursor |
| Windows nested external consumer optional-off | 1/1 | public-header/CMake consumption |

TSan was not repeated because M0.5 changes no concurrent production code; the
authoritative M0.2 ThreadPool/AsyncRecorder/PackStepper result remains 3/3.
Core-only checks establish the ordinary-C++ precursor but do not claim the
actual WASM/native parity reserved for M10. Hosted CI workflow counts and path
filters were updated for 54 tests and the structural gate, but no hosted run is
claimed before push.

Artifacts:

- `p9c1-pre-refactor-scalar-fixtures-2026-07-12.md`
- `p9c1-status-coverage-refresh-2026-07-12.md`
- `p9c1-status-coverage-refresh-2026-07-12.json`
