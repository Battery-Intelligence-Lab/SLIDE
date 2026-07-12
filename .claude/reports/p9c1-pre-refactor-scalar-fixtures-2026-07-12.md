# M0.5 / 9C-1 pre-refactor scalar-kernel fixtures (2026-07-12)

The initial CPU/Dual/CUDA exact-bit fixtures use production source state
`b170b14`, before any production physics source was changed. Each trace feeds every `double` representation bit,
in order, through two independent 64-bit recurrences and checks the exact value
count.  They are backend-local baselines: cross-backend equality is neither
assumed nor required.

Exact hash enforcement is deliberately opt-in through
`SLIDE_ENABLE_RECORDED_SCALAR_BITS=ON`; CMake restricts opt-in to x64 Windows
Clang 21.1.8. Successful reproduction additionally depends on a compatible CRT,
ISA, and, for CUDA, GPU because Release uses `-march=native` and hardware is not
configuration-guarded. These are
compiler/CRT/ISA/configuration provenance fixtures rather than portable
mathematical oracles. Other toolchains retain the value-count, analytic,
finite-difference, backend-parity, and structural gates.

## Named traces

| Path | Case | Configuration | Values | FNV-1a | Independent mix |
|---|---|---|---:|---:|---:|
| CPU | Chen2020, 4 heterogeneous lanes, NCH=5, initial observation + dt={1e-5,7.25,19}, signed currents | Clang 21.1.8 Debug/ThinLTO (`-O0`) | 944 | `f461aa18b8f29e9d` | `e99152e8e9b27262` |
| CPU | same, NCH=8 | Clang 21.1.8 Debug/ThinLTO (`-O0`) | 1136 | `81d095bfc708d897` | `de4eb9a14a355bc4` |
| CPU | same, NCH=12 | Clang 21.1.8 Debug/ThinLTO (`-O0`) | 1392 | `983fa84ba68fc867` | `dfc61b3e371e80f1` |
| CPU | same, NCH=5 | Clang 21.1.8 Release/fast-math | 944 | `cd87b2262663a983` | `7966eedf67f341b5` |
| CPU | same, NCH=8 | Clang 21.1.8 Release/fast-math | 1136 | `5c766159af330d3e` | `ad91f8f42bf550ce` |
| CPU | same, NCH=12 | Clang 21.1.8 Release/fast-math | 1392 | `21b0119f4f1ec474` | `1feb3cf7a32e5b35` |
| CPU | same, NCH=5 | Clang 21.1.8 Release/fast-math/ThinLTO | 944 | `00a78a61e6bd84d6` | `196a34e9b1d00094` |
| CPU | same, NCH=8 | Clang 21.1.8 Release/fast-math/ThinLTO | 1136 | `f4d3b8eb337172a6` | `4e06b19d061124f7` |
| CPU | same, NCH=12 | Clang 21.1.8 Release/fast-math/ThinLTO | 1392 | `8dca38c943fc058a` | `ac4b6059db321d83` |
| Dual | Complete P7-G2: 61 times + 61 voltages + 61x10 raw tangents | Clang 21.1.8 Debug/ThinLTO (`-O0`) | 732 | `6dd17b952b7e315a` | `c286fa543130b91b` |
| Dual | same | Clang 21.1.8 Release build (strict-FP producer) | 732 | `06dc78ab74d5e7a4` | `d10596df0d2c2c5a` |
| Dual | same | Clang 21.1.8 Release/ThinLTO build (strict-FP producer) | 732 | `56b17a9abb12f827` | `0ba390f9803cf187` |
| CUDA | Chen2020, 4 heterogeneous lanes, NCH=12, dt={1e-5,7.25,19}, signed currents | CUDA 13.0, sm_89 RTX 4000 Ada, driver 596.72 | 1388 | `a59c34d1685ddb19` | `2e88cc95b54eb713` |

The CPU trace includes all arena values and terminal voltage at every recorded
point.  The CUDA trace includes the initial arena plus every downloaded arena and
device terminal-voltage vector.  The Dual trace covers the complete committed
P7-G2 sensitivity output, not selected samples.

## Orthogonal pre-checks and validation

- Three independent read-only audits mapped the duplicated CPU/CUDA/Dual
  expression trees and agreed on the exact association that must be preserved.
- The existing long-double closed-form modal oracle and the three-step centred-FD
  sensitivity arbiters remain independent; neither uses the fixture fingerprint.
- Debug CPU fixture: 27/27 assertions.
- Release CPU fixture: 27/27 assertions.
- Debug P7-G2 including full FD arbiters: 3863/3863 assertions.
- Release P7-G2 including full FD arbiters: 3863/3863 assertions.
- CUDA backend-local fixture: 14/14 assertions.

Audit discovery: the Dual observation omitted the entropic OCV correction away
from reference temperature. P7-G2 is at the reference temperature, so its frozen
bits remained the refactor gate. Commit `8561eca` corrected the off-reference
defect test-first and documented its intentional output transition separately
from the bit-preserving extraction.

The ThinLTO rows were recovered from detached historical fixture commit `efcb7e1` with
`SLIDE_WITH_CUDA=ON`, `ENABLE_IPO=ON`, and the same Clang/CUDA/sm_89 toolchain.
They exactly match the post-refactor CUDA-on host paths, and therefore extend the
pre-refactor evidence rather than accepting post-refactor values as a baseline.

## Spectral-diffusion addendum

The later consumer audit found `SpectralDiffusion.hpp` outside the initial
CPU/CUDA/Dual map. Before changing that production source, commit `8a08514`
recorded its complete eight-lane heterogeneous trace (initial state plus 600
Euler steps, 48,080 values): Debug/ThinLTO `d2601e8e7d13249b` /
`5eabfd035f26594c`, Release/fast-math `8f7f609c10bb3d1a` /
`8332970bbbba9a67`, and Release/fast-math/ThinLTO `785b850afc76b768` /
`3a0476924499700f`. All three pairs match after migration.

WSL Ubuntu 24.04 / Clang 18.1.3 portable numerical, finite-difference, and
structural gates pass with exact hashes disabled. A historical Linux comparison
retained exact Dual and CPU NCH=8 hashes but changed CPU NCH=5/12 across the
code-generation boundary. Those values were not blessed as post-refactor
baselines; they demonstrate why the Windows capture provenance is explicit.
