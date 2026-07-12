# M0.5 / 9C-1 pre-refactor scalar-kernel fixtures (2026-07-12)

These exact-bit fixtures were recorded at commit `b170b14`, before any production
physics source was changed.  Each trace feeds every `double` representation bit,
in order, through two independent 64-bit recurrences and checks the exact value
count.  They are backend-local baselines: cross-backend equality is neither
assumed nor required.

## Named traces

| Path | Case | Configuration | Values | FNV-1a | Independent mix |
|---|---|---|---:|---:|---:|
| CPU | Chen2020, 4 heterogeneous lanes, NCH=5, initial observation + dt={1e-5,7.25,19}, signed currents | Clang 21.1.8 Debug | 944 | `f461aa18b8f29e9d` | `e99152e8e9b27262` |
| CPU | same, NCH=8 | Clang 21.1.8 Debug | 1136 | `81d095bfc708d897` | `de4eb9a14a355bc4` |
| CPU | same, NCH=12 | Clang 21.1.8 Debug | 1392 | `983fa84ba68fc867` | `dfc61b3e371e80f1` |
| CPU | same, NCH=5 | Clang 21.1.8 Release/fast-math | 944 | `cd87b2262663a983` | `7966eedf67f341b5` |
| CPU | same, NCH=8 | Clang 21.1.8 Release/fast-math | 1136 | `5c766159af330d3e` | `ad91f8f42bf550ce` |
| CPU | same, NCH=12 | Clang 21.1.8 Release/fast-math | 1392 | `21b0119f4f1ec474` | `1feb3cf7a32e5b35` |
| Dual | Complete P7-G2: 61 times + 61 voltages + 61x10 raw tangents | Clang 21.1.8 Debug | 732 | `6dd17b952b7e315a` | `c286fa543130b91b` |
| Dual | same | Clang 21.1.8 Release (strict-FP producer) | 732 | `06dc78ab74d5e7a4` | `d10596df0d2c2c5a` |
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

Audit discovery: the Dual observation omits the entropic OCV correction away
from reference temperature.  P7-G2 is at the reference temperature, so its frozen
bits remain the refactor gate.  The off-reference defect will be corrected in a
separate test-first change, with its intentional output transition documented
rather than hidden inside the refactor.
