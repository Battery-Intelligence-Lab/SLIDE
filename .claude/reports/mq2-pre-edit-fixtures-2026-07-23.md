# MQ.2 pre-edit fixtures (2026-07-23)

## Provenance

Captured after preregistration commit `c48f545a8490` and before any MQ.2
source/test edit. `git status --porcelain` was empty.

The executables are the successful MQ.1 trees built at report-only commit
`af898d1`; this is applicable to current source because:

```text
git diff --name-only a49224e..c48f545 -- src tests cmake benchmark
<empty>
```

Thus the code and tests in the binaries equal current pre-edit code. These runs
freeze fixtures; they do not claim that future source has passed.

## Touched-binary assertion/test-case floor

Each executable was launched directly and returned zero. Debug and fast-math
Release produced the same summaries:

| Catch2 binary | Debug assertions / cases | Release assertions / cases |
|---|---:|---:|
| `core_AgeingKernel` | 384 / 6 | 384 / 6 |
| `core_SpmFactory` | 58 / 7 | 58 / 7 |
| `core_SpmPipeline` | 18 / 1 | 18 / 1 |
| `core_SpmObservables` | 21 / 2 | 21 / 2 |
| `core_SpmElectrical` | 39 / 2 | 39 / 2 |
| `core_SpmStress` | 19 / 3 | 19 / 3 |
| `core_Sei` | 57 / 3 | 57 / 3 |
| `core_SurfaceCrack` | 82 / 3 | 82 / 3 |
| `core_ThermalLumped` | 27 / 2 | 27 / 2 |
| `core_StateArena` | 626 / 11 | 626 / 11 |
| `core_CompiledCurve` | 172 / 7 | 172 / 7 |
| `core_PackSolver` | 878 / 26 | 878 / 26 |
| `core_ModeC` | 52 / 6 | 52 / 6 |
| `core_PackStepper` | 214 / 8 | 214 / 8 |
| `core_PackTopology` | 148 / 9 | 148 / 9 |
| `core_NetlistCsv` | 821 / 9 | 821 / 9 |
| `core_Recorder` | 222 / 7 | 222 / 7 |
| `core_AsyncRecorder` | 471 / 12 | 471 / 12 |
| `core_ParameterSet` | 2,375 / 8 | 2,375 / 8 |
| `core_ParserAllocation` | 1,143 / 13 | 1,143 / 13 |
| `core_Experiment` | 436 / 14 | 436 / 14 |
| `core_Simulation` | 175 / 6 | 175 / 6 |
| `core_ThreadPool` | 605 / 9 | 605 / 9 |

New tests may increase these values; no pre-existing assertion or case may
disappear without a row-specific explanation.

## Embedded digit fixtures

These constants are already executable assertions in the pre-edit tests:

| Fixture | Strict Debug `(values, fnv1a, mixed)` | Fast-math Release | Fast-math IPO |
|---|---|---|---|
| all-mask ageing | `1077, 87119b1b6fa83b81, 55296fa07292792f` | `1077, dc8f92a59dd8d67f, 59085638e06725d3` | `1077, 0be580cf849e57e1, 4c6451971b74f789` |
| Experiment parser | `96, 1e404a3be7365df5, 94a28984eed0159b` | same | same |
| Experiment runner | `503, c779e41caac8339c, 32dec619a41e0324` | `503, e9616b6ce3131bb6, b72cecb529e7c3c0` | `503, 9d97787b3976978b, 4c7bcd0565164aa2` |
| Chen values | `a5d211a63dd77522, 8ccfefe72e03094e` | `e5c3dc16e932ad03, b4e4cfcce5167740` | same as Release |
| compiled SPM input | `134b5b4650709d1a, 1e63d406b3d98297` | `705555174f1b2e23, d9357b2ee6af55e4` | same as Release |
| BPX values | `04bdc30e3e6ca305, e75964cae8725037` | `d3f393daa406ad80, f7d10dea65f3c8f5` | same as Release |

The ParameterSet metadata fixtures are configuration-independent:

```text
Chen: entries=53 values=11537 strings=159 bytes=3491
      fnv1a=8ea717ef0459a7b7 mixed=745379e447872b1f
BPX:  entries=35 values=4135 strings=105 bytes=1906
      fnv1a=f4a1c008bf74be9a mixed=fa654b7fe7b9ab89
compiled-input values=11635
```

The CUDA CPU/device trace is:

```text
values=1388 fnv1a=a59c34d1685ddb19 mixed=2e88cc95b54eb713
```

MQ.1 additionally froze `core_CudaSpmBatch` at 433,671 assertions / 4 test
cases. MQ.2 intentionally changes only the Experiment runner and SurfaceCrack/
all-mask-ageing hashes; every other constant above is a no-op gate.

## Behavior-change starting points

- The passing US06 assertions pin the old lagged drive samples:
  `solution.current[9] == 1.0` and `solution.current.back() == -0.5`.
  MQ.2 replaces them with the registered right-endpoint sequence.
- `RecordingReader::open` currently rejects a correctly computed zero CRC via
  the explicit `expected_crc == 0` clause; the registered 3,584-byte fixture
  isolates that branch before its behavior fix.
- `SpmPipeline::advanceEuler` currently has an assertion-only output-size
  boundary and no all-current finite precheck. The registered invalid-input
  tests must fail against this baseline before the fix.
- The factory row-fill rewrite may not begin until its dedicated NCH=12,
  thermal+all-ageing, nine-lane padded-state digest is captured in an
  oracle-first test commit. That fixture is intentionally separate because the
  current test suite does not expose the complete arena bytes.

No performance, sanitizer, coverage, hosted-CI, package, or cross-platform
claim is made by this capture.
