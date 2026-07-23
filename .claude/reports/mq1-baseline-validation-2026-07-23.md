# MQ.1 three-lane baseline — VALIDATION (2026-07-23)

## Result

**PLAN-level PASS; one preregistered environment subclaim FALSIFIED and one provenance cadence
qualified.** Source commit `a49224e`, plus the report-only preregistration commit `af898d1`,
matches the M0.10 baseline in three fresh, independent native lanes:

| Lane | Configuration | Discovery | Full CTest | Delta from M0.10 |
|---|---|---:|---:|---:|
| Debug | Clang 21.1.8, Debug, ThinLTO, CUDA off | 58 | 58/58 | 0 |
| Release | Clang 21.1.8, Release, `-Ofast -ffast-math`, IPO off, CUDA off | 58 | 58/58 | 0 |
| CUDA | Clang 21.1.8 + nvcc 13.0.48/MSVC 19.50, Release, host-C++ ThinLTO, sm_89 | 58 | 58/58 | 0 |

The preregistered CUDA witness also passed with exactly **433,671 assertions in 4 test cases**.
Each lane's data-loading SEI binary passed when launched directly from the repository root:
**57 assertions in 3 test cases**, exit 0. That direct launch, not ordinary CTest's
`<build>/bin` working directory, is the discriminator for the landed CWD-independence fix.

No source or test file changed during the gate. Immediately after the final witness:

```text
HEAD = af898d1c0172
git status --porcelain = <empty>
git diff --name-only a49224e..HEAD =
.claude/reports/mq1-baseline-preregistration-2026-07-23.md
```

Clean-`af898d1` snapshots were captured immediately before Debug, Release, and the initial CUDA
configure, again in every lane's post-build discovery, and after the final witness. A separate
snapshot was **not** captured immediately before CUDA corrective configure attempts 2 and 3, so
the preregistration's literal “before each configure” observation cadence was not fully met.
However, the successful CUDA tree was observed clean before execution and an immediate rebuild at
that state returned `ninja: no work to do`; this strongly corroborates that the tested artifacts
were current for the final registered source tree without claiming the missing intermediate
snapshots occurred.

Preregistration:
`.claude/reports/mq1-baseline-preregistration-2026-07-23.md`.

## Lane provenance

All lanes used fresh tree-local output; no executable glob and no pre-existing or stale artifact
from another tree entered the evidence. Each successful build was followed by a no-op rebuild
(`ninja: no work to do.`). All successful caches record
`CMAKE_HOME_DIRECTORY:INTERNAL=C:/D/git/SLIDE`; the common tools were CMake 4.2.3 and Ninja 1.13.2.

### Debug

Cache/provenance:

```text
CMAKE_BUILD_TYPE:STRING=Debug
CMAKE_CXX_COMPILER:UNINITIALIZED=C:/Program Files/LLVM/bin/clang++.exe
ENABLE_IPO:BOOL=ON
SLIDE_ENABLE_RECORDED_SCALAR_BITS:BOOL=ON
SLIDE_WITH_CUDA:BOOL=OFF
build.ninja lines containing -flto=thin: 334
```

Discovery found 58 commands: 55 build-local executables, three intentional CMake-driven
structural tests, and zero foreign executable paths. The lane substitution was:

```text
53: Test command: C:\D\git\SLIDE\build-mq1-debug\bin\Debug\unit_test_core_CudaDisabled.exe
Test #53: unit_test_core_CudaDisabled
Total Tests: 58
```

The three structural commands named the current
`C:/D/git/SLIDE/tests/structural/{pc10_single_source,p9c_architecture,p9c5_header_selfcontained}.cmake`
scripts and the configured Clang compiler.

Full gate:

```text
100% tests passed, 0 tests failed out of 58
```

Direct-CWD witness:

```text
All tests passed (57 assertions in 3 test cases)
DIRECT_CWD_EXIT=0
```

### Plain fast-math Release

Cache/provenance:

```text
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_CXX_COMPILER:STRING=C:/Program Files/LLVM/bin/clang++.exe
ENABLE_IPO:BOOL=OFF
SLIDE_ENABLE_RECORDED_SCALAR_BITS:BOOL=ON
SLIDE_WITH_CUDA:BOOL=OFF
build.ninja lines containing -Ofast: 228
build.ninja lines containing -ffast-math: 228
build.ninja lines containing -flto=thin: 0
```

Discovery again found 58 commands, 55 build-local executables, three current-source structural
commands, zero foreign paths, and the expected `CudaDisabled` substitution:

```text
53: Test command: C:\D\git\SLIDE\build-mq1-release\bin\Release\unit_test_core_CudaDisabled.exe
Test #53: unit_test_core_CudaDisabled
Total Tests: 58
```

Full gate and direct-CWD witness:

```text
100% tests passed, 0 tests failed out of 58
All tests passed (57 assertions in 3 test cases)
DIRECT_CWD_EXIT=0
```

### CUDA Release / host-C++ ThinLTO

The successful evidence tree is `build-mq1-cuda-vsenv`. Cache/provenance:

```text
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_CUDA_ARCHITECTURES:UNINITIALIZED=89
CMAKE_CUDA_COMPILER:STRING=C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.0/bin/nvcc.exe
CMAKE_CUDA_FLAGS:STRING=--allow-unsupported-compiler
CMAKE_CXX_COMPILER:STRING=C:/Program Files/LLVM/bin/clang++.exe
CMAKE_RC_COMPILER:FILEPATH=C:/Program Files (x86)/Windows Kits/10/bin/10.0.26100.0/x64/rc.exe
ENABLE_IPO:BOOL=ON
SLIDE_ENABLE_RECORDED_SCALAR_BITS:BOOL=ON
SLIDE_WITH_CUDA:BOOL=ON
build.ninja lines containing -flto=thin: 336
```

The compile database contains exactly one `CudaSpmRuntime.cu` command. It invokes nvcc, contains
`--allow-unsupported-compiler`, generates compute/sm 89 code, and contains the new
`SLIDE_ROOT_DIR` definition. This confirms that the directory-wide definition reaches CUDA
compilation and that nvcc accepts its quoting.

Discovery found 58 commands, 55 build-local executables, three current-source structural commands,
zero foreign paths, and the CUDA implementation rather than the disabled shim:

```text
53: Test command: C:\D\git\SLIDE\build-mq1-cuda-vsenv\bin\Release\unit_test_core_CudaSpmBatch.exe
Test #53: unit_test_core_CudaSpmBatch
Total Tests: 58
```

Full gate and direct-CWD witness:

```text
100% tests passed, 0 tests failed out of 58
All tests passed (57 assertions in 3 test cases)
DIRECT_CWD_EXIT=0
```

The additive verbose CUDA witness reported:

```text
All tests passed (433671 assertions in 4 test cases)
100% tests passed, 0 tests failed out of 1
```

## CUDA configuration falsifiers

The original preregistration repeated the prior shorthand that `cl.exe` on `PATH` is sufficient
and required that literal configure command to exit zero. That registered environment
precondition is **FALSIFIED**. A **fresh** CMake 4.2 CUDA tree disproved the advice:

1. With only the VS host-compiler bin directory added, nvcc compiled the probe but CMake invoked
   `vs_link_exe` with `--rc=...\llvm-rc.exe --mt=""`; manifest embedding reported
   `MT: command "mt ..." failed` followed by `no such file or directory`.
2. Explicit SDK `rc.exe`/`mt.exe` resolved the manifest step, but changing those tools in the
   partial cache forced compiler rediscovery. CMake selected `cl.exe` without the complete VS
   library environment, and the link failed:
   `LINK : fatal error LNK1104: cannot open file 'kernel32.lib'`.
3. A new equivalent-configuration tree under the complete `vcvars64.bat` environment, while
   explicitly retaining Clang as `CMAKE_CXX_COMPILER`, configured, completed all 339 reported Ninja
   build steps, and passed the PLAN-level lane gate plus every runtime witness above.

Therefore the durable environment rule is: a fresh Windows CUDA configure requires the full VS
x64 developer environment (host compiler plus SDK manifest/resource tools and
`INCLUDE`/`LIB`/`LIBPATH`), not merely `cl.exe` on `PATH`. The two failed attempts reused one
retained partial tree (`build-mq1-cuda`); the successful equivalent configuration lives in
`build-mq1-cuda-vsenv`. No untracked files were deleted.

The distinctive failure signatures are preserved in this report. While the retained partial tree
exists, `build-mq1-cuda/CMakeFiles/CMakeConfigureLog.yaml` corroborates attempt 2's
`kernel32.lib` failure. Cache regeneration erased attempt 1 from that log, so its empty manifest
tool failure survives in the command transcript and this tracked report; the report remains the
durable record if MQ.3 later proposes deleting the untracked tree.

## Adversarial interpretation and scope

This result restores the named baseline; it does not prove the test suite is complete. The
quality backlog already owns the fast-math-unsafe `SpmObservables` range gate
(`fastmath-unsafe-stoichiometry-gate`), the assert-only `advanceEuler` span boundary
(`advance-euler-asserts-only`), and untested CUDA rejection/cadence arms
(`test-gap-untested-rejection-arms`). They remain MQ.2 inputs, not hidden behind the green count.

Two additional build observations also go to MQ.2 disposition:

- Clang 21 warns that `-Ofast` is deprecated; the project already passes `-ffast-math`, so the
  spelling is redundant policy debt rather than a gate failure. This is a compiler-driver warning
  on Release compilation, including v4 targets; the earlier “zero core warnings” source audit does
  not erase it.
- the CUDA build emits the retained legacy `benchmark_LP_cases` narrowing warning through
  `Deep_ptr.hpp`; that diagnostic is in the legacy/benchmark surface, not the v4 core. Neither
  observation is one of the 68 existing survivors, so MQ.2 must add and disposition both rather
  than silently dropping them.

No performance claim is made. No ASan/UBSan, TSan, coverage, hosted-CI, installed-package, Linux,
or macOS claim is made. The reported CTest durations are deliberately omitted as evidence.
