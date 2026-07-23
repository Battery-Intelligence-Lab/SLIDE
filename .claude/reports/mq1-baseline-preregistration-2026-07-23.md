# MQ.1 three-lane baseline — PREREGISTRATION (2026-07-23)

Written before configuring or running the MQ.1 gate. The source under test is clean
`Claude` commit `a49224e`; this report-only commit does not change build inputs.

## Claim and falsifier

**Claim.** Current `HEAD` still satisfies the M0.10 baseline in three independent native
configurations: Debug, fast-math Release, and CUDA Release.

**Sharpest falsifier.** A nominally green run can execute stale or cross-lane binaries. The
existing `build-release/bin` contains old `*_fast`/`*_ipo` executables outside its Ninja graph,
and prior sanitizer work showed that a shared runtime directory can make one build run another
build's tests. Therefore MQ.1 uses three new, tree-local build directories and judges only CTest's
registered commands, never a binary glob:

- `build-mq1-debug`
- `build-mq1-release`
- `build-mq1-cuda`

Before each suite runs, `ctest -N -V` must enumerate exactly 58 tests. Every command that runs a
built SLIDE test executable must resolve inside that lane's own build tree; the three intentional
structural tests must invoke CMake with the current source-tree scripts and current lane compiler.
A missing executable, a path from another build/source tree, or any count other than 58 fails the
lane before simulation.

## Fixed configurations

All three lanes use CMake 4.2.3, Ninja 1.13.2, and Windows Clang 21.1.8, targeting
`x86_64-pc-windows-msvc` on an Intel Core Ultra 9 285. Exact scalar fingerprints are enabled
because this is the recorded capture host/toolchain. Debug and CUDA retain the M0.10 ThinLTO
configuration; the plain fast-math Release lane has IPO off, matching its recorded cache.
Fresh-cache defaults keep every other optional dependency off. Builds and tests run sequentially
across lanes; test execution is serial within each lane.

```powershell
cmake -S . -B build-mq1-debug -G Ninja `
  -DCMAKE_BUILD_TYPE=Debug `
  -DCMAKE_CXX_COMPILER="C:/Program Files/LLVM/bin/clang++.exe" `
  -DENABLE_IPO=ON `
  -DSLIDE_WITH_CUDA=OFF `
  -DSLIDE_ENABLE_RECORDED_SCALAR_BITS=ON
cmake --build build-mq1-debug --parallel 2

cmake -S . -B build-mq1-release -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_CXX_COMPILER="C:/Program Files/LLVM/bin/clang++.exe" `
  -DENABLE_IPO=OFF `
  -DSLIDE_WITH_CUDA=OFF `
  -DSLIDE_ENABLE_RECORDED_SCALAR_BITS=ON
cmake --build build-mq1-release --parallel 2

$env:Path = "C:/Program Files/Microsoft Visual Studio/18/Community/VC/Tools/MSVC/14.50.35717/bin/Hostx64/x64;$env:Path"
cmake -S . -B build-mq1-cuda -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_CXX_COMPILER="C:/Program Files/LLVM/bin/clang++.exe" `
  -DCMAKE_CUDA_COMPILER="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.0/bin/nvcc.exe" `
  -DCMAKE_CUDA_ARCHITECTURES=89 `
  "-DCMAKE_CUDA_FLAGS=--allow-unsupported-compiler" `
  -DENABLE_IPO=ON `
  -DSLIDE_WITH_CUDA=ON `
  -DSLIDE_ENABLE_RECORDED_SCALAR_BITS=ON
cmake --build build-mq1-cuda --parallel 2
```

The CUDA host-compiler path exists at registration time. `nvidia-smi` identifies an
`NVIDIA RTX 4000 Ada Generation`, compute capability 8.9; nvcc is CUDA 13.0.

## Registered gates

The following bands are fixed before the runs:

1. **Source identity:** before each configure and after the final witness, the worktree is clean
   and `HEAD` remains this preregistration commit. Its only build-input difference from source
   commit `a49224e` is this Markdown report, so all compiled sources remain the registered tree.
2. **Configure/build:** every command above exits zero. Release compile commands contain the
   project's `-Ofast` and `-ffast-math` flags and no IPO flag; Debug and CUDA retain ThinLTO.
   The CUDA cache records `SLIDE_WITH_CUDA=ON`, architecture 89, and the unsupported-host opt-in,
   and its compile database contains the production `.cu` translation unit.
3. **Lane identity:** `ctest --test-dir <tree> -N -V` discovers exactly 58 tests, with every
   built-test command rooted in `<tree>/bin`; the three CMake-driven structural commands name
   `C:/D/git/SLIDE/tests/structural`, this source tree, and the lane's configured compiler.
   Debug and Release must register `unit_test_core_CudaDisabled` and must not register
   `unit_test_core_CudaSpmBatch`; CUDA must make the inverse substitution. The decisive discovery
   and provenance lines are quoted verbatim in the validation report before execution.
4. **Full suites:** the unfiltered command
   `ctest --test-dir <tree> --output-on-failure --no-tests=error -j1` reports exactly
   **58/58 passed** in Debug, Release, and CUDA. A skip, `Not Run`, missing executable, failed
   test, or count delta is a failure; no per-executable substitute is accepted.
5. **Direct-CWD witness:** CTest's intentional `<build>/bin` working directory cannot distinguish
   the old `../..` data fallback from the new source-root definition. After each full suite,
   run that lane's tree-local `unit_test_core_Sei.exe` directly while the process working directory
   is `C:/D/git/SLIDE`. Each of the three runs must exit zero with zero failed assertions/test cases
   and no data-file-open error. This is the shortest test that discriminates the landed CWD fix.
6. **CUDA witness:** the CUDA lane's `unit_test_core_CudaSpmBatch` must actually execute and report
   exactly **433,671 assertions in 4 test cases**, the frozen M0.10 witness. Merely compiling CUDA
   or passing a CPU fallback is insufficient. After the unfiltered suite passes,
   `ctest --test-dir build-mq1-cuda -V -R '^unit_test_core_CudaSpmBatch$'` captures this summary;
   it is additive evidence, not a substitute for the full suite.
7. **Failure handling:** any red result leaves MQ.1 unticked and is investigated against the first
   failing test. A count delta must be explained or ledgered before any MQ implementation lands.

## Adversarial scope check

Commits since M0.10 include build-policy/path changes, shared Arrhenius routing, new test oracles,
and a topology-complexity fix; they are not documentation-only and therefore require all three
lanes. The all-mask `AgeingKernel` recorded trace is the strongest existing Release oracle for the
Arrhenius refactor. A fresh CUDA compile is the only present check that the new directory-wide
`SLIDE_ROOT_DIR` string definition survives nvcc's host path.

The quality backlog already records a fast-math-unsafe `SpmObservables` validity gate and untested
CUDA rejection arms. MQ.1 establishes the existing 58-test baseline; it does **not** claim those
known coverage gaps are absent or resolved. MQ.2 owns their disposition. Concurrent work elsewhere
on this machine makes wall-clock numbers unusable, so this gate makes no timing or performance
claim. It also makes no ASan/UBSan, TSan, coverage, hosted-CI, cross-platform, or installed-package
claim; those are not MQ.1's named gate.
