# P9-G1 whole-project sanitizer evidence — 2026-07-12

## Claim and falsifier

P9-G1 requires the full Debug suite under ASan+UBSan and the concurrency suites under TSan. A lane counts only if production and test compile commands carry the requested flags, final executables contain the sanitizer runtime, and CTest resolves build-local binaries. Any report, leak, UB, race, crash, missing test, or ordinary assertion failure is red.

## Infrastructure

- `cmake/Sanitizers.cmake` applies homogeneous directory-wide compile/link instrumentation before any SLIDE target is created and rejects unsupported compiler, IPO, Windows TSan, or incompatible sanitizer combinations.
- Top-level runtime outputs live under `${CMAKE_BINARY_DIR}/bin/<Config>`; `PROJECT_IS_TOP_LEVEL` prevents an embedded SLIDE from rewriting its parent's policy.
- `.github/workflows/core-sanitizers.yml` reproduces the two Ubuntu lanes and proves flags, symbols, paths, and counts before CTest. The action is committed; no hosted-green claim is made before push.
- TSan label: `unit_test_core_ThreadPool`, `unit_test_core_AsyncRecorder`, `unit_test_core_PackStepper`, all serial.

## Local toolchain and structural proof

- WSL2 Ubuntu 24.04, Clang/Clang++ 18.1.3, CMake 3.31.10, Ninja 1.11.1.
- ASan+UBSan production `ThreadPool.cpp`, legacy `Cell_SPM.cpp`, and the ThreadPool test TU contain `-fsanitize=address,undefined`; both `libslide_core.a` and legacy `libsrc.a` expose `__asan_` and `__ubsan_` references; the test executable defines `__asan_init`; CTest discovers 51 executables under `build-m02-asan-clang/bin/Debug` and none under source-root `bin`.
- TSan production/test commands contain `-fsanitize=thread`; `libslide_core.a` exposes `__tsan_`; the test executable defines `__tsan_init`; CTest selects exactly three build-local tests.
- Negative configure probes pass: sanitizer+IPO and Windows TSan are rejected with explicit fatal diagnostics.

## Runs and findings

1. Initial ASan+UBSan full suite: **50/51**. `unit_test_core_NetlistCsv` aborted before parser entry with a six-byte global-buffer-overflow from a hard-coded 50-byte view over a 44-byte embedded-NUL literal.
2. P9-B38 fix: derive the view extent as `sizeof embedded_nul_csv - 1`. Focused instrumented NetlistCsv: **1/1**.
3. Final ASan+UBSan with `detect_leaks=1`, strict strings, abort-on-error, and halt-on-UB: **51/51**, zero sanitizer finding.
4. Initial four-target TSan link: the three concurrency binaries linked; `AsyncRecorderAllocation` failed structurally because its deliberate global new/delete replacements duplicate TSan's mandatory allocator interceptors. It is excluded explicitly from TSan and remains covered by ASan+UBSan.
5. Final TSan with halt-on-error: PackStepper, ThreadPool, AsyncRecorder **3/3**, zero race report.
6. Ordinary Windows Clang Debug and fast-math Release: **51/51** each. Core-only and nested external-consumer smokes: **1/1** each.

## Interpretation

P9-G1 is locally green. The sanitizer infrastructure is a real whole-project lane rather than a wrapper-only or shared-binary false green. P9-B38 is a confirmed and fixed test-harness defect. The TSan allocation-counter exclusion is structural and narrow: the ordinary AsyncRecorder suite exercises the same drain-thread lifecycle under TSan, so the excluded allocation suite contributes no unique concurrency path.
