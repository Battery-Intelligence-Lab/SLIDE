# Handoff — 2026-07-10 — P8-G0 optionality and portability

## Goal

Close the registered P8-G0 gate without relying on cached builds, zero-test CTest output, or optional SDKs.

## What changed

- Added `SLIDE_BUILD_CORE_TESTS` and a dependency-light smoke that compiles public core headers and checks StateArena alignment/geometry, ThreadPool exact-once execution, deterministic reduction, and explicit disabled CUDA/zstd/Arrow behavior.
- Added a nested-superproject consumer. This forced the public core include root from `CMAKE_SOURCE_DIR` to `PROJECT_SOURCE_DIR`, so `add_subdirectory(SLIDE)` works from a parent project.
- Made `SLIDE_WITH_CUDA` private to the implementation boundary; CUDA SDK types remain confined to the `.cu` file.
- Made Eigen discovery find-package-first while preserving the existing MKL/vectorization/MPL2 options and pinned CPM fallback.
- Added Linux/macOS/Windows core-only CI and made `cmake/**` trigger installed-wheel CI.
- Recorded D-25: “dependency-free” in P8-G0 means free of optional toolchain/runtime dependencies; Eigen remains the required cold-path linear-algebra dependency.

## What was tested

- Fresh root `SLIDE_CORE_ONLY=ON` configure/build with CUDA, MATLAB, zstd, Arrow, and Python disabled: smoke 1/1.
- Fresh nested-superproject configure/build through the public `slide_core` target: external consumer 1/1.
- Full optional-off Debug CPU suite: 49/49.
- Full optional-off Release CPU suite: 49/49.
- Rebuilt CPython 3.13 Windows wheel, installed outside the source tree: 9 passed, 2 expected optional skips (PyBOP and CUDA).
- `clang-format --dry-run --Werror` for the new C++ smoke and `git diff --check`.

## Key results and interpretation

The pre-change core-only configuration built but registered 0 tests, so P8-G0 was initially falsified. The new root and nested gates both pass, the full CPU suites remain green, and the isolated wheel works. No physics implementation changed, so this is a portability/build-boundary improvement rather than a numerical branch. Cross-platform CI definitions are committed but cannot be called green for this commit until the branch is pushed and those jobs execute.

## Recommended next step

Implement P8-G5: replace the v3-era docs with tested v4 installation and C++/Python/MATLAB quickstarts, extension guides, optional-dependency matrix, and declared PyBaMM gaps; extract and execute every quickstart in its available toolchain.
