---
layout: default
title: Dependencies and compatibility
nav_order: 5
---

# Dependencies and compatibility

## Dependency and capability matrix

| Capability | CMake/API switch | Default | Requirement | Behavior when absent |
|---|---|---:|---|---|
| C++ core | `SLIDE_CORE_ONLY=ON` | off | C++20, CMake 3.31, Eigen 3.4 | Installed Eigen is preferred; the pinned source fallback is fetched when needed |
| Python wheel | `SLIDE_BUILD_PYTHON=ON` | off | Python ≥3.10, nanobind, scikit-build-core, NumPy | No Python extension is built |
| MATLAB | `SLIDE_WITH_MATLAB=ON` | off | MATLAB with MEX compiler | Explicit `ON` fails configuration if MATLAB is unavailable |
| CUDA SPM batch | `SLIDE_WITH_CUDA=ON` | off | CUDA compiler and toolkit | CPU remains available; `available_devices()`/`availableDevices()` omit CUDA |
| zstd recording | `SLIDE_WITH_ZSTD=ON` | off | installed zstd or pinned source fallback | raw shuffled blocks remain available; requesting zstd is rejected during recorder configuration |
| Arrow/Parquet | `SLIDE_WITH_ARROW=ON` | off | installed Arrow and Parquet CMake packages | CSV and native mmap remain; Parquet returns `NotImplementedYet` |
| Plotting | Python extra `plot` | off | Matplotlib | numerical solve and data access remain available |
| MAT-file export | Python extra `matlab` | off | SciPy >=1.10 | solves still run; `Solution.save_data(..., to_format="mat")` is unavailable |
| PyBaMM reference tests | Python extra `test` | off | pinned PyBaMM 26.6.2.0 | hermetic C++ fixture comparisons still run |
| PyBOP fitting | Python extra `test` | off | pinned PyBOP 25.11 on its supported Python range | `slide.pybop` fitting helpers are unavailable; core sensitivities remain |

The core is **dependency-light**, not literally dependency-free: Eigen is required for cold-path spectral decomposition and sparse pack solving. CUDA, MATLAB, zstd, Arrow, TBB, and the v3 dependency stack are not required by a CPU core-only build.

## PyBaMM-shaped API: known gaps

“PyBaMM-compatible” means familiar setup and result access, not full object or numerical equivalence.

- v4.0 provides SPM only, with registered Chebyshev orders `nch={5,8,12}`. SPMe, DFN, arbitrary expression trees, custom spatial methods, and `var_pts` are not accepted.
- Chen2020 and the supported SPM subset of BPX 1.x are absorbed. Native curve overrides and state-dependent solid diffusivity are rejected rather than silently approximated.
- Isothermal and lumped-thermal compositions are available; arbitrary PyBaMM thermal submodels are not.
- Python custom controls/events receive time, local time, voltage, current, and power. Callbacks that require arbitrary internal model variables are outside the current surface.
- Signature-compatibility fields currently ignored are experiment-level `temperature` and `termination`, plus step-level `temperature`, `tags`, `description`, and `skip_ok` (and `direction` on parsed string steps). Step voltage/current/custom terminations remain active.
- `varied()` lanes currently cannot be combined with structured/custom experiment steps.
- CUDA dispatch supports the explicitly preflighted base-isothermal, fixed-duration constant-current single/`varied()` workload. Unsupported compositions fail before a solve.
- `simulateS1` sensitivities require one fixed-duration constant-current segment and cannot be combined with parameter sweeps.
- `Solution` exposes time, terminal-voltage/voltage, and current series; arbitrary PyBaMM model variables are not recorded.
- Cross-tool agreement is tolerance-based. Different spatial discretisations and interpolants mean SLIDE and PyBaMM are not expected to be bit-identical.

## Portability evidence

The committed CI matrices cover Linux, macOS, and Windows core-only consumers plus installed Python wheels on Python 3.10 and 3.13. The MATLAB MEX gate runs on an explicitly provisioned, licensed R2025b machine; public-project batch licensing is not treated as proof of the external C++ interface. CUDA validation likewise requires an explicitly provisioned GPU machine.
