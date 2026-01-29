# SLIDE — Claude Runbook (authoritative)

You are the coding agent responsible for improving SLIDE (Simulator for Lithium-Ion Degradation) with:
- Portability first, excellent performance and parallelisation second
- Clean, extensible, maintainable architecture for battery simulation
- Seamless Python + MATLAB interfaces (planned)
- Excellent docs, tests, CI, versioning, changelog discipline
- Keep develop/TODO.md up-to-date with structured milestones

This is a multi-year project. Keep develop/TODO.md current so we can pair program effectively.

## Project Overview

**SLIDE** is a C++20 library for fast lithium-ion battery degradation simulation featuring:
- Single Particle Model (SPM) with coupled bulk thermal model
- Multiple degradation mechanisms (SEI growth, LAM, surface cracks, Li-plating)
- StorageUnit -> Cell -> Module -> Battery hierarchy
- Pack-level simulation with thermal management systems
- Battery tester-like interface (CC, CV, CCCV, profiles)

## Architecture Layers

```
Layer 4: Battery          Complete system with cooling & converter
Layer 3: Module_s/p       Series/parallel configurations
Layer 2: Cell             Cell_SPM (physics-based), Cell_ECM (equivalent circuit)
Layer 1: StorageUnit      Abstract base defining electrical/thermal interface
Layer 0: State/Params     State_SPM, DEG_ID, SEIparam, LAMparam, etc.
```

### Key Components
| Component | Location | Purpose |
|-----------|----------|---------|
| StorageUnit | `src/StorageUnit.hpp` | Abstract base for all battery units |
| Cell_SPM | `src/cells/Cell_SPM/` | Single Particle Model cell (6 impl files) |
| Cell_ECM | `src/cells/Cell_ECM/` | Equivalent Circuit Model (1-3 RC pairs) |
| Module_s/Module_p | `src/modules/` | Series/parallel configurations |
| Battery | `src/system/Battery.hpp` | Complete battery system |
| Cycler | `src/procedures/Cycler.hpp` | CC, CV, CCCV test procedures |
| Procedure | `src/procedures/Procedure.hpp` | High-level aging protocols |
| CoolSystem | `src/cooling/` | Thermal management (HVAC, open-loop) |
| Deep_ptr | `src/types/Deep_ptr.hpp` | Smart pointer with deep copy semantics |

## Non-negotiables

1. Do NOT introduce runtime dependence on repo-relative paths
2. Every PR must:
   - Update CHANGELOG.md (Unreleased section) for user-visible changes
   - Add/adjust tests for changed behavior
   - Keep formatting/lint clean
3. Keep public API small and stable. Hide implementation details
4. Optional dependencies only (OpenMP, Armadillo, HiGHS, CUDA). Core must build without them
5. Make library detections robust across Windows, macOS, Linux

## Domain-Specific Conventions

### Sign Convention
- **Positive current = discharge** (current flowing out of cell)
- **Negative current = charge** (current flowing into cell)

### Units
| Quantity | Unit | Notes |
|----------|------|-------|
| Capacity | Ah | Ampere-hours |
| Voltage | V | Volts |
| Current | A | Amperes |
| Temperature | K | Kelvin internally (use `25.0_degC` literal for conversion) |
| Time | s | Seconds |
| Power | W | Watts |

### State Management
- `State_SPM` holds 29+ state variables (electrical, thermal, degradation, concentration)
- States are std::array-based for performance
- Use `getStates(span)` / `setStates(span)` for bulk operations
- Backup/restore pattern for rollback on constraint violations

### Degradation Models
- Selected via `DEG_ID` structure (flags for each model)
- Models: SEI growth, LAM (loss of active material), surface cracks, Li-plating
- Parameters in `src/cells/Cell_SPM/param/`: SEIparam, LAMparam, CSparam, StressParam

## Code Quality Standards

- **C++20** required (CMake 3.31+, concepts, ranges, std::span)
- No naked new/delete in core - use `Deep_ptr`, `unique_ptr`, standard containers
- Use `std::span` for array views in hot paths (zero-copy)
- Avoid allocations in inner loops; use scratch buffers or thread-local pools
- Provide deterministic RNG seeding options
- State classes inherit from `std::array` for performance

### Patterns to Follow
```cpp
// StorageUnit ownership via Deep_ptr (enables polymorphic deep copy)
Deep_ptr<StorageUnit> su = makeBattery(...);

// State access via span (zero-copy)
std::vector<double> states(su->getNstates());
su->getStates(states);

// Cycler for electrical operations
Cycler cycler(su.get(), "test_id");
cycler.CC(current, vlimit, tlimit, dt, ndt_data);

// Status codes for error handling
Status result = su->setCurrent(I, false);
if (result != Status::Success) { /* handle */ }
```

## Performance Guidelines

- Provide baseline microbenchmarks (see `benchmark/` folder)
- Optimize only when benchmarks show wins. Record numbers in `/benchmark/README.md`
- Prefer clear loops over clever meta-programming
- Known issues to watch:
  - `redistributeCurrent_new` can require 2500+ iterations (optimization target)
  - MSVC builds ~3x slower than Clang (vectorization differences)

## Bindings Strategy

### Python (Planned)
- Will use pybind11 or nanobind
- Expose numpy arrays without copies where possible
- Build wheels via CMake + scikit-build-core
- Run pytest in CI
- Target: PyBaMM-compatible Experiment interface

### MATLAB (Planned)
- Use MEX (MATLAB Executable) as primary route
- Provide `+slide` package with OO wrappers calling MEX
- Keep API symmetric with Python where reasonable
- Current MATLAB code is for post-processing only (reading CSV outputs)

## Documentation Requirements

docs/ should include:
- Installation (C++, Python planned, MATLAB planned)
- Quickstart examples for each interface
- API reference (Doxygen for C++; Python docstrings when available)
- "How to add a new cell type" guide
- "How to add a new degradation model" guide

## Release Discipline

- Use SemVer (MAJOR.MINOR.PATCH)
- CHANGELOG.md follows Keep a Changelog format
- Tag releases; generate GitHub Releases notes from changelog
- VERSION source of truth: CMakeLists.txt PROJECT_VERSION field
- Current version: 3.0.0 (merged slide-pack, C++20)

## Testing

- Framework: Catch2 v3
- Location: `tests/unit/` for unit tests, `tests/integration/` for integration
- Run: `cmake -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build && ctest --test-dir build`
- Use `REQUIRE()` not raw `assert()` for better error messages

## Working Style (Claude Code best practices)

- Always start by exploring and planning; do not jump to edits without understanding
- Make small, reviewable commits
- Prefer refactors behind feature flags/options when risk is high
- Update develop/TODO.md with progress and discoveries
- Use `.claude/discussions.md` to log significant design decisions

## Quick Reference: Adding New Features

### Adding a New Cell Type
1. Create `src/cells/Cell_NewType/Cell_NewType.hpp`
2. Inherit from `Cell` class
3. Implement: `V()`, `getOCV()`, `setCurrent()`, `timeStep_CC()`, state methods
4. Add to `src/cells/cells.hpp`
5. Create `tests/unit/Cell_NewType_test.cpp`

### Adding a Degradation Model
1. Add parameters to struct in `src/cells/Cell_SPM/param/`
2. Add case to switch in `Cell_SPM_degradation.cpp`
3. Update `DEG_ID` documentation
4. Add test case demonstrating the model
