# SLIDE Development Discussions

This document tracks design decisions, architecture evolution, and session notes for long-term pair programming collaboration on SLIDE.

---

## Active Design Decisions

### [2026-01-29] Initial Codebase Review & Collaboration Setup

**Status**: Decided

**Context**: First comprehensive review of SLIDE codebase to establish pair programming infrastructure for multi-year development.

**Key Findings**:
1. Architecture is well-structured: StorageUnit → Cell → Module → Battery hierarchy
2. C++20 codebase with modern practices (Deep_ptr, std::span, constexpr)
3. Good test infrastructure (Catch2) but some tests commented out
4. Python/MATLAB bindings planned but not implemented
5. TODO.md has 365 lines of unstructured items

**Decision**: Create comprehensive documentation and tracking infrastructure:
- Updated CLAUDE.md with SLIDE-specific domain knowledge
- Updated style guides (cpp-style.md, python-style.md) for SLIDE
- Create skills for common operations
- Reorganize TODO.md with priorities
- Use this discussions.md for ongoing decision tracking

---

### [Pending] Python Bindings Strategy

**Status**: Open

**Context**: Python bindings are a high-priority feature request for broader adoption.

**Options**:
1. **pybind11** - Mature, well-documented, widely used
   - Pros: Large ecosystem, good numpy integration
   - Cons: Compile times, header-only can increase binary size

2. **nanobind** - Modern successor to pybind11 by same author
   - Pros: Faster compile, smaller binaries, better C++17/20 support
   - Cons: Smaller community, fewer examples

3. **SWIG** - Language-agnostic wrapper generator
   - Pros: Can generate Python and MATLAB simultaneously
   - Cons: Less Pythonic output, steeper learning curve

**Decision**: TBD - Need to evaluate compile time impact and numpy integration quality

---

### [Pending] `redistributeCurrent_new` Performance

**Status**: Open

**Context**: The current redistribution algorithm for parallel modules can require up to 2500 iterations (TODO.md line 79), causing significant slowdowns.

**Options**:
1. **Newton-Raphson solver** - Quadratic convergence
2. **Direct linear algebra** - Solve coupled equations directly
3. **Adaptive tolerance** - Accept less precision when converging slowly
4. **Better initial guess** - Use previous solution as starting point

**Decision**: TBD - Need profiling data to identify root cause

---

## Architecture Evolution Log

### Version 3.x Roadmap

**Current**: v3.0.0 (slide-pack merged, C++20)

#### Planned Milestones
- [ ] **v3.1.0**: Python bindings (pybind11/nanobind)
- [ ] **v3.2.0**: MATLAB MEX interface
- [ ] **v3.3.0**: PyBaMM-compatible Experiment interface
- [ ] **v3.4.0**: SUNDIALS solver integration (optional)
- [ ] **v4.0.0**: GPU acceleration (CUDA, optional)

#### Completed Milestones
- [x] v3.0.0 (Current)
  - Merged SLIDE and slide-pack
  - Upgraded to C++20
  - Added Cell_ECM, Cell_Bucket
  - Converted shared_ptr → unique_ptr/Deep_ptr
  - Added std::span for state access
  - Status class for error handling
  - Model_SPM optimization: 34s → 12s (2.8x speedup)
  - File I/O elimination: 30s improvement

---

## Performance Benchmarks History

| Date | Version | Test | Time | Hardware | Notes |
|------|---------|------|------|----------|-------|
| 2026-01-29 | v3.0.0 | 5000 1C CC cycles | <1 min | - | Baseline from README |
| 2026-01-29 | v3.0.0 | CC+CV cycles | <2 min | - | Baseline from README |
| 2026-01-29 | v3.0.0 | EPFL battery 1hr CC | ~2 sec | - | From TODO.md |
| 2026-01-29 | v3.0.0 | EPFL battery + aging | ~3.5 sec | - | From TODO.md |

---

## Known Technical Debt

| Priority | Issue | File | Notes |
|----------|-------|------|-------|
| High | `redistributeCurrent_new` 2500+ iterations | Module_p_impl.cpp | Performance bottleneck |
| High | MSVC 3x slower than Clang | CMakeLists.txt | Vectorization flags |
| Medium | Raw `assert()` in tests | Cycler_test.cpp | Should use Catch2 REQUIRE() |
| Medium | Mixed `#define` and `constexpr` | settings.hpp:101 | DATASTORE_BATT macro |
| Medium | `StorageUnit::copy()` returns raw ptr | StorageUnit.hpp | Should return unique_ptr |
| Medium | C-array `double Tneighb[]` params | Module.hpp:284 | Use std::span |
| Low | Inconsistent Status codes | Various | Some int, some Status enum |
| Low | Global `settings::isParallel` | settings.hpp | Encapsulate thread control |
| Low | `therm.Qcontact` 6x overestimated | Module thermal | Calculation bug |

---

## Session Notes

### [2026-01-29] Initial Codebase Exploration

- **Focus**: Complete exploration of SLIDE library, documentation setup, collaboration infrastructure
- **Decisions**:
  - Style guides updated in-place (not new files) for SLIDE
  - Completed TODO items archived to separate COMPLETED.md
  - Session notes use summary bullets + rationale format
- **Next Steps**:
  - Create skill files for common operations
  - Restructure TODO.md into prioritized format
  - Run test suite to verify build still works
- **Blockers**: None

---

## 2026-07-07 — v4 architecture decided (Fable session)

- **Decision**: Full v4 refactor architecture fixed and recorded in [/PLAN.md](../PLAN.md) §3–§4 (14-entry
  decision log with rejected alternatives). Headlines: SoA state arena with archetype batches (D-01),
  compile-time model composition with per-batch dispatch (D-02), netlist compiled flat to sparse MNA with
  three solver modes — sparse Newton / Thomas ladder / waveform-relaxation+Baumgarte (D-03..05), exponential
  modal propagator replacing forward Euler on diffusion (D-07), Strang multirate (D-08), strangler migration
  with digit-diff parity gates (D-14).
- **Evidence**: three-scout audit 2026-07-07; findings condensed in PLAN.md §2 with file:line citations.
  Chebyshev nch≠5 root-caused (centre-node `-0.5` sign + Eigen inverse aliasing; already fixed on `Claude`).
- **Next Steps**: Phase 0 bug fixes (2 Opus agents, in flight); Phase 1 core data model after Volkan reviews
  PLAN.md §3/§4 and §7 open questions.
- **Blockers**: PLAN.md §7 Q1–Q7 await Volkan's answers (assumptions recorded, work proceeds on them).

---

## 2026-07-08 — Phase 1: Q8 Release re-confirm + production diffusion kernel

- **Q8 fully closed.** P1-G0 pilot re-run in Release/-O3 (`build-release`): decisive band rel ≤ 1e-12
  HOLDS (max_rel 5.63e-15); exact bit-identity (Debug: 0) FALSIFIED under -O3 cross-TU FMA contraction
  (max_abs 2.26e-17). `-ffp-contract=off` rejected — can't fix cross-TU on the test target alone, and
  forcing it globally recompiles legacy (§5.1 forbids). Bonus hypothesis falsified & recorded; decisive
  band untouched. Pilot's exact-zero CHECK scoped to Debug. See `handoff-2026-07-08-phase1-diffusion.md`.
- **Production `SpectralDiffusion<NCH>` landed** (`src/core/SpectralDiffusion.hpp`): vectorised SoA
  forward-Euler diffusion. Validated vs the legacy-shaped oracle — Debug bit-identical (proves math
  identity), Release rel 3.8e-15 ≤ 1e-12. Two-statement inner form; once-allocated scratch (PC-1).
- **Design decision surfaced (OPEN for Volkan/Fable):** ThermalLumped + ageing kernels are NOT
  self-contained — they need the observable-reconstruction layer (c_surf = C·z + D·flux, overpotential,
  dOCV, Rdc; §3.7/D-10) and the §3.11 BatchView/StepCtx interface, both DEFERRED this session (diffusion
  took plain spans as a stopgap). Recommend designing that layer + doing **P1-G3 (Chebyshev/Carslaw &
  Jaeger oracle)** — which exercises the same C/D output path — before porting thermal/ageing.
- **Review channels down:** advisor unavailable, Fable agent out of usage credits — no external review
  this session. Decisions rest on the Debug bit-identity arbiter + derived accumulation bounds.

## Quick Links

- [CLAUDE.md](CLAUDE.md) - Main runbook
- [cpp-style.md](cpp-style.md) - C++ conventions
- [python-style.md](python-style.md) - Python conventions
- [develop/TODO.md](../develop/TODO.md) - Active development items
- [develop/COMPLETED.md](../develop/COMPLETED.md) - Archived completed items
- [CHANGELOG.md](../CHANGELOG.md) - Release history
