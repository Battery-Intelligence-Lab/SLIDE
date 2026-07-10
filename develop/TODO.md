# SLIDE Development TODO

> **Note**: This is a living document. See [COMPLETED.md](COMPLETED.md) for archived completed items.
>
> **Disclaimer**: Some items are informal notes. Priority may shift based on user needs.
>
> **2026-07-10 — v4 refactor active.** The authoritative vision/roadmap is [/PLAN.md](../PLAN.md)
> (architecture, decision log, phased gates, open questions). Phases 0–8 and 9A are complete;
> Phase 9B systematic adversarial bug-hunt is active; its first Experiment
> validation/rollback, initial recorder, and pack-solver safety defects are fixed and recorded in the Phase-9B ledger. Historical "Critical Bugs"
> below may be superseded by PLAN.md. Update PLAN.md §8, not just this file.

---

## Immediate (Next 1-2 PRs)

### Critical Bugs
- [ ] `redistributeCurrent_new` requires up to 2500 iterations - major performance issue
- [ ] MSVC vs Clang ~3x performance gap - check vectorization/SSE2/AVX flags
- [ ] `therm.Qcontact` is 6x overestimated in thermal model
- [ ] T_MODEL == 2 causes thermal runaway in `test_CyclerVariations_high`
- [ ] `test_Cycler_CoolSystem` passes only when T_MODEL==2
- [ ] `i < getNSUs() - 1` is a bug - fix immediately!

### High-Priority Fixes
- [ ] CV is not doing intended thing for both SLIDE and slide-pack
- [ ] Series module voltage limit handling - what to do if one cell reaches max?
- [ ] CCCV for ageing CV should be done with remaining voltage
- [ ] Cycler CV controls voltage limit unnecessarily
- [ ] Bugfix: Module_s::getI() returns value if (getNSUs() >= 0) but should be >

---

## Short-Term (This Quarter)

### Code Quality
- [ ] Replace `assert()` with Catch2 `REQUIRE()` in tests
- [ ] Convert `#define DATASTORE_BATT` to constexpr (settings.hpp)
- [ ] `StorageUnit::copy()` should return `unique_ptr` not raw pointer
- [ ] `double Tneighb[]` in Module should use `std::span`
- [ ] Unified Status return codes (some return int, others Status enum)
- [ ] Check #CHECK and #TODO tags in the code
- [ ] Copy functions are commented out - review and fix
- [ ] Make more methods const: Vmin(), Vmax(), VMIN(), VMAX(), Cap()
- [ ] Improve const correctness throughout
- [ ] Change const string& to string_view

### Testing
- [ ] Snapshot testing framework
- [ ] Automated tests against PyBaMM
- [ ] Add static analysers: include-what-you-use, valgrind, clang-tidy
- [ ] Test cases testing shared_ptr logic should be removed (using unique_ptr now)

### Features
- [ ] Integrate NLopt for determineOCV optimisation
- [ ] Make modules writable to JSON files (nested structure)
- [ ] Add GITT (Galvanostatic Intermittent Titration Technique) function
- [ ] SOC/Temperature dependent RC pairs for Cell_ECM
- [ ] setVoltage function for better CV period

---

## Medium-Term (Next Release)

### Architecture
- [ ] Voltage should be inside the states (discrete algebraic state)
- [ ] Make overpotential (etap, etan) removable - need voltage model
- [ ] T_MODEL and T_ENV should not be constants
- [ ] Cell_SPM should not hold all ageing model parameters - use composition
- [ ] `settings::isParallel` global should be encapsulated
- [ ] Create Module_p_ApproxPI variant
- [ ] Better findCurrent algorithm for parallel modules
- [ ] EmptyStorageUnit to remove if(nullptr) parent checks
- [ ] More hierarchy for ECM and SPM models
- [ ] Factory methods for battery creation
- [ ] Create another class for thermal model
- [ ] CellData, CellDataStorage, CellDataWriter - need policy design

### Data & Storage
- [ ] Variable data storage with enums and deserializers
- [ ] Create enum for storable data
- [ ] Add generic state term to cover time, Ah, Wh
- [ ] Matio and parquet data types support
- [ ] Create file type to compactly save and retrieve data

### Bindings
- [ ] Python bindings via pybind11/nanobind
- [ ] MATLAB MEX interface
- [ ] PyBaMM-compatible Experiment interface (drop-in replacement)
- [ ] Julia wrapper consideration (Brady suggestion)

### Build & CI
- [ ] CPack installation improvements
- [ ] Fix Ccache not working
- [ ] Add package manager installation option
- [ ] Configure clang-format, cmake-format

---

## Long-Term (Future Releases)

### Architecture Evolution
- [ ] SUNDIALS solver option (Martin Robinson suggestion)
- [ ] GPU acceleration (optional)
- [ ] Higher order integration (Runge Kutta 4) and adaptive time stepping
- [ ] Header-only library portions for easy compilation

### Performance
- [ ] SmallVector and SmallArray optimizations
- [ ] Template-based free functions to replace state-dependent functions
- [ ] Use just one memory allocation for Model_SPM
- [ ] Memoize Cap() calculations
- [ ] Threadlocal vectors for dynamic-sized stack
- [ ] `double degState[CELL_NSTATE_MAX]` causes 700 kB unnecessary memory

### API Improvements
- [ ] Instead of taking unique_ptr or raw ptr, functions should take object references
- [ ] Begin and end functions for StorageUnit to traverse children
- [ ] Visitor pattern for hierarchical structures
- [ ] setStates to support r-values
- [ ] viewStates, viewVariations returns span

---

## Known Issues (Tracking)

### Numerical/Physics
- [x] Chebyshev discretisation only works for nch = 5 — fixed: zero eigenvalue detection used absolute threshold; now uses MATLAB-matching relative normalisation and per-electrode detection
- [x] Model.Dn and Dp created as nch elements but require nch+1 elements — confirmed not a bug: naming confusion between state-space D vector (size nch+1, correct) and state scalars Dp/Dn
- [ ] determineOCV is inefficient - can reduce search space
- [ ] Why does getOCV not include entropic coefficient?
- [ ] Cell_SPM::setSOC does not actually set SOC
- [ ] Should we include entropic effect in OCV or not?
- [ ] Qrev definition difference between old and new models

### Module Behavior
- [ ] Module requires number of cells to construct (unnecessary?)
- [ ] Battery class distributes current equally - consider Schimpe-style optimization
- [ ] redistributeCurrent() PI Control causes high current error
- [ ] Why do we need getRtot()? Maybe getThevenin better
- [ ] Cell_ECM: time step should be less than smallest tau or oscillation

### Testing Issues
- [ ] In `test_specificDeg` in `Procedure_test.cpp` there is a weak cell
- [ ] std::span<double>& for set states doesn't work when vector given

### Documentation
- [ ] Renew MATLAB scripts to read results
- [x] Documentation of v4 features (tested P8-G5 installation, quickstarts, compatibility, and extension guides)
- [ ] Automatise SOC -> OCV conversion

---

## Code Review Items (From Martin Robinson)

- [ ] Config file model: have a config file, model file, and options
- [ ] Class naming: standard format, Capital letter Camel case
- [ ] Template projects in OxRSE template-project-cpp
- [ ] gui_starter_library (Jason Turner)

---

## Comparing Against slide-pack

- [ ] paperCode::thermalModel does not work for slide-pack
- [ ] Should capacity check also contain CV phase?

---

## Procedure Module
- [ ] Reduce dynamic_pointer_cast in Procedure - let polymorphism work
- [ ] Capacity checking protocol in Procedure::CheckUp removed - review
- [ ] Writing functions should be outside of Procedure
- [ ] Markers to indicate end of sections deleted - review

---

## From SLIDE v2 (Legacy)

- [ ] Convert strings in function parameters to string references
- [ ] Inline small functions
- [ ] Why does getstates call setstates twice?
- [ ] BasicCycler.cpp is very large - consider splitting
- [ ] Do not use exceptions for normal system operation
- [ ] Cycler follows patterns from old BasicCycler

---

## JOSS Publication

- [ ] JOSS folder and GitHub workflow - complete submission

---

*See [.claude/discussions.md](../.claude/discussions.md) for design decision history*
*See [COMPLETED.md](COMPLETED.md) for archived completed items*
