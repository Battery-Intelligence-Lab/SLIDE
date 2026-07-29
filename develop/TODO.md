# SLIDE Development TODO

> **Note**: This is a living document. See [COMPLETED.md](COMPLETED.md) for archived completed items.
>
> **Disclaimer**: Some items are informal notes. Priority may shift based on user needs.
>
> **2026-07-12 — v5 refactor active.** The authoritative vision/roadmap is [/PLAN.md](../PLAN.md)
> (architecture, decision log, phased gates, open questions). Phases 0–8 and 9A are complete;
> Phase 9B systematic adversarial bug-hunt is complete; its Experiment transaction/parser,
> BPX parser/allocation, recorder, pack-solver/stepper/thermal, byte-shuffle, and compiled-curve index defects are fixed and recorded in the Phase-9B ledger. Historical "Critical Bugs"
> below may be superseded by PLAN.md. Update PLAN.md §8, not just this file.
> M0.7/9C-3 cold-file splitting is complete (see CHANGELOG and `.claude/reports/p9c3-cold-file-split-validation-2026-07-12.md`).
> M0.8/9C-4 is complete: `tests/support/CoreSpmTestHarness.hpp` is the single test-scaffolding idiom (MC-4) for successful build, unit-tagged observation, constant-current traces, and voltage-error metrics; nine binaries are migrated and construction/allocation/topology/backend subjects keep direct factory calls by allowlist. Native Debug, fast-math Release, and CUDA each pass 57/57; all 43 pre-existing binaries meet or exceed their frozen assertion floors; fifteen registered mutations turn a gate red. Recorded honestly: M0.8 did NOT reduce lines (see the MC-1 item below), and two preregistered mutations were falsified — the aggregate floor cannot see a single removed assertion (floors re-frozen), and the allocation-window rule lives in the allocation binaries, not the structural gate. Report: `.claude/reports/p9c4-shared-test-harness-validation-2026-07-13.md`.
> M0.9/9C-5 is complete: MC-5 is now enforced, not just labelled. `SpmFactory.hpp` no longer includes `SpmPipeline.hpp` — the parameter blocks and row layouts it stored by value became value-type headers (`{Sei,SurfaceCrack,Lam,LithiumPlating}Params.hpp`, `SpmBatchLayout.hpp`, `AgeingModelMask.hpp`), so a public-consumer translation unit fell from 26 core headers (eight of them kernels) to 21 with none. Every core header carries an `@surface` tier (api/support/internal) and `tests/structural/p9c5_public_surface.cmake` forbids any api or support header from including an internal one, restricts bindings/docs to api headers, and pins each api header's declaration count and public types. Naming: one name per concept (`n_lanes()`/`n_rows()`, `deviceAllocationCount`/`deviceWideSynchronizationCount`/`deviceArenaBytes`, `checked_lane_count`/`checked_shape`; `validate…`→Status, `valid…`/`is…`→bool, no `get` prefix). All 53 binaries digit-identical in Debug and in a like-for-like Release baseline; 57/57 Debug/Release/CUDA; 13 registered mutations turn a gate red (two real gate holes were found by the battery and fixed). Report: `.claude/reports/p9c5-public-surface-validation-2026-07-13.md`. M0.10/9C-6 dead-code and line-debt sweep is next.
>
> **2026-07-14 — PLAN extended to THE superior-stack goal.** Volkan's 16-point directive + the
> integrator-research follow-up are absorbed into PLAN.md: goals 15–23, validation contract §1.4
> (VC-1..VC-5), model-family tiers §3.23 (ECM core tier, semi-empirical/SimSES, lead-acid Schiffer —
> Q16 overturned by D-43 — and a storage-system tier), `slide.pybamm` drop-in shim §3.24, integrator &
> predictive-acceleration research §3.25, D-43..D-49, ladder extended M13–M19 with `v6.0.0` at M19.
> The M0.10 UNMET test-file line debt now has its owner: **M1.0** is the next unticked box. An
> independent orthogonal review of the revision (12 findings, all applied) is at
> `.claude/reports/plan-revision-orthogonal-review-2026-07-14.md`.
> Same day, the **Newman instrumentation & design wave** was added: §3.26 + D-50..D-52 + new M18
> (thermodynamic identity gates, EIS by analytic linearisation, Ragone/sizing optimisation with exact
> forward-sensitivity Jacobians, Jacobian service, ICA/DVA, (opt) PSD tier), each box shipping
> oracle-first unit tests per MC-2/MC-4; former M18/M19 renumbered M19/M20 — `v6.0.0` is now M20.
>
> **2026-07-21 — code-quality pass (Debug only).** A nine-way disjoint review of `src/core` with an
> adversarial verifier per finding raised 117, of which 68 survived. Eight landed (PC-10 Arrhenius
> duplication; `SLIDE_ROOT_DIR` never defined so binaries depended on their launch directory;
> `project_warnings` never linked to `slide_core` or the unit tests; a shadowed `derivative`; Boost as
> `SYSTEM`; `PathVar` globals `inline const`; plus real oracles for the byte-shuffle codec and
> `parseValue`'s exponent branch, both of which were previously untestable-by-construction).
> **60 verified findings remain unapplied and need owners** — see §5 of
> `.claude/reports/code-quality-pass-2026-07-21.md`; the full machine-readable set is
> `.claude/reports/code-quality-pass-2026-07-21-survivors.json`. Highlights: O(cells²) ladder detection,
> a line-for-line duplicated topology derivation, `PackStepper::substeps` not subdividing `dt`,
> `crc32` defined twice in the recording subsystem, a fast-math-unsafe FP validity gate.
> **Release and CUDA were not re-run for this pass.**
>
> **Stale build artifacts (found during M0.9, worth cleaning):** `build-release/bin/Release/` still holds `*_fast`/`*_ipo` probe binaries dated 2026-07-12 that are NOT in the current Ninja build graph. One of them fails an assertion. They are not evidence about current code — the live Release lane already compiles with `-ffast-math` and passes that test — but they will keep poisoning any glob-based measurement until the directory is reconfigured from scratch.
>
> **2026-07-23 — MQ quality wave inserted; AGENTS.md created.** PLAN.md now carries milestone **MQ**
> between M0 and M1 (Volkan's clean/re-derive/hunt directive): MQ.1 restore the Debug-AND-Release
> baseline, MQ.2 disposition all 68 quality-pass survivors (the "60 unapplied" above minus what landed
> since: linear ladder detection in 452ce0c; `substeps * dt` semantics are retained and documented,
> with the executable equivalence gate still pending — the report's landed list is authoritative),
> MQ.3 repo hygiene (the stale-build
> note above becomes a deletion PROPOSAL there), MQ.4 derivation inventory `docs/derivations/INDEX.md`,
> MQ.5 independent re-derivation of the mathematics, MQ.6 logic hunt, MQ.7 structural performance
> hunt, MQ.8 closeout. `AGENTS.md` at the repo root is the standing operating contract for autonomous
> Codex sessions; PLAN.md remains the single source of truth. Next unticked box: **MQ.1**.
>
> **2026-07-23 — MQ.1 baseline restored.** Fresh tree-local native Debug/ThinLTO,
> fast-math Release, and CUDA Release with host-C++ ThinLTO each discover and pass 58/58; CUDA executes
> 433,671 assertions / 4 cases and all three repo-root direct-CWD witnesses pass 57 / 3.
> No source or test changed. The prior `cl.exe`-on-`PATH` CUDA shorthand was falsified twice;
> fresh CUDA configuration requires the complete x64 `vcvars64.bat` environment. Known test gaps
> and the two observed build-policy warnings pass to MQ.2; no performance/sanitizer/coverage/CI/
> package/cross-platform claim is made. Report:
> `.claude/reports/mq1-baseline-validation-2026-07-23.md`. Next unticked box: **MQ.2**.
>
> **2026-07-24 — MQ.2 P0/P1/P2 reconciled and landed.** The four flagged commits were
> registration-only. The
> +313-line exact PackSolver oracle boundary landed at `6c7ab74`; Debug, fast-math Release, and
> host-C++/CUDA-tree binaries each pass 949/30, and four division-site plus two oracle mutations
> turn it red. P1's three pack-algebra findings are APPLIED at `8630207` under the unchanged
> three-configuration exact gate and five adversarial mutations. P2's Mode-C scratch finding is
> APPLIED at `bf185dd`: seven roles, one transactional publication, six unchanged hot fills, and
> six adversarial mutations. Current evidence census over 71 decisions: 21 APPLIED, 0 REFUTED,
> 8 named deferrals, 42 pending. The PackSolver test-file debt is now 1,510 lines and remains
> explicitly owned by M1.0.
>
> **2026-07-29 — MQ.2 S1 PackStepper ownership landed.** Commit `62f1587` applies the shared
> gather/scatter and prefix-uniqueness owners and pins `substeps=N` as N full-`dt` advances under
> one frozen electrical solve and thermal assembly. Focused Debug, fast-math Release, and
> host-C++/CUDA-tree gates pass PackStepper 230/9, PackSolver 949/30, allocation 14/2, and
> structural 1/1; all three restored trees also pass the full 58/58 suite. Nine mutations turn
> red and exact hashes were restored. Census: 24 APPLIED, 0 REFUTED, 8 named deferrals, 39 pending
> among the original 71.
>
> **2026-07-29 — MQ.2 S1.1 moved-owner hardening landed.** Commit `574cfaf` adds exhaustive
> no-throw, self-guarded ownership transfer for `SolverWorkspace`, `PackSolver`, and `PackStepper`.
> Focused floors PackTopology 148/9, PackSolver 984/32, PackStepper 333/13, allocation 14/2, and
> structural 1/1 pass in Debug, fast-math Release/IPO-off, and host-C++ CUDA; all three restored
> trees pass 58/58. Eight mutations turn red and exact hashes were restored. The original census
> remains 24 APPLIED / 39 pending / 8 deferrals. The post-baseline registry has four S1.1 IDs
> APPLIED and the Eigen warning DEFERRED to B1; combined: 28 APPLIED / 39 pending / 9 deferrals
> across 76 stable IDs.
>
> **2026-07-29 — MQ.2 S2 old-production boundary frozen.** The exact ten-assertion scripted case
> at `4a9d29a` plus capture `d94a0ba` reaches 9/10 against unchanged production and prints the
> abandoned hidden-step residual `2.5`; solution and warm-state restoration pass.
>
> **2026-07-29 — MQ.2 S2 source-step rollback landed.** Commit `7b3fd14` passes PackSolver
> 994/33 and all companion focused gates in Debug, fast-math Release/IPO-off, and host-C++ CUDA;
> every restored tree passes 58/58. Three independent restore mutations turn separate assertions
> red at 9/10. Original census: 26 APPLIED / 37 pending / 8 deferrals; combined 76-ID registry:
> 30 APPLIED / 37 pending / 9 deferrals.
>
> **2026-07-29 — MQ.2 T1 topology ownership landed.** Commit `ad6d785` plus structural
> hardening `30464e4` and formatting `479f55f` replaces three archetype maps with one sorted
> slot record and both branch-graph/BFS derivations with one guarded owner. PackTopology 164/11
> and all companion focused gates pass in Debug, fast-math Release/IPO-off, and host-C++ CUDA;
> all three final trees pass 58/58 and no-op rebuilds. Nine mutations turn red and exact hashes
> are recorded. Original census: 29 APPLIED / 34 pending / 8 deferrals; combined 76-ID registry:
> 33 APPLIED / 34 pending / 9 deferrals. Next is T2 thermal topology.

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
- [ ] **MC-1 test-file line debt — explicitly owned by M1.0.** M0.8 shared the test *mechanics* (MC-4) but grew the files; M0.10 swept the production side and did not touch the test side. Still oversized: `core_ParserAllocation_test.cpp` (1,229), `core_PackSolver_test.cpp` (1,780 after MQ.2 S2), `core_Experiment_test.cpp` (1,168), `core_AsyncRecorder_test.cpp` (1,118), `core_ParameterSet_test.cpp` (1,007). Any split must preserve assertion counts and exact-site coverage.
- [x] M0.10 production-side line debt: done. `SpmFactory.cpp` 713 → 494 + `SpmBatch.cpp`; `PackSolver.cpp` 910 → 610 + `PackSolverIterative.cpp`; `AsyncRecorder.cpp` 890 → 396 + `AsyncRecordingCodec.cpp` + `detail/AsyncRecordingFormat.hpp`; `Recorder.cpp` 721 → 293 + `RecordingFormat.cpp`; `SpmPipeline.hpp` 783 → 745 + `SpmDiffusionRhs.hpp`; the 1,281-line coverage reporter → the `slide_coverage` package. `checkedAdd`/`checkedMultiply` were duplicated and now have one definition. Only `SpmPipeline.hpp` (745) and `CyclerV2.cpp` (771) exceed 700, both with written justifications. Core lines grew 17,363 → 17,595 (+1.3%) — predicted before the work; there was no dead code to delete.
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
