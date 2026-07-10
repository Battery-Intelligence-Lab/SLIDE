# Changelog {#changelog}

[TOC]

This changelog contains a non-exhaustive list of new features and notable bug-fixes (not all bug-fixes will be listed). 

<br/><br/>
# Unreleased

## Added
* **Nonlinear SPM pack boundary and chord safeguards** (`src/core/SpmFactory.*`, `SpmPipeline.hpp`, `PackSolver.*`): every concrete SPM batch now exposes an allocation-free frozen-state Thevenin tangent that includes surface-flux, OCV, exchange-current, film, and ohmic current dependence. Mode A applies a true residual correction through its cached sparse Jacobian, refreshes only when contraction degrades, limits trial-current changes, and falls back to eight-stage source stepping with atomic solution rollback. A heterogeneous four-lane Kokam SPM pack converges to equal terminal voltage within `1e-10 V`; a 100-step evolving-state CC segment stays within the P2-G4 limit of ten numeric factorizations.
* **Phase-2 Thevenin pack solver foundation** (`src/core/PackSolver.*`): one indirect linearization call per archetype batch feeds a flat affine Thevenin system. Mode A stamps the compiled netlist into Eigen SparseLU with one symbolic analysis, persistent numeric-factorization reuse, warm currents, diagnostics, and explicit invalidation; detected series-of-parallel packs use allocation-free analytical Mode B layers. Affine P2-G2 current agreement is `≤1e-12 A`, nested packs solve in one global loop within two iterations, invalidated solves digit-match cold results, and 100-step heterogeneous 4p/nested segments retain one numeric factorization. Full Debug and Release suites are 34/34.
* **Phase-2 pack compile foundation** (`src/core/PackTopology.*`): value-semantic cell/series/parallel combinators flatten nested authoring trees into one electrical netlist with explicit cell/resistor branches, stable hierarchical paths, lexicographically stable archetype batch/lane addresses, canonical nodal sparsity, connectivity/index-1 metadata, and ladder eligibility. D-21 is implemented as an independent canonical cell/boundary thermal graph with sorted single-evaluation edges, CSR incident gathers, fixed-order allocation-free `q_ext` assembly, atomic cold validation, and exact pair-energy conservation. Full Debug and Release suites are 33/33.
* **D-21 compiled pack thermal-adjacency design** (`PLAN.md` §3.4/§4): closes Q9 before Phase 2 with an electrical-topology-independent cell/boundary graph, canonical sorted thermal pair list, fixed incident-order gather into arena-owned `q_ext`, single-evaluation edge fluxes, no atomics, explicit stage/rollback/snapshot semantics, and five registered thermal implementation gates.
* **PAY-1 reproducible payoff benchmark and base-batch fast path** (`benchmark/benchmark_PAY1_spm.cpp`): executes the registered 10,000-cell, 1C, one-simulated-hour workload against an equal-work legacy loop, alternates run order, reports median and conservative range-derived speedups, and refuses timing success if final state/voltage parity fails. The base isothermal archetype now avoids unused full-profile RHS observables, caches self-invalidating derived transport terms, assigns accurate ODE/algebraic row roles, fuses Euler diffusion/observation, compacts rollback copies, and exactly coalesces equal lanes after an every-step equality check with an automatic heterogeneous fallback. Full Release result across three 36,000,000-cell-step repetitions: `6.97×` median and `6.16×` conservative speedup (target `≥5×`), with `8.67e-19` maximum state error and `8.88e-16 V` voltage error. Debug and Release suites remain 32/32.
* **P1-G1 full single-cell trajectory parity gate** (`tests/parity/P1G1_spm_test.cpp`): runs the production factory, composed SPM pipeline, and `EulerLegacy` stepper in lockstep with legacy `Cell_SPM` for a 1200 s mid-SOC 1C discharge and a low-SOC steep-OCV-tail discharge. The adapter preserves the legacy Kokam reference temperature exactly (298.0 K, rather than silently normalising it to 298.15 K). Maximum terminal-voltage error is `4.44e-16 V` Debug / `8.88e-16 V` Release; every mapped physical and cumulative state satisfies the registered `1e-15 + 1e-12*scale` band. Full suites are 32/32.
* **P1-G4 bitwise restart gate** (`tests/unit/core_P1G4_restart_test.cpp`): a coupled thermal + SEI + stress-driven crack + stress-driven LAM + lithium-plating batch is compared against a split run that snapshots the full padded arena halfway, destroys all pipeline/stepper scratch, cold-builds replacements, restores only arena bytes, and resumes. Final arena state and terminal voltages are byte-for-byte identical in Debug and Release, proving no later-step physics depends on hidden non-arena state. Full suites are 31/31.
* **P1-G2 zero-allocation batch gate** (`tests/unit/core_P1G2_allocation_test.cpp`): replaces all global scalar/array/aligned allocation forms, cold-builds and warms a 10,000-lane SPM batch, then proves an accepted `EulerLegacy` step moves the allocation counter by exactly zero in Debug and Release. The same gate measures 240 bytes/cell of arena state, within the ≤300-byte PC-6 contract. Full Debug and Release suites are 30/30.
* **`EulerLegacy` and Phase-1 `Simulation` façade** (`src/core/EulerLegacy.hpp`, `src/core/Simulation.*`): the allocation-free parity stepper evaluates the concrete batch RHS, advances only rows marked `ode`, leaves algebraic/input/cumulative rows out of numerical integration, computes accepted-state terminal voltage through the shared observable path, and then updates arena-owned elapsed time, Ah, and Wh exactly once outside the RHS. A once-allocated full-state backup makes failed accepted-state validation atomic, and stress history is stored from the correct pre-step stress after acceptance. The minimal façade builds through the registry and solves constant-current experiments, including event-aligned partial final steps and epsilon-safe decimal step ratios, into sample-major multi-lane voltage traces. Full Debug suite is 29/29; affected Release tests pass.
* **D-02 SPM composition registry and cold factory** (`src/core/SpmFactory.*`): runtime `SpmModelOptions` now select one of 12 explicitly compiled batch archetypes—`nch={5,8,12}` × `{isothermal, thermal, isothermal+ageing, thermal+ageing}`—once at build time. `SpmBatch::rhs` performs one indirect call per batch into a fully concrete pipeline with no per-cell dispatch. Individual SEI/crack/LAM/plating selections remain batch-level masks and disabled stages are skipped before lane loops. The factory flattens `CellDesign`, curves, spectral coefficients, thermal constants, mechanism parameters, initial state, row roles, and restart-safe stress history; malformed geometry, curves, masks, modifiers, lane counts, or mechanism parameters return `Status` without replacing an existing valid batch. Registry selection and generic RHS boundaries pass Debug and Release.
* **Validated per-batch spectral model compiler** (`src/core/SpectralModel.hpp`): replaces the v4 path’s dependence on the global `Model_SPM<>` singleton and accepts physical particle radii for each batch. Cold builds now `Status`-fail on invalid geometry, Eigen solver failure, complex contamination, singular modal transforms, non-finite coefficients, or disagreement with the analytic `tan(μ)=μ` spectrum at the registered `nch={5,8,12}` bands. Default-geometry A/B/C/D coefficients, state transform, centre map, and stress integration matrix are exactly identical to legacy in Debug and Release; custom-radius scaling and atomic failure are covered. This closes audit items C1–C3.
* **Compile-time composed SPM RHS pipeline** (`src/core/SpmPipeline.hpp`): one fixed batch archetype now enforces the D-23 evaluation order—clear `ydot`, reconstruct concentration/electrical observables once, optionally reconstruct stress once, then accumulate diffusion, thermal, SEI, surface-crack, LAM, and lithium-plating contributions. Effective diffusivity and molar flux are shared between observable and diffusion stages instead of recomputed. Enabled scratch buffers are allocated in the cold constructor while disabled mechanisms allocate no workspace; all optional branches are forced through compilation, and the registered gate proves full derivative zeroing plus adaptive trial-vector rebinding in Debug and Release. Full Debug suite is 26/26.
* **Yang lithium-plating ageing** (`src/core/LithiumPlating.hpp`): scalar-generic Tafel current matches direct `Cell_SPM::LiPlating` calls within `1e-13` at charge and discharge points in Debug and Release. Its additive RHS updates negative diffusion modes, lost lithium, and the new plated-lithium thickness arena row. Together with SEI, cracking, stress, and LAM, every legacy SPM ageing mechanism now has a v4 kernel.
* **LAM ageing mechanisms 1–4** (`src/core/Lam.hpp`): ports Dai stress-driven electrode thinning, Delacourt–Safari flux-driven active-fraction loss, Kindermann NMC dissolution, and Narayanrao proportional active-area loss. The additive RHS composes direct area loss with `3ε/R`, and every raw positive/negative thickness, area, and active-fraction rate matches direct `Cell_SPM::LAM` calls within `1e-12` in Debug and Release. Full Debug suite is 24/24.
* **Surface-crack ageing mechanisms 1–5** (`src/core/SurfaceCrack.hpp`): ports Laresgoiti stress, Dai stress, Deshpande–Bernardi concentration-gradient, Barai throughput, and Ekstrom side-reaction models with batch-level masks. Optional graphite diffusivity loss preserves scalar sensitivities below its legacy cap. The additive RHS owns crack surface, extra crack-driven SEI flux/lithium loss, and negative diffusivity changes. Every mechanism and both diffusion settings match direct `Cell_SPM::CS` calls within `1e-12` in Debug and Release; full Debug suite is 23/23.
* **Shared SPM stress observables and restart-safe history** (`src/core/SpmStress.hpp`): allocation-free Dai hydrostatic maxima and Laresgoiti graphite stress are reconstructed from the same full concentration observable used by ageing/recording. Both match the public legacy stress functions within `1e-12` in Debug and Release. Previous Dai/Laresgoiti values and their interval are explicit algebraic arena rows, so adaptive integrators do not advance them and checkpoints no longer lose the hidden `s_*_prev` memory required by crack/LAM laws.
* **SEI ageing mechanisms 1–4** (`src/core/Sei.hpp`, PLAN.md Phase 1): scalar-generic batch kernels cover the legacy kinetics-limited, linear-diffusion-limited, Christensen–Newman, and fitted SEI variants, including optional Ashwin porosity loss. Batch-level model masks branch outside lane sweeps; a once-allocated scratch feeds an additive RHS that updates negative diffusion modes, SEI thickness, lost lithium, negative active fraction, and active area. All four mechanisms match direct `Cell_SPM::SEI` calls within `1e-13` relative error in Debug and Release. Lost-lithium and active-fraction state are now explicit arena rows.
* **`ThermalLumped` core RHS** (`src/core/ThermalLumped.hpp`, PLAN.md Phase 1): cold-compiles `ThermalDesign` to heat capacity and environmental conductance, then accumulates `dT/dt = (Q_internal + q_ext + hA(T_env-T))/C_th` across whole SoA batches without allocation. Internal heat comes from the shared electrical observable stage and `q_ext` is the reserved pack-coupling input. Generated heat energy and thermal elapsed time are arena rows rather than legacy hidden object members, making snapshots self-contained. Heterogeneous-lane analytic balance, additive-RHS, invalid-parameter/state, and Debug/Release tests are green.
* **Shared SPM electrical observable stage** (`src/core/SpmObservables.hpp`, PLAN.md §3.12/D-23): one scalar-generic, allocation-free pipeline now reconstructs particle concentrations and computes surface stoichiometry, Butler–Volmer exchange current/overpotential, electrode and cell OCV, entropic coefficients, evolving resistance, terminal voltage, and reversible/reaction/ohmic heat for every lane. Registered `{0 A @ 298.15 K, +20 A @ 298.15 K, -10 A @ 310 K}` cases match legacy `Cell_SPM` OCV and resistance exactly; terminal-voltage error is at most `4.441e-16` in Debug and zero in Release. Diffusivity, electrode thickness, active area, and resistance live in the canonical arena layout so degradation and arbitrary integrator trial vectors cannot observe stale frozen parameters. Invalid surface concentrations use the common `Status` channel. Full Debug suite 19/19 and affected Release tests 4/4 green.
* **Compiled parameter curves** (`src/core/CompiledCurve.hpp`, PLAN.md §3.11/D-16): `IndexedPiecewiseLinear` preserves nonuniform measured knots exactly while replacing per-evaluation binary search with a uniform O(1) segment-index accelerator; it matches legacy Kokam graphite OCV interpolation bit-for-bit at registered non-knot queries. `UniformLut` provides the planned 4096-point canonical form for smooth injected functions and rejects builds that miss the registered relative-error gate. Invalid parameters and numerical build failures now use the common `Status` channel.
* **First-class electrode/cell description hierarchy** (`src/core/CellDesign.hpp`, PLAN.md §3.3): value-semantic `ActiveMaterial`, `ElectrodeDesign`, `SeparatorDesign`, `ElectrolyteDesign`, `ThermalDesign`, and `CellDesign` types now model the physical battery on the cold path, with compact `ElectrodeParams` for compiled kernels. The scoped v4 `Domain` has the plan's canonical negative-first order and cannot be mixed implicitly with legacy's positive-first enum; diffusion and concentration bridges now map explicitly, preventing silent electrode swaps as the factory grows.
* **P1-G3 (part b) — Carslaw–Jaeger/Crank constant-flux transient oracle** (`tests/unit/core_ChebyshevTransient_test.cpp`): validates the Chebyshev A/B/C/D maps together at `nch={5,8,12}` against the analytic spherical-diffusion series, including the centre `Cc/cc_coeff` path that carried the historical `nch!=5` sign bug. The registered relative-error band `<=1e-6` holds at every surface/interior/centre node for both electrodes: worst errors at dimensionless time `τ=0.2` are `8.254e-7`, `2.356e-11`, and `3.949e-12`; at `τ=1.0` all are `<=4.661e-13` across Debug/Release. The test explicitly documents SLIDE's flux sign convention: positive molar flux depletes the particle (`D∂c/∂r=-j`). Together with part (a), P1-G3 is complete.
* **Phase-1 rebindable RHS/observable contract** (`src/core/BatchView.hpp`, PLAN.md §3.12/D-23): fixed batch geometry is now separated from per-evaluation trial-vector pointers through `BatchView`; `RhsViews` rebinds state/derivative vectors without copying or allocation and enforces the mandatory derivative-zero stage before accumulating component kernels. `StepCtx` carries segment-frozen applied-current inputs, while `BatchBuilder` records per-row ODE/algebraic/cumulative/input roles so generic integrators can exclude non-ODE rows. The shared scalar-generic SPM concentration kernel reconstructs surface, all interior nodes, and centre from arbitrary trial vectors. It matches legacy `Cell_SPM::getC` exactly in Debug and to `2.665e-15` relative error in Release, and recovers heterogeneous uniform profiles at every node to `9.027e-15`. Arena move semantics and slice/lane bounds were hardened alongside the view layer.
* **P1-G3 (part a) — analytic Chebyshev eigenvalue oracle** (`tests/unit/core_ChebyshevEigenvalues_test.cpp`): independent-mathematics validation of the `Model_SPM` spectral operator. The folded, surface-condensed dimensionless diffusion operator has eigenvalues `λ_k·R² = −μ_k²` where `μ_k` are the roots of `tan μ = μ` (μ₀ = 0 is the mass mode — why exactly one zero eigenvalue is forced). Tests the operator alone (no trajectory, no output map). Confirmed: `Model_SPM` reproduces the analytic spectrum (μ₁ to 3.6e-14 at nch=12), exactly one zero eigenvalue, all others real and negative. Convergence is the true Chebyshev spectral curve: fundamental rel err 2.0e-5 / 2.3e-10 / 3.6e-14 at nch = 5 / 8 / 12; ≈⌈nch/2⌉ modes resolved to 1e-3. Quantifies for the first time that at the default nch=5 the fundamental diffusion eigenvalue is accurate to only ~2e-5 (fine for voltage). Note: this falsified the plan's initially-assumed flat 1e-10 band; the test enforces the measured spectral-convergence bands instead.
* **v4 core keystone (`slide::core`, `src/core/`)**: `StateArena` (contiguous 64-byte-aligned variable-major SoA state storage; snapshot/restore = single `memcpy`) and `BatchBuilder` (build-time state layout from component `StateSpec` declarations; hands out integer `StateSlice` handles; reserves the cross-batch thermal-flux seam `q_ext`). Header-only, not linked into the legacy library (strangler migration, PLAN.md §5). Tests include allocation-counting groundwork for the zero-per-step-allocation gate (P1-G2).
* **Production diffusion kernel `slide::core::SpectralDiffusion<NCH>`** (`src/core/SpectralDiffusion.hpp`): vectorised-across-lanes forward-Euler solid diffusion for one `CellBatch` archetype, sweeping whole SoA `StateArena` rows (one variable across all lanes). Batch-shared model/geometry in `DiffusionParams<NCH>`; per-lane temperature and current density as inputs; a per-lane scratch (`D_eff`, `flux`) allocated once at construction (zero per-step heap allocation). Legacy-Euler stepping mode (the exponential modal propagator is deferred to Phase 3). Validated against the legacy-shaped oracle on a heterogeneous 8-lane batch: bit-identical in Debug (proving the vectorised sweep is the same math), rel-drift 3.8e-15 ≤ 1e-12 in Release.
* **P1-G0 parity-drift pilot**: `SpectralDiffusionLegacyKernel` (op-order-faithful replica of the legacy `Cell_SPM` forward-Euler diffusion update on `StateArena` rows) + lockstep pilot test. Debug/-O0: **bit-identical** (drift exactly 0 over 1200 steps). Release/-O3 re-confirm (`build-release`, clang-21): drift max_abs 2.26e-17, **max_rel 5.63e-15** — exact bit-identity is lost to cross-TU FMA contraction under `-Ofast`, but the registered decisive band (rel ≤ 1e-12) holds with ~3 orders of margin, so the v4 parity gates keep the 1e-12 digit-diff (PLAN.md §7 Q8, fully closed). The exact-zero check is scoped to Debug; the decisive relative band gates CI in both configs.

## Fixed
* **CMake object-source propagation**: legacy implementation `.cpp` files were declared `PUBLIC` on object libraries and those object targets were exported publicly through `src`, causing every test/application to compile the same Cell_SPM/module/cooling/procedure sources again. Implementations are now private, `src` archives each object exactly once (including the nested `Cell_SPM` target), and public include paths are explicit. Clean Debug and Release trees build and all 28 tests pass with the corrected link graph.
* **NaN/Inf validation remains effective under Release fast-math** (`src/core/Numeric.hpp`): clang can fold `std::isfinite` to true under the project’s `-Ofast` flags, which caused Release-only acceptance of invalid compiled parameters. Core validation now classifies IEEE-754 exponent bits behind a deliberate integer compiler barrier; bit-injected NaN regression tests cover thermal parameters and curve tables in both build types.
* **`FixedData.hpp` is self-contained**: it now includes `<iostream>` and `<iterator>` for the `std::cout` and iterator tags it uses, instead of compiling only when unrelated headers happened to provide them transitively.
* **Build restored under clang ≥ 21**: `fmt` bumped 11.0.2 → 11.2.0; 11.0.2's consteval format-string checking rejects its own internal format calls (`format-inl.h`/`os.cc`), so the tree did not compile at all. (P0-C1)
* **Default `Cell_ECM`/`Cell_Bucket` OCV is no longer −55 782 V.** The default `ocv_coefs` polynomial was a truncated 3-term fit that evaluates to ≈ −55 782 V at every SOC, so default ECM/Bucket cells (and any `Module_p` built from them) were unusable. `getOCV()` now interpolates the cell's SOC-voltage table by default (restoring the pre-Aug-2024, 2.7–4.2 V behaviour); the polynomial path remains available by supplying a fit via `set_ocv_coefs` (analytical parallel-module work). Unit-test voltage expectations updated 3.15 V → 3.45 V: the 3.15 V values dated from the deleted standalone `Cell_Bucket` class whose dummy OCV ramped 2.0–4.3 V; the alias `Cell_ECM<0>`'s table ramps 2.7–4.2 V, giving 3.45 V at SOC 0.5. (P0-C2)
* **`Cell_SPM` unit test expectations updated to the current electrode-thickness calibration** (`thickp` 86.87357 µm, `thickn` 74.883947 µm, and the crack-surface value derived from `thickn`). The values were recalibrated in Aug 2024 ("improved SOC estimation for SPM()", 147ee4c); the test still asserted the pre-2024 70/73.5 µm values. Test-only change. (P0-C3)
* **`Cycler::setCurrent`'s dead `v_now` out-parameter removed** (private API; it was never assigned, and its only caller ignored it). (P0-C5)
* **`Cycler::CV` throughput**: energy now integrates with the trapezoid rule `(v_before + v_after)/2` per step — consistent with the `Cycler::CC` fix (B3) — instead of the end-of-step voltage only, and `th.time()` is now accumulated (it was never incremented in CV, so CV/CCCV time throughput read 0 s). (P0-C6)
* **`Module::setStates` rollback** now restores each child to *its own* slice of the pre-set states. Previously the rollback loop re-applied `SUs[i]` for every index and never advanced the per-child offset, so when a child's state was rejected the already-set children were left corrupted. (bug A1)
* **Parallel-module current solver (`Module_p::setCurrent`) is now instance/thread-safe.** The solver's working vectors and Jacobian matrix were function-`static`, i.e. shared across all `Module_p` instances and threads. This produced wrong branch currents, out-of-bounds access (Eigen assertion / crash) when modules had different child counts, and data races under the parallel time-step fan-out. They are now per-call locals. (bug A2)
* **Parallel-solver Jacobian is refactorised every iteration.** The LU factorisation was computed once from the initial resistance estimate and never refreshed even though the secant estimates update each iteration, degrading Newton to a slow/failing chord iteration. The Jacobian is now rebuilt from the updated estimates and refactorised each iteration. (bug A3)
* **Contact-resistance heat (`Module_p` `Qcontact`) no longer double-counts.** The cumulative branch-current accumulator was declared outside the per-resistor loop and never reset, so the contact heat grew quadratically. Each contact resistor now correctly sees the summed current of the cells behind it. (bug A4)
* **`Module_p::setVoltage`** no longer keeps its solver vectors/matrix in function-`static` storage (same instance/thread-sharing hazard as A2). (bug A2)
* **Analytical parallel-branch solver** returns `Status::Invalid_SUs` instead of dereferencing an unchecked `dynamic_cast` when a child is not a `Cell_ECM<1>`. (bug A5)
* **Parallel-solver early-exit residual** uses a proper L∞ norm (`max`) instead of a malformed running sum of maxes. (bug A7)
* **Time-series cell data storage (`CellDataStorage<storeTimeData>`) now appends** each timestep instead of dropping history: the ill-formed `data.assign(data.end(), {...})` became `data.insert(data.end(), {...})`. The `#include "CellDataWriter.hpp"` was also moved out of `namespace slide` to global scope — the mid-namespace include pulled standard headers into `namespace slide` and nested `slide::slide::CellDataWriter`, so the header could not compile when included. (bug B1)
* **`Histogram::add()` no longer writes out of bounds** on a default-constructed histogram (empty bins): it returns early when there are no bins. Affects `slide::EmptyHistogram` and any histogram used before `initialise()`. (bug B2)
* **`Cycler::CC` energy throughput (`Wh`) is now actually accumulated**, integrating with the trapezoid rule `(v_before + v_after)/2` per step. Previously it multiplied by a `vi` sample that `Cycler::setCurrent` never assigned, so CC energy throughput was always 0 Wh. (bug B3)
* **SEI degradation-model documentation aligned with the code**: the unreachable-model message now says ids "0 to 4" (cases 0–4 exist); the inline literature references for SEI ids 1 and 2 were swapped relative to the formulas (id 1 = kinetics-limited/Ning & Popov, id 2 = kinetics + SEI-layer diffusion/Pinson & Bazant); `DEG_ID` now documents id 4. No behavioural change. (bug B4)
* **Degradation forward-Euler invariant documented and debug-guarded** in `Cell_SPM::timeStep_CC`: the loop integrates every state index (including the algebraic I/V slots) and is only correct because degradation leaves those derivatives at 0 — now asserted in debug builds. Digit-identical behaviour. (bug B5)

<br/><br/>
# SLIDE v3.0.0 (aka slide-pack merged into SLIDE)

## New features and important updates
* We now require `CMake 3.21` and `C++20` capable compiler.
* Chebyshev discretisation is moved from MATLAB to C++ code. 
* Unit testing is added in continuous development cycle. 
* SOC calculation in `Cell_SPM` is changed with lithium fractions instead of coloumb counting. 
* SLIDE and SLIDE-pack has different LAM parameters as they use different set of fitting data.
* `copy()` method is added to all classes to clone the class. 
* `determine_OCV` functions now use absolute error which gives a better estimation of initial lithium fractions. 
* `Electrode_SPM` class is created.
* As we have several dependencies now, we print the licenses of these packages in `bin/third_party.txt`. 

### Battery pack simulation support is added
* SLIDE is merged with slide-pack. Now, SLIDE is now now capable of simulation large battery packs.

### Performance updates
* `slide::Model_SPM` -> `slide::Model_SPM*`  and `makeModel()` speeded up the performance from 34 seconds to 12 seconds. 18900x326 double = 47 MB RAM is also saved.
* Removed file reading when individually creating cells! 30 seconds speed-up for creating EPFL system.
* `DEG_ID` variables are changed from `int` to `uint_fast8_t` -> 144 byte to 36 byte size reduction.

### Other updates
* `getDaiStress` is simplified by removing unnecessary R multiplication and division. 
* More Doxygen comments are added. 
* `paperCode.hpp` is deleted. `paperCode.cpp` is made an executable.
* Different battery cell types: Cell_Bucket, Cell_ECM and Cell_PbA are added. 
* [CHANGELOG.md](https://github.com/davidhowey/SLIDE/blob/master/CHANGELOG.md) is added. 
* `*.CMake` files are added for a better build. 
* CPack support for unit tests are added. 
* `slide::Clock` class is created. `<chrono>` for timing is adapted for correct timing on Mac. 
* `Cell_KokamNMC.cpp` is deleted and it became an header. 
* DEG_ID improved, zero is not counted anymore.
* Free functions to simplify other classes are added: `free::check_Cell_states`
* Catch outside `LiPlating` is removed since it is not needed to be handled since it only throws when id is wrong. 
* `slide.hpp` for including the library is added. 
* `NULL` -> `nullptr`
* `shared_ptr` classes are turned to `unique_ptr`. Copying `shared_ptr` is eliminated. 
* print arguments for some functions are removed (e.g., `V(bool print), getOCV(bool print)`)
* `cube` and `sqr` functions added for utility. 
* `getIndex` is removed
* remainder for integers are eliminated. 
* `i` in Cycler::rest is removed. 
* `CONST_PI` constant defined. 
* `double Rdc` is removed from `Cell.hpp`
* Free functions to call member functions under `free` namespace.
* `develop` folder is added for developer-related matters.
* Unnecessary printing statements with `prdet` is removed to reduce cluttering. 


### slide_pack changes: 
* `(verbose_gl > v_noncrit_gl)`  ->  `(settings::verbose >= printLevel::printNonCrit)`
* `(verbose_gl >= v_crit_gl)`   -> `(settings::verbose >= printLevel::printCrit)`
* Module functions are combined.
* `#if TIMING` is removed, profiler should be used if needed. 
* `StorateUnit` -> `StorageUnit`

### C++20 upgrades
* `std::span` for state assignments. `XYdata_ss` for span. 

## Notable Bug-fixes
* `if(abs(Ii < 0.1))` in Cycler. 
* Typo `bockDegAndTherm` is fixed. 
* `Cell_ECM` had `-` in the equation, corrected

### Small bug fixes:
* `dt / 3600.0;` should be `dti / 3600.0;` and converted. 

## API changes
* `TIME_INF` variable is created for no time limit in cycling. 
* Most functions now use smart pointers. 
* Template / auto parameter deduction is used. 
* `Cell_Bucket` class is removed now `Cell_ECM` class with `0` RC pairs can be used as `Cell_Bucket`.
* Literal operators for units such as `degC` are added so you can write `25.0_degC` instead of `Kelvin + 25.0`.
* `io` namespace is added. 
* `Cell_ECM::R_C_pair` and `Cell_ECM::R_Tau_pair` classes are added for passing parameters. 
* isCharging() -> I() < 0  and isDischarging() -> I() > 0 are added to enforce consistent current direction representation. 
* `Status` class is created to hold error codes. They have member functions like  `status.good()` to test if the operation successful. 
* `setI_iterative` is removed. 
* Unnecessary `get` is removed from getters: `getV()` ->` V()`,  `getT()` -> `T()`, `getVcheck()` -> `checkV()`.
* `Geometry_SPM` class is created. 
* `operator[]` is added to `Module` class to reach SUs. 



## Dependencies
* `Eigen` library is added for matrix operations. 
* `Catch2` library is added for testing and.
* We require at least CMake 3.21.
* Now a compiler with C++20 support is required. 

## Build system / developer changes
* `CPM` package manager is added.
* `OBJECT` libraries are being used for subfolders instead of `STATIC` libraries. 
* Added many warnings. 
* `CTest` for testing. 
* Subfolders now includes linker options. 
* Added `CPMLicenses.cmake` to print the licenses of third-party packages in `bin/third_party.txt`. 


<br/><br/>
# SLIDE v2

Modern version of slide. ([\#8](https://github.com/davidhowey/SLIDE/pull/8)) Deprecated but you can find a copy of Slide V2 in [this branch](https://github.com/davidhowey/SLIDE/tree/SLIDE_v2).

## New features / updates
- Most of the functions, especially `estimateOCVparameters()` should be faster now. Unnecessary file reading inside classes is removed. `fitAMnAndStartingPoints` was calling readOCVinput all the time. It is fixed. 
- SLIDE should be able to run on different platfroms/compilers through CMake.
- Name capitalisation bugs are fixed. 
- Documentation is updated. Creating wiki page from `*.pdf`s. 
- False position method for current finding is implemented for SLIDE. 
- For formatting code `clang-format` is adopted. 
- `std::endl` are replaced by `"\n"`. 
- Auxillary functions and data types are added (such as `fixed_data`, `XYdata`, `linspace`, `logspace`, `linstep`, etc.).
- Parallelisation is wrapped into `slide::run` function and `<thread>` library is started being used. 
- File structure is changed, `data`, `results`, and `docs` folders are created. 
- `loadCSV_1col` can take `std::array` or `std::vector` as argument now. It also reuses `loadCSV_mat`
- Binary search (with `std::lower_bound`) for `linInt` is implemented. Detecting if the data is fixed step. Remove linear search from interpolation! O(n) -> O(log(n)) or O(1) for fixed-step data. Nan initialization for interpolation is removed. 
- Used linspace, logspace when searching for parameters! 
- Remove empty destructors etc. Rule of zero. 
- changed pow to scientific notation  `7*pow(10,-10)` ->  `7e-10`. 
- Try-catch block based program flow is mostly removed due to performance penalty. 
- Model Cp, Cn, Vn, Vp, Q are 2D std::arrays, and shorter code. 
- OCVcurves is being moved to a struct. 
- State default constructor is removed since it sets all states to zero, thanks to uniform initialiser. `State() = default`
- State `setStates(State& s)`, `State::setIniStates` are removed. 
- name inputs are deleted from `fitAMnAndStartingPoints`.
- ValidOCV, calculateError, readOCVinput, discharge, fitAMnAndStartingPoints size input removed. 
- `linInt_noexcept`, `discharge_noexcept`, `writeCharacterisationParam`,  `linstep_fix and`, `logstep_fix`, `calcSurfaceConcentration`, `calcDiffConstant`, `calcMolarFlux`, `calcOverPotential` are added. 
- `findCVcurrent_bin` (false position method) is added. 
- `log(x + sqrt(1+x^2))` is changed with `asinh(x)`
- `slide::State` now inherits from `std::array`
- Interpolation is in OCVcurves now. Simplified.  
- New cost function for `determine_OCV`, it is now much faster. 
- RK4 is implemented and FWEuler is seperated. 
- Unit test to control results is written in MATLAB. But only folders, so make also for files. 
- overwriteCharacterisationStates, overwriteGeometricStates are moved to cell, use cell versions. 
- validState is moved to Cell. 
- void `checkModelparam();` is created to check MATLAB param. 
- void `initialise(slide::State &s_ini)` is removed.
- `linInt_matrix` is removed. 
- `oneLevelCharacterisationFit` is removed. Its contents moved to `hierarchicalCharacterisationFit`. 
- A struct created for stress parameters.
- Modern C++ features are now incorporated into the code. 
- Ahpi_vec and Ahni_vec (100k double heap vectors) are removed. `slide::fixed_data` type is used. 
- Reading into double and converting whole vector in integer -> reading into integer.
- Header guards use `#pragma once` now. 
- Change pointers to references whereever possible. 
- printlevels are converted to enum
- Removed many magical numbers:  273 -> `PhyConst::Kelvin`
- Same parameters are defined multiple places are cleaned. For example Kokam parameters.
- paths are moved to settings.
- Data recording system is converted to push_back.
  - Instead of increasing stack memory to hold large arrays, now we use `std::vector`. 
  - Most of the raw arrays are changed with `std::array`; hence, the need for passing number of array elements is eliminated. 
  - `#define` directives are changed with `constexpr` values.
  - `<filesystem>` library is used for handing file/directory operations (`mkdir` and other C-style functions are removed). 
  - TDM-GCC/compiler specific features are eliminated (e.g., dynamic array on stack). 
  - Release mode is added to CMake. 
  - C-based header files are replaced by C++ versions (e.g., `<stdio.h>` to `<cstdio>`)
  - Range-based for loops adopted whereever possible. 

## Notable Bug-fixes
- Using `#define ns` for number of states was causing problems with `<chrono>` library literal `ns`; therefore, it was not possible to compile SLIDE other than TDM-GCC with C++1 standard. Bug was fixed by removing preprocessor directives.  
- Now heap allocation is used instead of stack area expansion. \
- `using namespace std;` from global scope is removed. `std::` is added all necessary functions. 
- warning: variables `j` and `timeCheck` used in loop condition not modified in loop body. `cycler.cpp mode==1, for (int j = 0; j < timeCheck; i++)` They are changed accordingly. 
- Important bug on linux, only use `std::abs` otherwise it may call int `abs(int)`.
- Object slicing due to copying cells is removed. 
- `abs(Icell - I) < pow(10, -10)` -> here `pow(10,-10)` may be ZERO. So `pow(10.0, -10.0)` is better. And `1e-10` is even better. 
- Remove used-defined classes from std namespace (such as state). A namespace called "slide" is created. `State` and `Cell_user` are moved from `std` to `slide`. 
- Others:
  - Typo `CalendarAgeig` -> `CalendarAgeing` fixed. 

## Tried but not implemented: 
- Fixed time step `*.csv` files were tried but no measurable effect in -O0.


## Dependencies
  * Required C++ standard is upgraded from C++11 to C++17. 

<br/><br/>
# SLIDE v1

This is the initial release of SLIDE. Deprecated but you can find a copy in [this branch](https://github.com/davidhowey/SLIDE/tree/SLIDE_v1).

## Dependencies
  * A compiler with C++11 support.
