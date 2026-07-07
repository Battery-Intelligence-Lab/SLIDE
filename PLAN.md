# SLIDE v4 — Architecture Refactor Plan (living document)

> **Status:** ACTIVE. Last updated 2026-07-07 by Fable (architect). Implementation delegated to Opus 4.8 agents.
> **How to use this document:** This is the single source of truth for the v4 refactor. Any session (Fable, Opus, human)
> continuing this work must (1) read this file first, (2) execute the next unblocked item in §6, (3) update §8 status
> ledger and this header. Requirements originate in `.claude/FABLE.md`. Do not re-litigate decisions in §4 without new
> evidence that overturns the cited rationale.
> **Operating rules:** ≤2 parallel agents (Volkan, 2026-07-07 — quota tightening; was ≤3); small reviewable commits on branch `Claude`; every user-visible change updates
> `CHANGELOG.md` (Unreleased); wall-clock benchmarks on this machine are UNRELIABLE (user runs parallel jobs) — judge
> performance by structural arguments (allocation counts, complexity, vectorizability), not timings.

---

## 1. Goals (from `.claude/FABLE.md`, condensed)

1. **Scale:** simulate 10⁴–10⁵ cells in packs; microsecond-level per-cell-step cost; GPU path eventually.
2. **Zero-overhead abstraction:** model choices (discretisation, ageing, cost functions) resolved at compile time or
   once per run off the hot path — never per-cell virtual dispatch in inner loops.
3. **Modularity:** Electrode as first-class object (kill `_p`/`_n` twin variables); injectable ageing models;
   swappable discretisations (future: 3D, spectral variants, DFN, SPMe).
4. **Pack solving:** replace tree-recursive algebraic solving with a flat sparse formulation usable by off-the-shelf
   solvers; support nested-topology inputs but *compile* them flat.
5. **Interfaces:** consistent across C++/Python/MATLAB (CasADi-style); Python layer PyBaMM-compatible
   (Experiment / ParameterValues / Simulation / Solution) so PyBaMM users can drop in SLIDE like `jax.numpy` for `numpy`.
6. **Memory-aware recording:** recompute derived quantities from state instead of storing; optional mmap/Parquet sinks.
7. **Portability:** Windows/macOS(x86+ARM)/Linux; parallelism works out of the box (no silent fallback-to-sequential);
   `slide::test::parallelisation()` diagnostic. Optional deps stay optional (core builds with none).
8. **Correctness discipline:** registered pass/fail bands before decisive runs; parity gates old-vs-new; oracles
   validated on non-degenerate cases (per user CLAUDE.md).

### 1.1 Performance contract (enforceable invariants — every phase gate re-checks these)

The guiding split: **the user-facing API may be arbitrarily pleasant because it runs once at setup; after
`compile()`/`build()`, only spans, indices, and precompiled kernels exist.** Strings, maps, options dicts,
nested constructors, parameter lookup — all cold path, all gone before the first step.

Structurally verifiable invariants (no wall-clock needed):

| # | Invariant | Enforced by |
|---|-----------|-------------|
| PC-1 | ZERO heap allocations per accepted simulation step | counting-allocator assert in tests (P1-G2) |
| PC-2 | ≤1 indirect call per BATCH per step; zero virtual dispatch per cell | code review + registry design (D-02) |
| PC-3 | All state contiguous, 64-byte aligned, SoA (SIMD-sweepable) | StateArena is the only state owner (D-01) |
| PC-4 | No exceptions, no locks, no I/O on the hot path | review + `Status` returns (D-09) |
| PC-5 | Inner loops see concrete types (templates), never function pointers/std::function | composition design (§3.2) |
| PC-6 | Memory ≤ ~300 B/cell state; batch-shared params stored once | 10⁵ cells ≈ tens of MB (P4-G3) |
| PC-7 | Rollback/checkpoint = memcpy, never tree traversal | arena snapshots (§3.1) |
| PC-8 | Pack solve cost: O(n) modes B/C, O(nnz) sparse mode A — never O(n²) dense, never nested iteration | compile() flattening (D-03) |

Any proposed change that violates a PC-invariant needs a PLAN.md §4 decision-log entry overturning it first.

## 2. Evidence base — what the audit found (2026-07-07, three scout reports)

All claims below are **[confirmed]** by code inspection with the cited locations, on branch `Claude` at commit `9dcf1ad`.

### 2.1 Current architecture (why it cannot reach the goals)

- Heap tree of `Deep_ptr<StorageUnit>` nodes; `StorageUnit` has ~40 pure virtuals (`src/StorageUnit.hpp:36`). No
  contiguous cell storage. `sizeof(Cell_SPM)` ≈ 1.2–1.5 KB [inferred], including TWO full `State_SPM` copies
  (`st` + rarely-used `s_ini`, `src/cells/Cell_SPM/Cell_SPM.hpp:38`).
- `State_SPM = State<19+2*nch>` → 29 model states + time/Ah/Wh, 256 B (`src/cells/Cell_SPM/State_SPM.hpp:25-50`);
  `nch` is a compile-time global (`src/settings/settings.hpp:41`); one `static Model_SPM<>` shared by ALL cells
  (`src/cells/Cell_SPM/Model_SPM.hpp:203-206`) — blocks heterogeneous geometry.
- Hot path: `Cycler::CC` calls `su->setCurrent` every loop iteration → `Module_p::setCurrent_previous_impl`
  (`src/modules/Module_p_impl.cpp:454`): damped **chord** iteration (quasi-Newton with a Jacobian factorised ONCE into
  function-`static` Eigen objects and never refactorised, `:480,524-548`) → linear convergence, ≤50-iter cap,
  historical 2500+ iteration blowups (`.claude/discussions.md:58`). ~150 virtual calls per cell per CC step.
- **`static` solver state is shared across all `Module_p` instances and threads** (`Module_p_impl.cpp:480,524-548`;
  `Module_p.cpp:127,156`; `:297-335` analytical path) — wrong currents / OOB for heterogeneous or multiple modules;
  data race under the `std::thread` fan-out in `src/utility/parallelisation.hpp:19`. Fatal for pack scale.
- Nested `Module_p` multiplies iteration counts per level: cost ≈ O(∏ per-level iters × cells) voltage evaluations.
- Time integration: fixed-step forward Euler; diffusion every substep, thermal+degradation once per `dt·nstep` window
  (`src/cells/Cell_SPM/Cell_SPM_dstate.cpp:229-334`). **No stability check**: modal eigenvalues grow ~O(nch⁴), so
  raising nch with fixed dt silently destabilises Euler.
- Recording: `Cell::storeData/writeData` are empty no-ops for SPM (`src/cells/Cell.hpp:109-110`); only ECM records.
  Only `time/Ah/Wh` are path-dependent — V, OCV, SOC, T, R are all recomputable from the 29 states.
- Two error channels coexist: `Status` enum and `throw int` (codes 10,11,14,98,99,101,104,106,108).

### 2.2 Chebyshev "only works for nch=5" — root cause (SOLVED, verify remains)

Two independent defects, both already fixed on branch `Claude`:
(a) hardcoded centre-node coefficient `-0.5`, correct value `-0.5·(-1)^N` with `N=nch+1` — accidentally right only at
nch=5 (even N); wrong SIGN for even nch (fixed at `Model_SPM.hpp:199`, consumed `Cell_SPM.cpp:195`);
(b) Eigen fixed-size `.inverse()` in-place aliasing corrupts `V` for nch≤4 (fixed with `.eval()`, `Model_SPM.hpp:176-177`).
Among {3,4,5,6,10}, nch=5 was the only survivor — fully explains the observation.
**Remaining gap:** `tests/unit/Chebyshev_test.cpp` is self-consistency only (round-trips, null-space, one linear
profile); NO external oracle. A shared-convention error would pass silently. → Phase 1 gate adds an analytic
transient-diffusion oracle (see §6, P1-G3). Release builds also silently skip the `assert(zero_pos==zero_neg)` guard
(`Model_SPM.hpp:195`).

### 2.3 Existing assets to harvest

- `setCurrent_analytical_impl` (`Module_p_impl.cpp:288-452`, Nilsu Atlan 2024): direct Thomas-style elimination for
  parallel branch currents — currently `Cell_ECM<1>`-only (hard `dynamic_cast`), disabled. Generalises to ANY cell via
  the Thevenin linearization interface (§3.4).
- Cached analytical inverse of the current-distribution matrix + pack DAE assembly (`Module_p_impl.cpp:44-146`),
  validated against `A22.inverse()`.
- `setVoltage` (`Module_p.cpp:115-206`): perturbation Jacobian refactorised per call — the correct Newton template.
- PI-control redistribution (Jorn) was REMOVED from the code (`Module_p.cpp:108`, `COMPLETED.md:31`). Its proper
  literature form is §3.4 mode C.
- pouch-cell-spectral project (user's, 2026): exponential-propagator integrator gave 26–52× over stiff BDF at sub-mV
  error; two logged "reference-convergence traps" (the oracle was wrong, not the model) — adopt both (§3.5, §5).
- dtw-cpp (user's): mmap store with header-CRC + offset-table validation (llfio); Arrow/Parquet optional-dep CMake
  wiring (note: its Arrow CPM build FAILS on Windows+Clang); working PyTorch-style `device=cpu/gpu/hpc` Python
  dispatch template with lazy data handles and preflight diagnostics.

### 2.4 Confirmed bug list (Phase 0 targets)

| # | Bug | Location | Failure |
|---|-----|----------|---------|
| A1 | Rollback restores wrong SU: loop uses `SUs[i]` not `SUs[j]`, offset never advances | `src/modules/Module.cpp:233-234` | silent state corruption on any failed `setStates` |
| A2 | Function-`static` solver matrices shared across instances/threads | `Module_p_impl.cpp:480,524-548`; `Module_p.cpp:127,156`; `:297-335` | wrong currents, OOB, data race |
| A3 | Jacobian LU never refactorised though `r_est` updates each iteration | `Module_p_impl.cpp:548,577` | chord/linear convergence → iteration blowups |
| A4 | `Qcontact` inner current accumulator never reset across `i` | `Module_p.cpp:274-279` | contact heat grows quadratically, wrong thermal load |
| A5 | `dynamic_cast<Cell_ECM<1>*>` unchecked | `Module_p_impl.cpp:325,350,428` | nullptr deref for non-ECM children |
| A6 | `Module_p::V()` uses only `Rcontact[0]` vs solver's cumulative drops | `src/modules/Module_p.hpp:45` vs `Module_p.cpp:60-95` | terminal V disagrees with equalisation model |
| A7 | Residual accumulation `error += max(|b(i)|, error)` (sum of maxes) | `Module_p_impl.cpp:500` | malformed convergence test (early-exit path) |
| B1 | `data.assign(data.end(), {...})` — overwrites instead of appends; ill-formed form | `src/recording/CellDataStorage.hpp:112` | time-series storage loses history / won't compile if instantiated |
| B2 | `Histogram` default ctor leaves `bins` empty; `add()` clamps then writes `bins[0]` | `src/types/Histogram.hpp:38-56,93` | heap corruption on default-constructed histogram |
| B3 | Wh throughput uses end-of-step voltage, self-flagged should be trapezoid | `src/procedures/Cycler.cpp:249` | systematic energy bias, worse at high C-rate |
| B4 | SEI model ids: default-case doc says 0–3 but cases 0–4 exist; DEG_ID docs swap SEI 1/2 | `Cell_SPM_degradation.cpp:101-103`, `DEG_ID.hpp:47-62` | users select wrong physics |
| B5 | Degradation Euler integrates ALL 32 indices incl. I and V slots (works only because those d_st are 0) | `Cell_SPM_dstate.cpp:331` | fragile; any model writing them corrupts I/V |

Deferred (design changes, handled by v4 core, not Phase 0): missing current-limit check (`Cell.hpp:72-78`), dead
adaptive `nOnce` stepping (`Cycler.cpp:219-224`), global shared `Model_SPM`/limits, `throw int` channel, dead SPM
recording, Euler stability guard.

### 2.5 Phase-0 outcomes & newly discovered defects (2026-07-07, Agent A report)

- **A6 REFUTED — not a bug.** `Module_p::V() = SUs[0]->V() − I()·Rcontact[0]` is exactly `getVall()[0]`:
  `I()` (sum of branch currents) IS the cumulative current through `Rcontact[0]`. Derived + proved by test
  (`phase0_A6_V_consistent_with_getVall` passes pre- and post-fix; existing `test_contactR` already asserted it).
  Documenting comment added; no behavior change. Audit table above stands corrected.
- **A3 nuance:** for linear cells (ECM/Bucket) `getRtot()` equals the exact incremental resistance, so the old
  chord converged fine; the stale-Jacobian defect only hurt NONLINEAR (SPM) cells. The A3 test is a registered
  correctness gate (exact conductance split to 1e-6 A, ≤5 iters), not a pre-fix failure.
- **NEW P0-C1 (build blocker):** `fmt` 11.0.2 consteval format checking rejected by clang 21.1.8 — the tree did
  not compile at all pre-Phase-0. Local build-dir patch only; REQUIRED: bump fmt ≥11.1 (or force
  `FMT_USE_CONSTEVAL 0`) in `cmake/Dependencies.cmake`.
- **NEW P0-C2 (major, pre-existing):** default `ocv_coefs` polynomial in `src/cells/Cell_ECM/Cell_ECM.hpp:48`
  yields OCV ≈ −55 782 V at every SOC → `Cell_Bucket`/`Cell_ECM`/`Module_p` test binaries fail at baseline, and
  legacy `Module_p` solving cannot converge on default ECM cells (inner `setCurrent` hardcodes `checkV=true`).
- **NEW P0-C3:** `Cell_SPM` test expects `thickp()==70e-6` but code gives 8.687e-5 — param/expected mismatch,
  needs root-cause (test wrong vs param regression).
- **NEW P0-C4 (minor):** same static-sizing hazard in the disabled Boost-integrator path
  (`Module_p_impl.cpp:213,230,256` free functions `parallel_model*`); leave for v4 (path is dead code).

Agent B outcomes (all five B-fixes landed; two audit corrections):

- **B3 was WORSE than audited:** `vi` was never assigned by `Cycler::setCurrent` (dead out-param), so CC
  energy throughput was always **0 Wh**, not merely end-of-step-biased. Registered prediction (trapezoid on a
  linear 4.0→3.0 V ramp = 3.5 Wh ± 1e-6) hit exactly post-fix.
- **B4 direction corrected:** `DEG_ID.hpp` had the physics-CORRECT SEI 1/2 pairing; it was the `.cpp` inline
  comments that were swapped (and `DEG_ID` was missing id 4). Docs aligned to formulas; no behaviour change.
- **B1 uncovered a latent header bug:** `#include "CellDataWriter.hpp"` sat inside `namespace slide`,
  corrupting standard headers and nesting `slide::slide::` — the recording header never compiled when included.
- **NEW P0-C5:** `Cycler::setCurrent`'s `v_now` out-param is dead (never assigned) — remove or fix.
- **NEW P0-C6 (minor):** `Cycler::CV`/CCCV energy still end-of-step-voltage biased (there `vi` IS assigned);
  make trapezoid for consistency.

---

## 3. Target architecture

Namespace `slide::core` (header tree `src/core/`), built alongside the legacy code (strangler pattern, §5). Legacy
`src/` stays green until parity gates pass.

### 3.1 State storage: SoA arena with archetype batches

A **CellBatch** holds N cells of identical model composition. All state lives in one contiguous, 64-byte-aligned
arena, **variable-major** (SoA): value of state-row `r` for cell lane `c` is `arena[r*stride + c]`.

```cpp
struct StateSpec  { std::string_view name; int rows; /* e.g. "zn" has nch rows */ Unit unit; };
struct StateSlice { int row_begin, rows; };            // integer handle a model component holds
class  StateArena {                                    // owns doubles; snapshot/restore = row-range memcpy
  std::span<double>       row(int r);                  // one variable across all cells (SIMD sweep)
  std::span<const double> row(int r) const;
  void snapshot(StateSlice, double* dst) const;  void restore(StateSlice, const double* src);
};
```

- Components (electrode diffusion, thermal, each ageing model) declare `StateSpec`s at **build time**; a
  `BatchBuilder` performs layout and hands each component its `StateSlice`. Adding a model appends rows — exactly the
  user's `beg/end` integer-offset idea, realised per-variable-row instead of per-cell-span so that kernels vectorise
  across cells. Layout is fixed after `build()`; never recomputed during simulation.
- A "cell" is `(batch_id, lane)`. Heterogeneous packs = several batches (one per distinct model composition).
  Per-cell parameters are also SoA rows (cell-to-cell variation stays vectorisable); parameters shared by the whole
  batch are scalars in the batch header.
- Cumulative quantities (`time, Ah, Wh`) are ordinary rows. Everything else observable (V, OCV, SOC, T_surface, R) is
  **derived on demand** — the recording system stores state snapshots, not observables (evidence: §2.1 last bullet).
- Rollback: `Cycler`-level checkpoints become row-range memcpy of the arena (no per-cell traversal).

### 3.2 Zero-overhead model composition

Physics kernels are free-function templates over spans; a cell model is a compile-time composition:

```cpp
template <int NCH> struct SpectralDiffusion { /* StateSpec spec(); void rhs(BatchView, ...); */ };
struct ThermalLumped;                            // same shape
template <class... Mechanisms> struct AgingPack; // SEI_Kinetic, SEI_Diffusion, LAM_Stress, CS_..., LiPlating...
template <class Diffusion, class Thermal, class Aging> struct SPM; // the composed cell model

using Kokam_SPM = SPM<SpectralDiffusion<5>, ThermalLumped, AgingPack<SEI_Kinetic, LAM_Stress>>;
```

- The library **explicitly instantiates** a registry of common compositions; a runtime factory (string/options-keyed,
  fed by the PyBaMM options dict, §3.9) selects an instantiation ONCE at build time. Hot loop: zero virtual calls per
  cell; exactly one indirect call **per batch per step** (amortised over thousands of lanes).
- Runtime-flexible ageing injection (user requirement) without per-cell virtuals: an ageing mechanism compiled into
  the registry can be toggled per-batch by a mask evaluated per batch, not per cell; truly novel user models enter via
  recompilation (C++) or via selecting the closest registered composition (Python/MATLAB) — same philosophy as CasADi:
  the dynamic language chooses among precompiled fast paths.
- `nch` stops being a global: it is a template parameter of the composition; different batches may use different nch.
- **Extension axes (Volkan, 2026-07-07 — "COMSOL-lite": a fixed compiled menu, not a runtime PDE engine).** Each
  slot of the composition is an independent axis a future model plugs into without touching the others:
  *discretisation* (`SpectralDiffusion<NCH>` now; FVM/FDM particle grids, DFN/SPMe electrolyte discretisations later —
  same slot, same `CellDesign` input); *geometry* (spherical particle now; plate/cylinder via a geometry tag the
  kernel is templated on — shape factors are compile-time constants); *thermal / heat transfer* (ThermalLumped now;
  Thermal1D-through-thickness and pack-coupled cooling later; `h_conv` is a §3.11 `ParamFn` so constant vs
  correlation vs user-injected costs the same at run time). New physics = new type in a slot + registry entry —
  the hierarchy, arena, and solver layers are untouched by construction.

### 3.3 Physical description hierarchy — Electrode as first-class (kill `_p`/`_n`)

Two layers, deliberately separate (this is how "good hierarchy" and "zero overhead" coexist, §1.1):

**Layer 1 — description (cold, value-semantic OO, mirrors real battery anatomy).** What a user or a parameter
set constructs; arbitrarily pleasant; exists only until `build()`:

```cpp
enum class Domain : int { neg = 0, pos = 1 };

struct ActiveMaterial { OCVCurve ocv; double cs_max, x_0, x_100; Arrhenius D_s, k_ct; };
struct ElectrodeDesign {                       // one physical electrode
  ActiveMaterial am;
  double thickness, porosity, active_frac, particle_radius;
  StressParams stress;                          // per-electrode mechanics
  std::vector<AgingMechanismSpec> aging;        // SEI/plating attach to neg, LAM/CAM-loss to either —
};                                              // mechanisms live on the electrode they act on, not on the cell
struct SeparatorDesign  { double thickness, porosity; };
struct ElectrolyteDesign{ double c0, D_e, t_plus; };   // consumed by SPMe/DFN later; SPM ignores it
struct CellDesign {                             // the actual battery, as built
  std::array<ElectrodeDesign, 2> electrode;     // electrode[Domain::neg], electrode[Domain::pos]
  SeparatorDesign separator; ElectrolyteDesign electrolyte;
  ThermalDesign thermal;                        // mass, Cp, surface area, h_conv
  double capacity_Ah, area;                     // form factor
};
```

**Layer 2 — compiled (hot).** `build()` flattens `CellDesign` + model choice into batch parameter rows and
selects kernels from the registry (§3.2). Per-domain quantities become 2-row state/param entries indexed by
`Domain d` — the 15-line hand-written `thickp/thickn/...` blocks (`Cell_SPM_dstate.cpp:212-226`) and the
verbatim-duplicated stress integrals (`Cell_SPM_degradation.cpp:492-665`) collapse to domain loops.

Consequences the current design cannot express, and this one gets for free:

- **Ageing attaches to the electrode/interface it physically acts on** (SEI + plating on the negative
  particle surface, LAM per electrode, cracking per particle) — replaces the cell-global `DEG_ID` flag block.
- **Asymmetric chemistries are natural** (blended cathode later = `std::vector<ActiveMaterial>` per electrode).
- **Future models slot in without touching the hierarchy**: SPMe/DFN consume `SeparatorDesign`/
  `ElectrolyteDesign` that SPM ignores; a 3D/spectral discretisation is a different Layer-2 compilation of the
  SAME `CellDesign`.
- Migration: legacy accessors kept as forwarders during transition (`thickp()` → `thick(Domain::pos)`).

### 3.4 Pack topology: netlist compiled to sparse form

Users keep the ergonomic nested constructors (series/parallel combinators), but construction only records topology.
`Pack::compile()` flattens everything to a flat netlist — nodes, branches (cell branches + contact resistances) — and
precomputes the sparse structure. The tree disappears from the runtime.

The ONLY interface the pack solver needs from a cell is the **affine Thevenin linearization**
`V_cell(I) ≈ ocv − r_eff·I` (per batch, vectorised). This removes all `dynamic_cast`-on-cell-type coupling (§2.4 A5).

Three solver modes behind one interface, selected per pack (all operate on the same compiled netlist):

- **Mode A — exact sparse Newton (default for ≤10³ cells):** modified nodal analysis; sparse LU (Eigen SparseLU in
  core; SuiteSparse KLU as optional dep) with symbolic factorisation cached at compile(), numeric refactorisation each
  Newton iteration (fixes §2.4 A3 by design). Liionpack (JOSS 2022, 10.21105/joss.04051) proves this staggered
  MNA-current-solve + per-cell-stepping architecture at pack scale.
- **Mode B — ladder fast path:** pure parallel/series ladders get the Thomas-style O(n) elimination (generalisation of
  Nilsu's `setCurrent_analytical_impl` via the Thevenin interface). Detected automatically at compile().
  **Upgrade (SOTA check 2026-07-07):** the published form of this solution — Lone, **Atlan**, Fasolato, Raimondo,
  **Drummond**, arXiv:2508.14454 (2025) — goes further: analytical current distribution converts the parallel-pack
  **DAE into plain ODEs** (algebraic constraint eliminated, not just solved fast; ~44% faster than direct DAE at
  n=135, gains grow with n). Mode B adopts the ODE reformulation where its conditions hold (parallel topology,
  ohmic interconnects, known R_k); the Thomas elimination remains the fallback for mixed ladders. This RESOLVES
  §7 Q6: "Ross's analytical solution" = this paper; Nilsu is a co-author and her code is its implementation.
- **Mode C — relaxation advance (target for 10⁴–10⁵ cells, GPU):** Jorn's "PI" idea in its literature form:
  **waveform relaxation** (Miekkala & Nevanlinna, SIAM J. Sci. Stat. Comput. 1987, 10.1137/0908046) with **Baumgarte
  constraint stabilisation** (replace `g=0` by `ġ+2αg=0`) on the voltage-equality constraints. O(n) per step, no
  matrix solve, embarrassingly parallel; convergence conditional on an index-1 topological criterion — `compile()`
  checks it and refuses Mode C otherwise. Gain α tuned like a Baumgarte damping constant, with the exact solve
  (Mode A) as the arbiter in tests. **Honesty note (SOTA check):** WR for packs is published (J. Energy Storage
  2022, S2352152X21014304) and Baumgarte is standard, but the WR+Baumgarte COMBINATION has no published precedent
  found — it is our synthesis; treat as a research contribution, gated on the Mode-A arbiter (P4-G1/G2). Mode C is
  only needed where B's conditions fail (arbitrary series-parallel meshes at 10⁴⁺ cells).

NOT chosen: one monolithic 100k-state DAE handed to IDA/KLU — KLU is serial and a global Jacobian factorises poorly;
liionpack deliberately avoids it too (§4 D-06).

### 3.4.1 Solver memory: per-pack `SolverWorkspace` (the Jacobian cache done right)

**Origin (Volkan, 2026-07-07):** the legacy function-`static` Eigen objects in `Module_p_impl.cpp` were not an
accident — they were a deliberate quasi-Newton memory: keep the factorised Jacobian and secant resistance estimates
across calls so the solver does NOT pay a Jacobian rebuild + factorisation every time step. The Phase-0 fix
(per-call locals, refactorise every iteration — `Module_p_impl.cpp:481,527,546`) restored correctness but destroyed
that memory: every `setCurrent` now re-seeds `r_est` from `getRtot()` and refactorises per iteration. v4 keeps the
idea and fixes the storage class. The named literature form of the idea is the **chord/Shamanskii method** (Newton
with a frozen, periodically refreshed Jacobian — Kelley, *Iterative Methods for Linear and Nonlinear Equations*,
SIAM 1995, ch. 5).

**The object.** One `SolverWorkspace` per compiled `Pack`, owned as a member, allocated exactly once at
`Pack::compile()` (PC-1: zero per-step allocations):

```cpp
struct SolverWorkspace {
  // persistent across time steps — the "memory" the statics were emulating
  VectorXd        x;        // last converged branch currents + node voltages → warm start
  VectorXd        r_eff;    // per-branch secant Thevenin resistances (NOT re-seeded from getRtot() each call)
  Factorization   fact;     // numeric LU of J(r_eff), reused across iterations AND steps
  int             age = 0;  // steps since last numeric refactor (diagnostics)
  bool            valid = false;
  // fixed at compile()
  SymbolicPattern pattern;  // analyzePattern() once; every refactor reuses it
  // scratch, sized once
  VectorXd        b, dx;
};
```

**Freshness policy — refresh only when the cache stops paying (registered rule, not vibes):**

- Iterate with the *cached* factorisation: each iteration = residual build O(n) + triangular solve O(nnz), no
  factorisation. Monitor the contraction ratio `ρ_k = ‖b_{k+1}‖_∞ / ‖b_k‖_∞`.
- Chord-method theory: convergence is linear with rate ∝ ‖J_cached − J_true‖. Between consecutive time steps
  `r_eff` drifts at the %-level (SOC/T move slowly), so ρ ≪ 1 and 1–2 iterations/step is the expected regime —
  this is exactly the win the statics bought.
- **Refresh trigger:** `ρ_k > 0.5` or `k > 4` iterations → rebuild J from the current secant `r_eff` (O(nnz)
  fill, pattern unchanged), numeric refactorisation only, reset `age`. Worst case degenerates to today's
  refactor-per-iteration behaviour; typical case is one factorisation amortised over many steps.
- **Tier-0 shortcut:** a batch whose cells declare `linear == true` (ECM with fixed R) has a *constant* J →
  factor once at `compile()`, never refresh; the chord iteration is then exact Newton and converges in 1 step.
  Detected structurally at compile, not by runtime probing (fixes the A3 class of bug by construction).
- Mode B (Thomas ladder) skips `fact` entirely — an O(n) elimination is cheaper than any cache — but still uses
  `x`/`r_eff` warm start. Mode C (waveform relaxation) uses `x` as the initial waveform iterate.

**Invalidation contract — the part the statics never had (and why they broke: shared across instances/threads,
wrong sizes, stale across topology changes — §2.4 A2/A4):**

| Event | Action |
|-------|--------|
| `Pack::compile()` | build `pattern`; `valid = false` |
| Checkpoint **restore** (rollback) | `valid = false`; keep `pattern`; keep `r_eff` as initial guess only |
| Topology event (contactor open/close) | recompile symbolic; full reset |
| Temperature/SOC step beyond registered threshold | force refresh at next solve |

- The workspace is **ephemeral state**: excluded from `StateArena` snapshots (it is reconstructible; snapshots
  stay pure-state memcpy, PC-7). The restore path MUST call `ws.invalidate()` — parity test required: solve
  after rollback digit-matches a cold solve on the same state.
- **Threading:** one workspace per pack instance; parallelism is across packs/batches, so no sharing and no locks
  (PC-4). Never function-`static`, never `thread_local` (a pack that migrates threads under the pool would pick
  up a foreign cache).

**Deliverable placement:** Phase 2 (Mode A/B carry the workspace from day one); the rollback-invalidation parity
test joins P2-G1.

### 3.4.2 Pack description layer — intuitive generation, compiled once

The user-facing construction API is a **value-semantic description tree** (cheap, copyable, printable); nothing
about how the user *writes* the topology may influence how it is *solved* — `compile()` erases the authoring shape
(fixes the FABLE.md complaint "depending how it is created, solving becomes a mess"). Mirrors the two-layer split of
the parameter system (§3.11): friendly at description time, flat arrays at run time.

```cpp
using namespace slide;
auto cell = CellDesign{ .neg = ElectrodeDesign{...}, .pos = ElectrodeDesign{...}, ... };   // §3.11

auto brick  = parallel(5, cell,  Link{ .R_contact = 0.5_mOhm });        // 5p
auto string_= series (14, brick, Link{ .R_busbar  = 0.1_mOhm });        // 14s5p
auto pack   = Pack(string_, PackOptions{ .cooling = CoolingDesign{...} });

// cell-to-cell heterogeneity — composes with §3.11 varied():
auto spread = parallel(5, varied(cell, { .capacity = Spread::normal(0.02, /*seed*/ 42) }), Link{...});

// escape hatch for non-ladder topologies + liionpack/PyBaMM users:
auto custom = Pack(Netlist::from_csv("topology.csv"));   // liionpack netlist schema (node, node, R, type)

pack.compile();   // ONE pass: archetype batching → flat netlist → sparsity pattern →
                  // ladder detection (Mode B) / index-1 check (Mode C) → SolverWorkspace alloc
```

- `series(n, child, link)` / `parallel(n, child, link)` nest arbitrarily; repetition shares one description by
  value. `Link` packs carry connection physics (contact/busbar R, later: fuses, contactors) — designated
  initializers, no positional soup (D-15).
- **Addressing:** compile assigns hierarchical path IDs (`"s03.p2"`, matching legacy `getFullID` style) in a flat
  string table; results are indexable by path or `(batch, lane)`. The tree itself is discarded after compile
  (kept only if the user retains the design object — it is just data).
- **Round-trip:** any compiled pack can emit its flat netlist (`pack.netlist()`) — serialisable, diffable,
  re-importable; this is also the PyBaMM/liionpack interop surface (Phase 7).
- Contactor/switch events (later phase): declared in the description; each discrete configuration gets its own
  cached symbolic pattern at compile so a switch flip is a pattern swap, not a re-analysis.

### 3.5 Time integration: exponential modal propagator + multirate

- **Diffusion (fast, linear, modally decoupled):** the Chebyshev eigendecomposition already gives
  `dz_k/dt = D·λ_k·z_k + B_k·j` per mode. For piecewise-constant flux j over a substep h the EXACT update is
  `z_k(t+h) = e^{Dλ_k h}·z_k(t) + ((e^{Dλ_k h} − 1)/(D·λ_k))·B_k·j`, and for the λ_k=0 mode `z_k += B_k·j·h`.
  (Units: D[m²/s]·λ[1/m²] = 1/s; exponent dimensionless. ✓) This is unconditionally stable — it deletes the
  O(nch⁴) Euler stability cliff (§2.1) outright — and exact for CC segments, so substep counts collapse.
  pouch-cell-spectral measured 26–52× on the analogous structure.
- **Multirate (Strang):** electrical (µs–s scale) inner steps inside thermal/degradation (min–h scale) outer steps —
  formalising what the current code does ad hoc, but with O(h²) splitting error and an error-controlled outer step.
  Beware stiff-operator order reduction to O(h); acceptable since degradation rates are slow and weakly coupled.
- **Event-aligned stepping:** integrate segment-boundary-aligned; never step across a current/power breakpoint
  (pouch-cell lesson: sampling through a jump gives an O(1) error that substeps cannot remove).
- Variable outer step with a proportional error controller comes after parity (Phase 3); parity mode replicates
  legacy fixed-step Euler digit-for-digit first (§5).

### 3.6 Error handling

Single channel: `Status` (extended) + `[[nodiscard]]`; hot path never throws; `throw int` eliminated. Constraint
violations surface as `Status` from `step()`, and the Cycler restores the arena checkpoint (§3.1).

### 3.7 Recording & I/O

- `Recorder` subscribes to a batch: stores arena snapshots at a cadence (every k-th accepted step) + cumulative rows;
  derived observables computed lazily at read time through the same kernels (one code path, no drift).
- Sinks: CSV (compat), custom binary mmap format (adopt dtw-cpp's 64-byte header + magic + version + endian marker +
  header-CRC32 + validated offset table idiom; implement with plain OS mmap wrappers, do NOT take the llfio dep),
  optional Arrow/Parquet behind `SLIDE_WITH_ARROW` (warning: dtw-cpp's Arrow-via-CPM breaks on Windows+Clang; prefer
  find_package-first).
- Histograms remain as bounded lossy summaries; fixed (B2) and moved behind the Recorder.

### 3.8 Parallelism & device

- Own minimal thread pool over batches (std::thread; no TBB dependency, no silent sequential fallback);
  `slide::test::parallelisation()` reports cores, backend, and a measured speedup smoke check.
- `device={cpu,gpu}` selected per Simulation (PyTorch-style, template = dtw-cpp `_api.py` lazy-dispatch design).
  GPU phase: one-cell-per-thread batched stepping (DiffEqGPU, arXiv:2304.06835, shows 20–100× vs vmap approaches for
  exactly this shape); MNA/relaxation coupling stays on host between batch steps. CUDA first; the SoA arena (§3.1) is
  already the correct device memory layout.

### 3.9 Language interfaces (CasADi-style symmetry)

- **C++ core API:** `Simulation(PackSpec, ModelOptions, ParameterSet, Experiment) → Solution`.
- **Experiment grammar:** ONE C++ parser of the PyBaMM instruction strings, reused by all bindings:
  verb ∈ {Discharge, Charge, Hold, Rest} × rate (C-rate | A | W | V) × termination (`for <n> h/min/s`,
  `until <x> V|A|C/n`, joinable `or until`), cycles via tuples/`*`.
- **ParameterSet:** flat `map<string, double|Curve>` keyed by PyBaMM's exact strings
  (`"Negative electrode thickness [m]"` …) + absorption table → internal `ElectrodeParams`. Ship `Chen2020` first.
  **Verified against PyBaMM stable 26.6.2.0, CalVer (SOTA check 2026-07-07):** key on the NEW names with
  deprecation-aware aliases for old ones — `"Negative electrode diffusivity"`→`"Negative particle diffusivity"`
  (v25.6.0, PR #3624); `"Exchange-current density for lithium plating"`→`"…for lithium metal electrode"` (v25.4.0);
  `"1 + dlnf/dlnc"`→`"Thermodynamic factor"` (v23.3). Behavioural notes: `update()` is now update-insert
  (`check_already_exists` deprecated v25.12.0); default constants no longer auto-added on construction (v25.12.0).
  Experiment grammar additions to cover: custom steps (v25.4.0), custom terminations (v25.10.0), `start_time`
  (v23.9). Solution additions: `.yp`, `.observe()` (v25.12.0) — out of scope v4.0, document as such.
- **Python:** nanobind + scikit-build-core wheels. `slide.Experiment`, `slide.ParameterValues("Chen2020")`,
  `slide.Simulation(...).solve()`, dict-like `Solution["Terminal voltage [V]"]` (+ `.entries`, call-interpolation,
  `.plot()`, `save_data(..., to_format="csv"|"matlab")`). Options dict (`{"SEI": "solvent-diffusion limited", ...}`)
  maps onto the composition registry (§3.2). Out of scope, documented as such: symbolic expression trees, arbitrary
  `spatial_methods`/`var_pts`, symbolic events.
- **MATLAB:** MEX + `+slide` package mirroring the Python surface (Phase 8).

### 3.10 v4 coding conventions (Volkan, 2026-07-07)

- **Modern C++20 idioms**: ranges, `std::span`, concepts, structured bindings; standard algorithms over
  hand-rolled loops. Reuse and extend the existing utility helpers (`src/utility/slide_algorithms.hpp`,
  `src/utility/free_functions.hpp` — e.g. the summing helpers) rather than re-implementing them inline.
- **Value semantics over pointer soup**: containers of CONCRETE objects (`std::vector<CellBatch<Model>>`,
  SoA rows) wherever the composition is known — never `std::vector<StorageUnit*>`/`Deep_ptr` in hot
  structures. (Volkan verified experimentally that `vector<actual_class>` works alongside `vector<SU*>`;
  the concrete form is the v4 default, polymorphic indirection allowed only at cold boundaries.)
- **Dependencies**: small, well-justified deps are acceptable (explicit user approval 2026-07-07), but the
  core must still build with none (non-negotiable #4) — so new deps enter as optional with a fallback.

### 3.11 Parameter system — COMSOL-grade injection, table-speed execution

COMSOL's model, adapted: parameters live on the **physical entity they belong to** (material → electrode →
cell), any scalar can be **promoted to a function** of state, everything carries **units** — and none of it
survives to the hot path.

**Grouping (no `function(...15 inputs)` — API rule):** every construction goes through designated-initializer
parameter packs mirroring §3.3's hierarchy; any API function needing >4 parameters takes a pack instead:

```cpp
auto neg = ElectrodeDesign{
  .am        = slide::materials::Graphite_Chen2020(),   // named library pack, COMSOL material-library style
  .thickness = 85.2_um,                                  // Quantity: unit-checked at build, raw SI double after
  .porosity  = 0.25,
  .particle_radius = 5.86_um,
};
neg.am.D_s = ParamFn::arrhenius(3.3e-14_m2ps, 35.0_kJpmol);          // promote constant → function of T
neg.am.k_ct = ParamFn::table1D(c_data, k_data);                      // or any user table/expression/callable
```

- `Quantity` = value+unit at Layer 1, checked and converted to SI once at `build()`; UDLs extend the existing
  `src/utility/units.hpp` (`25.0_degC` already there). Zero runtime cost (PC-4/PC-5 untouched).
- Named packs form the **materials library** (`slide::materials::…` — Graphite/NMC/LFP + electrolytes),
  which is also where PyBaMM parameter-set absorption lands (Chen2020 → packs, §3.9).
- **Cell-to-cell variation = COMSOL parameter case, generalised:** any Layer-1 field can be marked
  `varied(span_or_distribution)` → compiles to a per-lane SoA row instead of a batch-shared scalar. This is
  the native mechanism for heterogeneous 10⁵-cell studies (manufacturing spread, ageing spread).

**The speed trick — canonicalisation at `build()` (D-16):** arbitrary injected functions would wreck PC-5
(opaque callables in inner loops). So `build()` compiles EVERY function-valued parameter into one of a fixed
set of hot-path forms:

1. batch-shared scalar (constant),
2. per-lane SoA row (varied constant),
3. **uniform-grid LUT** + linear interpolation — index = fma+floor, 2 loads, lerp: branchless, SIMD-able,
4. separable product of 1D LUTs for (c,T)-dependence — `D(c)·arrh(T)`, 2 lookups + multiply — with full 2D
   bilinear as fallback,
5. analytic Arrhenius (exp kept inline; it beats a LUT of log-spaced T for typical ranges — revisit if not).

User expressions/callables of ANY complexity are **tabulated once at build** over their declared domain with a
registered accuracy gate: refine grid until max rel. interp error vs the source < 1e-6 or 4096 points
(ASSUMED band; tighten per-parameter if physics demands). Hot-path cost is therefore INDEPENDENT of what the
user injected — COMSOL flexibility, table speed. Existing `OCVcurves`/`XYdata`/`interpolation.hpp` (which
binary-search nonuniform data) get resampled to uniform grids at build; the search disappears.

**Kernel-side rule:** kernels receive `BatchView` (arena rows + resolved param views) and a `StepCtx` — never
loose scalar arguments; adding a parameter never changes a kernel signature.

**Traceability (cheap, cold):** every pack supports `describe()` → name/value/unit/provenance dump (which
library pack, which overrides, which tabulation error) — the COMSOL report equivalent, and the §1 "trail
another scientist can trust".

## 4. Decision log (do not re-litigate without overturning the evidence)

| ID | Decision | Rationale | Rejected alternative |
|----|----------|-----------|---------------------|
| D-01 | SoA arena, variable-major, batch archetypes | SIMD/GPU across cells; snapshot=memcpy; matches user's integer-offset idea | AoS per-cell objects (today): 1.5 KB/cell, pointer-chased, unvectorisable |
| D-02 | Compile-time composition + instantiation registry; dispatch per BATCH | zero per-cell virtuals; dynamic languages select precompiled paths (CasADi precedent) | per-cell virtual `StorageUnit` (today: ~150 virtual calls/cell/step) |
| D-03 | Netlist compiled flat at `Pack::compile()` | nested trees multiply solver iterations per level (§2.1); off-the-shelf sparse solvers become usable | keep runtime tree + recursive equalisation |
| D-04 | Cell↔pack interface = affine Thevenin (ocv, r_eff) per batch | removes cell-type coupling (A5); enables Modes A/B/C uniformly | solver dynamic_casts per cell type (today) |
| D-05 | Three pack modes: sparse Newton / ladder Thomas / WR+Baumgarte | scale coverage 10⁰→10⁵ cells; Liionpack proves A; Nilsu's code seeds B; Miekkala-Nevanlinna + Baumgarte theory grounds C | single algorithm for all scales |
| D-06 | NO monolithic global DAE via IDA+KLU | KLU serial; global Jacobian factorisation dominates; liionpack independently rejected it | 100k-state IDA solve |
| D-07 | Exponential modal propagator for diffusion | exact for CC segments; unconditionally stable (kills O(nch⁴) Euler cliff); 26–52× precedent | keep forward Euler + add CFL guard |
| D-08 | Strang multirate electrical/thermal+degradation | formalises existing split with O(h²) + error control | MRI-GARK (reserve for later if splitting order limits) |
| D-09 | Single `Status` error channel, no hot-path throw | two channels today (§2.1); exceptions on hot path cost and corrupt partial state | keep dual channel |
| D-10 | Store state snapshots, derive observables lazily | only time/Ah/Wh are path-dependent (§2.1); memory ~256 B/cell/snapshot vs unbounded frames | store per-step observable frames (today's ECM path) |
| D-11 | Experiment grammar parsed in C++, shared by all bindings | single source of truth, CasADi-style symmetry | per-language parsers |
| D-12 | Own mmap writer (dtw-cpp idiom, no llfio dep); Arrow optional find_package-first | llfio is a heavy dep pinned to `develop`; dtw-cpp's Arrow CPM build breaks on Windows+Clang | adopt llfio / require Arrow |
| D-13 | Own thread pool over batches; no TBB | user requirement: parallelism must not silently fall back | oneTBB optional dep |
| D-14 | Strangler migration with digit-diff parity gates | big-bang rewrite of a research-validated code loses physics silently | in-place rewrite |
| D-15 | Parameters grouped in designated-initializer packs on physical entities; any API function needing >4 params takes a pack | COMSOL grouping precedent; kills `function(...15 inputs)` and the 4-positional-x0/x100 ctor style | flat global parameter bag / long positional signatures |
| D-16 | Every function-valued parameter canonicalised at `build()` to {scalar, SoA row, uniform LUT, separable LUT product, Arrhenius}; build-time accuracy gate | hot-path cost independent of injected expression complexity; preserves PC-5 | opaque std::function/callables in kernels |
| D-17 | Units checked at description layer (`Quantity` + UDLs), raw SI doubles after `build()` | dimensional safety with zero runtime cost | runtime unit objects (cost) or no checking (CLAUDE.md violation) |
| D-18 | Per-pack `SolverWorkspace`: warm start + chord/Shamanskii Jacobian reuse with contraction-monitored refresh + explicit invalidation (§3.4.1) | keeps the speed the legacy statics bought (Volkan's quasi-Newton memory) without their races/staleness; Kelley ch.5 grounds the refresh rule | (a) function-statics (races, cross-instance pollution — §2.4 A2/A4); (b) refactorise every iteration (current Phase-0 state: correct, memoryless, pays O(n³/nnz) per iteration) |
| D-19 | Pack construction = value-type combinator tree + `Netlist` escape hatch; `compile()` erases authoring shape (§3.4.2) | intuitive generation AND solver independence from nesting style; liionpack netlist schema = free PyBaMM interop | (a) runtime tree solved recursively (today — nesting multiplies iterations); (b) netlist-only API (hostile for the 99% ladder case) |
| D-20 | Mode B adopts the analytical DAE→ODE reformulation for parallel packs (Lone/Atlan/Fasolato/Raimondo/Drummond, arXiv:2508.14454) where conditions hold; Thomas elimination as fallback | removes the algebraic constraint entirely (exact, no iteration, ~44% faster than direct DAE at n=135, growing with n); resolves Q6 — this IS Ross's solution, Nilsu's code implements it | route pure-parallel packs through generic MNA (Mode A) or WR (Mode C) — both pay for a constraint this case doesn't need |

## 5. Migration strategy & verification discipline

1. **Legacy stays green.** New core grows in `src/core/` (`slide::core`); legacy `src/` untouched except Phase 0 bug
   fixes. CI runs both.
2. **Parity harness (the arbiter).** `tests/parity/` runs registered scenarios through legacy AND core and
   digit-diffs: single Kokam SPM 1C CC discharge; CCCV cycle; 3s2p ECM pack with contact R; SPM pack 2p. For
   digit-diff, core runs in **legacy-Euler mode** (same scheme, same dt). Registered band, written BEFORE the run:
   max |ΔV| ≤ 1e-12 V, max |Δstate| ≤ 1e-12 (rel) over the full trajectory.
3. **Scheme upgrades validated separately** (never against loose-tolerance references — pouch-cell trap, logged
   twice): exponential propagator vs a tolerance-CONVERGED reference (tighten until the reference moves < 0.01 mV);
   registered band: |V_expm − V_ref| < 0.1 mV over a 1C discharge with a current step.
4. **Oracles on non-degenerate cases:** parity scenarios include asymmetric electrodes, nonzero contact R,
   heterogeneous initial SOC — never only the symmetric/uniform case.
5. **Baselines recorded first.** Before each phase: run full ctest, record failing-test names + counts in §8; every
   commit re-runs; report deltas ("2 failing {a,b} → 3: +c, caused by me").
6. **No timing claims.** Machine runs concurrent jobs; performance evidence = allocation counts, iteration counts,
   complexity, vectorisation reports — never wall clock (until user provides a quiet machine).
7. **Payoff checkpoints — the refactor must prove itself before it is allowed to grow (Volkan, 2026-07-07).**
   Insurance first: the strangler strategy (item 1) means legacy stays green the whole time, so the worst case of an
   underperforming v4 core is deleting `src/core/` — nothing user-facing is ever bet on it. On top of that, each
   early phase ends with a **wall-clock checkpoint run by Volkan on a quiet machine** (the only trusted timing
   source), with the band registered here BEFORE the run and an explicit abort threshold:
   - **PAY-1 (Phase 1 exit):** 10⁴ identical Kokam SPM cells, 1C CC discharge, 1 h simulated — v4 batch vs a loop of
     legacy `Cell_SPM`. Hypothesis [inferred, from devirtualisation + SoA + expm substep collapse; pouch-cell
     precedent 26–52× on the integrator alone]: ≥5×. **Abort threshold: <2× → STOP; profile, find where the model
     was wrong, redesign before any Phase-2 work.** No new phase on top of an unproven core.
   - **PAY-2 (Phase 2 exit):** 16s4p heterogeneous pack, CC cycle — v4 compiled pack vs legacy `Module_s/Module_p`.
     Hypothesis: ≥10× (flat solve replaces nested iteration, D-03; workspace kills refactorisations, D-18).
     Abort threshold: <3× → stop and re-examine before Phases 3–4.
   - **PAY-3 (Phase 4 exit):** 10⁵-cell pack advances on the quiet machine within memory budget (<1 GB state) —
     feasibility, the original goal (§1). Fails → Mode C redesign before GPU work.
   Structural proxies (allocations, virtual calls, iteration/factorisation counts, bytes/cell) are tracked
   continuously as leading indicators; wall clock at the checkpoints is the confirming evidence. Rule: phase N+1
   does not start before phase N's checkpoint passes or Volkan explicitly waives it.

## 6. Phased roadmap

Effort tags are relative. Every phase ends: tests green, CHANGELOG updated, §8 ledger updated, small commits pushed.

### Phase 0 — Stabilise legacy (IN FLIGHT, 2 Opus agents, disjoint files)
Fix §2.4 bugs. Agent A (solver/state): A1–A7 in `src/modules/*` + regression tests + CHANGELOG. Agent B (data/misc):
B1–B5 in `src/recording/`, `src/types/Histogram.hpp`, `src/procedures/Cycler.cpp`, degradation docs + tests;
CHANGELOG lines to `.claude/changelog-phase0b.md` (avoid merge conflict; architect merges).
**Gate P0-G1:** full ctest delta vs recorded baseline shows only intended changes; each fix has a test that failed
before the fix.

### Phase 1 — Core data model + SPM kernels (the keystone)
Deliver: `StateArena`/`StateSpec`/`StateSlice`/`BatchBuilder`, `Domain`/`ElectrodeParams`, `SpectralDiffusion<NCH>`,
`ThermalLumped`, ageing kernels (SEI/LAM/CS/plating ported mechanism-by-mechanism), composition registry + factory,
legacy-Euler stepping mode. Single-cell `Simulation` façade.
**Gates:** P1-G1 parity single-cell scenarios (§5.2 band). P1-G2 one batch of 10⁴ identical cells steps with ZERO
per-step heap allocations (assert via allocation-counting new). P1-G3 Chebyshev external oracle: transient sphere
diffusion with constant-flux BC vs the analytic series solution (Carslaw & Jaeger form), registered band rel. err
< 1e-6 at nch=5,8,12; plus cross-check vs Howey Spectral_li-ion_SPM conventions.

### Phase 2 — Pack layer
Deliver: netlist combinators + `Pack::compile()` (flatten, sparsity, ladder detection, index-1 check), Mode A sparse
Newton (Eigen SparseLU; KLU optional), Mode B Thomas ladder, Thevenin batch interface, `SolverWorkspace` (§3.4.1)
with warm start, chord refresh policy, and invalidation contract.
**Gates:** P2-G1 parity vs legacy `Module_s`/`Module_p` on 3s2p (band §5.2) + rollback-invalidation test (solve after
restore digit-matches cold solve). P2-G2 Mode B ≡ Mode A on ladders to
1e-10 A. P2-G3 nested-constructed pack (p-in-p-in-s) compiles flat and solves in ONE Newton loop (no nested
iteration), Newton iterations ≤ 8 on the 4p heterogeneous-resistance case that historically blew up. P2-G4 workspace
efficacy: on a 100-step 4p CC segment, count of numeric factorisations ≤ 10 (vs 1 per iteration today) at identical
converged currents (1e-10 A) — an iteration/factorisation COUNT, not a timing (§5.6).

### Phase 3 — Integration upgrade
Deliver: exponential modal propagator, Strang multirate, event-aligned segmentation, arena checkpoints/rollback,
variable outer step with error controller.
**Gates:** P3-G1 expm vs converged reference band (§5.3). P3-G2 energy/charge conservation: |ΔAh_in − ΔAh_out −
ΔAh_stored| < 1e-9 Ah on a full cycle. P3-G3 nch=12 runs stable where legacy Euler diverges (registered demonstration).

### Phase 4 — Mode C relaxation solver (scale path)
Deliver: WR+Baumgarte advance, gain selection rule, index-1 topology check, Mode A as arbiter.
**Gates:** P4-G1 steady-state branch currents within 0.1% of Mode A on 16p heterogeneous pack. P4-G2 constraint drift
bounded per theory (registered α-dependent band). P4-G3 10⁵-cell pack advances (structural check: memory < 1 GB,
zero per-step allocations; no timing claim).

### Phase 5 — Experiment & Cycler v2
Deliver: C++ PyBaMM-grammar parser, Experiment→segment compiler, CC/CV/CCCV/power/rest + drive cycles on the new core,
termination conditions.
**Gate P5-G1:** grammar round-trip test suite (every documented string form parses; malformed strings produce
diagnostics); CCCV parity vs legacy Cycler.

### Phase 6 — Recording & I/O
Deliver: Recorder (snapshot cadence + lazy derived), CSV sink, mmap binary sink (header-CRC idiom), optional Parquet.
**Gate P6-G1:** derived-vs-stored equivalence test (V recomputed from snapshot == V recorded live to 1e-12);
mmap file survives the hardened-open validation tests (truncated/corrupt-header cases).

### Phase 7 — Python bindings + PyBaMM compat
Deliver: nanobind module, wheels (scikit-build-core), `Experiment/ParameterValues/Simulation/Solution`, Chen2020
absorption table, options→registry map, `device=` dispatch (dtw-cpp template), pytest in CI.
**Gate P7-G1:** the PyBaMM getting-started Tutorial-5 experiment script runs with `import slide as pybamm`-style swap
and produces a voltage curve within a registered band of PyBaMM's own SPM (band set after a converged-reference run,
expected ~10 mV model-difference scale — document, don't hide, the modelling differences).

### Phase 8 — MATLAB MEX, GPU, docs, release
MEX `+slide` package symmetric with Python; CUDA one-cell-per-thread batch stepping (host-side coupling); docs site
update (installation, quickstarts ×3 languages, "add a cell model", "add an ageing mechanism"); v4.0.0 SemVer release,
CHANGELOG consolidation. Gates defined when phase opens.

### Beyond v4.0 (recorded so design choices don't foreclose them)

- **WASM GUI (Volkan, 2026-07-07):** Emscripten build of the dependency-free core + a browser front-end. Costless to
  keep open: core already must build with zero deps (non-negotiable #4), no threads assumed outside the pool (§3.8),
  no filesystem dependence in the hot path. Only rule it adds NOW: no platform API in core without a portable seam.
- Thermal 1D/2D per-cell models; blended electrodes (`vector<ActiveMaterial>`); SYCL/HIP (Q5); f32 storage (Q1).

## 7. Open questions for Volkan (OPEN/ASSUMED ledger)

| ID | Question | Current assumption |
|----|----------|-------------------|
| Q1 | float32 state option for GPU/memory? | ASSUMED f64 everywhere in v4.0; f32 storage + f64 accumulate later (dtw-cpp precedent) |
| Q2 | Keep `Cell_ECM<N_RC>` template generality in v4 core? | ASSUMED yes as a composition (`ECM<NRC>`), it's cheap |
| Q3 | Legacy API: keep as façade over core after parity, or hard-break at v4.0? | ASSUMED façade through v4.x, delete in v5 |
| Q4 | KLU/SuiteSparse as optional dep acceptable? | ASSUMED yes (optional, Eigen SparseLU default) — matches non-negotiable #4 |
| Q5 | GPU: CUDA-only first? | ASSUMED yes; SYCL/HIP revisit after CUDA lands |
| Q6 | Ross's analytical parallel solution — is `setCurrent_analytical_impl` (Nilsu 2024) the code you meant, or is there a separate derivation to recover? | **RESOLVED 2026-07-07 [confirmed]**: arXiv:2508.14454 (Lone, Atlan, Fasolato, Raimondo, Drummond 2025) — Nilsu co-authored it; her code implements it. Adopted as Mode-B upgrade (D-20) |
| Q7 | PyBaMM version to target for the parameter absorption table? | ASSUMED latest stable at Phase 7 start; key-rename check mandatory |

## 8. Status ledger

| Date | Item | State |
|------|------|-------|
| 2026-07-07 | 3-scout audit (architecture/bugs, numerics, external) | DONE — findings in §2, full reports in session transcripts |
| 2026-07-07 | Chebyshev nch≠5 root cause | SOLVED pre-session on branch `Claude` (§2.2); oracle test still MISSING (P1-G3) |
| 2026-07-07 | PLAN.md v1 written | DONE (this file) |
| 2026-07-07 | Phase 0 (A1–A7, B1–B5) | DONE — 11 commits on `Claude` (7529039…91faaf9); A6 refuted (non-bug, documented); regression tests added (Module_p_phase0, CellDataStorage, Histogram, Cycler_energy); CHANGELOG consolidated. Baseline had 4/6 test binaries failing PRE-EXISTING (§2.5 P0-C2/C3); agents' tests all pass; no new failures introduced |
| 2026-07-07 | Phase 0 follow-up (P0-C1 fmt/clang21, P0-C2 ocv_coefs, P0-C3 thickp, P0-C5/C6 Cycler) | IN FLIGHT — goal: full ctest green |
| 2026-07-07 | Solver-memory design (§3.4.1, D-18, P2-G4) | DONE — Volkan clarified the legacy statics were intentional quasi-Newton Jacobian memory; design keeps the memory, adds invalidation + thread safety. Agent cap now ≤2 (header) |
| 2026-07-07 | Pack description layer (§3.4.2, D-19) | DONE — combinator tree + Netlist escape hatch; compile() erases authoring shape. Per Volkan: agents = Opus HIGH (not xhigh), no "ultrathink" in agent prompts |
| 2026-07-07 | SOTA verification (report: `.claude/reports/sota-verification-2026-07-07.md`) | DONE — C1/C3/C5/C6 CONFIRMED, C2/C4 NUANCED (liionpack maintenance-mode; WR+Baumgarte combo unpublished = our synthesis), C7 PyBaMM 26.6.2.0 keys verified. Q6 RESOLVED, D-20 added (arXiv:2508.14454 Mode-B upgrade). Description-layer zero-cost pattern confirmed (CasADi/Eigen/Halide precedent) provided erasure is total |
| — | Phases 1–8 | NOT STARTED — Phase 1 is next; do not start before Volkan reviews §3/§4 |
