# SLIDE v4 — Architecture Refactor Plan (living document)

> **Status:** ACTIVE. Last updated 2026-07-10. Implementation in progress.
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
  **CORRECTED 2026-07-09 [confirmed]:** the ageing/thermal port falsifies "only time/Ah/Wh" — legacy hides
  path-dependent state OUTSIDE `State_SPM`: Dai-stress previous-step memory `sparam.s_dai_p_prev`
  (`Cell_SPM_dstate.cpp:244`, CONSUMED by LAM via `|s_dai_p − s_dai_p_prev|/s_dt`, `Cell_SPM_degradation.cpp:377`)
  and the `Therm_Qgen`/`Therm_time` accumulators (`Cell_SPM.hpp:62-64`, `Cell_SPM_dstate.cpp:290-291`). v4 RULE:
  anything a later step reads MUST be an arena row (else D-10 lazy-derivation and checkpoint-restart break
  silently). Gate P1-G4 (bitwise restart) exists to catch exactly this class.
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
  batch are scalars in the batch header. **Dual use (named 2026-07-09):** lanes need not be electrically coupled —
  the same machinery runs 10⁵ INDEPENDENT single-cell variants as one batch (UQ/Monte-Carlo over `varied()` spread,
  fitting sweeps; DifferentialEquations.jl `EnsembleProblem` shape). One infrastructure, two products: packs and
  parameter studies.
- **Stated assumption (Fable review 2026-07-07):** the batch-amortisation wins (PC-2's one indirect call per batch,
  D-01's SIMD width) presume FEW compositions with MANY lanes each (realistic packs: 1–3 archetypes). Per-cell
  heterogeneity in *composition* (not parameters — those are `varied()` rows) fragments the pack toward
  single-lane batches ≈ the old per-cell dispatch. `compile()` emits a diagnostic ("N batches, min lanes = k")
  and warns below a registered lane floor; per-lane ageing-mask rows are the escape hatch before fragmenting.
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
  same slot, same `CellDesign` input. Concrete menu recorded 2026-07-09, audit §5.4: conservative FVM — PyBaMM's
  default, valuable for like-for-like P7-G1 parity; 3-parameter parabolic profile — ultra-cheap tier for 10⁵-cell
  screening; eigenfunction/Duhamel modal basis — analytic λ from tan μ = μ, kills the eigensolve-conditioning
  concern; Legendre–Galerkin — symmetric by construction); *geometry* (spherical particle now; plate/cylinder via a geometry tag the
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
  **Upgrade (SOTA check 2026-07-07, PROVISIONAL — Volkan: "never trust without measuring"):** the published form
  of this solution — Lone, **Atlan**, Fasolato, Raimondo, **Drummond**, arXiv:2508.14454 (2025) — goes further:
  analytical current distribution converts the parallel-pack **DAE into plain ODEs** (algebraic constraint
  eliminated, not just solved fast; ~44% faster than direct DAE at n=135). This RESOLVES §7 Q6: "Ross's analytical
  solution" = this paper; Nilsu is a co-author and her code is its implementation.
  **Regime caution — the paper's assumptions clash with our main use case:** derivation assumes linear-in-current
  ohmic drop with known/fixed R_k, and reported validity degrades with cell heterogeneity (report cites ~1e-4 σ
  cap for n>100) — while SLIDE packs are SPM cells (r_eff varies with SOC/T) with %-level `varied()` spread, at
  10⁴–10⁵ ≫ the paper's n=135 experimental validation. So: implement behind a compile()-verified regime check +
  runtime flag, NEVER the silent default, and admit it only through gate **P2-G5** (below) with Mode A as arbiter.
  Thomas elimination stays the unconditional ladder fast path.
- **Mode C — relaxation advance (target for 10⁴–10⁵ cells, GPU):** Jorn's "PI" idea in its literature form:
  **waveform relaxation** (Miekkala & Nevanlinna, SIAM J. Sci. Stat. Comput. 1987, 10.1137/0908046) with **Baumgarte
  constraint stabilisation** (replace `g=0` by `ġ+2αg=0`) on the voltage-equality constraints. O(n) per step, no
  matrix solve, embarrassingly parallel; convergence conditional on an index-1 topological criterion — `compile()`
  checks it and refuses Mode C otherwise. Gain α tuned like a Baumgarte damping constant, with the exact solve
  (Mode A) as the arbiter in tests. **Honesty note (SOTA check):** WR for packs is published (J. Energy Storage
  2022, S2352152X21014304) and Baumgarte is standard, but the WR+Baumgarte COMBINATION has no published precedent
  found — it is our synthesis; treat as a research contribution, gated on the Mode-A arbiter (P4-G1/G2). Mode C is
  only needed where B's conditions fail (arbitrary series-parallel meshes at 10⁴⁺ cells).

**Solver robustness (added 2026-07-09, orthogonal review; SPICE solved all of these decades ago):**

- **Consistent initialization:** `Pack::compile()` is followed by `Pack::initialize()` — a consistent-IC solve of the
  branch currents at t=0 (heterogeneous initial SOC, cold workspace, `valid=false`) BEFORE the first step; a
  checkpoint restore re-runs it through the `ws.invalidate()` path. (PyBaMM/IDA do the same consistent-IC solve.)
- **Newton safeguarding:** damped updates with per-iteration TRIAL LIMITING — clamp trial branch currents/states so
  OCV-table and Butler-Volmer evaluations stay inside their tabulated/asinh domains (SPICE `pnjlim` idiom). A trial
  iterate outside the OCV table must be limited, never fed to interp/asinh to manufacture NaNs that a successful
  factorisation then hides.
- **Rescue path:** on non-convergence, source-stepping homotopy (ramp the applied terminal current from a solvable
  problem; SPICE2 gmin/source stepping, Xyce continuation). If rescue fails: `Status` failure with DEFINED atomicity —
  the staggered scheme rolls back ALL batches of the pack to the pre-step checkpoint; a half-stepped pack is never
  observable.
- **Detection on the compiled netlist:** ladder (Mode B) and index-1 (Mode C) checks run on the COMPILED netlist, so a
  `Netlist::from_csv` input electrically equal to a ladder gets the same fast path as combinator input.

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
  **Divergence guard (2026-07-09):** ρ_k ≥ 1 on two consecutive iterations → force refresh; STILL ρ_k ≥ 1 →
  escalate to damped full Newton (§3.4 safeguarding) — never iterate a diverging chord to the cap. `age` feeds the
  temperature/SOC-threshold force-refresh row of the invalidation table (it was tracked-but-unused).
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
  **Numerical note (registered 2026-07-09 so Phase 3 doesn't rediscover it):** `(e^{Dλh}−1)/(Dλ)` cancels
  catastrophically for |Dλh| ≪ 1 (slow modes hit this every step) — implement as `expm1(x)/x` with a φ₁ Taylor
  branch below a small-|x| threshold (standard exponential-integrator practice), NOT `exp(x)−1`; otherwise the
  cancellation noise eats the §5.3 validation band.
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
- **Async recording pipeline (Volkan 2026-07-09: "we don't stop") — the simulation NEVER blocks on I/O:**
  snapshots go into a double/triple-buffered RING; the sim thread (or GPU) writes a snapshot and continues
  immediately; a dedicated CPU recorder thread pool drains the ring — compress, then sink. On GPU: snapshots
  leave on a side CUDA stream into PINNED host buffers, so the D2H copy overlaps compute. Compression: SoA rows
  with byte-shuffle then LZ4/zstd blocks (the Blosc2 pattern — shuffle exposes the near-constant high bytes of
  smoothly-evolving f64 states; delta-along-time optional on top); LOSSLESS default; blocks carry the D-12
  header-CRC idiom. Cache hygiene: the snapshot copy uses streaming/non-temporal stores so recording never
  evicts the hot arena from cache. **Backpressure is explicit, never silent:** sink can't keep up ⇒ policy
  ∈ {block, thin the cadence + LOG the thinning count}, chosen per run — silently dropped frames violate the
  §1 traceability standard.

### 3.8 Parallelism & device

- Own minimal thread pool over batches (std::thread; no TBB dependency, no silent sequential fallback);
  `slide::test::parallelisation()` reports cores, backend, and a measured speedup smoke check.
- **Determinism rule (2026-07-09):** every cross-batch reduction (pack current sums, thermal adjacency, residual
  norms) uses FIXED-ORDER accumulation (serial or static-tree) — never scheduling-order-dependent atomics.
  Run-to-run reproducibility at the 1e-12 parity scale depends on it (GROMACS/LAMMPS deterministic-reduction
  precedent). Corollary: per-lane adaptive dt is REJECTED (destroys the SoA sweep) — batch-uniform step, error
  controlled on the batch max-norm.
- `device={cpu,gpu}` selected per Simulation (PyTorch-style, template = dtw-cpp `_api.py` lazy-dispatch design).
  GPU phase: one-cell-per-thread batched stepping (DiffEqGPU, arXiv:2304.06835, shows 20–100× vs vmap approaches for
  exactly this shape); MNA/relaxation coupling stays on host between batch steps. CUDA first; the SoA arena (§3.1) is
  already the correct device memory layout.
- **GPU fit at pack scale (Volkan question 2026-07-09: 2·10⁴-cell packs on GPU?) — yes; memory is a non-issue.**
  PC-6 gives 2·10⁴ × ~300 B ≈ 6 MB state; with params + ydot + observables scratch, tens of MB against GBs of VRAM.
  The design is device-shaped by construction: SoA variable-major rows = coalesced access (thread c ↔ lane c); zero
  virtuals/allocations/exceptions (PC-1/2/4/5) means kernels port as `__device__` free functions — virtual dispatch
  is not even EXPRESSIBLE efficiently on device, so D-02 is a GPU prerequisite, not just a CPU optimisation; the
  arena is ONE `cudaMalloc` at build. The four REAL concerns, recorded now so the CUDA phase starts honest:
  (1) *occupancy* — 20k lanes alone under-fill a modern GPU (~10⁵⁺ threads to saturate): parallelise element-wise
  over (row × lane) — the diffusion sweep is 29 × 20k ≈ 5.8·10⁵ elements — and/or run ensembles (§3.1 UQ framing);
  (2) *coupling* — the staggered host solve needs one (ocv, r_eff) exchange per batch step (~320 KB at 20k lanes,
  ~0.1 ms over PCIe): amortise with several device substeps per exchange (the multirate shape already does this), or
  move Mode B's scan / Mode C's relaxation on-device and kill the transfer entirely (Mode C is the GPU-native mode —
  O(n), Jacobi-style, no factorisation); (3) *launch overhead* — capture the per-step kernel sequence in a CUDA
  Graph (µs-scale launch latency otherwise dominates µs-scale steps); (3b) *kernel fusion* (SOTA-kernel lever,
  2026-07-09) — the §3.12 pipeline runs as separate row sweeps on CPU, but on device {observables → addRhs chain
  → advance} FUSES into ONE kernel per batch step: arithmetic intensity is low, so every extra global-memory pass
  costs ~a full step; the observables scratch lives in registers/shared memory inside the fused kernel (CUDA
  Graphs is the launch lever, §3.7's async ring is the I/O-overlap lever); (4) *f64 throughput* — 1/32–1/64 of f32 on
  consumer GPUs: the f32-storage instantiation (Q1 insurance + §3.12 scalar-generic kernels) is the consumer-GPU
  path; datacenter parts (A100/H100 class) run f64 natively fine.
  **Cache residency [inferred, roofline — NOT a timing claim; measured only at PAY checkpoints]:** the workload is
  memory-bound (arithmetic intensity ~a few flops/byte), and a 2·10⁴-cell working set (state + ydot + obs ≈ 20–30 MB)
  fits INSIDE the L2 of current GPUs (H100 50 MB, RTX 4090 72 MB) — stepping runs at cache bandwidth, not HBM.
  Traffic ≈ 20k × ~32 rows × 8 B × 2 ≈ 10 MB/step ⇒ µs-scale per step at L2/HBM bandwidths. "Very quick very long"
  simulations (Volkan 2026-07-09) are the PRODUCT of three multipliers, not GPU alone: exponential propagator (exact
  on CC segments — substep count collapses, D-07) × multirate outer steps (thermal/ageing at min-scale, D-08) × GPU
  batch stepping (per-cell cost collapses). Order-of-magnitude: 10 simulated years at ~10 s effective outer steps
  ≈ 3·10⁷ steps × µs-scale ⇒ device time in minutes [inferred — every factor here is a stated assumption; the
  hypothesis is registered as the GPU-phase PAY target, to be banded properly when that phase opens].

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
  **Same-parameters ⇒ same-results contract (Volkan 2026-07-09) — two unit-test layers, gate P7-G3:**
  (1) *absorption round-trip*: every absorbed key (Chen2020 table, BPX file) goes PyBaMM-named value → internal
  packs → `describe()` and must return the IDENTICAL SI value — exact equality, no simulation, runs in plain C++
  CI. (2) *behavioural parity fixtures*: registered scenarios (C/50 quasi-static — essentially OCV replay — and
  1C transient) compared against COMMITTED PyBaMM reference traces, generated once by a pinned-version script
  and checked in (C++ CI stays hermetic, no Python dependency; an optional pytest job regenerates fixtures and
  flags upstream drift). Bands registered per scenario BEFORE comparison. Honest limit, stated up front:
  identical parameters CANNOT give bit-identical results across tools — Chebyshev-spectral vs FVM discretisation,
  different OCV interpolants — so the band tests CONVERGED-model agreement (both discretisations → the same
  continuum SPM as they refine), tightest at low C-rate [inferred: sub-mV at C/50; set after the converged run].
- **Python:** nanobind + scikit-build-core wheels. `slide.Experiment`, `slide.ParameterValues("Chen2020")`,
  `slide.Simulation(...).solve()`, dict-like `Solution["Terminal voltage [V]"]` (+ `.entries`, call-interpolation,
  `.plot()`, `save_data(..., to_format="csv"|"matlab")`). Options dict (`{"SEI": "solvent-diffusion limited", ...}`)
  maps onto the composition registry (§3.2). Out of scope, documented as such: symbolic expression trees, arbitrary
  `spatial_methods`/`var_pts`, symbolic events.
- **MATLAB:** MEX + `+slide` package mirroring the Python surface (Phase 8).
- **Sensitivities / gradients — PyBOP entry (Volkan directive 2026-07-09: "we need to enter the PyBOP ecosystem").**
  Deliverable: forward sensitivities `∂(V, observables)/∂θ` for a registered fitting-parameter set (Q10), surfaced
  PyBaMM-style (`Solution` sensitivities; the `simulateS1`-shaped call PyBOP's gradient optimisers consume). Three
  routes, ranked (D-24): (1) **dual-number forward mode** through the scalar-generic kernels (§3.12 — the mechanism
  is already reserved; exact to roundoff, cost ≈ (1+n_θ)× per swept param, right for the few-parameter fits PyBOP
  runs); (2) **analytic modal sensitivities** where structure permits — the diffusion update is LINEAR in z, so
  ż_θ = D·A·z_θ + (∂_θ of coefficients) rides the SAME exponential propagator exactly, and ∂V/∂(z,T,δ,…) is
  closed-form for SPM (chain rule through §3.12 observables) — cheapest and exact for CC segments; (3) adjoint for
  many-parameter gradients — beyond v4.0 (recorded, not started). Arbiter for any route: central finite differences
  on the same trajectory at registered step sizes.

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
3. **exact nonuniform piecewise-linear curve + uniform segment-index accelerator** — preserves committed measured
   knots and P7-G3 absorption round trips bit-for-bit while replacing binary search with O(1) fma/floor + lookup;
   added 2026-07-10 after the Kokam OCV tables confirmed nonuniform knots,
4. **uniform-grid LUT** + linear interpolation — index = fma+floor, 2 loads, lerp: branchless, SIMD-able,
5. separable product of 1D LUTs for (c,T)-dependence — `D(c)·arrh(T)`, 2 lookups + multiply — with full 2D
   bilinear as fallback,
6. analytic Arrhenius (exp kept inline; it beats a LUT of log-spaced T for typical ranges — revisit if not).

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

### 3.12 Kernel/integrator contract — RHS-form kernels over structs-of-spans, no framework (D-22)

**Design principle (Volkan, 2026-07-09):** flexible, expandable classes with NO unnecessary abstraction; extreme
speed; off-the-shelf ODE integrators/methods usable; intuitive, sparse structure. This section is the §6 Phase-1
"observable layer + BatchView/StepCtx" blocker, resolved under that principle. **REVISED 2026-07-09 after orthogonal
review (D-23) — three critical defects fixed before implementation.** The hot layer is a short list of plain
constructs — no inheritance anywhere:

1. **`BatchView` / `StepCtx` — structs of spans, REBINDABLE per eval.** Layout (slices, strides, param views) is
   fixed at `compile()`; the (y, ydot) BASE POINTERS rebind per RHS call. Reason (D1, critical): an adaptive
   integrator evaluates f at ITS OWN trial vectors — CVODE clones internal N_Vectors and never hands the arena back —
   so views bound once to arena rows would silently integrate stale state. Rebinding = two pointer swaps on a fixed
   layout (hot-safe). Zero-copy holds for the IC load and the accepted-step writeback, NOT for trial evals.
2. **RHS-form kernels behind a FIXED eval pipeline.** Every RHS eval runs, in order: (i) **zero ydot** (one memset —
   stepper's job); (ii) ONE shared **observables stage** `computeObservables(y, obs)` into a once-allocated per-batch
   scratch — c_surf `= C·z + D·flux` + centre node, Butler-Volmer η, OCV/entropic dOCV, Rdc, V (SPICE "device load"
   precedent: thermal AND every ageing mechanism read the same obs, computed once, never per-component); (iii) the
   component `addRhs(y, ydot, obs, ctx)` chain. Accumulation is REQUIRED semantics — legacy SEI/plating add into the
   same zn slots diffusion owns (`Cell_SPM_dstate.cpp:210`) — which is exactly why a missing zero pass is a
   guaranteed wrong answer, not a style issue. Observables stay the ONE code path shared with the Recorder (D-10).
   (`SpectralDiffusion::stepEuler` survives ONLY as the §5.2 parity special case, pinned to legacy op order.)
3. **ODE-row mask.** `build()` marks which arena rows are integrable ODE states. Legacy's algebraic I/V slots and
   kinked cumulative rows (d|Ah|/dt = |I|, nonsmooth at I=0) re-create the B5 bug class and make adaptive steppers
   chatter — algebraic/cumulative rows are advanced by the stepper OUTSIDE the RHS, never handed to an integrator.
4. **Steppers own time:** (a) `EulerLegacy` (parity, §5.2); (b) `ExponentialModal` (Phase 3, D-07 — exact, uses the
   modal structure directly, not the generic RHS); (c) **generic RHS hook** — a plain `rhs(t, span y, span ydot)`
   callable for ANY user integrator; (d) **CVODE adapter as ARBITER, not production path** (`SLIDE_WITH_SUNDIALS`,
   optional; core ships (a)+(b)+(c) with no deps — non-negotiable #4). Scope (D5): small lane counts only — BDF needs
   a linear solver, the batch Jacobian is block-diagonal (~30×30 per lane), and stock serial N_Vector/SUNLinSol
   cannot exploit that at L=10⁴; a custom SUNLinSol is REJECTED bloat (if CVODE can't run it off-the-shelf, it stays
   an arbiter). A Boost.odeint adapter is REJECTED (its state-type algebra requirements force resize/copy machinery
   for no arbiter value beyond CVODE; hook (c) covers odeint users). Integrator calls are SEGMENT-SCOPED: never
   across a current/power breakpoint (§3.5) — `CVodeSetStopTime` + re-init per segment, else BDF history across the
   discontinuity collapses the order.
5. **Events are first-class.** Termination/protection conditions (V reaches Vlim, T limit, plating onset η ≤ 0) are
   declared zero-crossing functions g(y) = 0 located by root-finding (CVodeRootInit; Modelica state events) — not
   the legacy per-step threshold checks, which overshoot by up to one dt. Feeds Phase 5 (Cycler v2 termination +
   machine-readable termination REASON in `Solution`, needed for PyBaMM compat anyway).

Two zero-cost insurance notes: kernels are `template<class Real>` (scalar-generic — Stan/CoDiPack/CppAD pattern);
they are templates anyway, and NOT hardcoding `double` in kernel bodies keeps dual-number forward sensitivities
(parameter fitting, PyBOP-style) and an f32 GPU instantiation open without retrofit — `real_t` is the default
binding, not a commitment. And the {rhs, observables, event indicators} triple is deliberately congruent with the
FMI Model-Exchange shape — keep it congruent (costs nothing), implement FMU export never-in-v4.

Requirement check: *flexibility* — a new mechanism is one `StateSpec` + one free function, nothing else touched;
*speed* — same vectorised row sweeps, zero virtuals, one ydot arena extra (PC-1/2/3/5 hold; per-batch dispatch
unchanged); *off-the-shelf* — any integrator drives one RHS callback per batch per eval; *sparse/intuitive* — cold
layer stays the §3.3 battery anatomy; hot layer is the constructs above.

**Arbiter bonus (§5.3), with its stated limit (D4):** CVODE at rtol ~1e-12 on the same RHS is an independent
converged reference for the exponential propagator — different mathematics, same model. Caveat: with `StepCtx`
inputs frozen per outer segment, CVODE converges to the operator-SPLIT trajectory; it arbitrates INTEGRATOR error,
not splitting error (splitting error is owned by D-08's outer-step controller and validated separately).

**Boundary [stated]:** off-the-shelf *ODE* (not DAE) applies per batch because pack algebra (voltage-equality
constraints) is solved by Modes A/B/C BETWEEN batch steps (staggered, liionpack-style, D-05); the monolithic-DAE
hand-off stays rejected (D-06). Multirate (D-08) composes: outer and inner splits each expose their own RHS.

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
| D-16 | Every function-valued parameter canonicalised at `build()` to {scalar, SoA row, exact indexed piecewise-linear curve, uniform LUT, separable LUT product, Arrhenius}; build-time accuracy gate | hot-path cost independent of injected expression complexity; preserves PC-5. Exact indexed form added for nonuniform measured tables because P7-G3 requires exact absorption round trips; uniform resampling would perturb committed knots | opaque std::function/callables in kernels |
| D-17 | Units checked at description layer (`Quantity` + UDLs), raw SI doubles after `build()` | dimensional safety with zero runtime cost | runtime unit objects (cost) or no checking (CLAUDE.md violation) |
| D-18 | Per-pack `SolverWorkspace`: warm start + chord/Shamanskii Jacobian reuse with contraction-monitored refresh + explicit invalidation (§3.4.1) | keeps the speed the legacy statics bought (Volkan's quasi-Newton memory) without their races/staleness; Kelley ch.5 grounds the refresh rule | (a) function-statics (races, cross-instance pollution — §2.4 A2/A4); (b) refactorise every iteration (current Phase-0 state: correct, memoryless, pays O(n³/nnz) per iteration) |
| D-19 | Pack construction = value-type combinator tree + `Netlist` escape hatch; `compile()` erases authoring shape (§3.4.2) | intuitive generation AND solver independence from nesting style; liionpack netlist schema = free PyBaMM interop | (a) runtime tree solved recursively (today — nesting multiplies iterations); (b) netlist-only API (hostile for the 99% ladder case) |
| D-20 | Mode B adopts the analytical DAE→ODE reformulation for parallel packs (Lone/Atlan/Fasolato/Raimondo/Drummond, arXiv:2508.14454) **PROVISIONALLY**: behind a compile()-verified regime check + flag, admitted only via measurement gate P2-G5 vs Mode A (Volkan 2026-07-07: never trust without measuring; paper's linearity/heterogeneity/known-R_k assumptions clash with SPM + `varied()` packs) | removes the algebraic constraint entirely — exact, no iteration — WHERE VALID; resolves Q6 (this IS Ross's solution, Nilsu's code implements it) | (a) unconditional adoption (regime unproven at our heterogeneity/scale/nonlinearity); (b) ignoring it (pays for a constraint the pure-parallel linear case doesn't need) |
| D-22 | Kernels are RHS-form free functions over `BatchView`/`StepCtx` structs-of-spans; observables = shared free functions (one code path for RHS internals AND recording); steppers own time; generic RHS adapter exposes `arena.raw()` zero-copy to CVODE/odeint/user integrators as OPTIONAL deps (§3.12) | user requirements 2026-07-09 (flexible + expandable, no unnecessary abstraction, off-the-shelf integrators, sparse intuitive structure); arena contiguity makes flat-y interop free; PC-1/2/3/5 preserved; CVODE-at-tight-rtol doubles as the §5.3 converged-reference arbiter | (a) virtual Component hierarchy (per-cell dispatch — PC-2 violation, abstraction without need); (b) integration hard-wired inside kernels (today's Cell_SPM: blocks adaptive/implicit methods entirely); (c) required SUNDIALS dep (violates non-negotiable #4) |
| D-24 | Forward sensitivities are a NAMED v4 surface (PyBOP entry, §3.9): dual-number forward mode through scalar-generic kernels as the general route; analytic modal sensitivities (same exponential propagator, linear-in-z structure; closed-form ∂V/∂states) where structure permits; central-FD as arbiter; adjoint deferred beyond v4.0; parameter set = Q10 | Volkan 2026-07-09: PyBOP ecosystem entry required; gradient optimisers need `simulateS1`-shaped sensitivities; §3.12 scalar-generic kernels make forward mode near-free to add; PyBaMM pays IDAS-sensitivity cost for the same surface | (a) FD-only "gradients" (noise floor wrecks optimiser line searches); (b) adjoint-first (right for n_θ ≫ 10, wrong for typical 3–8-parameter cell fits, much higher implementation risk); (c) full runtime AD dependency (CoDiPack/Enzyme as REQUIRED dep — violates non-negotiable #4) |
| D-23 | §3.12 REVISED after orthogonal review (2026-07-09): rebindable BatchView (adaptive integrators evaluate f at THEIR trial vectors — CVODE clones internals, zero-copy holds only at IC/writeback); fixed eval pipeline {zero ydot → one shared observables stage → addRhs chain}; ODE-row mask (algebraic I/V + kinked cumulative rows never handed to an integrator); events as g(y)=0 root-finding; CVODE scoped to small-N ARBITER; kernels scalar-generic `template<class Real>` | D1/D2/D3 were guaranteed-wrong-answer defects if implemented as first drafted; legacy ageing accumulates into shared rows (`Cell_SPM_dstate.cpp:210`) and hides `_prev`/accumulator state (§2.1 correction); B5 precedent for algebraic rows; SPICE load-stage/limiting, Modelica events, SUNDIALS rootfinding, Stan/CoDiPack scalar-generic precedents | (a) view bound once to the arena (silently integrates stale state); (b) per-component observable recomputation (2–3× redundant asinh/interp per eval) or ad-hoc hidden scratch; (c) Boost.odeint adapter (state-type algebra forces copy machinery; the plain rhs hook covers odeint users); (d) custom SUNLinSol to run CVODE at pack scale (bloat — arbiter role only); (e) per-lane adaptive dt (destroys the SoA sweep) |

## 5. Migration strategy & verification discipline

1. **Legacy stays green.** New core grows in `src/core/` (`slide::core`); legacy `src/` untouched except Phase 0 bug
   fixes. CI runs both.
2. **Parity harness (the arbiter).** `tests/parity/` runs registered scenarios through legacy AND core and
   digit-diffs: single Kokam SPM 1C CC discharge; CCCV cycle; 3s2p ECM pack with contact R; SPM pack 2p. For
   digit-diff, core runs in **legacy-Euler mode** (same scheme, same dt). Registered band, written BEFORE the run:
   max |ΔV| ≤ 1e-12 V, max |Δstate| ≤ 1e-12 (rel, with an abs floor of 1e-15 for states crossing zero — rel is
   undefined at sign changes, and z-modes DO cross zero; floor ASSUMED 2026-07-09, tighten if a gate trips on it)
   over the full trajectory.
3. **Scheme upgrades validated separately** (never against loose-tolerance references — pouch-cell trap, logged
   twice): exponential propagator vs a tolerance-CONVERGED reference (tighten until the reference moves < 0.01 mV);
   registered band: |V_expm − V_ref| < 0.1 mV over a 1C discharge with a current step.
4. **Oracles on non-degenerate cases:** parity scenarios include asymmetric electrodes, nonzero contact R,
   heterogeneous initial SOC — never only the symmetric/uniform case. For COMPOSED kernels (diffusion + thermal +
   ageing coupled) where no analytic series exists, the oracle is the **Method of Manufactured Solutions** (pick a
   solution, derive the source term that makes it exact — Roache; Salari & Knupp): it covers the coupling terms
   that P1-G3's single-physics series cannot, and is mathematics independent of both parity and CVODE references.
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
   - **PAY-4 (Phase 7 exit; cross-tool positioning, quiet machine — registered TARGETS, not abort gates; Volkan
     2026-07-09: "1000s× faster than liionpack, faster than PyBaMM at single cell"):**
     (a) single-cell SPM, 1C discharge + CCCV cycle, vs PyBaMM 26.x on its FASTEST path (IDAKLU) — SLIDE faster
     both per-solve AND including setup (PyBaMM pays seconds of model build; a compiled SLIDE composition builds
     in µs). Hypothesis [inferred: expm exact on CC vs adaptive BDF, pouch-cell 26–52× precedent; no DAE
     machinery]: ≥10× per-solve. (b) pack: identical netlist + SPM cells, 16s4p and ~100p, vs liionpack —
     hypothesis [inferred: compiled SoA batches + cached sparse workspace vs Python orchestration + per-cell
     casadi solves + per-step MNA in numpy/scipy]: ≥10³×. Protocol honesty: identical model + documented
     tolerances; setup and solve timed separately; liionpack is maintenance-mode (§8 SOTA row) so this is
     positioning, not a moving target; failure ⇒ profile and record, not stop.
   Structural proxies (allocations, virtual calls, iteration/factorisation counts, bytes/cell) are tracked
   continuously as leading indicators; wall clock at the checkpoints is the confirming evidence. Rule: phase N+1
   does not start before phase N's checkpoint passes or Volkan explicitly waives it (PAY-4 exempt from the
   phase-blocking rule — it is a positioning measurement).

## 6. Phased roadmap

Effort tags are relative. Every phase ends: tests green, CHANGELOG updated, §8 ledger updated, small commits pushed.

### Phase 0 — Stabilise legacy (DONE 2026-07-07, incl. follow-ups P0-C1..C6 — §8)
Fix §2.4 bugs. Agent A (solver/state): A1–A7 in `src/modules/*` + regression tests + CHANGELOG. Agent B (data/misc):
B1–B5 in `src/recording/`, `src/types/Histogram.hpp`, `src/procedures/Cycler.cpp`, degradation docs + tests;
CHANGELOG lines to `.claude/changelog-phase0b.md` (avoid merge conflict; architect merges).
**Gate P0-G1:** full ctest delta vs recorded baseline shows only intended changes; each fix has a test that failed
before the fix.

### Phase 1 — Core data model + SPM kernels (the keystone)
Deliver: `StateArena`/`StateSpec`/`StateSlice`/`BatchBuilder`, `Domain`/`ElectrodeParams`, `SpectralDiffusion<NCH>`,
`ThermalLumped`, ageing kernels (SEI/LAM/CS/plating ported mechanism-by-mechanism), composition registry + factory,
legacy-Euler stepping mode. Single-cell `Simulation` façade. Arena scalar behind `using real_t = double;` (Q1);
`BatchBuilder` reserves the cross-batch thermal-flux seam — declaration only, D-21 designs it before Phase 2 (Q9).

**Progress + critical path (2026-07-09 review, [confirmed] by code inspection of `src/core/` + `tests/unit/core_*`):**
DONE: `StateArena`/`BatchBuilder` (+`q_ext` seam, tests), `SpectralDiffusionLegacy` parity kernel (P1-G0 both configs,
Q8 closed), production `SpectralDiffusion<NCH>` (heterogeneous 8-lane oracle test), rebindable `BatchView`/`StepCtx`
and ODE-row roles (§3.12/D-23), the shared full particle-concentration observable (surface/interior/centre: Debug
legacy parity exact, Release 2.665e-15, and 9.027e-15 analytic round-trip), all of P1-G3, and the §3.3
`Domain`/physical-description/`ElectrodeParams` layer, and the shared kinetics/OCV/resistance/voltage/heat
observable stage, `ThermalLumped` with arena-owned heat/time accumulators, SEI ageing mechanisms 1–4, surface-crack
mechanisms 1–5, LAM mechanisms 1–4, porosity/diffusivity coupling, and shared Dai/Laresgoiti stress observables with explicit
previous-step arena state, plus the compile-time fixed composed SPM RHS pipeline (mandatory zeroing, shared
observable/stress stages, rebind-safe trial-vector evaluation). NOT STARTED: the rest of the Deliver list — registry + factory,
`Simulation` façade — and
gates P1-G1/G2/G4.
**DEPENDENCY (surfaced 2026-07-08 handoff): the observable-reconstruction layer.** Thermal and ageing kernels are not
self-contained state→state maps: they need `c_surf = C·z + D·flux` (+ centre node, the §2.2 output path), Butler-
Volmer overpotentials, OCV/entropic-coefficient interp, Rdc — i.e. the derived-observables layer (D-10, §3.7) plus
the `BatchView`/`StepCtx` kernel interface (§3.11) that `SpectralDiffusion` deferred (plain spans as stopgap).
**Ordering:** (1) observable layer + `BatchView`/`StepCtx` — **DONE 2026-07-10: rebindable views, row roles,
full concentration, kinetics, OCV, resistance, voltage, and heat** → (2) P1-G3 Chebyshev
oracle FIRST (it validates exactly the C/D surface-concentration path the new layer exposes) → (3) `ThermalLumped`
(as `addRhs`, §3.12) — **DONE 2026-07-10** → (4) ageing kernels + fixed composed pipeline (same form; promote legacy `_prev`/accumulator members to arena rows —
§2.1 correction) — **DONE 2026-07-10** → (5) registry/factory + façade + `EulerLegacy` stepper → P1-G1/G2/G4 → PAY-1. The three
2026-07-09 core review marks and Model_SPM build audit C1–C3 were resolved on 2026-07-10 (see §8 rows).
**Gates:** P1-G0 parity-drift pilot (Q8): 1-cell legacy-Euler 1C CC, measure |ΔV|/|Δstate| drift legacy vs v4 kernel;
outcome closes Q8 (keep 1e-12 band, or pin op-order/`-ffp-contract=off`, or loosen with ulp argument) BEFORE P1-G1
runs. P1-G1 parity single-cell scenarios (§5.2 band as resolved by P1-G0). P1-G2 one batch of 10⁴ identical cells
steps with ZERO per-step heap allocations (assert via allocation-counting new). P1-G3 Chebyshev external oracle:
transient sphere diffusion with constant-flux BC vs the analytic series solution (Carslaw & Jaeger form), registered
band rel. err < 1e-6 at nch=5,8,12; plus cross-check vs Howey Spectral_li-ion_SPM conventions. **Extended
2026-07-09 (math audit `.claude/reports/chebyshev-math-audit-2026-07-09.md`):** (a) analytic EIGENVALUE oracle —
the folded operator's nonzero eigenvalues are known in closed form: λ_k·R² = −μ_k² with μ_k the roots of
tan μ = μ (μ₁ ≈ 4.4934, μ₂ ≈ 7.7253; μ₀ = 0 is the mass mode, which is WHY forcing one zero eigenvalue is
correct). Registered: rel err < 1e-10 for k ≤ nch/2 at nch = 5, 8, 12; all eigenvalues real (ratio < 1e-12) and
negative. Validates the operator independent of ANY trajectory — different mathematics from both parity and the
C&J series (which shares the same μ_k but tests A,B,C,D jointly).
**IMPLEMENTED + band FALSIFIED/corrected 2026-07-09 (Opus)** — `tests/unit/core_ChebyshevEigenvalues_test.cpp`,
14th test, green. The operator IS exactly the tan μ = μ spectrum (μ₁ matches to 3.6e-14 at nch=12), R1/R2/R3
confirmed. But the ASSUMED "1e-10 @ k ≤ ⌈nch/2⌉" is **FALSIFIED** — accuracy is the honest Chebyshev
spectral-convergence curve, not a flat floor. MEASURED (registered pre-run, recorded per CLAUDE.md §3): fundamental
μ₁ rel err 1.98e-5 (nch=5) → 2.25e-10 (nch=8) → 3.63e-14 (nch=12), ~1.5 digits/node; modes resolved to rel < 1e-3
= 1 / 3 / 6 (≈ ⌈nch/2⌉ — Fable's mode COUNT was right, the TOLERANCE should have been ~1e-3, not 1e-10). Notable
[confirmed]: at the production default **nch=5 even the fundamental diffusion eigenvalue is accurate to only ~2e-5**
(sub-mV on voltage, fine, but now quantified — a possible reason nch=5 was the historical "known-good" value).
Re-registered bands (test-enforced, ~3× margin): fundamental < {5e-5, 1e-9, 1e-12}; resolved-to-1e-3 count ≥
{1,3,6}; μ₁ rel err strictly monotone-decreasing in nch (spectral-convergence signature). **P1-G3(b) IMPLEMENTED
2026-07-10:** `tests/unit/core_ChebyshevTransient_test.cpp` checks the Carslaw–Jaeger/Crank constant-flux transient
at every surface/interior/centre node for both electrodes and nch={5,8,12}. The registered rel≤1e-6 band holds;
worst at τ=0.2 is {8.254e-7, 2.356e-11, 3.949e-12}; τ=1.0 all ≤4.661e-13. This validates A/B/C/D jointly and
closes centre-path audit C4. It also records the SLIDE convention D∂c/∂r=-j (positive j depletes). **P1-G3 COMPLETE.** (c) v4
`build()` **IMPLEMENTED 2026-07-10** in `src/core/SpectralModel.hpp`: adopts (a) as a Status-failing gate, replacing
the Release-silent assert + unchecked `EigenSolver` in the v4 path (audit C1–C3 closed). The remaining extension is an MMS check on the
composed RHS (§5.4). **P1-G4 (added 2026-07-09) bitwise restart:**
run 0→2T equals run 0→T + arena snapshot/restore + T→2T DIGIT-FOR-DIGIT — the cheapest test that catches hidden
non-arena state (the §2.1-correction class: `s_dai_p_prev`, `Therm_Qgen`/`Therm_time`) and workspace-invalidation
bugs; HPC checkpoint-restart discipline.

### Phase 2 — Pack layer
Deliver: netlist combinators + `Pack::compile()` (flatten, sparsity, ladder detection, index-1 check), Mode A sparse
Newton (Eigen SparseLU; KLU optional), Mode B Thomas ladder, Thevenin batch interface, `SolverWorkspace` (§3.4.1)
with warm start, chord refresh policy, and invalidation contract.
**Gates:** P2-G1 parity vs legacy `Module_s`/`Module_p` on 3s2p (band §5.2) + rollback-invalidation test (solve after
restore digit-matches cold solve). P2-G2 Mode B ≡ Mode A on ladders to
1e-10 A. P2-G3 nested-constructed pack (p-in-p-in-s) compiles flat and solves in ONE Newton loop (no nested
iteration), Newton iterations ≤ 8 on the 4p heterogeneous-resistance case that historically blew up. P2-G4 workspace
efficacy: on a 100-step 4p CC segment, count of numeric factorisations ≤ 10 (vs 1 per iteration today) at identical
converged currents (1e-10 A) — an iteration/factorisation COUNT, not a timing (§5.6). P2-G5 Mode B-ODE (D-20)
admission — registered BEFORE implementation: branch currents vs Mode A arbiter on (a) linear ECM 16p, R spread 1%:
max |ΔI| ≤ 1e-8 A (paper's home regime — must be exact); (b) SPM 16p, `varied()` 2% capacity + 5% resistance
spread, 1C CCCV cycle: max |ΔI|/I_branch ≤ 1e-3 over trajectory; (c) same at 256p. Failure of (b) or (c) does NOT
kill the fast path — it CONFINES it: compile() regime check then restricts Mode B-ODE to the measured-valid
envelope and logs why; falsification is a deliverable, record the numbers either way.

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
diagnostics); CCCV parity vs legacy Cycler. Terminations located by event root-finding (g(y)=0, §3.12 item 5), not
post-step threshold checks; `Solution` carries a machine-readable termination REASON (event/limit/error/final-time —
PyBaMM `solution.termination` equivalent, required for Phase-7 compat anyway).

### Phase 6 — Recording & I/O
Deliver: Recorder (snapshot cadence + lazy derived), CSV sink, mmap binary sink (header-CRC idiom), optional Parquet.
**Gate P6-G1:** derived-vs-stored equivalence test (V recomputed from snapshot == V recorded live to 1e-12);
mmap file survives the hardened-open validation tests (truncated/corrupt-header cases).

### Phase 7 — Python bindings + PyBaMM compat
Deliver: nanobind module, wheels (scikit-build-core), `Experiment/ParameterValues/Simulation/Solution`, Chen2020
absorption table, options→registry map, `device=` dispatch (dtw-cpp template), pytest in CI. **BPX JSON reader**
(Faraday Institution standard) as a second parameter-absorption source next to the PyBaMM-name table — the
ecosystem-neutral interop PyBaMM and BattMo already speak; one cold-path parser. Also surface the §3.1 ensemble
framing in the Python API: `varied()` lanes ARE a UQ/parameter-sweep engine (10⁵ independent single-cell variants =
one batch — DifferentialEquations.jl `EnsembleProblem` shape), not only manufacturing spread in packs.
**Gate P7-G1:** the PyBaMM getting-started Tutorial-5 experiment script runs with `import slide as pybamm`-style swap
and produces a voltage curve within a registered band of PyBaMM's own SPM (band set after a converged-reference run,
expected ~10 mV model-difference scale — document, don't hide, the modelling differences).
**Gate P7-G2 (PyBOP, added 2026-07-09):** (a) forward sensitivities vs central-finite-difference arbiter on the same
trajectory — registered band set per parameter BEFORE the run (FD step chosen by the standard √ε·scale rule, checked
for FD-noise floor); (b) one end-to-end PyBOP fitting example (e.g. GITT-style D_s + R identification on synthetic
SLIDE data with known truth) recovers the truth within a registered tolerance using a GRADIENT-based optimiser — this
is the "entered the ecosystem" proof, not an API checkbox.
**Gate P7-G3 (parameter fidelity, added 2026-07-09):** the §3.9 same-parameters⇒same-results contract — absorption
round-trip EXACT on every Chen2020/BPX key; fixture parity on the registered scenarios (bands per §3.9).

### Phase 8 — MATLAB MEX, GPU, docs, release
MEX `+slide` package symmetric with Python; CUDA one-cell-per-thread batch stepping (host-side coupling); docs site
update (installation, quickstarts ×3 languages, "add a cell model", "add an ageing mechanism"); v4.0.0 SemVer release,
CHANGELOG consolidation. Gates defined when phase opens.

### Beyond v4.0 (recorded so design choices don't foreclose them)

- **WASM GUI (Volkan, 2026-07-07):** Emscripten build of the dependency-free core + a browser front-end. Costless to
  keep open: core already must build with zero deps (non-negotiable #4), no threads assumed outside the pool (§3.8),
  no filesystem dependence in the hot path. Only rule it adds NOW: no platform API in core without a portable seam.
- Thermal 1D/2D per-cell models; blended electrodes (`vector<ActiveMaterial>`); SYCL/HIP (Q5); f32 storage (Q1).
- **FMU export (FMI Model-Exchange):** §3.12's {rhs, observables, event indicators} is already the FMI ME shape —
  export is a thin wrapper later; only rule NOW: keep the contract FMI-congruent (costs nothing).
- ~~Forward sensitivities / dual-number AD~~ **PROMOTED into v4 scope 2026-07-09** (Volkan: PyBOP entry required) —
  now §3.9 sensitivities + D-24 + P7-G2 + Q10. Only the ADJOINT (many-parameter gradients) stays beyond v4.0.

## 7. Open questions for Volkan (OPEN/ASSUMED ledger)

| ID | Question | Current assumption |
|----|----------|-------------------|
| Q1 | float32 state option for GPU/memory? | **DECIDED 2026-07-07 (Volkan)**: f64 everywhere in v4.0 — forced by the §5.2 parity band (≤1e-12 rel is unreachable in f32, ~1e-7 ulp); memory fine (10⁵ × ~300 B ≈ 30 MB ≪ PAY-3's 1 GB). Insurance: arena scalar behind a single `using real_t = double;` alias; mmap header type field width-aware — f32-storage/f64-accumulate later is a new arena instantiation, not a rewrite |
| Q2 | Keep `Cell_ECM<N_RC>` template generality in v4 core? | **DECIDED 2026-07-07 (Volkan)**: yes — entailed by registered gates (P2-G5a linear-ECM arbiter, §5.2 3s2p ECM parity case, §3.4.1 Tier-0 constant-Jacobian shortcut). Cost = one registry entry (D-02) |
| Q3 | Legacy API: keep as façade over core after parity, or hard-break at v4.0? | **DECIDED 2026-07-07 (Volkan)**: façade through v4.x, delete in v5 — the parity harness is every phase gate's arbiter (D-14, §5.2) and needs legacy compiled+runnable through Phase 8 |
| Q4 | KLU/SuiteSparse as optional dep acceptable? | **DECIDED 2026-07-07 (Volkan)**: yes (optional, Eigen SparseLU default) — matches non-negotiable #4. Condition: Phase-2 CMake detection degrades silently to Eigen (no configure failure) on all 3 platforms (§1 goal 7) |
| Q5 | GPU: CUDA-only first? | **DECIDED 2026-07-07 (Volkan)**: yes; SYCL/HIP revisit after CUDA lands. Only present-cost rule: no platform API in core without a portable seam (already required for WASM, §6 Beyond-v4) |
| Q6 | Ross's analytical parallel solution — is `setCurrent_analytical_impl` (Nilsu 2024) the code you meant, or is there a separate derivation to recover? | **RESOLVED 2026-07-07 [confirmed]**: arXiv:2508.14454 (Lone, Atlan, Fasolato, Raimondo, Drummond 2025) — Nilsu co-authored it; her code implements it. Adopted as Mode-B upgrade (D-20) |
| Q7 | PyBaMM version to target for the parameter absorption table? | **DECIDED 2026-07-07 (Volkan)**: latest stable at Phase 7 start; key-rename check mandatory. PyBaMM is CalVer — pinning today buys nothing; absorption table already keyed to verified 26.6.2.0 names + deprecation aliases (§3.9) |
| Q8 | Parity band (§5.2, ≤1e-12 rel) vs SoA kernel floating-point reassociation: different summation order + FMA contraction give O(1 ulp)/step, amplified over 10³–10⁴ Euler steps; 1e-12 rel ≈ 4 ulps. Pin legacy op-order + `-ffp-contract=off` in parity mode, or loosen the band? | **RESOLVED 2026-07-07 [confirmed] — band KEPT at 1e-12; parity mode uses an op-order-pinned kernel.** P1-G0 pilot ran (`tests/unit/core_P1G0_pilot_test.cpp`, kernel `src/core/SpectralDiffusionLegacy.hpp`): 1200×1 s lockstep 1C steps, legacy `Cell_SPM` Euler vs op-order-replica core kernel → **max_abs = 0, max_rel = 0, bit-identical** (H0 confirmed, registered pre-run; ctest 12/12). Why cheap: the modal update `dz_k = D·A_k·z_k + B_k·j` is diagonal — NO dot products, so op-order pinning costs nothing. Standing conditions: (1) §5.2 parity runs use the legacy-shaped kernel (production vectorised kernels are validated via §5.3 converged-reference bands, NOT the 1e-12 digit-diff); (2) ~~pilot ran in Debug/-O0 — re-confirm drift==0 in the Release config before P1-G1 sign-off~~ **DISCHARGED 2026-07-08 [confirmed]**: Release/-O3 (`build-release`, clang-21) re-run gives **max_abs = 2.26e-17, max_rel = 5.63e-15** → decisive 1e-12 gate HOLDS (~3 orders margin); H0 (exact bit-identity) FALSIFIED under -O3 cross-TU FMA contraction (legacy update lives in the prebuilt `src` lib, kernel is header-only in the test TU; `-Ofast`/`-ffast-math` contracts them differently). `-ffp-contract=off` NOT applied: on the test target alone it cannot reach 0 (legacy `src` side stays contracted), and forcing it globally would recompile legacy — PLAN §5.1 forbids that. Drift is bounded (modal Euler map non-expansive for stable modes; mean-mode accumulation ~1e-15 over a 10-cycle run, ≥3 orders below the gate). Pilot's `max_abs==0` CHECK scoped to Debug (where it holds); decisive `rel≤1e-12` REQUIRE is the CI gate in both configs. P1-G1 UNBLOCKED — register its scenarios mid-SOC or check the steep-OCV tail (dV/dcs amplification watch-point) |
| Q10 | Which parameters get first-class forward sensitivities in v4 (P7-G2)? Cost is per-parameter (dual-number sweep ≈ +1× per θ), so the set should be the fitting-relevant one, not everything | ASSUMED (2026-07-09, from typical PyBOP/GITT practice): solid diffusivities D_s(p,n), reaction rate constants k_ct(p,n), film/ohmic resistance, electrode capacities/stoichiometry limits (x_0, x_100), thermal h_conv — ~10 params. Volkan to confirm/trim when Phase 7 opens |
| Q9 | Pack-level thermal coupling has no compiled representation (`Pack::compile()` emits an electrical netlist only), but legacy modules exchange heat between children + `CoolSystem`, and T-states sit inside P2-G1's parity band | **DECIDED 2026-07-07 (Volkan)**: reserve the seam now, design later — Phase 1 `BatchBuilder` reserves a cross-batch thermal-flux interface (declaration only, no implementation); a full D-21 (thermal adjacency compiled alongside the netlist) must be written and logged in §4 BEFORE Phase 2 opens. P2-G1 is unpassable until D-21 exists |

## 8. Status ledger

| Date | Item | State |
|------|------|-------|
| 2026-07-07 | 3-scout audit (architecture/bugs, numerics, external) | DONE — findings in §2, full reports in session transcripts |
| 2026-07-07 | Chebyshev nch≠5 root cause | SOLVED pre-session on branch `Claude` (§2.2); oracle test still MISSING (P1-G3) |
| 2026-07-07 | PLAN.md v1 written | DONE (this file) |
| 2026-07-07 | Phase 0 (A1–A7, B1–B5) | DONE — 11 commits on `Claude` (7529039…91faaf9); A6 refuted (non-bug, documented); regression tests added (Module_p_phase0, CellDataStorage, Histogram, Cycler_energy); CHANGELOG consolidated. Baseline had 4/6 test binaries failing PRE-EXISTING (§2.5 P0-C2/C3); agents' tests all pass; no new failures introduced |
| 2026-07-07 | Phase 0 follow-up (P0-C1 fmt/clang21, P0-C2 ocv_coefs, P0-C3 thickp, P0-C5/C6 Cycler) | DONE — full ctest 10/10 GREEN (baseline this session: 4/10 failing {Cell_Bucket, Cell_ECM, Cell_SPM, Module_p}, prediction hit exactly). P0-C1: fmt 11.2.0, full 1144-target clang-21 build green. P0-C2: default `ocv_coefs` was a truncated 3-term fit (author's own 8-term fit found in b0a1c82's integration test; even that spans only 2.53–3.46 V, cannot span the 2.7–4.2 V window) → default now empty ⇒ `getOCV()` falls back to OCV-table interp; poly stays opt-in via `set_ocv_coefs`; test expectations 3.15→3.45 V (3.15 dated from the DELETED standalone Cell_Bucket with 2.0–4.3 V ramp). P0-C3: thickp/thickn recalibrated intentionally in 147ee4c ⇒ TEST was stale, updated (incl. CS and CSurf values derived from x_init·Cmax: 35562.14/14694.92, closed-form matched actuals to all printed digits). P0-C5: dead `v_now` removed (Cycler::setCurrent is private; single internal caller). P0-C6: CV energy now trapezoid + `th.time()` accumulated in CV (was never incremented ⇒ CV/CCCV time throughput 0 s); registered test failed pre-fix exactly as predicted (time 0, end-of-step Riemann 3.49986 vs 3.5±1e-6 band). Also fixed en route: Module_p test SPM-SOC tolerance 1e-15→5e-5 (SPM SOC is Li-fraction-estimated, not coulomb counting — intentional since v3) |
| 2026-07-07 | Solver-memory design (§3.4.1, D-18, P2-G4) | DONE — Volkan clarified the legacy statics were intentional quasi-Newton Jacobian memory; design keeps the memory, adds invalidation + thread safety. Agent cap now ≤2 (header) |
| 2026-07-07 | Pack description layer (§3.4.2, D-19) | DONE — combinator tree + Netlist escape hatch; compile() erases authoring shape. Per Volkan: agents = Opus HIGH (not xhigh), no "ultrathink" in agent prompts |
| 2026-07-07 | SOTA verification (report: `.claude/reports/sota-verification-2026-07-07.md`) | DONE — C1/C3/C5/C6 CONFIRMED, C2/C4 NUANCED (liionpack maintenance-mode; WR+Baumgarte combo unpublished = our synthesis), C7 PyBaMM 26.6.2.0 keys verified. Q6 RESOLVED, D-20 added (arXiv:2508.14454 Mode-B upgrade). Description-layer zero-cost pattern confirmed (CasADi/Eigen/Halide precedent) provided erasure is total |
| 2026-07-07 | §3/§4/§7 review + Q1–Q7 sign-off | DONE — Fable architecture review (agent acc785905bea7cbdf): all six ASSUMED defaults survive adversarial reading (Q1/Q2/Q3 entailed by registered gates); Volkan accepted all six. Two gaps found OUTSIDE the Q-ledger, logged as Q8 (parity band vs FP reassociation — OPEN, pilot-first rule, gates P1-G1) and Q9 (pack thermal coupling missing from compile() — seam reserved Phase 1, D-21 due before Phase 2). R3 (archetype fragmentation) closed with stated assumption + compile() diagnostic in §3.1. Phase 1 UNGATED |
| 2026-07-07 | P1-G0 parity-drift pilot (Q8) | **PASSED — drift exactly 0 (bit-identical)** over 1200×1 s 1C lockstep steps, legacy `Cell_SPM` Euler vs `SpectralDiffusionLegacyKernel` on `StateArena`. H0 (==0) registered pre-run and CONFIRMED; Q8 closed, §5.2 band stays 1e-12. ctest 12/12 (baseline 11/11 + pilot). Residual: re-confirm in Release config before P1-G1 sign-off |
| 2026-07-08 | P1-G0 Release re-confirm (Q8 standing condition 2) | **DISCHARGED.** `build-release` (clang-21, -O3/-Ofast). Release: max_abs=2.26e-17, max_rel=**5.63e-15** → decisive rel≤1e-12 HOLDS (~3 orders margin). H0 exact bit-identity FALSIFIED under -O3 cross-TU FMA (registered falsification, numbers recorded); Debug/-O0 still max_abs=0. Legacy NOT recompiled (-ffp-contract=off would need global legacy change, §5.1 forbids; test-target-only can't reach 0). Pilot's `max_abs==0` CHECK scoped `#ifndef NDEBUG`; both configs green (Release 1/1, Debug 1/1). **Q8 FULLY CLOSED; P1-G1 UNBLOCKED.** Watch-point logged: register P1-G1 mid-SOC or check steep-OCV tail |
| 2026-07-08 | Production `SpectralDiffusion<NCH>` kernel (§3.2/§3.5) | DONE — vectorised-across-lanes forward-Euler diffusion on SoA `StateArena` rows; `DiffusionParams<NCH>` (batch-shared A/B/D0/D_T/a/thick/sgn), per-lane T/i_app spans, once-allocated D_eff/flux scratch (PC-1). Validated vs legacy-shaped oracle on a heterogeneous 8-lane batch (600 steps): **Debug max_abs=0** (H_math arbiter — vectorised sweep == legacy math exactly), **Release max_rel=3.8e-15** (Q8 rel≤1e-12 HOLDS, ~3 orders margin; sub-ulp vectorisation reassoc, expected per Q8). ctest Debug 13/13 (+1, no regressions). DEFERRED (needs review): §3.11 BatchView/StepCtx bundling — T/i_app passed as spans for now |
| 2026-07-09 | Progress review vs plan (Fable; advisor unavailable, no external review) | DONE — Phase-1 state written into §6 Phase 1 "Progress + critical path". Three defects were marked for follow-up: diffusion lane-count OOB, foreign-slice/moved-arena safety, and a stale Release result comment. **RESOLVED 2026-07-10** by the rebindable-view increment and regression tests. Remaining cold-path design notes: `BatchBuilder` invalid input becomes a Status failure when the factory lands; the factory also replaces tests' hand-copied `elec_surf` geometry literal. |
| 2026-07-09 | Kernel/integrator contract designed (§3.12, D-22) | DONE — Volkan's requirements (flexible/expandable, no unnecessary abstraction, extreme speed, OFF-THE-SHELF integrators, sparse intuitive structure) resolved the Phase-1 observable-layer blocker: RHS-form kernels (`addRhs` into a once-allocated ydot arena) over `BatchView`/`StepCtx` structs-of-spans; observables = free functions shared by RHS + Recorder (one code path, D-10); steppers own time — EulerLegacy (parity) / ExponentialModal (Phase 3) / generic RHS adapter (`arena.raw()` zero-copy to CVODE `N_VMake_Serial`, Boost.odeint, or user callable; optional deps, core builds with none). CVODE at tight rtol doubles as the §5.3 converged-reference arbiter for expm validation. NOT reviewed externally (advisor down; Fable-only design) — Opus should read §3.12 critically before implementing |
| 2026-07-09 | Orthogonal review of the plan (independent Fable critic agent a4c48623f018bf968 + cross-simulator harvest) | DONE — 13 defects, 3 critical, ALL in the fresh §3.12 (caught BEFORE implementation): D1 CVODE clones internal vectors → BatchView must REBIND per eval (the "zero-copy N_VMake_Serial" claim as first written was WRONG — recorded, corrected); D2 ydot zeroing contract was unstated (accumulate semantics + legacy shared-row adds `Cell_SPM_dstate.cpp:210` = guaranteed wrong answer); D3 [confirmed by grep] legacy hides path-dependent state outside `State_SPM` (`s_dai_p_prev` dstate:244 → LAM degradation:377; `Therm_Qgen`/`Therm_time` Cell_SPM.hpp:62-64) falsifying §2.1's "only time/Ah/Wh" — §2.1 corrected, v4 rule added (later-step reads ⇒ arena row). Fixes folded: §3.12 rewritten (D-23: rebindable view, eval pipeline, ODE-row mask, events-as-rootfinding, CVODE=arbiter-only, odeint REJECTED, scalar-generic kernels), §3.4 solver robustness (consistent init, SPICE-style trial limiting, source-stepping homotopy rescue, rollback atomicity, compiled-netlist detection), §3.4.1 divergence guard, §3.5 expm1/φ₁ cancellation note, §3.8 fixed-order reductions + GPU-fit analysis (2·10⁴ cells ≈ 6 MB state; occupancy/coupling/launch/f64 concerns recorded — answers Volkan's GPU question), §5.2 band abs floor (ASSUMED 1e-15), §5.4 MMS oracle for composed kernels, P1-G4 bitwise-restart gate, P5-G1 events + termination reason, Phase-7 BPX + ensemble API, §3.1 dual-use (packs AND UQ sweeps), Beyond-v4 FMI/AD notes |
| 2026-07-09 | P1-G3 part (a) IMPLEMENTED (Opus) — Chebyshev eigenvalue oracle | DONE — `tests/unit/core_ChebyshevEigenvalues_test.cpp` (14th test, full suite 14/14 green, no regressions vs the 13/13 baseline). Independent-math oracle: `Model_SPM` A-eigenvalues vs analytic roots of tan μ = μ. Operator CONFIRMED (μ₁ to 3.6e-14 @ nch=12; one zero mass mode; rest real+negative). Fable's ASSUMED "1e-10 @ k ≤ ⌈nch/2⌉" band FALSIFIED with data (recorded, CLAUDE.md §3): true spectral convergence, fundamental 2.0e-5/2.3e-10/3.6e-14 @ nch 5/8/12, ~⌈nch/2⌉ modes to 1e-3. nch=5 fundamental only ~2e-5 accurate [confirmed]. Bands re-registered from the run. STILL OPEN: P1-G3 part (b) C&J transient series incl. centre node. Marks C1–C3 (Model_SPM build() hardening) still deferred to the v4 model build |
| 2026-07-09 | Chebyshev math audit + speed targets + async recording + parameter fidelity (Volkan directives) | DONE — (1) **Chebyshev audit** `.claude/reports/chebyshev-math-audit-2026-07-09.md`: full derivation trail (u=rc → heat eq → odd folding → surface-node Schur condensation → nonsymmetric eigensolve → modal form) [confirmed vs Model_SPM.hpp]; KEY FIND: operator eigenvalues are analytic — λ_k·R² = −μ_k², tan μ_k = μ_k (mass mode = μ₀ = 0 explains the forced zero eigenvalue) → new P1-G3 eigenvalue oracle (rel 1e-10, k ≤ nch/2) + centre-node c(0,t) coverage (cc_coeff path had NO external oracle) + v4 build() Status-gate replacing Release-silent assert/unchecked EigenSolver (marks C1–C3 in Model_SPM.hpp, grep `REVIEW-MARK(2026-07-09`); discretisation menu added §3.2 (FVM/parabolic/Duhamel/Legendre). (2) **PAY-4** cross-tool targets (≥10³× liionpack, ≥10× per-solve vs PyBaMM IDAKLU single cell — hypotheses [inferred], positioning not abort gates). (3) **§3.7 async recording ring** (pinned buffers + side stream, Blosc-style shuffle+zstd, explicit backpressure, non-temporal snapshot copies) + §3.8(3b) GPU kernel fusion lever. (4) **§3.9 same-params⇒same-results contract** + P7-G3 (absorption round-trip exact; committed PyBaMM fixture parity; converged-model band honesty). Ageing-embedding question answered: declare()/StateSlice + addRhs IS the embedding design (§3.1+§3.12); ageing port itself still not started |
| 2026-07-09 | PyBOP ecosystem entry (Volkan directive) + GPU speed/memory analysis | DONE — forward sensitivities promoted from beyond-v4 into v4 scope: §3.9 sensitivities surface (dual-number forward mode via §3.12 scalar-generic kernels; analytic modal sensitivities on the same exponential propagator where linear structure permits; central-FD arbiter), D-24 (adjoint deferred; FD-only and required-AD-dep rejected), P7-G2 gate (sensitivity band vs FD arbiter + end-to-end PyBOP gradient-fit recovering known truth on synthetic data), Q10 (fitting-parameter set, ASSUMED ~10 params, Volkan confirms at Phase 7). §3.8 GPU: cache-residency roofline recorded [inferred, no timing claim] — 20k-cell working set fits GPU L2; "very quick very long" = expm × multirate × GPU product, hypothesis registered for GPU-phase PAY banding |
| 2026-07-10 | Phase-1 rebindable view + first observable | DONE — `BatchView`/`RhsViews` rebind arbitrary integrator trial vectors without copies; mandatory ydot-zero stage and ODE/algebraic/cumulative/input row roles implemented (D-23). Concentration reconstruction is one scalar-generic free function shared by future RHS/Recorder paths; hardened arena moves/slice bounds and diffusion lane bound. Initial full Debug suite 15/15 green. |
| 2026-07-10 | P1-G3(b) analytic transient + centre observable | DONE — full surface/interior/centre concentration reconstruction matches legacy exactly in Debug and to 2.665e-15 in Release, with a 9.027e-15 all-node uniform round trip. Independent constant-flux spherical series validates A/B/C/D and `Cc/cc_coeff` jointly at nch={5,8,12}; registered rel≤1e-6 holds (worst 8.254e-7 at nch=5, τ=0.2; nch=8/12 ≤3.949e-12). P1-G3 is complete; full Debug suite 16/16 and affected Release tests 2/2 green. |
| 2026-07-10 | Physical description hierarchy + canonical Domain | DONE — added value-semantic ActiveMaterial/Electrode/Separator/Electrolyte/Thermal/Cell designs and hot `ElectrodeParams`. Resolved a plan-vs-legacy ordering hazard: scoped v4 Domain is negative-first; legacy is positive-first; all bridges map explicitly and tests prevent silent electrode swaps. Full Debug suite 17/17 and affected Release tests 4/4 green. |
| 2026-07-10 | Compiled parameter-curve forms (D-16) | DONE — exact nonuniform piecewise-linear tables now use an O(1) uniform segment-index accelerator and reproduce legacy Kokam OCV interpolation bit-for-bit; smooth injected curves compile to a validated ≤4096-point uniform LUT. Malformed/tolerance-failing builds return the common Status channel and invalidate prior data atomically. Also fixed missing self-contained includes in `FixedData.hpp`. Full Debug suite 18/18 and Release curve target green. |
| 2026-07-10 | Shared SPM electrical observables | DONE — canonical negative-first arena layout now includes all evolving inputs needed by reconstruction (`D`, thickness, active area, SEI/electrode/current-collector resistance). One scalar-generic stage computes concentrations, stoichiometry, exchange current, overpotentials, electrode/cell OCV, resistance, terminal voltage, and heat without hot allocations. Legacy parity at zero/discharge/charge and two temperatures: OCV/R exact, max voltage error 4.441e-16 Debug and zero Release. Invalid concentrations return Status. Full Debug suite 19/19; affected Release 4/4. |
| 2026-07-10 | `ThermalLumped` + fast-math-safe validation | DONE — scalar-generic additive RHS consumes shared internal heat plus `q_ext` and convection; rho*Cp*V and h*A are cold-compiled. Generated energy and thermal elapsed time replace hidden legacy members with snapshot-safe arena rows. Analytic heterogeneous-lane balance and invalid-input tests pass Debug/Release. Release testing exposed that `std::isfinite` is folded away by `-Ofast`; core now uses an IEEE exponent-bit check with an integer compiler barrier, also applied to compiled curves and covered by bit-injected NaN tests. Full Debug suite 20/20; affected Release 2/2. |
| 2026-07-10 | SEI ageing mechanisms 1–4 | DONE — scalar-generic masked batch stage ports kinetics-limited, linear-diffusion-limited, Christensen–Newman, and fitted variants plus optional Ashwin porosity loss. Additive RHS updates negative z modes, SEI thickness, lost lithium, active fraction, and active area; lost-lithium/active-fraction values are canonical arena state. Direct legacy `Cell_SPM::SEI` parity holds at rel≤1e-13 for every mechanism in Debug/Release. Full Debug 21/21; affected Release 2/2. |
| 2026-07-10 | Shared Dai/Laresgoiti stress stage | DONE — stress is reconstructed from the full concentration observable with no heap work; Dai maxima for both electrodes and the Laresgoiti negative stress match public legacy functions at rel≤1e-12 in Debug/Release. The previously hidden `s_dai_*_prev`, `s_lares_n_prev`, and interval are algebraic arena rows: checkpointed but excluded from ODE integration. Full Debug 22/22; affected Release green. |
| 2026-07-10 | Surface-crack mechanisms 1–5 | DONE — masked scalar-generic stage ports Laresgoiti, Dai, Deshpande–Bernardi, Barai, and Ekstrom laws plus optional negative-diffusivity loss. Additive RHS covers crack surface, crack-driven extra SEI flux/lithium loss, and `Dn`; all mechanisms × both diffusion settings match direct `Cell_SPM::CS` at rel≤1e-12 in Debug/Release. Sensitivity type is preserved below the legacy diffusion-rate cap. Full Debug 23/23. |
| 2026-07-10 | LAM mechanisms 1–4 | DONE — ports Dai stress thinning, Delacourt–Safari flux loss, Kindermann dissolution, and Narayanrao active-area loss. Additive RHS composes direct area and `3ε/R` contributions. All six raw geometry rates per mechanism match direct `Cell_SPM::LAM` at rel≤1e-12 in Debug/Release. Full Debug 24/24. |
| 2026-07-10 | Yang lithium plating | DONE — scalar-generic side-current matches direct `Cell_SPM::LiPlating` at rel≤1e-13 for charge/discharge in Debug/Release; additive RHS updates negative z, LLI, and explicit plated-layer thickness. All individual legacy SPM ageing mechanisms now have v4 kernels. Full Debug 25/25. |
| 2026-07-10 | Fixed composed SPM RHS pipeline | DONE — compile-time composition enforces zero-derivative → one shared observable pass → optional one shared stress pass → additive diffusion/thermal/SEI/crack/LAM/plating order. Diffusion consumes the observable stage's exact effective diffusivity and molar flux. Registered test proves mandatory whole-arena zeroing and trial-vector rebinding; all optional branches compile. Debug/Release pipeline gates green; full Debug 26/26. |
| 2026-07-10 | Validated per-batch spectral compiler (audit C1–C3) | DONE — independent v4 cold build supports registered nch={5,8,12}, physical per-domain radii, and Status-fails invalid geometry, eigensolver failure, complex contamination, non-invertible transforms, non-finite output, or analytic-spectrum mismatch. Default A/B/C/D, transforms, centre map, and integration matrix are exactly legacy-identical in Debug/Release; custom scaling and atomic failure covered. Full Debug 27/27. |
| — | Phases 1–8 | Phase 1 IN PROGRESS (DONE: StateArena/BatchBuilder + tests, P1-G0, production `SpectralDiffusion<NCH>`, rebindable BatchView/StepCtx + row roles, physical description/Domain/ElectrodeParams, shared concentration/electrical/heat/stress observable stages, `ThermalLumped`, all SEI/crack/LAM/plating ageing mechanisms, fixed composed SPM pipeline, P1-G3. **NEXT: registry/factory** → façade + EulerLegacy stepper → P1-G1/G2/G4 → PAY-1); Phases 2–8 not started; D-21 (pack thermal — adopt LAMMPS-style static adjacency pair list with fixed-order accumulation, per §3.8 determinism rule) still due before Phase 2 |
