# SLIDE v4 — Architecture Refactor Plan (living document)

> **Status:** ACTIVE. Last updated 2026-07-10 (Phase 9A complete). Phases 0–8 and 9A COMPLETE;
> remaining work: Phase 9B/9C (bug-hunt + simplification), 10 (PyBOP integration tests),
> and 11 (v4.0.0 release).
> **How to use this document:** single source of truth for the v4 refactor. Any session (Fable, Codex, Opus, human)
> continuing this work must (1) read this file first, (2) execute the next unblocked item in §6, (3) update §8 and
> this header. Requirements originate in `.claude/FABLE.md`. Do not re-litigate decisions in §4 without new evidence.
> Full pre-compression history (original §2 evidence base, verbose phase text, 60-row ledger) is archived VERBATIM at
> `.claude/summaries/plan-archive-2026-07-10-phases0-8.md` and in git (`PLAN.md` @ `1f18d8e`).
> **Operating rules:** ≤2 parallel agents (Volkan, 2026-07-07); small reviewable commits on branch `Claude`; every
> user-visible change updates `CHANGELOG.md` (Unreleased); wall-clock benchmarks on this machine are UNRELIABLE (user
> runs parallel jobs) — judge performance by structural arguments (allocation counts, complexity, vectorizability).
> **NEW (Volkan, 2026-07-10): simulations cost time — keep every test/gate simulation SHORT** (≤ a few hundred steps,
> seconds of wall-clock, small lane counts unless the gate is specifically about scale). Prefer analytic oracles,
> derivations, and structural counters over long trajectories. Unlimited thinking, rationed simulating.

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

## 2. Evidence base — audits (compressed 2026-07-10; full detail in the archive file)

> NOTE: citations of the form "§2.1/§2.2/§2.4 A5/§2.5" inside §3–§5 refer to the ORIGINAL evidence-base
> numbering, preserved verbatim in the archive file — not to the subsections below.

### 2.1 Legacy audit (2026-07-07) — DISCHARGED

The three-scout audit (architecture limits, bug tables A1–A7/B1–B5 + P0-C1..C6, Chebyshev nch≠5 root cause,
assets to harvest) drove Phases 0–1 and is fully discharged; complete text in the archive. Standing corrections
that must NOT be re-litigated (killed ideas stay killed):

- **A6 REFUTED — not a bug.** `Module_p::V()` consistency was derived and proved by test.
- **"Only time/Ah/Wh are path-dependent" FALSIFIED:** legacy hides `s_dai_p_prev` and `Therm_Qgen`/`Therm_time`
  outside `State_SPM`. v4 RULE: anything a later step reads MUST be an arena row. Enforced by P1-G4 (bitwise restart).
- **Chebyshev "only nch=5 works":** two defects (centre-node sign `-0.5·(-1)^N`, Eigen in-place `.inverse()`
  aliasing), both fixed; the operator is now validated by the analytic tan μ = μ eigenvalue oracle and the
  Carslaw–Jaeger constant-flux transient (P1-G3). At nch=5 the fundamental eigenvalue is only ~2e-5 accurate
  [confirmed] — sub-mV on voltage, but quantified.
- Legacy solver statics were intentional quasi-Newton memory (Volkan) — kept properly as `SolverWorkspace` (D-18).

### 2.2 Implementation audit — the 67 Codex commits (2026-07-10, Fable + review agent, artifact-checked)

**Verdict: plan followed; no fabrication or gate-rigging found.** Sampled gate tests exist and enforce the
registered numeric bands verbatim (P1-G1/G2, P2-G5, P3-G1..G3, P4-G1/G3, P7-G1/G3, P8-G1/G2 checked against
test source, file:line evidence in the audit report). Falsified hypotheses are recorded prominently, not buried
(PAY-2 10×, Chebyshev 1e-10 band, Q8 bit-identity, PAY-1's initial 0.43× and PAY-2's initial 0.69× failures).
Register-before-run held in 4/4 commit-order spot-checks (with one exception below). Non-negotiables hold:
optional deps default OFF, CHANGELOG discipline kept, no committed binaries, no repo-relative runtime paths in
core, no reverts/fixups. **Full ctest re-run this session: Debug 49/49, Release 49/49 [confirmed].**

Deviations found — each is a Phase-9A work item:

| # | Finding | Severity |
|---|---------|----------|
| AUD-1 | Resolved: P2-G1 voltage/current assertions now enforce 1e-12. Current Debug/Release maxima are 5.33e-15 V and 3.70e-13 A; the current limit has only 2.71× worst-case headroom. | resolved |
| AUD-2 | Resolved: commit `09e8ec5` introduced implementation and bands together. The 0.2 µV / 20 µA / 0.2 µAh and 2e-12 V values are explicitly classified as post-hoc regression envelopes; derivation and current reproduction are in the Phase-9A report. | resolved |
| AUD-3 | Resolved by D-27/Q11: retroactive quiet-host/operator waiver for v4.0; PAY-1/2/4 remain qualified development-host evidence, never portable guarantees. | resolved |
| AUD-4 | Resolved by D-26: mandatory CVODE arbiter waived for the exact diagonal subflow after making the analytic oracle independent in branch and operation order; original full-trajectory CVODE coverage is explicitly not claimed. | resolved |
| AUD-5 | `docs/` was still v3-era; resolved by the tested v4 user/contributor guide and docs CI in P8-G5 | resolved |

### 2.3 Architecture quality assessment (2026-07-10, Fable, [confirmed] by direct inspection)

**Not bloated, not spaghetti.** The entire v4 core is **13,151 lines across 48 files** for all of Phases 1–8.
Hot-path discipline holds structurally: zero `throw` in `src/core`, zero TODO/FIXME/HACK markers, legacy coupling
is a single `Status.hpp` include, and the only three `std::function` uses are cold-path (Experiment custom-step
callbacks, recorder drain hook). 704 REQUIRE assertions across the core unit tests.

Named debt (the Phase-9C target list — evidence, not vibes):

- **Physics triplicated:** the exact modal update (`exp(x)·z + dt·φ₁(x)·B·j`, `expm1(x)/x` with Taylor branch)
  exists independently at `SpmPipeline.hpp:472` (CPU), `CudaSpmRuntime.cu:91-92` (GPU), and through `Dual.hpp`
  (sensitivities). One scalar-generic source must feed all three.
- `ParameterSet.cpp` (1,293 lines) mixes PyBaMM absorption table + BPX JSON reader + expression AST — split.
- `Experiment.cpp` (883) mixes parser and runner; `SpmPipeline.hpp` (721) and `SpmFactory.cpp` (714) are at the
  review-size threshold.
- SEI/CS/LAM/plating kernels (`Sei.hpp` 277, `SurfaceCrack.hpp` 253, `Lam.hpp` 233, `LithiumPlating.hpp` 126)
  repeat the same mask/scratch/lane-sweep scaffolding — one idiom, four physics bodies.

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

#### D-21 compiled pack thermal adjacency — DECIDED 2026-07-10

`Pack::compile()` emits a thermal graph alongside (but independent from) the electrical netlist. Electrical nodes
and thermal endpoints are not conflated: a cell lane has one stable `CellIndex`, resolved after archetype batching
to `{batch_index, lane}`, while a thermal endpoint is either that cell index or a named boundary/coolant node.
The description layer supplies undirected cell↔cell links and cell↔boundary links with SI conductance `G [W/K]`;
geometry (`k·A/L`) is reduced to `G` on the cold path. Compile rejects self-links, duplicate endpoint pairs,
non-finite/negative conductance, unknown endpoints, and thermal links to non-thermal cells. Zero conductance is
removed. Canonical endpoint ordering plus lexicographic edge sorting makes authoring-tree shape irrelevant.

The compiled hot representation is a static pair list plus a fixed incident-list gather:

```cpp
struct ThermalEdge { std::uint32_t low, high; double conductance; }; // canonical low < high
struct ThermalIncident { std::uint32_t edge; std::int8_t sign; };    // +flux enters this cell
struct CompiledThermalGraph {
  std::vector<ThermalEdge> edges;          // lexicographically sorted, immutable
  std::vector<std::uint32_t> offsets;      // CSR offsets, one range per cell/boundary
  std::vector<ThermalIncident> incidents;  // edge-index order within every range
  std::vector<double> edge_flux;           // compile-sized ephemeral scratch
  std::vector<double> boundary_heat;       // compile-sized diagnostics/sink coupling
};
```

For each RHS stage, edge `e=(low,high)` is evaluated exactly once as
`edge_flux[e] = G[e]·(T_high - T_low)`. Edge evaluation may run in parallel because it writes disjoint slots.
Every cell then gathers its incident edge fluxes in increasing edge-index order into its arena-owned `q_ext` row;
no atomics, scheduling-order reductions, or per-step allocations are permitted. Boundary heat uses the same fixed
gather. This is the LAMMPS-style static-neighbour representation chosen by Q9, satisfies §3.8's deterministic
reduction rule, is GPU-shaped, and preserves pairwise energy exchange by reusing one stored flux with opposite signs.

`q_ext` is a derived stage input, not an integrated state: pack stepping overwrites every thermal batch's complete
`q_ext` row before any coupled RHS evaluation. Standalone batches may still supply it directly. The explicit Phase-2
stage order is `{checkpoint all arenas → solve electrical currents at current T → assemble thermal edge flux from
the same trial/accepted temperature vector → evaluate/advance every batch → commit}`. There is no hidden
electrical↔thermal fixed-point iteration. Adaptive/multistage pack steppers must reassemble from each stage's trial
temperatures; a per-batch off-the-shelf integrator instead sees `q_ext` frozen over its declared outer split, whose
splitting error belongs to D-08. On any failure all arenas (including `q_ext`) roll back and the electrical workspace
is invalidated; thermal edge scratch is reconstructible, excluded from snapshots, and simply overwritten next stage.

Mandatory Phase-2 thermal gates: (1) two-cell analytic exchange conserves pair energy to roundoff; (2) a cross-batch
edge digit-matches the same two lanes placed in one batch; (3) equivalent combinator/CSV descriptions compile to
byte-identical sorted edge/incident arrays and fixed-order `q_ext`; (4) rollback+workspace invalidation then solve
digit-matches a cold solve; (5) allocation counter remains unchanged during assemble+accepted step. P2-G1 includes
temperature and `q_ext` parity against legacy module/cooling behavior, so electrical-only parity cannot close it.

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
  Experiment additions **implemented 2026-07-10**: custom explicit/algebraic/differential steps (v25.4.0), named
  custom terminations (v25.10.0), and datetime-normalised `start_time` cut/rest scheduling (v23.9), all through
  the shared C++ runner. Solution additions: `.yp`, `.observe()` (v25.12.0) — out of scope v4.0, document as such.
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
  maps onto the composition registry (§3.2). Out of scope, documented as such: arbitrary model expression trees,
  `spatial_methods`/`var_pts`, and callbacks requiring model variables beyond the exposed time/V/I/power subset.
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
   callable for ANY user integrator; (d) a **future optional CVODE adapter**, never a production pack path
   (`SLIDE_WITH_SUNDIALS`; core ships (a)+(b)+(c) with no deps — non-negotiable #4). D-26 waives the adapter's
   mandatory v4 arbiter role. If later shipped, its scope stays small lane counts only — BDF needs
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

**Historical arbiter proposal (§5.3), superseded by D-26:** CVODE at tight tolerance could still provide a
different numerical path for the same RHS, but the mandatory v4 run was not implemented. The exact diagonal
subflow is instead checked against an independently evaluated closed form. A future CVODE adapter would arbitrate
generic-RHS interoperability and full-trajectory numerical integration; it would still not arbitrate splitting
error, which is owned by D-08's outer-step controller and separate tests.

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
| D-21 | `Pack::compile()` emits an independent, canonical static thermal pair list + fixed incident-list gather; each edge flux is evaluated once, then accumulated into arena `q_ext` in fixed edge order with no atomics (§3.4) | closes Q9 before Phase 2; preserves cross-batch heat exchange, deterministic reductions, GPU-shaped adjacency, snapshot/rollback discipline, and separation of electrical vs thermal topology | (a) hide heat exchange in recursive modules (authoring shape changes physics/ordering); (b) atomically scatter edge flux from parallel tasks (schedule-dependent digits); (c) fold thermal links into electrical MNA nodes (wrong topology/units) |
| D-22 | Kernels are RHS-form free functions over `BatchView`/`StepCtx` structs-of-spans; observables = shared free functions (one code path for RHS internals AND recording); steppers own time; generic RHS adapter exposes `arena.raw()` zero-copy to CVODE/odeint/user integrators as OPTIONAL deps (§3.12) | user requirements 2026-07-09 (flexible + expandable, no unnecessary abstraction, off-the-shelf integrators, sparse intuitive structure); arena contiguity makes flat-y interop free; PC-1/2/3/5 preserved. The later D-26 waives CVODE's mandatory arbiter role without closing this optional interoperability seam. | (a) virtual Component hierarchy (per-cell dispatch — PC-2 violation, abstraction without need); (b) integration hard-wired inside kernels (today's Cell_SPM: blocks adaptive/implicit methods entirely); (c) required SUNDIALS dep (violates non-negotiable #4) |
| D-24 | Forward sensitivities are a NAMED v4 surface (PyBOP entry, §3.9): dual-number forward mode through scalar-generic kernels as the general route; analytic modal sensitivities (same exponential propagator, linear-in-z structure; closed-form ∂V/∂states) where structure permits; central-FD as arbiter; adjoint deferred beyond v4.0; parameter set = Q10 | Volkan 2026-07-09: PyBOP ecosystem entry required; gradient optimisers need `simulateS1`-shaped sensitivities; §3.12 scalar-generic kernels make forward mode near-free to add; PyBaMM pays IDAS-sensitivity cost for the same surface | (a) FD-only "gradients" (noise floor wrecks optimiser line searches); (b) adjoint-first (right for n_θ ≫ 10, wrong for typical 3–8-parameter cell fits, much higher implementation risk); (c) full runtime AD dependency (CoDiPack/Enzyme as REQUIRED dep — violates non-negotiable #4) |
| D-23 | §3.12 REVISED after orthogonal review (2026-07-09): rebindable BatchView (adaptive integrators evaluate f at THEIR trial vectors — CVODE clones internals, zero-copy holds only at IC/writeback); fixed eval pipeline {zero ydot → one shared observables stage → addRhs chain}; ODE-row mask (algebraic I/V + kinked cumulative rows never handed to an integrator); events as g(y)=0 root-finding; external integrators scoped to small-N use; kernels scalar-generic `template<class Real>`. D-26 later waives the mandatory CVODE arbiter for v4.0. | D1/D2/D3 were guaranteed-wrong-answer defects if implemented as first drafted; legacy ageing accumulates into shared rows (`Cell_SPM_dstate.cpp:210`) and hides `_prev`/accumulator state (§2.1 correction); B5 precedent for algebraic rows; SPICE load-stage/limiting, Modelica events, SUNDIALS rootfinding, Stan/CoDiPack scalar-generic precedents | (a) view bound once to the arena (silently integrates stale state); (b) per-component observable recomputation (2–3× redundant asinh/interp per eval) or ad-hoc hidden scratch; (c) Boost.odeint adapter (state-type algebra forces copy machinery; the plain rhs hook covers odeint users); (d) custom SUNLinSol to run CVODE at pack scale (bloat — arbiter role only); (e) per-lane adaptive dt (destroys the SoA sweep) |
| D-25 | P8-G0's registered "dependency-free" wording means **free of optional runtime/toolchain dependencies**, not literally stdlib-only: Eigen 3.4 remains the required cold-path eigensolver/sparse-linear-algebra dependency. Prefer installed Eigen; use the pinned CPM source fallback otherwise. | The implemented core has always required Eigen for validated spectral compilation and sparse Mode A. Replacing both solvers during a portability verification gate would be a new numerical architecture with much higher correctness risk. The gate's actual subject is absence of CUDA/MATLAB/zstd/Arrow/legacy leakage. | (a) pretend the fallback download is zero-dependency; (b) post-hoc Eigen removal without independent spectral/sparse arbiters; (c) require a preinstalled Eigen package and break clean first builds. |
| D-26 | Waive CVODE as a mandatory v4 Phase-3 arbiter; use the strengthened analytic diagonal-subflow oracle. The test evaluates `exp(rate·h)·z₀ + expm1(rate·h)·forcing/rate` through an independent long-double branch/order and exact zero-rate limit. | For frozen coefficients and piecewise-constant flux the modal ODE is diagonal and has an exact closed form, with no reference-integration error. P3-G1 covers both electrodes and `nch={5,8,12}`; P1-G3 independently covers spectral compilation. **Limit:** the originally registered full 1C/current-step voltage comparison did not run. This waiver does not validate a generic RHS adapter, nonlinear full-cell integration, slow splitting, or SUNDIALS interoperability. | (a) add SUNDIALS solely to numerically approximate an analytically closed subflow; (b) claim that the former copied-φ₁ test was algorithmically independent; (c) imply the unrun full-trajectory CVODE gate passed. |
| D-27 | Retroactively waive the §5.7 quiet-host/Volkan-operator condition for the already completed v4.0 PAY-1, PAY-2, and PAY-4 checkpoints. Keep every number permanently labelled qualified development-host evidence. | Equal-work/correctness and structural gates passed; PAY-1/2 alternate order and publish conservative within-run ratios. The literal quiet/operator condition was not met, and busy-host interference can favour either side. Any unqualified speed claim requires a future named-host run with committed raw output, repeated timings for both tools, and load/thermal-stability criteria. | (a) block the completed architecture on a rerun that this session cannot make quiet or human-operated; (b) pretend the busy host satisfied the protocol; (c) assert that quiet conditions can only improve ratios. |

## 5. Migration strategy & verification discipline

1. **Legacy stays green.** New core grows in `src/core/` (`slide::core`); legacy `src/` untouched except Phase 0 bug
   fixes. CI runs both.
2. **Parity harness (the arbiter).** `tests/parity/` runs registered scenarios through legacy AND core and
   digit-diffs: single Kokam SPM 1C CC discharge; CCCV cycle; 3s2p ECM pack with contact R; SPM pack 2p. For
   digit-diff, core runs in **legacy-Euler mode** (same scheme, same dt). Registered band, written BEFORE the run:
   max |ΔV| ≤ 1e-12 V, max |Δstate| ≤ 1e-12 (rel, with an abs floor of 1e-15 for states crossing zero — rel is
   undefined at sign changes, and z-modes DO cross zero; floor ASSUMED 2026-07-09, tighten if a gate trips on it)
   over the full trajectory. Phase 9A additionally adopts the same 1e-12 digit-parity scale for P2-G1 branch
   current; §5.2 had not originally named a current band, so this is not mislabelled as a prior registration.
3. **Scheme upgrades validated separately** (never against loose-tolerance references — pouch-cell trap, logged
   twice). D-26 replaces the planned tolerance-converged CVODE run with an independently evaluated exact
   diagonal-subflow oracle at `nch={5,8,12}`. The originally registered `|V_expm − V_ref| < 0.1 mV` full
   1C/current-step voltage run did **not** occur and is not claimed; see the explicit scope limits in D-26.
4. **Oracles on non-degenerate cases:** parity scenarios include asymmetric electrodes, nonzero contact R,
   heterogeneous initial SOC — never only the symmetric/uniform case. For COMPOSED kernels (diffusion + thermal +
   ageing coupled) where no analytic series exists, the oracle is the **Method of Manufactured Solutions** (pick a
   solution, derive the source term that makes it exact — Roache; Salari & Knupp): it covers the coupling terms
   that P1-G3's single-physics series cannot, and is mathematics independent of both parity and CVODE references.
5. **Baselines recorded first.** Before each phase: run full ctest, record failing-test names + counts in §8; every
   commit re-runs; report deltas ("2 failing {a,b} → 3: +c, caused by me").
6. **No unqualified timing claims.** Machine runs concurrent jobs; performance evidence = allocation counts,
   iteration counts, complexity, and vectorisation reports. Existing wall-clock PAY results remain permanently
   labelled qualified development-host evidence under D-27.
7. **Payoff checkpoints — the refactor must prove itself before it is allowed to grow (Volkan, 2026-07-07).**
   Insurance first: the strangler strategy (item 1) means legacy stays green the whole time, so the worst case of an
   underperforming v4 core is deleting `src/core/` — nothing user-facing is ever bet on it. On top of that, each
   early phase was intended to end with a **wall-clock checkpoint run by Volkan on a quiet machine**. That literal
   condition was not met for PAY-1/2/4; D-27 waives it retroactively for v4.0 while permanently qualifying the
   results. Future unqualified claims still require the named quiet-host protocol. Registered thresholds remain:
   - **PAY-1 (Phase 1 exit) — QUALIFIED PASS 2026-07-10:** 10⁴ identical Kokam SPM cells, 1C CC discharge, 1 h simulated — v4 batch vs a loop of
     legacy `Cell_SPM`. Hypothesis [inferred, from devirtualisation + SoA + expm substep collapse; pouch-cell
     precedent 26–52× on the integrator alone]: ≥5×. **Abort threshold: <2× → STOP; profile, find where the model
     was wrong, redesign before any Phase-2 work.** No new phase on top of an unproven core. The reproducible Release
     harness alternates core/legacy order and validates equal final work. Three full 36·10⁶-cell-step repetitions:
     core 0.436–0.471 s (median 0.460), legacy 2.900–3.466 s (median 3.209), median speedup 6.97× and conservative
     `min(legacy)/max(core)` 6.16×; max mapped-state error 8.67e-19 and max |ΔV| 8.88e-16 V. Target and abort gate
     both clear even under the conservative range calculation.
   - **PAY-2 (Phase 2 exit) — QUALIFIED ABORT-GATE PASS 2026-07-10; 10× HYPOTHESIS FALSIFIED:** heterogeneous 16s4p
     pack, 1,800 s CC cycle — v4 compiled pack vs legacy `Module_s/Module_p`. Three alternating Release repetitions:
     compiled median 0.003089 s vs legacy 0.015049 s = 4.87×; conservative `min(legacy)/max(core)` = 4.83×.
     The ≥10× hypothesis did not hold, but the <3× abort threshold clears. The initial 0.69× measurement forced
     analytic SPM tangents, compiled periodic-brick proof/coalescing, compressed periodic rollback/cumulatives, and
     outer CC substepping before rerun. Final legacy-arbiter differences: state 2.72e-8, current 2.70e-6 A, voltage
     2.92e-8 V (legacy parallel solve itself uses a 1e-6 A residual). Timing remains development-host-qualified.
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
     **MEASURED 2026-07-10 (development machine; qualified, reproducible harness):** PyBaMM 26.6.2.0 IDAKLU
     (`rtol=1e-6`, `atol=1e-8`) warm solve median 12.373 ms vs SLIDE 0.682 ms = **18.13×**; cold
     setup-plus-first-solve ratio **71.4×**. liionpack 0.3.12/CasadiManager (its pinned PyBaMM 24.9.0) vs SLIDE:
     16s4p median **1,363×**, conservative comparator-solve/SLIDE-max **1,165×**; 1s100p **1,160×** median,
     **1,111×** conservative. Both targets clear. Common Chen2020/SOC/current/timestep and 100 µΩ cell resistance
     were used; liionpack's required 0.1 nΩ connectors, older pinned PyBaMM, and mandatory output recording are
     explicit qualifications. Final pack voltage difference is 0.90 mV/cell. See `benchmark/PAY4.md`.
   Structural proxies (allocations, virtual calls, iteration/factorisation counts, bytes/cell) are tracked
   continuously as leading indicators; wall clock at the checkpoints is the confirming evidence. Rule: phase N+1
   does not start before phase N's checkpoint passes or Volkan explicitly waives it (PAY-4 exempt from the
   phase-blocking rule — it is a positioning measurement).

## 6. Phased roadmap

Every phase ends: tests green (Debug AND Release), CHANGELOG updated, §8 updated, small commits pushed.
Gate simulations obey the header rule: SHORT registered scenarios; think first, simulate last.

### Phases 0–7 — COMPLETE (2026-07-07 → 2026-07-10; verbose gate text in the archive, outcomes in §8)

- **Phase 0** legacy stabilisation: A1–A7/B1–B5 + P0-C1..C6 fixed with regression tests; full ctest 10/10.
- **Phase 1** core data model + SPM kernels: P1-G0..G4 passed (trajectory parity 8.88e-16 V, zero per-step
  allocations at 10⁴ lanes, 240 B/cell, bitwise checkpoint-restart); **qualified PAY-1 6.97× median (target ≥5×)**.
- **Phase 2** pack layer (compile/netlist/D-21 thermal graph, Modes A/B, `SolverWorkspace`): P2-G1..G5 passed
  (3s2p parity 5.33e-15 V; Mode B vs A ≤3.1e-13 at 256p); **PAY-2 4.87× — 10× hypothesis FALSIFIED, ≥3× abort
  gate cleared** after periodic-brick + analytic-tangent redesign; timing is qualified development-host evidence.
- **Phase 3** exponential/modal integration + multirate + step control: P3-G1 ≤2e-12 vs closed form, P3-G2
  ≤1e-9 Ah cycle balance, P3-G3 stable at nch=12/dt=1000 s where Euler diverges.
- **Phase 4** Mode C waveform relaxation: P4-G1 ≤0.1% vs Mode A, P4-G2 drift within registered bound,
  **P4-G3/PAY-3 10⁵-cell pack in 24 MB, zero measured-step allocations/factorisations**.
- **Phase 5** Experiment grammar + Cycler v2 + event root-finding: P5-G1 passed (band provenance flagged — AUD-2).
- **Phase 6** Recorder + hardened I/O (CSV, CRC mmap, optional Parquet): P6-G1 passed, derived==live exactly.
- **Phase 7** Python wheel + PyBaMM compat + sensitivities + PyBOP entry: P7-G1 0.778/0.221 mV (band 15/8);
  P7-G2 dual-vs-FD worst 53 nV, PyBOP L-BFGS-B recovery ~1e-9 relative; P7-G3 C/50 0.134/0.034 mV, 1C
  11.577/2.028 mV; **PAY-4 18.13×/71.4× vs PyBaMM IDAKLU, ~1,165×/1,111× conservative vs liionpack (qualified,
  dev machine)**; installed-wheel CI (3 OS × Python 3.10/3.13) + pinned-fixture drift job.

### Phase 8 — MATLAB, CUDA, recording, runtime (COMPLETE; G6 moved to Phase 11)

DONE 2026-07-10: **P8-G0** optional-toolchain-free root and nested-superproject core builds now run a real
external-consumer smoke (1/1 on Windows) with CUDA/zstd/Arrow disabled-capability checks; public CPU headers contain
no optional SDK types, CUDA metadata is private, installed Eigen is preferred with a pinned fallback, and a
Linux/macOS/Windows core-only CI matrix plus complete CMake-trigger coverage is committed. Full optional-off CPU
suites pass Debug 49/49 and Release 49/49. Per D-25, "dependency-free" never meant removing required cold-path Eigen.
The rebuilt CPython 3.13 Windows wheel installs outside the source tree and passes 9 tests with 2 expected optional
skips (PyBOP and CUDA).
**P8-G1** MATLAB `+slide` over one stateless MEX (Tutorial-5 0.778/0.221 mV; 1.15e-14 V vs the
Python wheel; stable `slide:*` errors; 50 lifecycle repetitions). **P8-G2** CUDA fused SPM backend (10,003-lane
heterogeneous gate; byte-exact device rollback; PAY-5 23.85× vs CPU exact batch on RTX 4000 Ada) + truthful
device dispatch through Python. **P8-G3** async compressed recording, CPU ring + pinned CUDA side-stream, explicit
block/thin backpressure, bitwise round-trips. **P8-G4** persistent thread pool, fixed-order reductions
bit-repeatable across worker counts.

DONE 2026-07-10: **P8-G5 docs** adds v4 installation, C++/Python/MATLAB quickstarts, cell-model and ageing
extension guides, an optional-dependency matrix, and explicit PyBaMM-shaped API gaps while retaining clearly
labelled v3 pages. One standard-library extractor checks links/front matter and runs the exact fenced snippets:
C++ external consumer 1/1, installed-wheel Python, and licensed MATLAB R2025b all produce seven finite samples
(Python/MATLAB final voltage 3.879196 V). Doxygen 1.14 builds a non-blank v4 API landing page with zero generator
errors (123 pre-existing `src/` comment warnings remain); production Jekyll/GitHub Pages builds eight themed v4
pages with the pinned theme, `/SLIDE` base URL, and correct edit links. CI runs C++/Python fences and scopes Pages
write/OIDC permission to deployment; MATLAB remains a licensed-machine gate.

### Phase 9 — Adversarial bug-hunt + code-quality hardening (NEW 2026-07-10 — release gate-keeper)

Rationale: 13k lines written in 3 days pass every registered gate, but gates only test what was anticipated.
This phase hunts what was NOT anticipated, then simplifies without behaviour change. Unlimited thinking,
rationed simulating: every new test is a SHORT scenario (≤ a few hundred steps) or no simulation at all.

**9A — audit debt: COMPLETE 2026-07-10.**

1. **AUD-1 resolved:** P2-G1 voltage/current assertions enforce 1e-12. Debug/Release short gates pass; worst
   current drift is 3.70e-13 A (36.95% of the band), so the removed slack was material.
2. **AUD-2 resolved:** `09e8ec5` introduced implementation and thresholds together. The 0.2 µV/20 µA/0.2 µAh
   and 2e-12 V values remain useful regression sentinels but are explicitly post-hoc, not independent accuracy
   guarantees. The charge scale and 45-bisection event scale are derived in the Phase-9A report.
3. **AUD-4 resolved by D-26:** mandatory CVODE is waived for the exact frozen diagonal subflow. P3-G1 now uses
   an independent long-double `exp`/`expm1(rate·h)/rate` evaluation with no copied production Taylor branch.
   The unrun full 1C/current-step voltage comparison is explicitly not claimed.
4. **AUD-3/Q11 resolved by D-27:** the literal quiet-host/operator requirement is waived retroactively for v4.0;
   PAY-1/2/4 stay qualified development-host evidence. Raw PAY-4 JSON is committed; no weak benchmark rerun was
   performed on the same busy host.

Full derivations and validation evidence: `.claude/reports/p9a-audit-debt-2026-07-10.md`.

**9B — systematic bug-hunt.** Per-subsystem adversarial passes; every confirmed bug gets a registered SHORT
failing test BEFORE its fix; every refuted candidate is recorded refuted (killed ideas stay killed). Keep a bug
ledger in §8. Subsystems: StateArena/BatchView/row-roles; PackTopology+compile; PackSolver A/B/C + workspace
invalidation; steppers (EulerLegacy, ExponentialModal, step-doubling rollback); PackStepper transaction;
Experiment parser + CyclerV2 events; Recorder/AsyncRecorder/mmap reader; ParameterSet/BPX/expression AST;
ForwardSensitivity/Dual; CUDA mirror; Python/MATLAB boundaries.

Hunt checklist (mechanised, not vibes):
- lane/row indexing at degenerate counts {1, 2, 3, non-SIMD-multiple, coalescing-boundary};
- rollback completeness: every arena row + workspace + warm state restored on EVERY failure path (inject failures
  mid-transaction);
- ignored `[[nodiscard]]` Status returns (grep + compiler flag);
- fast-math folding hazards: audit every finite/NaN check against the `-Ofast` build (one instance already bitten
  and fixed — `std::isfinite` folded away; the IEEE exponent-bit idiom must be used everywhere);
- uninitialised scratch on first use after cold build and after restore;
- accumulator double-counting across substeps/rejected trials;
- event bisection corners: event true at segment start, two events in one step, event exactly at a breakpoint;
- degenerate packs: 1s1p, single lane, zero-R link, single-cell "pack";
- compiled-curve domain edges: query at/beyond first/last knot, one-bin tables (the BPX 1/65,536 class);
- snapshot cadence boundaries (first/last step, cadence > run length);
- CUDA tail lanes (non-multiple of block), device/host divergence under rejected steps.

Tooling (CI lanes; core itself stays dependency-free):
- ASan+UBSan lane on the full Debug suite; TSan lane on ThreadPool/AsyncRecorder/PackStepper tests;
- bounded fuzz drivers for the three parsers (Experiment grammar, BPX JSON, netlist CSV): malformed input must
  Status-fail atomically — never crash, never UB, never partial state; minutes of fuzzing per CI run, with a
  committed regression corpus;
- allocation-counter coverage extended to the PackStepper thermal path and recorder enqueue;
- error-branch coverage measured (llvm-cov) — every `return Status::…` failure branch in `src/core` exercised
  by at least one test.

**Gates:** **P9-G1** ASan/UBSan/TSan lanes green. **P9-G2** parser fuzz: bounded campaign, zero crash/UB/leak,
every rejection atomic; corpus committed. **P9-G3** error-branch coverage 100% on `src/core` Status-failure
branches (measured, exceptions listed and justified). **P9-G4** bug ledger complete in §8: each bug has
{failing test first, fix, ledger row}; refutations recorded.

**9C — advanced simplification (the Linus test: good taste removes special cases; guards are a smell).**
Rule for EVERY change: no-op verified — digit-identical outputs on recorded cases (user CLAUDE.md §4: no-op
claims need digit-identical outputs, not code inspection), full Debug+Release suites green, PC-1..8 preserved,
no new required deps. Targets (from §2.3):
1. **One physics source:** extract the modal update + shared observable scalar kernels into a single
   scalar-generic header consumed by `SpmPipeline` (CPU), `CudaSpmRuntime.cu` (`__host__ __device__`), and the
   `Dual` instantiation. Gate: all three paths digit-identical to their pre-refactor outputs on a recorded case.
2. **One ageing-kernel idiom:** unify the SEI/CS/LAM/plating mask/scratch/lane-sweep scaffolding; four physics
   bodies remain, one pattern.
3. **Split oversized cold files:** `ParameterSet.cpp` → absorption / BPX reader / expression AST TUs;
   `Experiment.cpp` → parser vs runner.
4. **Shared test harness:** factor the repeated build-batch/run/compare scaffolding in `tests/unit/core_*` into
   one helper; suite counts must not drop.
5. **Public-surface audit:** classify each `src/core` header API vs implementation detail; move detail out of
   public reach; naming-consistency pass (one verb per concept across files).
6. **Dead-code sweep** (`SpectralDiffusionLegacy` is deliberate — §5.2 parity kernel, KEEP).
**Gate P9-G5:** every simplification lands digit-identical; core line count recorded before/after in §8 (expect
net reduction; growth requires written justification).

### Phase 10 — PyBOP integration test suite (NEW; Volkan 2026-07-10)

All fits use SHORT synthetic data (≤600 s simulated, dt ≥ 10 s, ≤4 lanes; reuse P7-G2 problem shapes). Derive
before running: identifiability and noise floors are computed analytically first, and bands registered from the
derivation, not from a trial fit.

- **10A** gradient regression (exists — P7-G2 L-BFGS-B via `simulateS1`): keep as the pinned smoke test.
- **10B** gradient-FREE optimiser through PyBOP (CMA-ES or XNES): same 2-parameter synthetic recovery with a
  registered looser band — proves the surface works without sensitivities.
- **10C** multi-parameter GITT-style fit: ≥4 parameters {D_s,n, D_s,p, R_contact, Q} on a synthetic pulse train
  (~10 pulses × 60 s); recovery bands registered from a pre-run identifiability analysis (normalised-sensitivity
  Gram matrix rank + conditioning — pure derivation, no simulation).
- **10D** noise robustness: registered Gaussian noise (e.g. σ = 1 mV) on the synthetic trace; recovery bands
  derived from the CRLB estimate BEFORE the one decisive run.
- **10E** version drift: pinned PyBOP 25.11 job blocking; a latest-PyBOP job non-blocking that flags API/numeric
  drift (mirror of the PyBaMM fixture-drift pattern).
- **10F** sensitivity CI regression: `simulateS1` vs central-FD arbiter (short horizon, 3 parameters) as an
  installed-wheel pytest — catches silent dual-path regressions.

**Gates:** **P10-G1** 10B/10C/10D recover truth within pre-registered bands in the installed-wheel CI matrix.
**P10-G2** the identifiability + CRLB derivations are committed docs (the trail another scientist can trust).

### Phase 11 — v4.0.0 release (was P8-G6; runs AFTER Phases 9–10)

Registered gate text unchanged (archive): CMake/package/MATLAB versions agree on `4.0.0`; consolidated CHANGELOG
section with migration/known-limit notes; dependency licenses; reproducible artifacts (CPU wheel, CUDA wheel or
build instructions, MEX package); final Debug/Release CPU, installed-wheel Python, MATLAB batch, CUDA, async
recording, and release-consistency checks green before the annotated `v4.0.0` tag. Legacy v3 façades stay
compiled/runnable (Q3). ADDITION: release notes state Phase-9 outcomes (bugs found/fixed/refuted counts) and the
qualified-timing caveats (AUD-3/Q11).

### Beyond v4.0 (recorded so design choices don't foreclose them)

- WASM GUI (Emscripten build of the dep-free core); thermal 1D/2D per-cell; blended electrodes
  (`vector<ActiveMaterial>`); SYCL/HIP (Q5); f32 storage instantiation (Q1); FMU export (keep the §3.12 triple
  FMI-congruent); ADJOINT sensitivities (many-parameter fits).
- NEW candidates (2026-07-10): structural-counter regression CI (allocations/factorisations/iterations asserted
  per commit — no timing); MSVC CI lane (portability goal 7; known ~3× slower codegen is a measurement, not a
  blocker); clang-tidy + include-what-you-use profile; SPMe/DFN discretisation slots (§3.2 menu);
  liionpack-netlist CSV import round-trip example in docs.

## 7. Open questions for Volkan (OPEN/ASSUMED ledger)

Q1–Q10 are DECIDED/RESOLVED — one-line records below; full reasoning in the archive and §4.

| ID | Question | Outcome |
|----|----------|---------|
| Q1 | f32 state? | f64 everywhere in v4.0; `real_t` alias keeps f32 open (2026-07-07) |
| Q2 | Keep `Cell_ECM<N_RC>`? | Yes — entailed by P2-G5a/parity/Tier-0 gates (2026-07-07) |
| Q3 | Legacy API fate? | Façade through v4.x, delete in v5 (2026-07-07) |
| Q4 | KLU optional dep? | Yes, optional; Eigen SparseLU default, silent degrade (2026-07-07) |
| Q5 | CUDA-only first? | Yes; SYCL/HIP later; no platform API without a portable seam (2026-07-07) |
| Q6 | "Ross's analytical solution"? | arXiv:2508.14454; adopted as Mode-B upgrade D-20 (2026-07-07) |
| Q7 | PyBaMM target version? | Latest stable at Phase-7 start = 26.6.2.0, keyed + aliased (2026-07-07) |
| Q8 | Parity band vs FP reassociation? | Band kept 1e-12; op-order-pinned parity kernel; Release drift 5.63e-15 (2026-07-08) |
| Q9 | Pack thermal coupling? | D-21 compiled thermal graph, designed before Phase 2, implemented (2026-07-10) |
| Q10 | First-class sensitivity parameters? | Ten-parameter fitting set; `h_conv` deferred to thermal composition (2026-07-10) |
| **Q11** | **AUD-3: PAY-1/2/4 quiet-machine confirmation runs — rerun on a quiet machine, or waive the §5.7 protocol retroactively?** | **RESOLVED by D-27:** waive retroactive quiet-host confirmation for v4.0; retain PAY-1/2/4 numbers as qualified development-host evidence. Busy-host interference can favour either side, so no directional claim is made. |
| **Q12** | **Release ordering: v4.0.0 after Phase 9 only, or after 9 AND 10?** | **ASSUMED after both** (Phase-11 ordering) — PyBOP integration is cheap relative to a post-release API fix; overturn by one word if speed matters more |

## 8. Status ledger (compact; the 60-row per-gate history is archived verbatim)

| Date | Item | State |
|------|------|-------|
| 2026-07-07 | Phase 0 + legacy audit + plan v1 + SOTA verification + Q1–Q7 | DONE (archive) |
| 2026-07-08 | P1-G0 Release re-confirm; production `SpectralDiffusion<NCH>` | DONE (archive) |
| 2026-07-09 | §3.12/D-22/D-23 kernel contract; orthogonal review (13 defects pre-implementation); Chebyshev math audit; PyBOP promotion; async-recording design | DONE (archive) |
| 2026-07-10 | Phases 1–7 all gates + PAY-1 (6.97×), PAY-2 (4.87×, 10× falsified), PAY-3 (10⁵ cells/24 MB), PAY-4 (18×/71× PyBaMM, ~1.1–1.2·10³× liionpack, qualified) | COMPLETE (archive rows 2026-07-10) |
| 2026-07-10 | Phase 8 G1 MATLAB, G2 CUDA (+PAY-5 23.85×), G3 async recording, G4 thread pool | COMPLETE (archive) |
| 2026-07-10 | Implementation audit of the 67 Codex commits (agent + Fable, artifact-checked); architecture quality assessment; full ctest Debug 49/49 + Release 49/49 | DONE — verdict §2.2/§2.3; AUD-1..5 opened as Phase-9A items |
| 2026-07-10 | PLAN.md compressed + Phases 9/10/11 added (bug-hunt+simplification, PyBOP integration, release); short-simulation operating rule added; P8-G6 → Phase 11 | DONE (this revision; archive created) |
| 2026-07-10 | P8-G0 optionality/portability | PASSED — core-only and nested external-consumer smoke 1/1 on Windows; Debug/Release optional-off CPU 49/49; rebuilt installed CPython 3.13 wheel 9 passed/2 expected skips; private CUDA metadata; no optional SDK header leaks; installed-Eigen-first/pinned fallback; 3-OS core and installed-wheel CI matrices cover CMake changes. Cross-platform jobs are committed but cannot run until pushed. |
| 2026-07-10 | P8-G5 tested v4 documentation | PASSED — exact C++/Python/MATLAB fences run in available toolchains (7 finite samples each); 20 local links/front matter checked; Doxygen 0 generator errors with a rendered main page; production Jekyll build emits 8 themed v4 pages; docs workflow linted and least-privilege. AUD-5 resolved. |
| 2026-07-10 | Phase 9A audit-debt closure | PASSED — AUD-1 enforces 1e-12 P2-G1 V/I bands (worst ΔI 3.70e-13 A); AUD-2 classifies same-commit Phase-5 thresholds as post-hoc regression envelopes with derivation; D-26 waives CVODE narrowly after strengthening the exact oracle; D-27 resolves Q11 while preserving qualified timing labels. Targeted Debug/Release 3/3, final full Debug/Release 49/49; raw PAY-4 JSON preserved. |
| — | NEXT | Phase 9B systematic bug-hunt → 9C simplification → Phase 10 → Phase 11 release |
