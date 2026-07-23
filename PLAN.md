# SLIDE — THE GOAL: the superior battery-simulation stack (living contract)

> **Status:** ACTIVE. Rewritten 2026-07-11 (Fable) around a SINGLE GOAL with one milestone ladder (§6);
> extended 2026-07-14 (Fable) with the breadth wave (goals 15–24, M13–M20) from Volkan's 16-point directive
> and the same-day Newman-instrumentation directive; extended 2026-07-23 (Fable) with the **MQ quality
> wave** (clean the repository / re-derive the mathematics / hunt logical + performance mistakes —
> inserted between M0 and M1 in ladder order) and the split-out standing operating contract
> **`AGENTS.md`** for autonomous Codex sessions.
> The v4 core is COMPLETE through Phase 9A; Phase 9B was in progress at rewrite time and is carried into M0.
> Everything before this rewrite is archived VERBATIM at
> `.claude/summaries/plan-archive-2026-07-11-v4-phase9b.md` (and the older
> `.claude/summaries/plan-archive-2026-07-10-phases0-8.md`) and in git history.
> **How to use this document:** single source of truth. (`AGENTS.md` at the repo root compresses the
> operating rules for agents that load it automatically; it points here and never overrides this file.)
> Any session (Codex, Fable, Opus, human) must
> (1) read this file first, (2) take the FIRST unticked box in §6 whose dependencies are ticked,
> (3) tick it, update §8, commit, and take the next box — the non-stop protocol is §0.3. Do not re-litigate
> §4 decisions without new evidence overturning them.
> **Operating rules:** small reviewable commits on branch `Claude`; every user-visible change updates
> `CHANGELOG.md` (Unreleased); wall-clock numbers on this machine are UNRELIABLE (parallel jobs) — judge
> performance by structural counters (allocations, iterations, factorisations, bytes/cell); **simulations
> cost time — every test/gate simulation stays SHORT** (≤ a few hundred steps, seconds of wall-clock, small
> lane counts unless the gate is specifically about scale). Unlimited thinking, rationed simulating.
> Python deps via `uv` only; never write outside the repo root; measured data are read-only.

## 0. Who you are — the engineer this plan assumes

You are a principal-level scientific programmer. Your C++ is Stroustrup-grade: value semantics, concepts,
ranges, zero-overhead abstraction as a discipline, not a slogan. Your GPU numerics are Giles-grade: you
think in memory traffic and occupancy before you think in FLOPs. You know the battery-modelling literature
(the Newman school, PyBaMM's model zoo, spectral and finite-volume discretisations) and you check it —
via the web or the reference repos named in this plan — before designing from memory. And you carry the
lazy-genius creed with pride:

- **Laziness is a virtue when it deletes code.** Never write the same fact twice (PC-10). The best special
  case is no special case. A net-negative diff that preserves behaviour is a triumph, and boredom is a
  signal to simplify, never to stop.
- **Think first — your scarce resource is wall-clock, not tokens.** Builds and simulations cost minutes;
  derivation costs nothing. Derive the expected answer, register the pass/fail band, THEN run the short
  decisive test. Prefer analytic oracles, conservation identities, and structural counters over long
  trajectories. When a milestone names reference repositories with concrete paths, go read them — hours
  of reading beat days of rediscovery.
- **Quality bar (non-negotiable, from the standing CLAUDE.md discipline):** every claim traceable to an
  artifact; outputs quoted verbatim; every equation unit-checked; every approximation names its regime and
  leading-order error; killed ideas stay killed (§2, §4); "no-op" claims require digit-identical outputs on
  recorded cases, never code inspection; a falsified registered hypothesis is a deliverable, recorded
  prominently.

### 0.3 The non-stop protocol (one goal, many boxes, no stopping)

1. **Always take the FIRST unticked box in §6 whose dependencies are ticked.** Never idle; never pause at a
   milestone boundary — ticking the last box of M(k) flows directly into the first box of M(k+1).
2. **A box is ticked only when:** its gates are green (Debug AND Release where applicable), CHANGELOG is
   updated, a §8 row is added, small commits are made, and the checkbox in this file is edited to `[x]`.
3. **Gate failure:** make three genuinely different attempts (different hypotheses, not retries of one idea).
   Still red → record FALSIFIED/BLOCKED in §8 with the numbers, leave the box unticked, and move to the next
   unblocked box. Return later with new evidence only.
4. **Design decisions inside the invariants are YOURS.** Before each milestone marked "design note first",
   write a short note in `.claude/designs/` (problem, options, choice, why, invariants touched), choose the
   strongest option yourself, add a §4 D-entry marked ASSUMED, and continue. This plan states WHAT and the
   invariants; the HOW within a milestone is yours — a better design that preserves the invariants is
   welcome and expected. You are trusted to think on the spur of the moment; record what you decide.
5. **Only three things wait for Volkan:** overturning a PC/EC/MC-contract row (needs a §4 entry he approves),
   irreversible or outward-facing actions (annotated release tags, publishing artifacts, deleting data),
   and writing outside this repository. Everything else: decide, record, proceed. Boxes waiting on Volkan
   are EXCLUDED from the rule-1 selection: leave them unticked with a §8 "awaiting go" row, continue past.
6. **Session end = context exhaustion only.** Then: update box states and §8, commit, and leave a ≤10-line
   handoff in `.claude/summaries/` pointing at the next box.

## 1. THE GOAL

**Ship the superior battery-simulation stack: PyBaMM's model breadth and validation culture at
compiled-SoA speed; PyTorch's one-word device choice; COMSOL's parameter/material/discretisation
flexibility; a PyBaMM drop-in Python subset; Python, MATLAB, and browser reach; SimSES's system-level
scope — SPM/SPMe/DFN/ECM, lead-acid Schiffer, and semi-empirical tiers on one zero-overhead
architecture, every tier cross-validated against PyBaMM and fit-tested through PyBOP; enjoyable for the
developer (std-library composability, helper one-liners, clean seams) and for the user (Studio, docs
site, honest advisors).** Three release checkpoints mark the road: `v4.0.0` (M2), `v5.0.0` (M12 — the
original v5 scope, goals 9–14), `v6.0.0` (M20 — the 2026-07-14 breadth wave, goals 15–24).

Definition of done = every MANDATORY §6 box through M20 ticked or recorded FALSIFIED/waived in §8
(optional boxes skipped only with recorded justification), and the `v6.0.0` tag gates green. The v4 goals (scale
10⁴–10⁵ cells, zero-overhead composition, electrode-first modularity, flat compiled pack solving, CasADi-style
interface symmetry, memory-aware recording, portability, correctness discipline) remain in force — they are
the floor, not the target. v5 adds:

9. **Model breadth:** SPMe and DFN tiers beside SPM; FVM discretisation slot beside Chebyshev-spectral;
   every tier batch-compiled, SoA, sensitivity-capable — never a slow "reference-only" path.
10. **Thermal fidelity:** 2+1D/3D spectral pouch thermal (Robin × cosine tensor bases, exact modal
    propagation) ported from `C:\D\git\pouch-cell-spectral`.
11. **Model hierarchy with honesty:** a regime advisor computes the dimensionless groups (Biot, diffusion
    time ratios, conductivity ratios) at build and tells the user which simplification is valid and what its
    leading-order error is; assumption violations are reported, never silent.
12. **Reach:** WASM build of the dependency-free core + SLIDE Studio, a browser design-simulate-analyse
    workspace (visual pack building, project save/load, plotting, citations panel).
13. **Scholarship built in:** every model, discretisation, solver algorithm, and parameter set carries its
    citations; `print_citations` reproduces the exact bibliography of what a run actually used
    (PyBaMM-precedent, done at compile-time cost only).
14. **Expressiveness:** packs and batches compose with the standard library — row spans, lazy ranges,
    named helpers (`total_current`, `max_temperature`, `weakest_cell`) that make developers enjoy building
    on SLIDE.

**The breadth wave (Volkan directives 2026-07-14; lands in M13–M20):**

15. **Model-family breadth beyond Li-ion porous-electrode:** a first-class ECM batch tier in the core
    (legacy `Cell_ECM` parity), a semi-empirical degradation tier (SimSES-class stress-factor calendar +
    cycle models), the lead-acid Schiffer weighted-Ah-throughput model (Q16 OVERTURNED for lead-acid —
    D-43), and a storage-system tier (power conversion, EMS, application profiles) — §3.23, D-48.
16. **PyBaMM drop-in subset:** `import slide.pybamm as pybamm` runs a defined, honestly documented subset
    of PyBaMM scripts unchanged — the JAX-for-NumPy positioning (§3.24, D-44).
17. **Device as UX:** `device="cuda"` is a one-word choice on every GPU-instantiated tier (recorded
    list, minimum SPM + ECM — M16.1), PyTorch-style; every other tier reports the capability gap
    EXPLICITLY, never silently; CPU stays the reference; per-tier CPU↔GPU digit gates (D-45; §3.8
    promoted to a user contract).
18. **Cross-validation as a standing contract:** every tier with a PyBaMM counterpart ships registered
    parity fixtures; every fit-relevant surface is exercised through PyBOP (§1.4 VC table).
19. **Gradients everywhere they pay:** forward Dual sensitivities on all tiers (PC-10 makes them
    near-free); Jacobians served to PyBOP; adjoint/Enzyme evaluated once with a recorded decision (M16.5).
20. **Toolchain breadth:** Clang + GCC + MSVC across Linux/Windows/macOS CI, all mandatory (M19.1).
21. **A real documentation website:** task-first tutorials per interface, theory pages linking the
    derivation docs, Doxygen API reference, benchmark page under D-27 discipline (M19.2).
22. **Speed research with registered hypotheses:** specialised integrators exploiting THIS architecture
    (exponential/Rosenbrock, splitting theory, projective/cycle-extrapolation predictors,
    Anderson-accelerated relaxation) admitted only through derive-first gates; a falsified candidate is a
    deliverable (§3.25, D-49; M17).
23. **Hardening as ladder work, not culture:** recurring verification/bug-hunt/performance/simplification
    passes are explicit boxes (D-46; M12.0, M20.1).
24. **Newman instrumentation:** thermodynamic identities as runtime gates (entropy production, full
    Bernardi energy balance incl. heat of mixing), polarization accounting that closes to roundoff,
    virtual reference electrode + plating-onset events, EIS from analytic linearisation of the compiled
    system, gradient-based cell design with exact analytical Jacobians, ICA/DVA diagnostics
    (§3.26, D-50..D-52; M18).

**Directive concerns that are ALREADY contract rows — do not re-derive, point here:** 10⁵-cell storage =
D-01 SoA arena + PC-3/PC-6 budgets; "don't compute voltage when not needed" = D-10 (snapshots stored,
observables derived lazily on demand); "cache variables/matrices/Jacobians" = D-18 chord/Shamanskii
`SolverWorkspace` + D-16 compiled-curve canonicalisation; "PI-control-like faster simulation" = the
error-controlled outer step (§3.5, D-08 — the literal proportional–integral step controller is M17.2)
and Mode C relaxation advance (§3.4); "helpers for standard algorithms / elegant code" = §3.13,
EC-2/EC-3; "small libraries that do the job without bloat (Highway, Thrust/CUB)" = D-47 measured-in
optional-dep policy; multi-rate integrators = D-08 (in force since v4); "different ageing formulas" =
the §3.2 mechanism slots + per-batch masks + the goal-15 semi-empirical tier; "easy battery
construction" = §3.4.2 value-combinators + Studio's pack designer (§3.22).

**Plan flexibility (restating §0.3(4), because it governs the new wave too):** boxes state WHAT and the
gates; the HOW belongs to the implementing agent. An intelligent agent is expected to IMPROVE this
ladder — split boxes, add gates, redesign within the invariants, or overturn a decision with evidence —
by writing the design note, adding the §4 D-entry (ASSUMED where §0.3(4) allows), and recording the
change in §8. The fixed points are the contracts (§1.1–§1.4), the evidence discipline (§5), and the
non-stop protocol (§0.3) — never the wording of a box.

### 1.1 Performance contract (enforceable invariants — every milestone gate re-checks these)

The guiding split: **the user-facing API may be arbitrarily pleasant because it runs once at setup; after
`compile()`/`build()`, only spans, indices, and precompiled kernels exist.** Strings, maps, options dicts,
nested constructors, parameter lookup — all cold path, all gone before the first step.

| # | Invariant | Enforced by |
|---|-----------|-------------|
| PC-1 | ZERO heap allocations per accepted simulation step | counting-allocator asserts |
| PC-2 | ≤1 indirect call per BATCH per step; zero virtual dispatch per cell | registry design (D-02) |
| PC-3 | All state contiguous, 64-byte aligned, SoA (SIMD-sweepable) | StateArena is the only state owner (D-01) |
| PC-4 | No exceptions, no locks, no I/O on the hot path | review + `Status` returns (D-09) |
| PC-5 | Inner loops see concrete types (templates), never function pointers/std::function | composition design (§3.2) |
| PC-6 | Memory budget per model tier: SPM ≤ ~300 B/cell (proven 240 B); SPMe ≤ ~1 KB; DFN ≤ ~16 KB (Q15); batch-shared params stored once | per-tier gate |
| PC-7 | Rollback/checkpoint = memcpy, never tree traversal | arena snapshots (§3.1) |
| PC-8 | Pack solve cost: O(n) modes B/C, O(nnz) sparse mode A — never O(n²) dense, never nested iteration | compile() flattening (D-03) |
| PC-9 | WASM-clean core: buildable single-threaded with no filesystem, no exceptions, no thread requirement; missing capabilities degrade EXPLICITLY (Status/report), never silently | M10 build + capability tests |
| PC-10 | ONE physics source: every formula exists in exactly one scalar-generic definition, instantiated for CPU, CUDA, Dual (sensitivities), and WASM | M0.5 header + review rule |

Any proposed change that violates a PC-invariant needs a §4 decision-log entry overturning it first.

### 1.2 Expressiveness contract (what "best developer tools" means, testably)

| # | Invariant | Enforced by |
|---|-----------|-------------|
| EC-1 | User-facing views/ranges over packs never own memory and never allocate per element | allocation-counter test on range iteration |
| EC-2 | Every per-cell quantity (I, V, T, SOC, any state row) is reachable as a row span or lazy range consumable by `std::ranges` / `std::` algorithms | compile-time `static_assert`s + example file in CI |
| EC-3 | Common questions are named one-liners built on EC-2 (`total_current(pack)`, `max_temperature(pack)`, `weakest_cell(pack)`, `soc_spread(pack)`) | unit tests + docs snippets |
| EC-4 | Every public API lands with a compiling, running docs snippet (the P8-G5 fence extractor runs them in CI) | docs CI |
| EC-5 | The description layer (CellDesign, PackDesign, Experiment, ParameterSet) serialises to versioned JSON with digit-exact round-trip — one schema shared by C++, Python, MATLAB, and Studio | C++ schema v1 built in M4.8; round-trip gates M4.8 + M11.8 |

### 1.3 Maintainability contract (nicely packed, easy to test, understand, maintain)

| # | Invariant | Enforced by |
|---|-----------|-------------|
| MC-1 | One concept per file; a file crossing ~700 lines is a review smell requiring a recorded justification in §8 | line-count check in review |
| MC-2 | Feature cohesion: model/kernel X ships as {header(s), `tests/unit/core_X_test.cpp`, docs snippet, citation entry} with matching names — a reviewer finds all four with one glob | review rule + MC-2 checklist per milestone |
| MC-3 | Every `src/core` header opens with a 3–6 line contract comment: what it owns, which PLAN.md § it implements, hot/cold classification | review rule (pattern already established) |
| MC-4 | Tests mirror source names; the shared test harness (M0.8) is the ONLY test scaffolding idiom | test-layout review |
| MC-5 | Public surface minimal and intentional; implementation detail lives in `detail::` or private headers | M0.9 surface audit, kept current |

### 1.4 Validation contract (cross-tool honesty, testably)

| # | Invariant | Enforced by |
|---|-----------|-------------|
| VC-1 | Every model tier with a PyBaMM counterpart (SPM, SPMe, DFN, ECM/Thevenin — counterpart names verified at each box, never from memory) ships committed parity fixtures against the pinned PyBaMM, bands registered BEFORE comparison (P7-G3 pattern); tiers with no counterpart record that fact, never silently skip | per-tier parity boxes (M5.6, M6.6, M7.6, M13.2) + drift jobs |
| VC-2 | Every fit-relevant surface (fitting parameters, `simulateS1` sensitivities) is exercised through PyBOP in installed-wheel CI | M1 suite + per-tier extensions (M13.5, M13.6, M16.3; Schiffer fit-relevance decided and recorded at M14.2) |
| VC-3 | Every tier that claims GPU support digit-matches its CPU instantiation on registered fixtures | existing CUDA gates + M16.2 |
| VC-4 | The `slide.pybamm` shim's supported subset is defined by a committed coverage matrix and proven by running pinned PyBaMM's own examples under it; everything else raises with a pointer, never silently diverges | G-SHIM (M15.3) |
| VC-5 | Cross-tool numbers (PyBaMM, liionpack, SimSES) are qualified per D-27 and never claimed from memory | §5.6/§5.7 discipline |

## 2. Where we stand (v4 evidence, compressed 2026-07-11; full audit text in the archives)

> NOTE: citations of the form "§2.1/§2.2/§2.4 A5/§2.5/§5.x/§8 rows" inside the kept-verbatim §3–§5 refer to
> the ORIGINAL v4 evidence-base numbering, preserved in the archive files — not to the subsections below.
> Likewise all `P*-G*` gate IDs and "Phase N" mentions inside §3–§5 are v4 history: every cited gate PASSED
> (unless stated otherwise) and every phase is complete; details live in the archives.

### 2.1 v4 state — COMPLETE through Phase 9A [confirmed by full ctest + artifact-checked audits]

- Phases 0–8 and 9A done in 82 commits (2026-07-07 → 2026-07-10); full ctest Debug 49/49, Release 49/49.
- Payoff record (ALL wall-clock numbers qualified development-host evidence under D-27):
  PAY-1 batch-vs-legacy 6.97×; PAY-2 pack 4.87× (**10× hypothesis FALSIFIED**, ≥3× abort gate cleared);
  PAY-3 10⁵-cell pack in 24 MB, zero measured-step allocations; PAY-4 18.13× warm / 71.4× cold vs PyBaMM
  26.6.2.0 IDAKLU and ~1,165×/1,111× conservative vs liionpack (`benchmark/PAY4.md`); PAY-5 CUDA 23.85×
  vs CPU exact batch (RTX 4000 Ada).
- The 67-commit implementation audit (2026-07-10) found the plan followed, no fabrication or gate-rigging;
  its five deviations AUD-1..5 were all closed by Phase 9A (D-25..D-27 record the waivers honestly).
- Architecture quality [confirmed by direct inspection 2026-07-10]: 13,151 lines / 48 files; zero `throw`,
  zero TODO/FIXME markers in `src/core`; legacy coupling = one `Status.hpp` include; three cold-path
  `std::function` only. Named debt became the M0 simplification boxes (modal-update triplication,
  `ParameterSet.cpp` 1,293 lines, `Experiment.cpp` parser+runner mix, four ageing kernels sharing one
  unfactored idiom).
- Phase 9B was IN PROGRESS at rewrite time: bug ledger P9-B01..B36
  (`.claude/reports/p9b-bug-ledger-2026-07-10.md`), P9-G2 parser fuzzing PASSED (D-29), P9-G1/G3/G4 open —
  all carried into M0. Working tree held uncommitted WIP in `src/core/PackSolver.cpp` +
  `tests/unit/core_ParserAllocation_test.cpp` (M0.1 lands or reverts it, recorded either way).

### 2.2 Standing corrections that must NOT be re-litigated (killed ideas stay killed)

- **A6 REFUTED — not a bug.** Legacy `Module_p::V()` consistency was derived and proved by test.
- **"Only time/Ah/Wh are path-dependent" FALSIFIED:** legacy hid state outside `State_SPM`. v4 RULE:
  anything a later step reads MUST be an arena row (enforced by bitwise-restart gate P1-G4).
- **Chebyshev nch≠5 defects** (centre-node sign, Eigen in-place `.inverse()` aliasing) fixed and validated
  by the tan μ = μ eigenvalue oracle + Carslaw–Jaeger transient; nch=5 fundamental eigenvalue ~2e-5 accurate.
- Legacy solver statics were deliberate quasi-Newton memory — kept properly as `SolverWorkspace` (D-18).
- PAY-2's 10× and the initial 0.43×/0.69× results, the Q8 bit-identity attempt, and the pouch project's
  retracted +1.47 K "ceiling" are recorded falsifications — do not re-open without new evidence.

### 2.3 Reference repositories on this machine (read them; do not re-derive what they already learned)

**`C:\D\git\pouch-cell-spectral`** — Volkan's MATLAB research code for coupled electro-thermal spectral
pouch simulation, benchmarked against COMSOL; the source for M8. Key entry points (scout-verified 2026-07-11;
digest at `.claude/references/pouch-cell-spectral-survey-2026-07-11.md`):
- `whitepaper/methodology.tex` — canonical method write-up; `src/+solver/pouch_cell_phase1.m` — core solver.
- `src/+spectral/robin_eigenvalues.m`, `robin_eigenfunctions.m`, `robin_normalization.m`,
  `robin_tab_projection.m`, `heat2d_eigenfunction.m` — Robin (convective-BC) Sturm–Liouville bases.
- `src/+spectral/cosine_modes.m` — Neumann cosine modes + exact tab projection; `cheby_multidomain.m` —
  4-domain through-thickness collocation; `solid_diffusion_operator.m` — radial Chebyshev.
- `docs/spectral.md`, `docs/3d_thermal_design.md`, `docs/transmission_line_erf.md`,
  `docs/spherical_diffusion_erf.md`, `mittag/theory/tab_heat_erfc.md` — erf/erfc short-time machinery.
- `mittag/` — method sandbox; its conclusion: AAA rational approximation supersedes erfc-image/residue/Padé
  for production kernel inversion (`mittag/matlab_1d/benchmark.m` compares six methods).
- `docs/howie_modal_temperature_coupling.md`, `docs/kappa_per_mode_derivation.md` — the OPEN off-diagonal
  cosine-coupling problem; **quarantined** (see D-36): the Howie electro-thermal closure was reopened
  (EX780–802; +0.875 K bias, 57% transfer-NRMS unresolved) — port the validated thermal/spectral machinery,
  NOT the contested closure.
- Validated headline (whitepaper, canonical COMSOL params): 6/6 cases |T_err| ≤ 0.55 K, V_RMSE 6.9–11.4 mV
  at 4C–8C; accepted bar "COMSOL-grade" = 22.2 mV / 0.128 K (COMSOL-vs-data residual, EX662).

**`C:\D\git\unibatt`** — Volkan's Rust/WASM battery diagnostics app (GLiDE); the pattern source for M10/M11.
Digest at `.claude/references/unibatt-survey-2026-07-11.md`. Key entry points:
- `glide-wasm/web/src/styles/main.css` — the Studio colour tokens (§3.22 quotes them).
- `glide-wasm/web/src/components/chart.js` — Wong colourblind-safe data palette + themed plot layouts.
- `scripts/build-single-html.js` — single self-contained offline HTML build (esbuild + base64 wasm),
  the shareable-artifact pattern M11.10 imitates; `glide-wasm/web/vite.config.js` — Vite + wasm + ES
  workers config and the COOP/COEP-vs-CDN caveat (workers cancelled via `terminate()`).
- `glide-wasm/CLAUDE.md`, `glide-wasm/CHECKLIST.md`, and
  `glide-wasm/docs/superpowers/specs/2026-03-31-glide-wasm-migration-design.md` — living invariants + design.
- Architecture facts to reuse: vanilla ES modules (no framework), worker-offloaded compute with progress,
  dark/light theme via CSS variables + localStorage, wasm loaded with graceful mock fallback.
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
  runtime flag, NEVER the silent default, and admit it only through gate **P2-G5** (v4 gate, PASSED — archive) with Mode A as arbiter.
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
- `device={"cpu","cuda"}` selected per Simulation (PyTorch-style, template = dtw-cpp `_api.py` lazy-dispatch design).
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
speed; off-the-shelf ODE integrators/methods usable; intuitive, sparse structure. This section resolved the
v4 Phase-1 "observable layer + BatchView/StepCtx" blocker (archive §6) under that principle. **REVISED 2026-07-09 after orthogonal
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

### 3.13 Expressiveness layer — packs compose with the standard library (M3; EC-1..4)

The SoA arena stays the ONLY owner of state (D-01); expressiveness is a zero-cost VIEW layer on top,
never a return to AoS objects (D-32):

```cpp
// Row spans — direct, allocation-free, std-algorithm-ready:
std::span<const double> I = pack.currents();          // one SoA row (or a stitched per-batch range)
auto  T_max  = std::ranges::max(pack.temperatures());
auto  I_tot  = std::reduce(I.begin(), I.end());        // parallel-pack terminal current

// Lazy per-cell proxies for cold-path convenience (iteration = index arithmetic, no ownership):
for (CellRef c : pack.cells())                         // random-access range of {batch*, lane}
    if (c.V() < 3.0) report(c.path());                 // "s03.p2"-style addressing kept

// Named one-liners built on the ranges (each is 1–3 lines, documented, tested):
slide::total_current(pack);  slide::max_temperature(pack);  slide::weakest_cell(pack); // argmin V
slide::soc_spread(pack);     slide::imbalance(pack);        // max-min / σ of per-cell SOC
```

- Multi-batch packs expose stitched ranges (`views::join` over per-batch row spans) — still lazy, EC-1.
- Library developers compose CONCRETE containers directly (`std::vector<CellBatch<Model>>` is the storage
  idiom; document the pattern with a worked example) — `vector<StorageUnit*>` stays banned in hot structures
  (§3.10). Observable ranges (V, SOC) evaluate through the SAME shared observable kernels (D-10) — lazily,
  cold-path, one code path with the Recorder.
- Determinism note: user-facing reductions are cold-path; anything feeding physics keeps §3.8 fixed-order.

### 3.14 Discretisation menu — FVM joins Chebyshev-spectral (M5)

The §3.2 discretisation slot gets its second citizen: **conservative finite volume** on the sphere
(PyBaMM's default method — enables like-for-like parity), with the recorded menu (parabolic-profile tier,
eigenfunction/Duhamel modal, Legendre–Galerkin) remaining future slot entries.

- `FvmDiffusion<NR>`: NR shells, nonuniform radii allowed; face fluxes with harmonic-mean D_s at faces;
  strict conservation by construction (telescoping flux sum). Same `StateSpec`/RHS-pipeline machinery
  (§3.12); scalar-generic (PC-10); rows = shell concentrations (SoA across lanes).
- Gates are order-of-convergence (MMS, observed order ≥ 1.9) + mass conservation to roundoff + converged
  cross-check against `SpectralDiffusion` (both → the same continuum; registered band) — all SHORT or
  simulation-free.
- Why it earns its place: PyBaMM-default parity fixtures stop carrying a discretisation-difference term
  (today's P7-G1 0.778 mV band is mostly that term); and DFN (§3.16) reuses the same FVM operators in x.

### 3.15 SPMe — electrolyte physics on the reserved slot (M6)

Model: the asymptotic SPMe of Marquis et al. 2019 (J. Electrochem. Soc. 166 A3693) — SPM plus an
electrolyte concentration field and correction overpotentials. `ElectrolyteDesign` (§3.3) stops being
ignored:

- **New state:** c_e(x) over three domains (neg | sep | pos), discretised on the §3.14 slot (FVM default —
  conservative; spectral variant later if profiled worthwhile). Porosity ε_k and Bruggeman transport
  corrections ε_k^b enter as compiled per-domain parameter rows (D-16 forms).
- **Physics shape:** ε_k ∂c_e/∂t = ∂_x(D_e,eff(c_e,T) ∂_x c_e) ± (1−t⁺) a_k j_k / F with domain-interface
  continuity; terminal voltage gains electrolyte concentration + ohmic correction terms and
  electrode-averaged reaction overpotentials. The FULL equation set with every prefactor is DERIVED and
  unit-checked in the M6.1 design note (`docs/derivations/spme.md`) against the paper + PyBaMM's
  implementation — never transcribed from memory.
- **Pipeline fit:** c_e rows are ODE rows (D-23 mask); corrections happen in the shared observables stage
  (§3.12) so Recorder/RHS/pack-Thevenin linearisation stay one code path; `linearizeThevenin` gains the
  electrolyte-resistance contribution.
- **Conservation identity (free gate):** total electrolyte lithium ∫ ε c_e dx is exactly conserved —
  roundoff-level check every SPMe test rides along.
- **Sensitivities:** the Dual instantiation (D-24) must compile for SPMe kernels from day one (PC-10 makes
  this automatic).

### 3.16 DFN — full Newman P2D, batch-compiled, no global DAE (M7; design note FIRST)

The hardest v5 physics milestone; D-06 (no monolithic DAE) is UPHELD. Structure (D-34):

- **Differential state per lane:** c_s at each macro node × radial shell (both electrodes) + c_e(x) + T +
  ageing rows. This bursts the SPM memory budget by design — PC-6 becomes per-tier (DFN ≤ ~16 KB/cell, Q15);
  10⁴ DFN cells ≈ 160 MB, still fine.
- **Algebraic state per lane:** (φ_s(x), φ_e(x), j(x)) satisfy charge conservation
  (∂_x(σ_eff ∂_x φ_s) = a F j; ∂_x(κ_eff ∂_x φ_e + diffusional-conductivity term) = −a F j) with
  Butler–Volmer closing j — an elliptic block-banded system. Solved INSIDE the eval pipeline as a
  per-lane damped Newton on a block-tridiagonal Jacobian, batched across lanes (SoA sweeps per Newton
  stage), workspace preallocated at build (PC-1), SPICE-style trial limiting (§3.4 idiom), warm-started
  from the previous step. This is the §3.12 "algebraic rows solved by the stepper outside the integrator"
  pattern generalised — CVODE/IDA per cell stays rejected.
- **Reuse, don't rewrite:** x-operators come from §3.14 FVM; particle models are the EXISTING slot kernels
  (spectral or FVM per node); the exponential propagator applies per-node to particle modes where the modal
  structure holds.
- **Limit ladder (the honesty gates):** DFN → SPMe (transport-limited corrections vanish) and DFN → SPM
  (fast electrolyte) within bands DERIVED from the asymptotics before running; plus nodewise charge
  conservation to roundoff, bounded Newton iteration counts, zero per-step allocations, and SHORT PyBaMM
  DFN parity fixtures.
- The M7.1 design note (state layout, Jacobian structure, failure/rollback semantics, GPU shape) gets an
  independent orthogonal-critique pass recorded in the note before implementation begins.

### 3.17 Thermal 2+1D/3D — the pouch-cell-spectral port (M8)

Port the VALIDATED spectral thermal machinery from `C:\D\git\pouch-cell-spectral` (§2.3 paths; D-36).
Scope: a `Thermal2D<NY,NZ>` (and `Thermal3D` by adding the third direction) composition slot for pouch
in-plane temperature fields, with through-thickness lumping justified by Bi_x ≪ 1 (checked by §3.18).

- **Bases:** per-direction Sturm–Liouville Robin eigenfunctions Y_n(ξ) = (Bi/b_n)·sin(b_nξ) + cos(b_nξ),
  roots of tan(b) = 2·Bi·b/(b²−Bi²) (Bi→0 degenerates to Neumann nπ — free limit gate); tensor products
  give 2D/3D; in-plane eigenvalue λ_nm = (k_y/ρC_p L_y²)β_n² + (k_z/ρC_p L_z²)γ_m². Direction-agnostic 1D
  machinery, called once per axis (the pouch repo proves this shape: `robin_*.m` + `docs/3d_thermal_design.md`).
- **Time stepping:** modal ODEs dc_nm/dt = −λ_nm c_nm + Q_nm are EXACTLY the diagonal form our exponential
  propagator (D-07) already handles — the φ₁ kernel comes from the ONE physics source (PC-10, M0.5). No new
  integrator.
- **Heat injection:** lumped cell heat Q from the electrochemical model, plus exact orthogonal projections
  for tab/edge sources (`cosine_modes.m` / `robin_tab_projection.m` idiom). Optional erfc short-time tab
  overlay (Gibbs suppression during pulse onset, `mittag/theory/tab_heat_erfc.md`) behind a flag, default
  off — port only with its registered Gibbs-kill gate.
- **Pack fit:** a Thermal2D cell exposes its surface-mean/boundary values to the D-21 compiled thermal
  graph unchanged; per-mode state rows are ordinary arena rows (SoA across lanes — a batch of pouch cells
  sweeps modes × lanes).
- **Quarantine (D-36):** the Howie distributed electro-thermal CLOSURE (off-diagonal cosine coupling,
  κ(T)/γ(T) modal products) is contested upstream (EX780–802) — NOT ported in v5. Distributed
  electrochemistry over the pouch plane is recorded as v6 candidate work, waiting on upstream closure.
- **Utilities worth porting as tools:** AAA rational approximation (Nakatsukasa–Sète–Trefethen, SIAM J.
  Sci. Comput. 2018) as an OFFLINE surrogate/kernel-compression utility (the pouch repo's own conclusion
  over erfc-image/residue/Padé); the erf toolkit for spherical-particle short-time behaviour
  (`docs/spherical_diffusion_erf.md`, verified to 3.2e-6 upstream) as analytic-oracle material for our
  particle gates.

### 3.18 Model hierarchy + regime advisor — simplifications with stated validity (M9)

Every simplification in the stack (lumped vs 2D thermal; SPM vs SPMe vs DFN; parabolic-profile tier) gets
its validity told to the user instead of assumed silently (D-37):

- **At `build()`/`compile()` (cold):** compute the dimensionless groups — Bi per direction (hL/k),
  solid-diffusion time ratio τ_d/τ_cycle = (R_p²/D_s)·(I_1C·rate/Q) shape, electrolyte transport group
  (the SPMe small parameter), σ_eff/κ_eff ratio, tab-injection nonuniformity fraction. Report via
  `pack.advise()` and `describe()`: {group, value, threshold, suggested tier, leading-order error term}.
- **Monitors (optional, cold cadence):** groups that move with state (D_s(c), h(T)) are re-checked at
  recording cadence; a violated assumption raises an explicit report/Status flag — never silent (§1
  traceability). No auto-switching of models mid-run (rejected in D-37): the advisor advises, the user
  chooses.
- **Validation:** each advisor rule carries ONE registered pair (full vs simplified SHORT run) confirming
  the predicted error order — e.g. Bi = 0.05 lumped-vs-2D difference matches the O(Bi) estimate band.
- The derivation doc (`docs/derivations/regime-map.md`) is a deliverable: the Π-group map another scientist
  can re-derive (the pouch repo's unstarted `nondimensional_analysis.md` intent, realised here).

### 3.19 Parameter & chemistry library — all of PyBaMM's sets, with provenance (M4)

Beyond Chen2020 (D-38): absorb EVERY Li-ion parameter set shipped by pinned PyBaMM 26.6.2.0 —
enumerated PROGRAMMATICALLY from the pinned install (Marquis2019, Ecker2015, Mohtat2020, ORegan2022,
Prada2013/LFP, Ai2020, OKane2022, Xu2019, … — the authoritative list comes from the enumeration script,
never hand-copied), plus BPX JSON (already shipped).

- **Pipeline:** one pinned-version Python dump script (P7-G3 pattern) emits committed JSON per set —
  SI values, OCP/OCV knot tables (nonuniform knots preserved bit-exactly, D-16 form 3), entropic
  coefficients, temperature functions, chemistry tag, **and the set's citations** (PyBaMM attaches them
  per set). C++ CI stays hermetic; a non-blocking latest-PyBaMM job flags upstream drift.
- **Round-trip gate per set:** absorbed key → internal packs → `describe()` returns the IDENTICAL SI value
  (exact equality, no simulation) — the P7-G3 discipline, now × all sets.
- **Chemistry breadth:** NMC, NCA, LFP, Si-Ox blends arrive free with the sets; sodium-ion rides the same
  SPM/SPMe mathematics where PyBaMM ships Na-ion sets (SOTA-check box in M4). Lead-acid (different model
  family) is explicitly OUT of v5 (Q16). Hysteresis OCP data lands where sets ship it (optional box,
  SOTA-check first).
- **License hygiene:** parameter VALUES are data (PyBaMM is BSD-3); every set records source + license note
  in its provenance block.

### 3.20 Citation registry — `print_citations` for exactly what a run used (M4; Volkan 2026-07-11)

PyBaMM-precedent (`pybamm.print_citations()`), done at zero hot-path cost (D-39):

```cpp
// cold path, build/compile time only:
struct Citation { std::string_view key; std::string_view doi; std::string_view bibtex; };
citations::cite(key::Chen2020);            // components mark what they use as they compile
// user surface (C++ / Python / MATLAB / Studio):
solution.citations();                      // deduplicated, EXACTLY what this run used
slide::print_citations();                  // session-aggregate convenience, BibTeX or plain text
```

- Every model tier, discretisation, solver algorithm (chord/Shamanskii, WR+Baumgarte, exponential
  integrator, AAA), parameter set, and ported method (Chu 2020, Lin 2022, Marquis 2019, DFN 1993, …)
  registers a Citation at its definition site (MC-2 makes this part of "shipping a feature").
- Gates: a Chen2020+SPM(spectral)+Mode-A run prints EXACTLY the registered expected list (no more, no
  fewer — asserted in a unit test); every registry entry has a DOI or stable key + syntactically valid
  BibTeX (format-checked test); SLIDE's own canonical citation(s) verified from `CITATION.cff`/README
  (M4.6 — verify, never guess).

### 3.21 WASM target — the dependency-free core in a browser tab (M10)

- **Toolchain:** Emscripten; single-threaded build (PC-9); exceptions already off the hot path; recorder
  degrades to in-memory sink (no filesystem) — EXPLICIT capability report, never silent.
- **Bindings:** embind by default (Q14) exposing the SAME cold surface as Python/MATLAB (§3.9 symmetry):
  designs, parameter sets by name, Experiment strings, run-with-progress-callback, Solution accessors as
  typed-array VIEWS (zero-copy into the wasm heap). Progress callback cadence keeps the UI thread honest —
  Studio runs the whole module inside a Web Worker anyway (cancel = terminate, the unibatt lesson).
- **Determinism gate:** wasm vs native on committed fixtures agree ≤ 1e-12 rel (IEEE doubles; libm
  differences are the only expected source — if the band trips, pin the offending function to our own
  scalar-generic implementation, which PC-10 makes trivial).
- **Size:** register the budget BEFORE the first successful build (ASSUMED ≤ 5 MB uncompressed core module;
  measure, record, tighten); `-Os`, no RTTI leakage, WASM SIMD flag evaluated structurally (v128 in the
  disassembly of the row sweeps), wall-clock claims stay qualified.
- **Deliverable:** npm-shaped package (`@slide/core`: .wasm + ESM loader + .d.ts) + node-based CI test +
  worker usage example. No CDN dependence anywhere (offline constraint, §3.22).

### 3.22 SLIDE Studio — visual design/simulate/analyse workspace in the browser (M11)

Not a COMSOL clone — an original, small, fast tool over our own schema. Architecture (D-41, unibatt-proven):

- **Stack:** vanilla JS ES modules + Vite; NO framework; hand-written CSS custom properties; compute in Web
  Workers over the M10 wasm module; Playwright headless E2E. Two build outputs: the Vite site AND a single
  self-contained offline HTML file (esbuild + base64-embedded wasm — imitate
  `C:\D\git\unibatt\scripts\build-single-html.js`).
- **Palette (from `glide-wasm/web/src/styles/main.css`, reuse verbatim as tokens):** dark
  `--bg-body:#0f0f23`, `--bg-panel:#1a1a2e`, `--bg-input/plot:#16213e`, text `#e0e0e0/#a0a0b8/#8888a8`,
  border `#2a2a4a`; light theme overrides; **data series = Wong colourblind-safe palette**
  `['#0072B2','#D55E00','#009E73','#F0E442','#CC79A7','#56B4E9','#E69F00']`; Oxford blue `#002147` for
  chrome ONLY, never for data (the unibatt rule). Dark/light toggle persisted in localStorage.
- **Views (left-nav workspace, one project open):**
  1. *Cell designer* — materials library (M4 sets), parameter forms with units + validation, OCP curve
     plot/editor; advisor chips (§3.18) live-update.
  2. *Pack designer* — SVG canvas: series/parallel groups, link resistances, drag/duplicate/delete; live
     compiled-netlist preview + advise() panel. Visual-programming feel without node-graph bloat.
  3. *Experiment editor* — step list ⇄ PyBaMM instruction strings, validated by the SAME C++ parser via
     wasm (one source of truth, D-11).
  4. *Run* — worker execution, progress bar, cancel, memory guard, run history.
  5. *Results* — uPlot (bundled; Q13) overlays: V/I/T/SOC, per-cell traces, pack aggregates; CSV/JSON export.
  6. *Project* — save/load `.slide.json` (versioned schema = the C++ description layer serialised, EC-5;
     round-trip gate against C++ digit-exact).
  7. *Citations & help* — the §3.20 panel: "print citations" for the current project's runs; quickstart doc.
- **Constraints:** fully offline-capable (no CDN at runtime); assert-data-not-pixels testing; accessibility
  basics (keyboard nav, contrast — the Wong palette already carries this).

### 3.23 Model-family tiers — ECM, semi-empirical, lead-acid Schiffer, storage-system (M13–M14)

Breadth beyond Li-ion porous-electrode models, all on the SAME arena/batch/registry machinery — a tier is
a composition (§3.2), never a parallel architecture:

- **ECM batch tier (M13).** 1–3 RC pairs on the §3.12 pipeline: rows {SOC, V_RC1..n, T}, OCV via D-16
  compiled curves, linear-when-fixed-R detected structurally → the Tier-0 constant-Jacobian fast path
  (§3.4.1) at pack scale. Legacy `Cell_ECM` (already validated) is the digit-parity oracle via the §5.2
  harness. This tier is also the natural home of grid-scale 10⁵-lane studies where SPM fidelity is
  unnecessary.
- **Semi-empirical degradation tier (M13).** SimSES-class stress-factor models: capacity/resistance change
  = calendar(t, T, SOC) + cycle(throughput, DOD, C-rate) superposition, integrated as slow SoA rows over
  macro steps. Formulas and coefficients enumerated from a PINNED SimSES release and their source papers
  (D-38 provenance discipline — programmatic extraction, committed JSON, citations web-verified at the
  box). Composable with ECM (the SimSES pairing) or SPM tiers.
- **Lead-acid Schiffer tier (M14).** The weighted-Ah-throughput lifetime model (Schiffer et al., J. Power
  Sources 2007 — citation web-verified at M14.1, never trusted from memory): SOC/acid-stratification/
  current-weighted throughput, corrosion-layer growth, gassing, capacity loss. A DIFFERENT model family:
  it gets its own composition tier and derivation doc, never a forced fit into particle-diffusion slots.
  Validation = the paper's own registered behaviours + limit checks. Whether VC-1 binds is decided at
  M14.1, separately for the ageing and electrical halves: "PyBaMM ships no Schiffer-ageing counterpart"
  is a from-memory claim to WEB-VERIFY there — PyBaMM's origins include lead-acid porous-electrode
  models, so the electrical half may have a counterpart even where the ageing half does not.
- **Storage-system tier (M14).** SimSES-scope system simulation as COLD orchestration over batch tiers:
  power-conversion efficiency curves (legacy `src/power_conversion` absorbed or replaced), energy-
  management strategies, application profiles (FCR, peak shaving, self-consumption) driving
  Experiment-grade current/power schedules. Hot path untouched — the system layer only schedules segments
  and reads observables (D-10). Every new file/profile parser inherits the D-29 fuzz discipline.

### 3.24 `slide.pybamm` — the drop-in subset shim (M15; D-44)

The JAX-for-NumPy positioning: `import slide.pybamm as pybamm` runs a DEFINED subset of PyBaMM scripts
unchanged, at SLIDE speed. Honesty is the design centre:

- **Coverage matrix first:** the pinned PyBaMM public API is enumerated programmatically; every symbol is
  classified {supported, mappable-later, out-of-scope} in a COMMITTED matrix. Out-of-scope forever (and
  said so): arbitrary expression-tree model surgery, custom spatial methods, solver internals.
- **The shim translates, never simulates:** a thin Python layer mapping PyBaMM's classes (`Simulation`,
  `Experiment`, `ParameterValues`, model classes for supported tiers, `Solution`'s dict surface) onto the
  native bindings (§3.9). Zero physics in the shim (D-11 one-source rule).
- **Unsupported ⇒ loud:** `NotImplementedError` naming the nearest native path — never a silent
  approximation of PyBaMM behaviour (§1 traceability applied to API semantics).
- **The gate is PyBaMM's own code:** a committed subset of the pinned release's example scripts runs
  unchanged under the shim within registered bands (G-SHIM); a non-blocking latest-PyBaMM job flags drift.

### 3.25 Integrator & predictive-acceleration research programme (M17; D-49)

Volkan's directive (2026-07-14): specialised integrators for THIS architecture (Leimkuhler-school
structure exploitation) and data-driven extrapolate-then-correct acceleration. This is RESEARCH — every
candidate enters through a derivation + registered-hypothesis gate, and a falsified candidate recorded
with its numbers is a full deliverable (§0.3(3)). Candidates, ranked by expected value:

1. **Cycle-extrapolation / projective integration for ageing (the highest-value target).** Degradation
   rows evolve smoothly over thousands of cycles: simulate k full cycles, fit the slow-row trajectory,
   extrapolate N cycles ahead (the Gear–Kevrekidis projective-integration shape), re-simulate, correct;
   accept/reject on a registered predictor-error test — extrapolation is NEVER silently trusted. Battery
   cycle-jumping precedents exist — SOTA-check at the box before designing. Hypothesis to register:
   order-of-magnitude fewer simulated cycles at bounded, measured degradation-state error, ON a scenario
   short enough that its full run is the oracle; the years-long claim is recorded separately as
   qualified extrapolation, never as the gate.
2. **PI step-size controller on the multirate outer step** (the Hairer–Wanner proportional–integral
   idiom; D-08's "error-controlled outer step" made literal): smoother step sequences and fewer
   rejections than a pure proportional rule; batch-uniform dt preserved (§3.8 corollary stands).
3. **Exponential Rosenbrock / Lawson methods** (the Hochbruck–Ostermann exponential-integrator family)
   for the nonlinearly coupled thermal + ageing outer system, reusing the ONE φ₁ source (PC-10): the
   linear-dissipative part is already exact (D-07); the open question is whether stiff nonlinear coupling
   ever limits the outer step in practice — derive the regime first, implement only if it does.
4. **Anderson acceleration on Mode C waveform relaxation** (Walker–Ni): cheap memory over past iterates,
   potentially large iteration-count cuts at 10⁴⁺ lanes; gate on iteration counters (structural, not
   wall-clock).
5. **Splitting / backward-error analysis (Leimkuhler & Reich):** our system is DISSIPATIVE, not
   Hamiltonian — symplectic integrators are the wrong tool, named here so nobody chases them; what that
   school DOES give us is splitting-error theory and the modified-equation view to justify/refine the
   D-08 Strang split.
6. **Parareal / parallel-in-time (Lions–Maday–Turinici):** likely FALSIFIED for us — batch lanes already
   saturate the hardware with spatial parallelism, and parareal's coarse-propagator overhead competes
   with an exact exponential fast path. Derive the speedup bound first; expect to record the kill,
   cheaply and simulation-free.

All citations in this section are from-memory pointers to be WEB-VERIFIED at M17.1 (M4.7 discipline)
before any implementation cites them. Each adopted method registers its Citation (§3.20) and its gate;
each rejected one gets a §8 FALSIFIED row with the derivation or measurement that killed it.

### 3.26 Newman instrumentation — thermodynamic gates, EIS, design optimisation (M18; D-50..D-52)

The Newman school's discipline applied to SLIDE: thermodynamic consistency as runtime invariants,
current distribution as the master diagnostic, closed-form limits as oracles, and the model as a DESIGN
instrument. All citations in this section are from-memory pointers — WEB-VERIFY at their boxes before
implementation cites them (M4.7 discipline).

- **Identities as gates (D-50).** (i) Entropy production σ = Σ flux·conjugate-force ≥ 0, summed in the
  shared observables stage — a sign error in transport coupling that no parity fixture targets turns
  this red; the derivation doc unit-checks every term. (ii) The FULL Bernardi energy balance
  (Bernardi–Pawlikowski–Newman 1985): Q = I(V − U_avg) + I·T·dU/dT + heat of mixing + side-reaction
  heat; gate = electrical work in − enthalpy change − heat out closes to roundoff over a registered
  cycle. Heat of mixing matters exactly at the high C-rates SLIDE claims; most tools drop it (record at
  the box what pinned PyBaMM does, for the parity band). (iii) Voltage-loss decomposition:
  V = OCV − Ση with every component from the observables stage and the sum closing to roundoff against
  the independently computed terminal V — polarization accounting as an every-step invariant AND a user
  observable.
- **Reference-electrode honesty.** φ_neg vs Li/Li⁺ at the separator interface as a first-class
  observable; plating onset (φ ≤ 0) is a §3.12 zero-crossing event, never a post-hoc threshold sweep.
- **Current-distribution advisor groups.** Wagner number, reaction penetration depth δ/L, κ_eff/σ_eff
  join the M9 Π-group table — they answer "when does uniform-reaction SPM lie", which Biot-class groups
  cannot.
- **EIS by analytic linearisation (D-51) — the sleeper capability.** The kernels are scalar-generic and
  the DFN algebraic layer is block-tridiagonal: linearise the compiled system at an operating point
  (Dual-number columns or the assembled Newton Jacobian) and solve (iωE − J)x̂ = b per registered
  frequency — Z(ω) with NO time stepping and no FFT noise. Requires optional double-layer capacitance
  rows per electrode (quasi-steady Butler–Volmer alone has no semicircle). Oracles: the closed-form
  R_ct∥C_dl semicircle, the porous-electrode transmission-line (de Levie) limit, and a Kramers–Kronig
  consistency check. Product surface: `impedance()` in C++/Python; PyBOP impedance-fit example.
- **Design mode (D-52) — the model as instrument.** `slide::design`: parameter sweeps ride the §3.1
  ensemble lanes (thickness, porosity, loading as `varied()` rows — one batch run sweeps 10⁵ designs);
  Ragone surfaces from registered power/energy protocols; gradient-based sizing (max energy subject to
  power and geometric bounds) driven by EXACT forward sensitivities — analytical Jacobians only, finite
  differences never (noise wrecks line searches; D-24's machinery is the source). Every optimisation
  benchmark registers a synthetic problem whose optimum is DERIVED analytically first — the optimiser
  must find what the derivation says exists.
- **Jacobian service.** dObservables/dθ exposed as arrays through the `simulateS1`-shaped surface in
  C++/Python for ANY consumer (PyBOP today, user optimisers tomorrow); central-FD arbiter per tier.
- **ICA/DVA.** dQ/dV and dV/dQ computed cold from Recorder snapshots; LLI/LAM signature example on a
  synthetic aged cell — degradation diagnostics the way experimentalists actually read cells.
- **PSD/MPM tier (optional).** Particle-size-distribution bins per node (Darling–Newman precedent;
  PyBaMM's MPM) on the same composition-slot machinery; gates: bin-moment conservation + collapse to
  the single-size limit, digit-exact.
- **D_s convention.** One derivation-doc section fixing chemical vs tracer diffusivity (where the
  thermodynamic factor lives) plus per-set absorption notes — documentation trap-removal, no code.

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
| D-10 | Store state snapshots, derive observables lazily | anything a later step reads is an arena row (§2.2 rule; the original "only time/Ah/Wh are path-dependent" claim was FALSIFIED and replaced by that rule), so snapshots capture everything; memory ~256 B/cell/snapshot vs unbounded frames | store per-step observable frames (v3 ECM path) |
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
| D-28 | Mode C success requires both maximum cell-current update and maximum KCL residual over every non-reference node to meet the caller tolerance. `constraint_bound` is only the maximum of that acceptance tolerance and a floating-point roundoff estimate, not a universal contraction theorem. | A tiny relaxation gain made current updates stagnate while KCL remained violated; terminal-only KCL also missed a constructed internal-node residual. The former general `(1-α)^k` claim is valid only for P4-G2's one-unknown parallel graph, where the test now derives it independently. | (a) current-delta-only convergence; (b) terminal-only KCL; (c) present the one-node Jacobi contraction as a bound for arbitrary coupled graphs. |
| D-29 | Parser fuzzing compiles an isolated `slide_core_fuzz` copy with ASan+UBSan, uses poisoned prior outputs plus independent success invariants, and runs both Debug and fast-math Release. Small committed corpora are complemented by generated 65,537-byte/4 MiB-plus boundary replays. | Instrumenting the shared production target leaked sanitizer/Windows CRT ABI choices into unrelated consumers. Deterministic double parsing alone could miss append/merge publication bugs, and calling the production netlist validator as the oracle could reproduce the same defect on both sides. | (a) mutate `slide_core` in place; (b) default/empty sentinels; (c) production-validator-only success checks; (d) cap every campaign at 65,536 bytes without separately reaching documented size limits. |
| D-30 | Production CPU parallelism uses one movable `BatchExecutor` per configured solver. It creates at most one persistent pool, parallelizes only identity-distinct archetype linearization/advance work, and keeps current gathers, result scatters, thermal assembly, reductions, and substep barriers in canonical serial order. Requested workers are capped to batch count and reported; one selected worker executes inline. | D-13 requires a pool over batches, but the original P8-G4 implementation had no production caller. A heap-owned pool keeps the executor nothrow-movable through `PackSolver`/`PackStepper` candidate publication. Cold identity checks prevent two tasks from aliasing one arena/cache. Serial fixed-order boundaries preserve reproducibility and rollback; configuration owns all thread creation before the hot path. | (a) standalone diagnostic-only pool; (b) embed a nonmovable `ThreadPool` directly in the solver; (c) scheduling-order atomics/reductions; (d) per-step or per-call thread creation. |
| D-31 | Release v4.0.0 (M2) BEFORE any v5 physics work — satisfied by M2.2 (all release gates green, artifacts built); the literal tag/publish (M2.3) is timed by Volkan and does NOT block M3+ | a released, gate-green rollback point under the new scope; release gates were already registered (archive P8-G6 text) | (a) big-bang v5 with v4 unreleased (loses the safety net); (b) idling the ladder on a tag formality (violates §0.3) |
| D-32 | Expressiveness = non-owning views over the SoA arena: `CellRef` proxies, stitched row ranges, named helper one-liners; user-facing reductions cold-path (or fixed-order when feeding physics) | std-library composability (Volkan 2026-07-11) without touching PC-1/2/3; observable ranges reuse the shared kernels (D-10) | (a) AoS cell objects / `vector<StorageUnit*>` back in hot structures; (b) helper zoo duplicating range logic (MC violation) |
| D-33 | SPMe = Marquis et al. 2019 asymptotic form; c_e on the §3.14 discretisation slot (FVM default); corrections in the shared observables stage; full equation set derived + unit-checked in `docs/derivations/spme.md` before coding | published, PyBaMM-compatible reduction; conservative discretisation gives the ∫εc_e identity as a free gate; observables stage keeps one code path (D-10) | (a) DFN-only electrolyte (pays P2D cost for every user); (b) ad-hoc empirical electrolyte resistance; (c) equations transcribed from memory |
| D-34 | DFN = semi-explicit index-1: per-lane batched block-tridiagonal damped Newton for (φ_s, φ_e, j) inside the eval pipeline, warm-started, trial-limited, preallocated; PC-6 becomes per-tier budget (Q15); D-06 upheld | keeps the batch/SoA/no-global-DAE architecture at P2D scale; reuses §3.4 solver idioms (limiting, workspace, rollback) and §3.14 operators | (a) monolithic DAE via IDA/KLU (re-rejected); (b) per-cell CVODE (unbatchable, dep-bound); (c) nested electrolyte-solid fixed-point without Newton (convergence unproven) |
| D-35 | FVM discretisation slot: conservative, harmonic-mean face diffusivities, nonuniform grids, scalar-generic; validated by MMS order + conservation + spectral cross-check | PyBaMM-default parity without a discretisation-difference term; DFN x-operators reuse it; conservation by construction | (a) FDM (non-conservative); (b) another spectral variant first (no parity payoff); (c) skipping the second discretisation (slot machinery stays unproven) |
| D-36 | Port pouch-cell-spectral's VALIDATED thermal machinery (Robin×cosine tensor bases, exact modal/Duhamel update via the ONE φ₁ source, tab projections, optional erfc overlay); QUARANTINE the contested Howie electro-thermal closure (EX780–802) until upstream closes it | the thermal/spectral parts have analytic oracles + whitepaper validation; the closure has a +0.875 K bias and 57% transfer-NRMS open upstream — porting it would import an open research problem as if settled | (a) port everything including the contested closure; (b) 3D FV thermal grid per cell (loses exactness + SoA modal cheapness); (c) waiting for upstream before porting the validated parts |
| D-37 | Regime advisor: Π-groups computed at build; `advise()`/`describe()` report {group, threshold, suggested tier, leading-order error}; optional cold-cadence monitors; NEVER auto-switch models silently | simplifications with stated validity is the §1 honesty standard applied to model choice; monitors keep long runs honest | (a) silent auto model switching (irreproducible, hides regime exits); (b) no advisor (users guess; support burden); (c) hot-path monitoring (PC-4 violation) |
| D-38 | Parameter library = ALL pinned-PyBaMM-26.6.2.0 Li-ion sets, enumerated programmatically by one dump script into committed JSON (values + OCP knots bit-exact + citations + license notes); BPX stays; per-set exact round-trip gates | multi-chemistry breadth for free, with provenance (Volkan 2026-07-11); P7-G3 discipline already proved the pattern | (a) hand-copied tables (transcription errors, no provenance); (b) runtime PyBaMM dependency (violates hermetic C++ CI); (c) Chen2020-only (single-chemistry library) |
| D-39 | Citation registry: components `cite()` at build/compile; `Solution::citations()` returns EXACTLY what the run used; `print_citations` in C++/Python/MATLAB/Studio; BibTeX validated by test; registering a citation is part of MC-2 "shipping a feature" | scholarship built in (Volkan 2026-07-11); PyBaMM precedent; zero hot-path cost | (a) static BIBLIOGRAPHY file only (can't answer "what did MY run use"); (b) citations in docstrings only (not machine-aggregatable) |
| D-40 | WASM via Emscripten: single-thread default (PC-9), embind surface mirroring §3.9, no-CDN/offline constraint, registered size budget, wasm-vs-native 1e-12 determinism gate | smallest portable step to the browser; unibatt evidence that COOP/COEP (threads) breaks offline/CDN-free deployment — so threads are NOT the v5 default | (a) pthreads-first build; (b) rewriting the core in JS/Rust; (c) server-side simulation (kills offline + privacy) |
| D-42 (ASSUMED, M0.9) | MC-5 is enforced as a three-tier header classification — `@surface api` (16 headers; the only ones a binding or doc may include), `support` (12; the value/vocabulary types the API's signatures are written in), `internal` (everything else, incl. all of `src/core/detail/`) — with the load-bearing rule that **no api or support header may include an internal one**. Value types are split out of kernel headers (`*Params.hpp`, `SpmBatchLayout.hpp`, `AgeingModelMask.hpp`) so the public factory can name them without dragging the physics in. `SpmScalarKernels.hpp` is `support` and DOES reach a user TU: PC-10 (one physics source) outranks MC-5 there, since the alternative is duplicating interpolation inside `CompiledCurve` | classification alone would have labelled the violation instead of removing it: a user TU was compiling eight kernel headers. R3 + R1 make the property inductive rather than merely measured, and `clang -MM` confirms it (26 headers/8 kernels → 21/0) | (a) tag-only audit (labels the leak, keeps it); (b) pimpl `SpmFactoryInput` (kills the aggregate, adds a cold allocation/indirection for no gain); (c) duplicate the curve interpolation to free `CompiledCurve` from `SpmScalarKernels` (violates PC-10, the stronger contract) |
| D-41 | Studio = vanilla ES modules + Vite + Web Workers over the M10 module; unibatt palette (Wong data colours, Oxford-blue chrome only); versioned `.slide.json` project schema = the C++ description layer (EC-5); single-file offline artifact; Playwright data-assert E2E | proven stack in `C:\D\git\unibatt` (shipped, offline-capable, ~2.3 MB single file); no framework = no churn; one schema = no drift between Studio and core | (a) React/Vue/framework stack (dependency churn, larger surface); (b) CDN-dependent runtime (breaks offline artifact); (c) pixel-golden screenshot tests (flaky, theme-dependent) |
| D-43 | Q16 OVERTURNED for lead-acid (Volkan directive 2026-07-14): the Schiffer weighted-Ah-throughput model enters as its own composition tier (M14) with its own derivation doc and validation gates | the user directive is the overturning evidence §4 requires; the model family is industry-standard for off-grid/renewable lead-acid lifetime prediction (SimSES ships it — verify at the M13.3 pin) | (a) keeping lead-acid out (directive says in); (b) forcing lead-acid into particle-diffusion slots (wrong physics shape) |
| D-44 | PyBaMM drop-in = `slide.pybamm` shim: committed programmatically-enumerated coverage matrix; translate-only layer over the native bindings; NotImplementedError-with-pointer for everything else; gated by pinned PyBaMM's own examples running unchanged (G-SHIM, §3.24) | JAX-for-NumPy positioning (Volkan 2026-07-14) with the §1 honesty standard applied to API semantics; D-11 keeps one physics/parser source | (a) claiming full compatibility (expression trees/spatial methods make it a lie); (b) forking PyBaMM internals (unmaintainable); (c) silent behavioural approximation of unsupported paths |
| D-45 | Device selection is a USER contract: `device={"cpu","cuda"}` per Simulation/batch, cold path only; device change = rebuild arenas, never hot-path branching; CPU is the reference implementation; every GPU-claiming tier carries a CPU↔GPU digit gate (VC-3) | PyTorch-style ergonomics (Volkan 2026-07-14) on the §3.8 design that already exists; digit gates keep both instantiations honest through the ONE physics source (PC-10) | (a) per-call device dispatch (hot-path branching); (b) GPU-only tiers (kills the reference/arbiter); (c) implicit device inference (surprising, irreproducible) |
| D-46 | Recurring hardening passes are LADDER BOXES: H1 (M12.0) and H2 (M20.1), each = {adversarial orthogonal review, bug ledger failing-test-first, mutation battery over new gates, sanitizer + fuzz lanes over new parsers, performance-counter re-baseline, simplification sweep with digit-identical gates} | Volkan 2026-07-14: verification/bugfix/perf/simplification passes must be scheduled work, not culture; M0/9B/9C proved the pattern pays (50 closed bugs, 20/20 mutation batteries) | (a) trusting per-milestone gates alone (misses cross-cutting interactions); (b) "continuous" background hardening with no box (never happens under §0.3(1)) |
| D-47 | Optional acceleration deps (Highway SIMD, CUB/Thrust device primitives; KLU already optional) enter ONLY behind a CMake option with a core fallback, admitted by a registered measurement gate on structural counters (+ qualified wall clock); a dep that loses its gate is recorded FALSIFIED and stays out | small well-justified deps allowed (§3.10) but the core must build with none (non-negotiable #4); measured-in beats fashionable-in (Volkan 2026-07-14: "does the job very well, unlike Boost") | (a) required deps; (b) Boost-class monoliths; (c) adopting on reputation without the gate |
| D-48 | Semi-empirical + storage-system tiers (SimSES scope): stress-factor calendar/cycle models as slow SoA rows; the system layer (converters, EMS, application profiles) is COLD orchestration over batch tiers reading D-10 observables; formulas/coefficients enumerated from a pinned SimSES release with D-38 provenance | grid-scale system studies are the SimSES use case Volkan wants absorbed (2026-07-14); cold orchestration keeps every PC row intact | (a) hot-path system coupling (violates PC-4 for no physics reason); (b) hand-copied coefficients (D-38 already rejected that); (c) a separate system simulator beside SLIDE (two sources of truth) |
| D-49 | Speed research (M17) is hypothesis-gated: every integrator/predictor candidate needs {derivation of the expected win, registered band, SHORT decisive test}; extrapolation predictors always carry a runtime accept/reject error test — never silently trusted; falsified candidates are deliverables | Volkan 2026-07-14 (Leimkuhler-school integrators, data-driven extrapolate-then-correct); the §3 registered-prediction discipline is what separates research from vibes | (a) implementing integrators before deriving their regime; (b) trusting extrapolation without a runtime rejection test (silent wrong answers); (c) treating a falsified candidate as wasted work |
| D-50 | Thermodynamic identities are RUNTIME GATES riding the shared observables stage: entropy production σ ≥ 0; full Bernardi energy balance (incl. heat of mixing) closing work − ΔH − heat to roundoff over a registered cycle; voltage-loss decomposition summing to the independently computed terminal V | Newman-school consistency (Volkan 2026-07-14): identities catch sign/coupling errors that no parity fixture targets; the observables stage (D-10) makes them one-source and cheap; unit tests ship per identity (MC-2) | (a) validating only at milestone gates (bugs live between them); (b) separate diagnostic code paths (drift vs the physics, violates D-10); (c) silently dropping mixing heat like most tools (any omission must be recorded with its regime) |
| D-51 | EIS = analytic linearisation of the compiled system at an operating point — solve (iωE − J)x̂ = b per registered frequency, J from Dual columns or the assembled Newton Jacobian, optional C_dl rows per electrode; NEVER time-domain sinusoid + FFT | exact to roundoff at the linearisation, no transient truncation/windowing noise; block-tridiagonal structure + scalar-generic kernels (PC-10) make it near-free; opens impedance fitting through PyBOP (Volkan: EIS is good) | (a) time-domain EIS (orders slower, noisy, tolerance-limited); (b) quasi-steady-only impedance (no semicircle — misleading); (c) a separate small-signal model (second physics source, violates PC-10) |
| D-52 | Design mode is a COLD layer over existing machinery: sweeps = ensemble lanes (§3.1 `varied()` rows), objectives/constraints read D-10 observables, gradients = EXACT forward sensitivities (D-24) — finite-difference gradients are FORBIDDEN in shipped optimisers; every optimisation benchmark registers an analytically derived optimum before the run | Volkan 2026-07-14: "more optimisation and analytical Jacobians"; FD noise wrecks line searches (D-24 already records this); ensemble lanes make a 10⁵-design sweep one batch run | (a) FD-gradient optimisers (noise floor); (b) a bespoke hot-path design engine (nothing about sizing needs the hot path); (c) unregistered benchmark problems (the optimiser "wins" nothing checkable) |

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
   continuously as leading indicators; wall clock at the checkpoints is the confirming evidence.
   **Rule (v5 restatement):** positioning PAYs (PAY-4 and all v5 PAYs 6–9) NEVER block the ladder — a missed
   target is recorded FALSIFIED with its numbers and work continues (§0.3(3)); only registered ABORT gates
   (the PAY-1/2/3 class, all already passed) can stop a milestone, and doing so requires the §8 row naming
   the abort threshold that tripped.
8. **Port hygiene (v5).** Every ported equation/algorithm cites its source repo file AND its paper at the
   definition site (§3.20 registry). Contested or unvalidated upstream results are QUARANTINED by name
   (D-36 lists them) — validated machinery only. A port is complete only when its oracle gates pass in OUR
   codebase; upstream validation is evidence for choosing what to port, never a substitute for our gates.
9. **UI & WASM verification.** Playwright headless E2E asserts DATA and DOM state, never pixels; the
   offline single-file artifact must run from `file://` with network access disabled (a test enforces it);
   wasm-vs-native determinism band 1e-12 rel on committed fixtures; Studio round-trip gates are digit-exact
   against the C++ schema serialisation (EC-5).
10. **Docs snippets are tests.** The P8-G5 fence extractor runs every new public-API snippet in CI (EC-4);
    a feature without its snippet is unfinished (MC-2).
11. **Derivation docs are deliverables.** New physics (SPMe, DFN, thermal 2D, regime map) lands with a
    `docs/derivations/*.md` another scientist can re-derive from scratch — assumptions where they enter,
    every equation unit-checked, regime + leading-order error named (user CLAUDE.md §2).

## 6. THE LADDER — one goal, ticked box by box (§0.3 protocol governs)

Every box: gates green (Debug AND Release where applicable; "where applicable" = Release skippable ONLY for
non-code boxes — docs, design notes — recorded in §8) → CHANGELOG → §8 row → small commits → tick.
Gate simulations obey the header rule: SHORT registered scenarios; think first, simulate last. Boxes marked
**(design note first)** require a `.claude/designs/` note before code. Boxes marked *(opt)* may be skipped
with a one-line §8 justification; all others are mandatory. Milestone order is dependency order — do not
reorder without a D-entry.
**Dependency rule (makes §0.3(1)/(3) precise):** within a milestone, a box depends only on lower-numbered
boxes it cites or builds on — otherwise boxes are independent and may be taken out of order when an earlier
box is BLOCKED; milestone M(k+1) opens when every MANDATORY M(k) box is ticked or recorded FALSIFIED/BLOCKED
in §8. **Artifact rule:** a box with no named gate is ticked only against a reviewer-checkable artifact
(test name, doc path, committed file) recorded in its §8 row.

### M0 — v4 closeout: finish the bug-hunt, land the simplifications (carries Phase 9B/9C)

The 9B hunt checklist, subsystem list, and tooling text are preserved verbatim in the archive §Phase 9;
the bug ledger lives at `.claude/reports/p9b-bug-ledger-2026-07-10.md` (P9-B01..B36 so far; refutations
recorded there — killed ideas stay killed).

- [x] M0.1 Land or revert the in-flight WIP (`src/core/PackSolver.cpp`,
      `tests/unit/core_ParserAllocation_test.cpp`): finish, test, commit — or revert with a recorded reason.
- [x] M0.2 P9-G1: ASan+UBSan lane green on the full Debug suite; TSan lane green on
      ThreadPool/AsyncRecorder/PackStepper tests. TSan runs via the Linux lane/WSL (unsupported on Windows);
      locally-run evidence recorded, hosted CI noted separately.
- [x] M0.3 P9-G3: error-branch coverage measured (llvm-cov) — 100% of `return Status::…` failure branches
      in `src/core` exercised; exceptions capped at ≤10, each naming a structural reason class
      (unreachable-by-construction, platform-specific, defensive-only) in the report.
- [x] M0.4 P9-G4: bug ledger complete — every bug has {failing SHORT test first, fix, ledger row};
      every refuted candidate recorded refuted.
- [x] M0.5 9C-1 ONE physics source (PC-10): extract the modal update (`exp(x)·z + dt·φ₁(x)·B·j`, expm1
      Taylor-branch φ₁) + shared observable scalar kernels into a single scalar-generic header consumed by
      `SpmPipeline` (CPU), `CudaSpmRuntime.cu` (`__host__ __device__`), and the Dual instantiation.
      Gate: all three paths digit-identical to pre-refactor outputs on named recorded cases — CPU/CUDA
      batch fixtures + the P7-G2 sensitivity traces for the Dual path (record the fixtures BEFORE
      refactoring).
- [x] M0.6 9C-2 one ageing-kernel idiom: unify SEI/CS/LAM/plating mask/scratch/lane-sweep scaffolding —
      four physics bodies, one pattern. Digit-identical gate.
- [x] M0.7 9C-3 split oversized cold files: `ParameterSet.cpp` → absorption / BPX reader / expression-AST
      TUs; `Experiment.cpp` → parser vs runner. Digit-identical + suite-count-preserved.
- [x] M0.8 9C-4 shared test harness: factor the repeated build-batch/run/compare scaffolding in
      `tests/unit/core_*` into one helper (MC-4); assertion counts must not drop.
- [x] M0.9 9C-5 public-surface audit: API vs detail classified per header (MC-5); naming-consistency pass
      (one verb per concept).
- [x] M0.10 9C-6 dead-code/line-debt sweep (`SpectralDiffusionLegacy` is deliberate — parity kernel, KEEP):
      split the coverage reporter and oversized recorder/solver/factory translation units plus
      `SpmPipeline.hpp` where a real concept boundary exists, or record a per-file cohesion justification.
      **Gate P9-G5:** every simplification digit-identical; core line count recorded before/after in §8
      (expect net reduction; growth requires written justification).

### MQ — Quality wave: clean, re-derive, hunt (2026-07-23 directive; sits between M0 and M1 in ladder order)

Volkan's directive: before the ladder resumes, spend one wave making what EXISTS top-quality science —
clean the repository, independently re-derive the mathematics, and hunt the remaining logical and
performance mistakes. This wave is deliberately LOOSER than most milestones: boxes state WHAT and the
gate; scope and method inside each box are yours (§0.3(4) governs — improvise, record). Two fixed
points: the evidence discipline (§5) applies in full — a re-derivation that finds nothing is a PASS
with the derivation committed, and a suspected bug becomes a failing SHORT test before it becomes a
fix — and the wave must not balloon: when a finding is real but large, give it an owner box in the
ladder (record where) instead of dissolving the wave into it. Testing depth, the GUI, and the rest of
the grand goal are NOT re-planned here — they are the existing ladder (tests ride every box per
MC-2/MC-4; the GUI is M11 SLIDE Studio; releases are M2/M12/M20).

- [x] MQ.1 Restore the Debug-AND-Release standard first (the wave's baseline): the 2026-07-21..23
      landings were validated Debug-only (recorded in §8). Run the full suite on native Debug,
      fast-math Release, AND CUDA at current HEAD; record counts against the M0.10 baseline
      (58/58 expected per lane); any delta is explained or ledgered before MQ code changes land
      on top.
- [ ] MQ.2 Quality-backlog owner box: triage EVERY remaining verified finding in
      `.claude/reports/code-quality-pass-2026-07-21-survivors.json` (68 survived its adversarial
      verification; the ladder-detection O(cells²) scan, the Sei/Lam/LithiumPlating Arrhenius PC-10
      duplication, two blind-test oracles, the warning-policy hookup, and the CWD dependence are
      already landed — the report's landed list is authoritative). Each finding ends in exactly one
      recorded state: applied (digit-identical gate, or a registered behaviour-change gate where the
      behaviour WAS the bug), REFUTED (with the evidence), or deferred-with-owner (a named ladder
      box). No silent drops — the report gains a disposition column covering all 68. Named landmines
      that may not be silently deferred: the line-for-line duplicated topology derivation, the
      thrice-reused Mode-C relaxation buffer, `crc32` defined twice, the fast-math-unsafe FP gate,
      the `substeps`-repeats-dt-rather-than-subdividing semantics (documented 2026-07-21; DECIDE the
      intended semantics here — a behaviour change needs its own registered gate and CHANGELOG line),
      and the SurfaceCrack Arrhenius association (the recorded remaining PC-10 debt — same decision
      shape: route through the one source and accept recorded bit changes, or record why not).
      MQ.1 additionally registered two build-policy findings outside the 68: Clang 21 deprecates
      the redundant `-Ofast` spelling, and the CUDA build exposes the retained legacy
      `benchmark_LP_cases` narrowing warning through `Deep_ptr.hpp`; add and disposition both here
      rather than silently dropping them.
- [ ] MQ.3 Repo hygiene sweep: audit the tracked tree (670 files at this writing) — stale or
      contradictory docs (develop/TODO.md header vs PLAN §8, README claims, CONTRIBUTING,
      LESSONS.md), dead scripts, anything tracked that should be ignored; fix `.gitignore` gaps.
      Close the repo-wide MC-3 header-contract debt identified by MQ.2
      (`mc3-contract-comments-missing` and `mc3-contract-comments-below-bar`): census every
      `src/core` header and add a content-aware structural gate requiring a substantive owner,
      exact PLAN section, and hot/cold classification in its leading 3–6 contract lines; a
      presence-only comment-count gate is insufficient.
      For UNTRACKED debris (the ~40 stale `build-*` trees, stray `build-*.log` files, the
      Windows-hazard file literally named `nul` at repo root) produce a deletion-PROPOSAL list in the
      hygiene report — do NOT delete untracked files or data yourself (§0.3(5): disk deletion is
      Volkan's one word; `data/` and `results/` stay read-only regardless). Gate: hygiene report
      committed; develop/TODO.md and README truthful against PLAN §8.
- [ ] MQ.4 Derivation-inventory matrix: enumerate every physics kernel, numerical method, and solver
      algorithm in `src/core` (+ the CUDA mirrors) and map each to its derivation artifact
      (`docs/derivations/*`, a plan §, a cited paper, or NONE). Commit the matrix as
      `docs/derivations/INDEX.md`. Every NONE becomes an MQ.5 work item; the matrix is the wave's
      coverage gate — MQ.5 cannot tick while a NONE row is unowned.
- [ ] MQ.5 Independent re-derivation sweep (the heart of the wave). For each inventory row: re-derive
      the result from first principles WITHOUT reading the implementation first, unit-check every
      equation, state each approximation's regime and leading-order error term, THEN diff against the
      code and any existing doc. Agreement → commit/refresh a derivation doc a reader can re-derive
      from scratch. Disagreement → failing SHORT test first, then fix or FALSIFY your own derivation;
      ledger row either way (a wrong derivation honestly recorded is also a deliverable). Priority
      order, highest consequence first: the modal update + φ₁ Taylor branch (M0.5's one source), the
      four ageing kernels (SEI/CS/LAM/plating — units and regime of every rate prefactor), Thevenin
      linearisation + chord/Shamanskii update, Mode B ladder elimination, the Mode C WR+Baumgarte
      convergence criterion (our own synthesis — deserves the hardest look), D-21 thermal-graph energy
      bookkeeping, Dual-number algebra at its boundary cases (pow at zero, saturated branches —
      M0.6 already caught two), OCV/entropic interpolation. Gate: `INDEX.md` shows zero unowned NONE
      rows; every discrepancy has a ledger row with its test.
- [ ] MQ.6 Fresh logic bug-hunt round (D-46 template, loop-until-dry): adversarial pass over
      `src/core`, the bindings, AND the test gates themselves (M0.8 falsified two of its own
      preregistered mutations — gates are code too); every candidate gets a failing SHORT test or a
      recorded refutation in the ledger; the round repeats until a full pass yields nothing new.
- [ ] MQ.7 Performance-mistake hunt, structural evidence only (D-27: wall-clock on this machine is
      not evidence): complexity audit of every `compile()`/`build()`/per-step path for super-linear
      scans (the O(cells²) ladder detection just fixed may have siblings), allocation audit against
      PC-1 on paths the counting allocator does not yet watch, PC-6 memory-budget re-check,
      counter-regression assertions extended wherever the missing assertion WAS the finding. Wins and
      non-findings both recorded with counters before/after.
- [ ] MQ.8 Wave closeout: full Debug + Release + CUDA suites green; CHANGELOG consolidated; a §8 row
      with the wave's numbers (findings applied/refuted/deferred, derivations added, bugs found and
      killed); develop/TODO.md synced. Then the ladder resumes at M1.0 — the grand goal is unchanged.

### M1 — PyBOP integration suite (carries Phase 10; ALL fits on SHORT synthetic data: ≤600 s simulated, dt ≥ 10 s, ≤4 lanes)

Derive before running: identifiability and noise floors computed analytically FIRST, bands registered from
the derivation, not from a trial fit.

- [ ] M1.0 Pay the MC-1 TEST-file debt owed by M0.8/M0.10 (explicit owner box — do not inherit again):
      split `core_ParserAllocation_test.cpp` (1,229), `core_PackSolver_test.cpp` (1,198),
      `core_Experiment_test.cpp` (1,168), `core_AsyncRecorder_test.cpp` (1,118),
      `core_ParameterSet_test.cpp` (1,007) along fixture/subsystem boundaries; the SUM of
      successor-binary assertion/test-case counts must be ≥ each predecessor's frozen floor, and floors
      are re-frozen per successor binary (M0.8 discipline); harness idiom (MC-4) preserved.
- [ ] M1.1 (10A) gradient regression kept pinned (P7-G2 L-BFGS-B via `simulateS1`).
- [ ] M1.2 (10B) gradient-FREE optimiser (CMA-ES or XNES) through PyBOP: 2-parameter synthetic recovery,
      registered looser band.
- [ ] M1.3 (10C) multi-parameter GITT-style fit: ≥4 params {D_s,n, D_s,p, R_contact, Q} on a synthetic
      pulse train (~10 × 60 s); bands from a pre-run identifiability analysis (normalised-sensitivity Gram
      rank + conditioning — pure derivation).
- [ ] M1.4 (10D) noise robustness: registered σ = 1 mV Gaussian; recovery bands from the CRLB estimate
      BEFORE the one decisive run.
- [ ] M1.5 (10E) version drift: pinned PyBOP job blocking; latest-PyBOP job non-blocking.
- [ ] M1.6 (10F) `simulateS1` vs central-FD arbiter as installed-wheel pytest (short horizon, 3 params).
      **Gates P10-G1** (M1.2–M1.4 recover truth within pre-registered bands in installed-wheel CI) and
      **P10-G2** (identifiability + CRLB derivations committed as docs).

### M2 — v4.0.0 release (was P8-G6; D-31)

- [ ] M2.1 Version sync (CMake/package/MATLAB = 4.0.0); consolidated CHANGELOG with migration/known-limit
      notes; dependency licenses; qualified-timing caveats (D-27) and Phase-9 outcome counts in the notes.
- [ ] M2.2 Reproducible artifacts: CPU wheel, CUDA wheel or build instructions, MEX package; final
      Debug/Release CPU, installed-wheel Python, MATLAB batch, CUDA, async recording, release-consistency
      checks green. Legacy v3 façades stay compiled/runnable (Q3).
- [ ] M2.3 Tag + publish: prepare everything; the annotated `v4.0.0` tag/publish is the §0.3(5) exception —
      one-word go from Volkan. DO NOT idle waiting: proceed to M3 immediately after M2.2, and M5+ physics
      proceeds even if the go is still pending (D-31 is satisfied by M2.2's release-ready state; this box
      stays excluded from the §0.3(1) selection until the go arrives).

### M3 — Expressiveness layer (§3.13; the developer-joy milestone)

- [ ] M3.1 **(design note first)** view/proxy API sketch: `CellRef`, `pack.cells()`, stitched row ranges,
      helper list; EC-contract mapping.
- [ ] M3.2 Row ranges + observable ranges implemented; `static_assert` std::ranges conformance; multi-batch
      stitching lazy (EC-1/EC-2).
- [ ] M3.3 Named helpers (`total_current`, `max_temperature`, `weakest_cell`, `soc_spread`, `imbalance`)
      with docs snippets run in CI (EC-3/EC-4).
- [ ] M3.4 Concrete-container developer pattern documented with a worked example
      (`std::vector<CellBatch<Model>>`; §3.10 rules restated where developers will read them).
- [ ] M3.5 Structural-counter regression CI: allocations/factorisations/iterations asserted per commit
      (the "beyond v4.0" candidate, landed).
- [ ] M3.6 *(opt)* MSVC CI lane (portability goal; known slower codegen is a measurement, not a blocker).
- [ ] M3.7 **Gate G-EXPR (mandatory):** an example file exercising ≥6 std::ranges/std algorithms on a pack
      compiles and runs in CI; range iteration shows ZERO allocations (counter test); docs snippets green.

### M4 — Parameter & chemistry library + citation registry (§3.19–3.20; before new physics, which needs both)

- [ ] M4.1 Citation registry core: `Citation`/`cite()`/`Solution::citations()`/`print_citations` in C++,
      Python, MATLAB; BibTeX syntax test; used-exactly semantics test (G-CITE: Chen2020+SPM+Mode-A run
      prints exactly the registered expected list).
- [ ] M4.2 Pinned-PyBaMM enumeration + dump script → committed JSON per Li-ion set (values, OCP knots
      bit-exact, entropic terms, temperature functions, chemistry tag, citations, license note).
- [ ] M4.3 Absorption + exact round-trip gate PER SET (P7-G3 pattern × all sets); deprecation-alias table
      extended as the enumeration reveals older keys.
- [ ] M4.4 SOTA-check box: Na-ion sets present in 26.6.2.0? hysteresis OCP data available per set? Land
      what exists behind the same pipeline; record what does not (one §8 line).
- [ ] M4.5 Chemistry tags surfaced in `describe()` (Studio consumes them in M11.3).
- [ ] M4.6 SLIDE's own canonical citation(s): read `CITATION.cff`/README, VERIFY, register — never guess.
- [ ] M4.7 Retro-register citations for everything already shipped (Chebyshev spectral, chord/Shamanskii,
      WR+Baumgarte, exponential integrator, liionpack netlist schema, Chen2020, BPX, PyBaMM compat).
      Every retro-registered DOI/reference is WEB-VERIFIED to resolve to the named paper (checked and
      recorded — plan-text references like the Miekkala–Nevanlinna DOI and the J. Energy Storage PII are
      unverified until this box).
- [ ] M4.8 C++ description-layer JSON serialisation (`.slide.json` schema v1, EC-5): CellDesign / PackDesign
      / Experiment / ParameterSet-reference round-trip digit-exact in C++; schema doc committed; reuses the
      hardened BPX JSON reader infrastructure (D-29 fuzz discipline applies). Python/MATLAB save/load ride
      this same implementation. (Studio M11.1/M11.8 EXTEND schema v1 — they do not invent it.)

### M5 — FVM discretisation slot (§3.14)

- [ ] M5.1 **(design note first)** grid/faces/BC/StateSpec layout; harmonic-mean justification;
      tail-lane and SIMD notes.
- [ ] M5.2 `FvmDiffusion<NR>` scalar-generic kernel on the §3.12 pipeline; registry entries for
      SPM<FVM, …> compositions.
- [ ] M5.3 MMS order gate: observed order ≥ 1.9 (manufactured solution; simulation-free or one short
      transient).
- [ ] M5.4 Conservation gate: particle mass to roundoff over a short registered cycle.
- [ ] M5.5 Cross-discretisation gate: FVM(refined) vs Spectral(refined) converge to the same continuum
      within a registered band on a SHORT scenario.
- [ ] M5.6 PyBaMM-default parity fixture re-run under FVM — register the hypothesis BEFORE the run
      (expected: band tightens vs P7-G1's 0.778 mV because the discretisation-difference term drops).
- [ ] M5.7 *(opt)* CUDA FVM sweep port (PC-10 header extension; digit-gate vs CPU).

### M6 — SPMe (§3.15)

- [ ] M6.1 **(design note first)** full derivation doc `docs/derivations/spme.md`: Marquis 2019 equation
      set, every prefactor unit-checked, small parameter + regime stated, PyBaMM-implementation
      cross-check; SOTA-check the citation.
- [ ] M6.2 Electrolyte rows (3-domain c_e on the FVM slot), Bruggeman/porosity parameter rows,
      interface handling.
- [ ] M6.3 RHS kernel + observables corrections (electrolyte ohmic + concentration overpotentials;
      Thevenin linearisation extended); scalar-generic; Dual compiles (sensitivities ride along).
- [ ] M6.4 Conservation gate: ∫εc_e exactly conserved (roundoff) — rides every SPMe test.
- [ ] M6.5 Limit gate: D_e → large ⇒ SPMe → SPM within a band DERIVED from the asymptotics.
- [ ] M6.6 PyBaMM SPMe parity fixtures (C/50 + 1C, SHORT), bands registered from the
      discretisation-difference derivation before comparison.
- [ ] M6.7 MMS gate on the coupled electrolyte-particle system.
- [ ] M6.8 Experiment/Cycler/Recorder/pack integration (SPMe cells in packs — Thevenin interface already
      abstracts the tier); events unchanged.
- [ ] M6.9 *(opt)* CUDA SPMe.
- [ ] M6.10 PAY-6 positioning vs PyBaMM SPMe (qualified dev-host; hypothesis registered before the run).

### M7 — DFN (§3.16; the summit — design rigour over speed of landing)

- [ ] M7.1 **(design note first, with a recorded orthogonal self-critique pass)**: state layout
      (N_x × N_r × 2 + c_e + T), algebraic-system structure, batched block-tridiagonal Newton, workspace,
      failure/rollback semantics, PC-6 tier budget (Q15), GPU shape. The note also records the Newman
      BAND(J) lineage of the batched block-tridiagonal Newton (citation registered §3.20, web-verified)
      and DECIDES internal nondimensionalisation for Jacobian conditioning (Newman-school practice;
      feeds M18's EIS solve).
- [ ] M7.2 x-operators from the FVM slot; per-node particle models from existing slot kernels.
- [ ] M7.3 Algebraic layer: batched damped Newton with trial limiting + warm start; zero per-step
      allocations (counter gate).
- [ ] M7.4 Charge-conservation gate: nodewise KCL to roundoff on every accepted step (cheap invariant).
- [ ] M7.5 Limit-ladder gates: DFN → SPMe and DFN → SPM within DERIVED bands (registered before runs).
- [ ] M7.6 PyBaMM DFN parity fixtures (SHORT: C/10 + 1C segments), bands registered first.
- [ ] M7.7 Newton-iteration bound gate (typical ≤ 5, hard cap with rescue path per §3.4 idioms).
- [ ] M7.8 Experiment/events/Recorder/pack integration; sensitivities forward-mode for the Q10 set.
- [ ] M7.9 *(opt)* CUDA batched-tridiagonal DFN.
- [ ] M7.10 PAY-7 positioning vs PyBaMM DFN IDAKLU (qualified; hypothesis registered first).

### M8 — Thermal 2+1D/3D port (§3.17; read the pouch repo first — paths in §2.3)

- [ ] M8.1 **(design note first)** port scope + citation set (Chu 2020, Lin 2022, pouch whitepaper —
      SOTA-check: titles/DOIs verified against the actual papers, never trusted from digests);
      quarantine list restated (D-36); mapping from MATLAB files to C++ headers.
- [ ] M8.2 Robin 1D machinery: root solver for tan(b) = 2·Bi·b/(b²−Bi²), eigenfunctions, normalisation,
      projections. Oracle gate: roots vs independent bisection at randomised Bi; Bi→0 → nπ limit.
- [ ] M8.3 `Thermal2D<NY,NZ>` tensor-product slot; modal state rows; exact Duhamel update through the
      ONE φ₁ source (M0.5). Gate: digit-match vs closed-form single-mode solutions.
- [ ] M8.4 Heat injection: lumped-Q uniform projection + tab/edge exact projections. Gate: projection
      coefficients digit-match committed reference outputs from the MATLAB `robin_tab_projection.m` /
      `cosine_modes.m` on a recorded case.
- [ ] M8.5 Analytic gates: Carslaw–Jaeger separable transients; energy conservation to roundoff;
      Bi→0 digit-match to `ThermalLumped`.
- [ ] M8.6 *(opt)* erfc short-time tab overlay behind a flag with its registered Gibbs-kill gate.
- [ ] M8.7 Thermal3D by the third direction call; emit Bi_x at build for later M9 consumption (no forward
      dependency — M9 reads what M8 already computed).
- [ ] M8.8 *(opt)* AAA rational surrogate as an offline utility (cold tooling, cite NST 2018).
- [ ] M8.9 Pack coupling: Thermal2D surface values ↔ D-21 thermal graph consistency gate (two-cell SHORT
      exchange conserves energy).

### M9 — Regime advisor + model-hierarchy guide (§3.18)

- [ ] M9.1 Derivation doc `docs/derivations/regime-map.md`: Π-groups, thresholds, leading-order error
      terms — re-derivable from scratch.
- [ ] M9.2 Build-time computation + `advise()`/`describe()` surfaces (C++/Python/MATLAB).
- [ ] M9.3 Cold-cadence monitors with explicit violation reports (never hot path, never silent).
- [ ] M9.4 Per-rule validation gates: one registered full-vs-simplified SHORT pair per rule confirming the
      predicted error order (lumped-vs-2D at Bi = 0.05; SPM-vs-SPMe at low/high C-rate; SPMe-vs-DFN at one
      registered stress case).
- [ ] M9.5 Docs: "which model when" user guide — the hierarchy, honestly stated.

### M10 — WASM core (§3.21)

- [ ] M10.1 Emscripten toolchain file + single-thread core build; capability report (no fs → in-memory
      sink; no threads → sequential batch loop) EXPLICIT (PC-9 gate).
- [ ] M10.2 embind surface (Q14): designs, parameter sets, Experiment strings, run-with-progress,
      Solution typed-array views; worker usage example.
- [ ] M10.3 Determinism gate: wasm vs native ≤ 1e-12 rel on committed fixtures for SPM, SPMe, AND DFN
      (goal 12 says the FULL core runs in a browser tab — all tiers exercised).
- [ ] M10.4 Size budget: registered BEFORE first successful build (ASSUMED ≤ 5 MB uncompressed; record,
      tighten); WASM SIMD evaluated structurally (v128 present in row-sweep disassembly).
- [ ] M10.5 npm-shaped package + node CI job.
- [ ] M10.6 PAY-8: wasm ≤ 2× native per-step on the fixture (qualified; register first).

### M11 — SLIDE Studio (§3.22; the user-joy milestone)

- [ ] M11.1 **(design note first)** information architecture (7 views), `.slide.json` schema EXTENSIONS
      over the M4.8 schema v1 (C++ round-trip stays the source of truth — extend, don't fork), palette
      tokens (from §3.22), offline constraint, worker protocol.
- [ ] M11.2 Scaffold: Vite + vanilla ES modules + workers + theme (dark/light, unibatt tokens).
- [ ] M11.3 Cell designer (materials library from M4, unit-validated forms, OCP plot/editor,
      advisor chips).
- [ ] M11.4 Pack designer (SVG canvas: series/parallel groups, link R, drag/duplicate; live netlist
      preview + advise panel).
- [ ] M11.5 Experiment editor (step list ⇄ PyBaMM strings via the wasm-hosted C++ parser — one source of
      truth, D-11).
- [ ] M11.6 Run panel (worker, progress, cancel = terminate, memory guard, run history).
- [ ] M11.7 Results explorer (uPlot bundled, Q13; overlays, per-cell traces, CSV/JSON export; Wong palette
      for data ONLY).
- [ ] M11.8 Project save/load; EC-5 round-trip gate digit-exact against C++ serialisation.
- [ ] M11.9 Citations & help view (print_citations for the project's runs — §3.20 surfaced to end users).
- [ ] M11.10 Single-file offline build (esbuild + base64 wasm; imitate
      `C:\D\git\unibatt\scripts\build-single-html.js`); `file://` + no-network test.
- [ ] M11.11 Playwright headless E2E: design → save → load → run (tiny sim) → plot → citations.
- [ ] M11.12 *(opt, v2)* parameter sweeps via batch lanes (the §3.1 ensemble framing in the browser),
      comparison runs, BPX import, netlist CSV import.
- [ ] M11.13 PAY-9: offline file loads and completes its first tiny sim; judged on STRUCTURAL budget
      (wasm instantiation + registered step count complete without error) plus a generous qualified
      wall-clock ceiling registered before the run (dev-host, ≤10 s ASSUMED — record, tighten).

### M12 — v5.0.0 release

- [ ] M12.0 Hardening pass H1 (D-46 template) over everything landed since M0: adversarial orthogonal
      review; bug ledger failing-test-first; mutation battery over gates added in M3–M11; sanitizer +
      fuzz lanes over parsers added since M0 (schema JSON, netlist, Experiment extensions); performance-
      counter re-baseline; simplification sweep (digit-identical gates); MC-1 line-debt check.
- [ ] M12.1 Docs sweep: every milestone's derivation docs + guides linked from the docs site; API reference
      regenerated; changelog consolidated with migration notes.
- [ ] M12.2 Version sync + artifacts: wheels (CPU/CUDA), MEX, npm wasm package, Studio site build + offline
      HTML; all CI matrices green. (No public deploy in this box.)
- [ ] M12.3 Release notes: model-tier matrix, qualified-timing caveats, citation instructions, known limits
      (distributed electro-thermal closure quarantined D-36 — unscheduled, beyond v6.0.0; adjoint = an
      M16.5 evaluation only; `.yp`/`.observe()` and FMU export unscheduled beyond v6.0.0; lead-acid
      Schiffer ARRIVES in v6 at M14 per D-43).
- [ ] M12.4 Annotated `v5.0.0` tag/publish + Studio Pages deploy — §0.3(5) exception, one-word go from
      Volkan (all outward publishing lives HERE, not in M12.2).

### M13 — ECM core tier + semi-empirical degradation tier (§3.23; D-48; goal 15)

- [ ] M13.1 **(design note first)** ECM batch tier: 1–3 RC pairs on the §3.12 pipeline (StateSpec rows,
      scalar-generic, registry entries); Tier-0 constant-Jacobian fast path (§3.4.1) detected
      structurally; legacy `Cell_ECM` digit-parity gates via the §5.2 harness across 1/2/3-RC fixtures.
- [ ] M13.2 ECM in packs + Experiment/Recorder/events; PyBaMM Thevenin-model parity fixture (VC-1, band
      registered first); PC-6 ECM budget registered (ASSUMED ≤ 64 B/cell — measure, tighten).
- [ ] M13.3 **(design note first)** semi-empirical tier: calendar + cycle stress-factor superposition as
      slow SoA rows; formulas/coefficients enumerated from a PINNED SimSES release + source papers (D-38
      provenance; citations web-verified at this box); Dual instantiation compiles.
- [ ] M13.4 Validation: closed-form throughput/time cases (pure derivation, registered bands); one
      registered scenario cross-checked against the pinned SimSES run (qualified, VC-5).
- [ ] M13.5 PyBOP fit of semi-empirical parameters on synthetic data (M1 discipline: identifiability +
      CRLB derivation first, bands registered before the one decisive run) (VC-2).
- [ ] M13.6 PyBOP fit of ECM parameters (R0, R_i/C_i, OCV scaling) on synthetic data — the most common
      PyBOP use case (M1 discipline; VC-2).

### M14 — Lead-acid Schiffer tier + storage-system tier (§3.23; D-43, D-48; goal 15)

- [ ] M14.1 **(design note first)** derivation doc `docs/derivations/schiffer.md`: weighted Ah throughput
      (SOC/stratification/current weighting), corrosion-layer growth, gassing, capacity loss — every
      equation unit-checked, regime + deliberate omissions stated; the Schiffer 2007 citation
      WEB-VERIFIED here.
- [ ] M14.2 Schiffer tier implementation (own composition tier, NOT particle-diffusion slots); SHORT
      validation cases reproducing the paper's registered behaviours; Recorder/Experiment/series-string
      pack integration. VC-1 binding decided at M14.1 per half (ageing vs electrical) and recorded;
      Schiffer VC-2 fit-relevance decided and recorded HERE (a fit box or a justified waiver).
- [ ] M14.3 Storage-system tier: converter efficiency curves (absorb or replace legacy
      `src/power_conversion`), EMS + application profiles (FCR, peak shaving, self-consumption) as cold
      orchestration over batch tiers; any new profile/file parser inherits the D-29 fuzz discipline.
- [ ] M14.4 System validation: energy conservation across the converter chain to roundoff; one registered
      application-profile scenario vs the pinned SimSES run (qualified, VC-5); docs quickstart.

### M15 — `slide.pybamm` drop-in shim (§3.24; D-44; goal 16)

- [ ] M15.1 **(design note first)** coverage matrix: pinned PyBaMM public API enumerated
      programmatically; every symbol classified {supported, mappable-later, out-of-scope}; matrix
      committed; NotImplementedError-with-pointer policy defined.
- [ ] M15.2 Shim implementation over the native bindings (translate-only; zero physics; D-11 symmetry).
- [ ] M15.3 **Gate G-SHIM (mandatory):** a committed subset of the pinned release's own example scripts
      runs UNCHANGED under `import slide.pybamm as pybamm`, with named `Solution` variables per script
      compared against committed pinned-PyBaMM outputs, band registered per script (VC-4 — script
      completion alone never ticks this); pass-list committed; non-blocking latest-PyBaMM drift job.
- [ ] M15.4 Docs page "shim coverage, honestly": the matrix rendered; unsupported paths named with their
      native alternatives.

### M16 — Device UX, gradients, measured-in acceleration (D-45, D-47; goals 17, 19)

- [ ] M16.1 `device=` API uniform across ALL tiers; the GPU-instantiated list is RECORDED in this box,
      minimum {SPM (exists), ECM (implement here — the cheap tier, Tier-0 structure)}; SPMe/DFN join if
      M6.9/M7.9 landed (each skip already carries its §8 line); any tier without a GPU instantiation
      reports the capability gap explicitly (PC-9 idiom — degrade EXPLICITLY, never silently).
- [ ] M16.2 CPU↔GPU digit gates for every tier on the M16.1 recorded list, on registered fixtures
      (VC-3) — the gate names the list, so it cannot pass vacuously.
- [ ] M16.3 PyBOP per-tier fit smoke: SPMe and DFN synthetic recoveries through PyBOP (M1 discipline,
      SHORT) — VC-2 closure across tiers.
- [ ] M16.4 *(opt)* Highway (SIMD row sweeps) / CUB–Thrust (device reductions) measured in per D-47; a
      dep that loses its gate is recorded FALSIFIED and stays out.
- [ ] M16.5 *(opt)* adjoint/Enzyme SOTA-check: derive the n_θ crossover vs forward Dual FIRST; prototype
      only if the derivation says it pays; decision recorded either way (extends D-24); if adopted,
      Jacobians reach PyBOP through the same `simulateS1` surface.

### M17 — Integrator & predictive-acceleration research (§3.25; D-49; goal 22)

- [ ] M17.1 SOTA + citation verification pass: cycle-jumping/projective-integration literature, PI
      step controllers, the exponential-Rosenbrock family, Anderson acceleration, parareal — every §3.25
      from-memory pointer WEB-VERIFIED or corrected; candidate ranking updated in the design note.
- [ ] M17.2 PI step-size controller on the multirate outer step (Hairer–Wanner idiom): registered
      step-rejection/accuracy gates on SHORT scenarios; batch-uniform dt preserved (§3.8).
- [ ] M17.3 **(design note first)** cycle-extrapolation predictor–corrector for ageing: simulate k
      cycles, fit slow rows, extrapolate N, re-simulate, correct; runtime accept/reject error test
      MANDATORY (D-49); registered hypothesis: ≥10× fewer simulated cycles at bounded, measured
      degradation-state error on the SHORT registered scenario whose full run is the oracle; the
      long-horizon claim is a separate qualified, non-blocking record.
- [ ] M17.4 *(opt)* exponential-Rosenbrock outer stage: derive FIRST whether stiff thermal/ageing
      coupling ever limits the outer step; implement only if the derivation says yes; otherwise record
      the kill.
- [ ] M17.5 *(opt)* Anderson acceleration on Mode C: gate on iteration counters at 10⁴ lanes.
- [ ] M17.6 Parareal speedup-bound derivation (expect FALSIFIED; record the kill with numbers — cheap,
      simulation-free).

### M18 — Newman instrumentation: thermodynamic gates, EIS, design optimisation (§3.26; D-50..D-52; goal 24)

Every box here lands with its named unit tests per MC-2/MC-4 (`tests/unit/core_X_test.cpp`, harness
idiom, oracle-first: closed forms and identities before fixtures) — a Newman feature without its
analytic-oracle test is unfinished.

- [ ] M18.1 **(design note first)** thermodynamic identity gates (D-50): entropy-production observable
      with σ ≥ 0 gating registered SPMe/DFN scenarios; full Bernardi energy balance incl. heat of mixing
      (Bernardi–Pawlikowski–Newman 1985, citation WEB-VERIFIED here); work − ΔH − heat closes to
      roundoff over a registered cycle; derivation doc `docs/derivations/energy-balance.md`, every term
      unit-checked; record what pinned PyBaMM includes, for the parity band.
- [ ] M18.2 Voltage-loss decomposition identity in the shared observables stage: components sum to the
      independently computed V to roundoff on every accepted step of registered scenarios; surfaced in
      Recorder/C++/Python with a docs snippet (EC-4).
- [ ] M18.3 Virtual reference electrode: φ_neg vs Li/Li⁺ observable; plating onset as a §3.12
      zero-crossing event, with a registered fixture proving the event fires at the analytic crossing
      (not the legacy per-step threshold overshoot).
- [ ] M18.4 Newman current-distribution groups (Wagner number, penetration depth δ/L, κ_eff/σ_eff) added
      to the M9 advisor + regime-map doc, each with one registered full-vs-simplified validation pair
      (M9.4 pattern).
- [ ] M18.5 **(design note first)** EIS (D-51): optional C_dl rows per electrode; (iωE − J) solve on the
      compiled structure per registered frequency grid; oracles = closed-form R_ct∥C_dl semicircle,
      transmission-line (de Levie) limit, Kramers–Kronig consistency check; `impedance()` in C++/Python;
      PyBOP impedance-fit example (VC-2).
- [ ] M18.6 **(design note first)** design mode (D-52): ensemble-lane sweeps (thickness/porosity/
      loading) → Ragone surface on a registered protocol; gradient-based sizing (max energy s.t. power +
      geometric bounds) with EXACT forward-sensitivity gradients; benchmark = a synthetic problem whose
      optimum is DERIVED analytically and registered before the run.
- [ ] M18.7 Jacobian service: dObservables/dθ arrays through the `simulateS1`-shaped surface
      (C++/Python), documented; central-FD arbiter per tier (extends M1.6 to SPMe/DFN/ECM).
- [ ] M18.8 ICA/DVA tool: dQ/dV and dV/dQ from Recorder snapshots (cold); LLI/LAM signature example on a
      synthetic aged cell; docs page.
- [ ] M18.9 *(opt)* PSD/MPM tier: particle-size bins per node on the composition-slot machinery
      (Darling–Newman / PyBaMM-MPM precedents, verified at this box); gates: bin-moment conservation +
      digit-exact collapse to the single-size limit.
- [ ] M18.10 D_s convention section in the derivation docs (chemical vs tracer diffusivity, where the
      thermodynamic factor lives) + per-set absorption-table audit note; no code.

### M19 — Toolchain matrix + documentation website (goals 20, 21)

- [ ] M19.1 Mandatory compiler/OS matrix over the APPLICABLE pairs (MSVC is Windows-only): Clang
      (Linux/Windows/macOS), GCC (Linux + macOS where runners exist), MSVC (Windows) — core, tests,
      wheels; the matrix itself is COMMITTED (MSVC promoted from M3.6 *(opt)*; its slower codegen is a
      measurement, never a gate).
- [ ] M19.2 Docs website expansion: task-first tutorials per interface (C++/Python/MATLAB/Studio), theory
      section linking every `docs/derivations/*.md`, Doxygen API reference integrated, benchmark page
      (D-27-qualified numbers only), "which model when" (M9.5) surfaced, EIS + design-mode tutorials
      (M18).
- [ ] M19.3 Benchmark refresh: PAY-class numbers re-recorded on current code under the D-27 protocol;
      structural counters asserted in CI for the new tiers (M3.5 extended).

### M20 — Hardening pass H2 + v6.0.0 release (D-46; goal 23)

- [ ] M20.1 Hardening pass H2 over M13–M19 (D-46 template): adversarial orthogonal review; bug ledger
      failing-test-first; sanitizer + fuzz lanes over every parser added since M12 (shim inputs, profile
      readers, Schiffer/semi-empirical parameter files); mutation battery; performance-counter
      re-baseline; simplification sweep (digit-identical); MC-1 line-debt sweep INCLUDING test files.
- [ ] M20.2 Version sync + artifacts (wheels CPU/CUDA, MEX, npm wasm package, Studio builds, shim
      package); consolidated CHANGELOG; release notes with the tier matrix, shim coverage matrix, and
      qualified-timing caveats.
- [ ] M20.3 Annotated `v6.0.0` tag/publish — §0.3(5) exception, one-word go from Volkan.

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
| **Q13** | Studio plotting library? | **ASSUMED uPlot, bundled** (small, fast, offline-friendly); switch to bundled-Plotly only if uPlot provably can't serve a view — record in the M11.1 note |
| **Q14** | WASM binding surface: embind or flat C API? | **ASSUMED embind** (fast to build, mirrors §3.9); revisit only if module size busts the M10.4 budget — a C API is the recorded fallback |
| **Q15** | DFN memory budget (PC-6 tier)? | **ASSUMED ≤ 16 KB/cell** (10⁴ DFN cells ≈ 160 MB); tighten in the M7.1 design note from the real state layout |
| **Q16** | Chemistry scope: sodium-ion? lead-acid? | **ASSUMED:** Na-ion YES where PyBaMM 26.6.2.0 ships sets (same SPM/SPMe maths); lead-acid NO for v5. **OVERTURNED for lead-acid 2026-07-14 (Volkan directive): Schiffer tier IN scope at M14 (D-43); Na-ion answer unchanged; still out of v5.0.0 itself** |
| **Q17** | Studio deployment? | **ASSUMED GitHub Pages + downloadable offline single-file HTML artifact per release** |
| **Q18** | Distributed (per-plane) electro-thermal coupling in v5? | **NO — quarantined (D-36)** until the pouch project closes EX780–802; v5 ships lumped-Q + tab projections only |

## 8. Status ledger (compact; the 60-row per-gate history is archived verbatim)

| Date | Item | State |
|------|------|-------|
| 2026-07-07 | Phase 0 + legacy audit + plan v1 + SOTA verification + Q1–Q7 | DONE (archive) |
| 2026-07-08 | P1-G0 Release re-confirm; production `SpectralDiffusion<NCH>` | DONE (archive) |
| 2026-07-09 | §3.12/D-22/D-23 kernel contract; orthogonal review (13 defects pre-implementation); Chebyshev math audit; PyBOP promotion; async-recording design | DONE (archive) |
| 2026-07-10 | Phases 1–7 all gates + PAY-1 (6.97×), PAY-2 (4.87×, 10× falsified), PAY-3 (10⁵ cells/24 MB), PAY-4 (18×/71× PyBaMM, ~1.1–1.2·10³× liionpack, qualified) | COMPLETE (archive rows 2026-07-10) |
| 2026-07-10 | Phase 8 G1 MATLAB, G2 CUDA (+PAY-5 23.85×), G3 async recording, G4 thread pool | COMPLETE (archive) |
| 2026-07-10 | Implementation audit of the 67 Codex commits (agent + Fable, artifact-checked); architecture quality assessment; full ctest Debug 49/49 + Release 49/49 | DONE — verdict now §2.1/archive; AUD-1..5 opened as Phase-9A items |
| 2026-07-10 | PLAN.md compressed + Phases 9/10/11 added (bug-hunt+simplification, PyBOP integration, release); short-simulation operating rule added; P8-G6 → Phase 11 | DONE (this revision; archive created) |
| 2026-07-10 | P8-G0 optionality/portability | PASSED — core-only and nested external-consumer smoke 1/1 on Windows; Debug/Release optional-off CPU 49/49; rebuilt installed CPython 3.13 wheel 9 passed/2 expected skips; private CUDA metadata; no optional SDK header leaks; installed-Eigen-first/pinned fallback; 3-OS core and installed-wheel CI matrices cover CMake changes. Cross-platform jobs are committed but cannot run until pushed. |
| 2026-07-10 | P8-G5 tested v4 documentation | PASSED — exact C++/Python/MATLAB fences run in available toolchains (7 finite samples each); 20 local links/front matter checked; Doxygen 0 generator errors with a rendered main page; production Jekyll build emits 8 themed v4 pages; docs workflow linted and least-privilege. AUD-5 resolved. |
| 2026-07-10 | Phase 9A audit-debt closure | PASSED — AUD-1 enforces 1e-12 P2-G1 V/I bands (worst ΔI 3.70e-13 A); AUD-2 classifies same-commit Phase-5 thresholds as post-hoc regression envelopes with derivation; D-26 waives CVODE narrowly after strengthening the exact oracle; D-27 resolves Q11 while preserving qualified timing labels. Targeted Debug/Release 3/3, final full Debug/Release 49/49; raw PAY-4 JSON preserved. |
| 2026-07-10 | Phase 9B bug ledger | IN PROGRESS — P9-B01..B36 additionally close ThreadPool production integration, deterministic failure selection, const/null callback support, cursor bounds, strict-FP reduction order, explicit thread linkage, batch-alias races, the legacy zero-worker/exception boundary, and CMake's no-module Clang probe portability. D-30 confines concurrency to identity-distinct batch work. **P9-G2 is PASSED:** isolated Debug/Release ASan+UBSan fuzz targets, independent/poisoned atomicity oracles, corpora/dictionaries, limit replays, mutation checks, and the two-job CI workflow are complete; Linux 60-second campaign counts are 149,399/73,981/94,933 Debug and 291,223/546,881/556,946 Release. P9-G1/G3/G4 remain open. Current Debug/Release subsystem gates: Experiment 203, ParameterSet 1,495, netlist CSV 728, parser allocation 913, Recorder 98, AsyncRecorder 238, async allocation 15, PackTopology 66, PackSolver 704, PackStepper 208, P2-G1 allocation 14, ThreadPool 600, ModeC 52, CompiledCurve 163, ThermalLumped 25 assertions. Fresh Linux Clang 18 core-only configure/build/smoke passes 1/1. Ledger and refutations: `.claude/reports/p9b-bug-ledger-2026-07-10.md`. |
| 2026-07-11 | PLAN.md rewritten as the v5 single-goal ladder (Fable): §0 persona + non-stop protocol; PC-9/10, EC, MC contracts; §3.13–3.22 (expressiveness, FVM, SPMe, DFN, thermal 2+1D/3D port, regime advisor, parameter+citation library, WASM, Studio); D-31..D-41; ladder M0–M12 (Phases 9B/9C/10/11 carried into M0–M2); orthogonal critique applied (32 findings: 2 blockers + 10 majors fixed — dependency/artifact rules, §5.7 PAY-blocking restated, M2.3 tag non-blocking, M3.7 gate mandatory, M4.8 EC-5 schema box added, DOI web-verification required, M10.3 covers DFN); reference digests committed under `.claude/references/` | DONE (this revision; archive `.claude/summaries/plan-archive-2026-07-11-v4-phase9b.md`) |
| 2026-07-12 | M0.1 PackSolver allocation-atomic reconfiguration (P9-B37) | PASSED — all fallible candidate construction precedes statically no-throw member publication; deterministic late OOM preserves/reuses the prior solver. Debug/Release focused gates: ParserAllocation 946, PackSolver 704, PackStepper 208, P2-G1 allocation 14 assertions; final full Debug/Release 51/51. Artifacts: P9B ledger B37 and `.claude/summaries/handoff-2026-07-12-m0-1-packsolver-atomicity.md`. |
| 2026-07-12 | M0.2 P9-G1 whole-project sanitizer lanes + P9-B38 | PASSED locally — Clang 18 ASan+UBSan flags and runtime symbols proven on production/test binaries, full suite 51/51; TSan symbols proven and ThreadPool/AsyncRecorder/PackStepper 3/3. AsyncRecorderAllocation is structurally excluded from TSan because its intentional global new/delete overrides collide with TSan interceptors; it passes ASan+UBSan. Native Debug/Release 51/51; core-only/nested consumer 1/1. Hosted workflow committed, not claimed run. Report: `.claude/reports/p9g1-sanitizers-2026-07-12.md`. |
| 2026-07-12 | M0.3 P9-G3 exact Status-failure coverage | PASSED locally — Linux Clang/LLVM 18.1.3 optional-off resolved all 336 active arms from 366 lexical failure arms: 329 measured, 7/10 hash-pinned structural exceptions (5 defensive-only, 2 platform-specific), 30 inactive optional branches, and 0 uncovered/unmapped. A fresh coverage session passed 52/52; native Debug and fast-math Release each passed 52/52; refreshed ASan+UBSan passed all 51 short tests and TSan passed 3/3. The gate binds each per-test profile to exact source/compile/binary identities, exports tests separately, supplies a whole-archive production anchor, and rejects stale or tool-mismatched evidence. Hosted workflow committed, not claimed run. Artifacts: `.claude/reports/p9g3-status-coverage-2026-07-12.{md,json}` and `tests/coverage/status_failure_exceptions.json`. MC-1 debt: exact proof closure intentionally kept component-local regressions and one auditable reporter while site hashes were frozen; oversized files are `scripts/status_failure_coverage.py` (1,281), `core_PackSolver_test.cpp` (1,128), `core_AsyncRecorder_test.cpp` (1,075), `core_Experiment_test.cpp` (861), `AsyncRecorder.cpp` (890), `Recorder.cpp` (721), `SpmFactory.cpp` (712), `ParameterSet.cpp` (1,496), `Experiment.cpp` (1,063), and `PackSolver.cpp` (911). Split `ParameterSet`/Experiment boundaries at M0.7, deduplicate fixtures at M0.8, and resolve remaining line debt at M0.10 without weakening exact-site evidence. |
| 2026-07-12 | M0.4 P9-G4 adversarial bug ledger | PASSED locally — the live ledger closes P9-B01..B50 with a failing SHORT regression or controlled mutation, an implemented fix, post-fix evidence, and a row for every defect; 23 refuted candidates remain recorded with proofs/tests so defensive guards and harness mistakes are not relabelled as product bugs. Final tranche B40..B50 closes schedule non-progress, extent conversion/storage overflow, solution and cycler allocation atomicity, exact LLVM tool selection, fast-math/FTZ-elided conductance validation, sparse resistor publication, duplicate spectral modes, 750,000× compressed-read allocation amplification, public recording allocation translation, and source-step overflow. Final native Debug/Release pass 52/52; Linux Clang 18 Release PackSolver passes 1/1; focused ASan+UBSan passes PackSolver and AsyncRecorderAllocation 2/2; the fresh exact coverage session passes 52/52 and 329 measured + 7/10 pinned exceptions, 0 failed, 30 inactive. TSan was not repeated because no concurrency production code changed; M0.3's authoritative 3/3 remains applicable. Hosted CI is committed, not claimed run. Artifacts: `.claude/reports/p9b-bug-ledger-2026-07-10.md`, `.claude/reports/p9g4-status-coverage-refresh-2026-07-12.{md,json}`. |
| 2026-07-12 | M0.5 / 9C-1 one-source SPM physics (PC-10) | PASSED locally — `SpmScalarKernels.hpp` supplies the production modal, observable, interpolation, and diffusion scalar algebra to CPU, CUDA, Dual, and `SpectralDiffusion`; the deliberately raw legacy diffusion oracle remains independent. Capture-host exact fixtures are explicit opt-in and retain all pre-migration bits across Debug/ThinLTO, Release with and without ThinLTO, and CUDA, including the 48,080-value spectral trace. Portable WSL Clang 18 analytic/FD/parity/structural gates pass without adopting post-refactor hashes. Adversarial review fixed the omitted off-reference Dual entropic tangent, ThinLTO-elided plating-scale overflow failures, and negative spectral lane allocation; controlled mutations prove the structural and independent numerical oracles turn red. Native Debug/Release/CUDA each pass 54/54; full WSL Clang 18 ASan+UBSan passes 54/54 without findings; exact Status coverage closes 339 active arms as 332 measured + 7 structural exceptions with zero uncovered/unmapped; core-only and nested consumers pass 1/1. TSan was not repeated because concurrency code did not change; M0.2's 3/3 remains applicable. Hosted CI is committed, not claimed run. Artifacts: `.claude/reports/p9c1-single-source-validation-2026-07-12.md`, `.claude/reports/p9c1-pre-refactor-scalar-fixtures-2026-07-12.md`, and `.claude/reports/p9c1-status-coverage-refresh-2026-07-12.{md,json}`. |
| 2026-07-12 | M0.6 / 9C-2 one ageing-kernel idiom | PASSED locally — `AgeingKernel.hpp` supplies checked field-major scratch, clearing for accumulated outputs, ordered model/lane traversal, and atomic pipeline-stage dispatch while SEI/crack/LAM/plating retain named scalar-generic physics and the raw legacy oracle remains independent. The pre-refactor 1,077-value all-mask A/B/A trace is exact in Debug/ThinLTO, fast-math Release, and fast-math Release/ThinLTO. Adversarial review fixed saturated-branch Dual tangents, the `pow(Dual,0)` NaN boundary, and signed extent/padding defects in ageing/stress/observable/transport/pipeline/arena scratch; it also added full-pipeline failure atomicity and exact-coverage source-shape boundaries. Mutation checks turn red for private loops, traversal/stage reorder, field aliasing, wrong Dual derivatives, zero-exponent NaNs, and merely-inline fast-math Status loss. Final native Debug/Release/CUDA pass 56/56; WSL Clang 18 ASan+UBSan passes 56/56 without findings; exact coverage closes 339 active arms as 332 measured + 7 pinned exceptions, 30 inactive, zero uncovered/unmapped; core-only and nested consumers pass 1/1 after rebuild/no-op checks. TSan was not rerun and no new TSan claim is made; hosted CI is committed, not claimed run. `SpmPipeline.hpp` is 746 nonblank lines after explicit construction proofs, recorded for M0.10 review under MC-1. Artifacts: `.claude/reports/p9c2-ageing-kernel-validation-2026-07-12.md`, `.claude/reports/p9c2-pre-refactor-ageing-fixture-2026-07-12.md`, and `.claude/reports/p9c2-status-coverage-refresh-2026-07-12.{md,json}`. |
| 2026-07-12 | M0.7 / 9C-3 cold-file split | PASSED locally — the 1,496-line fixture-baseline `ParameterSet.cpp` is separated into core absorption/SPM compilation, BPX mapping/files, stack-local expression AST, strict JSON, and narrow curve/expression/DOM seams; the 1,063-line fixture-baseline `Experiment.cpp` is separated into parser, cohesive Cycler transaction, and shared semantics. All pre-split ParameterSet metadata/value/factory and independent Experiment parser/runner fingerprints remain exact in Debug/ThinLTO, fast-math Release, and Release/ThinLTO; CTest stays 56. Adversarial review caught missing runner source ownership, non-atomic reused JSON (P9-B54), a structural seam false-negative, stale scratch after a valid NCH replacement, and five duplicate allocation-return coverage sites. Mutations turn red for the stale-scratch guard and both shared Status mappers. Final native Debug/Release/CUDA pass 56/56; WSL ASan+UBSan passes 56/56 with zero findings; exact coverage passes 369 lexical / 339 active = 333 measured + 6 exceptions, 30 inactive, zero uncovered/unmapped after retaining the initial 380/350/337+7+6 FAIL; core-only/nested consumers pass 1/1 and all rebuilds are no-op. Archive audit finds 18 members, six split objects once, exactly three intended strong seams, and no exported parser internals. `CyclerV2.cpp` is a reviewed 771-line MC-1 cohesion exception to preserve the frozen non-IPO runner transaction. Artifacts: `.claude/reports/p9c3-pre-refactor-cold-fixtures-2026-07-12.md`, `.claude/reports/p9c3-cold-file-split-validation-2026-07-12.md`, and `.claude/reports/p9c3-status-coverage-refresh-2026-07-12.{md,json}`. |
| 2026-07-13 | M0.8 / 9C-4 shared test harness | PASSED locally — `tests/support/CoreSpmTestHarness.hpp` is the single test-scaffolding idiom (MC-4): caller-owned successful build, unit-tagged (`CurrentA` / `CurrentDensityApm2`) terminal-voltage observation, constant-current `ExponentialModal` traces on a caller-supplied grid including the initial sample and every partial final interval, and a transparent max-abs/RMS voltage metric — with no default options, chemistry, tolerance, reference, file I/O, stored batch/stepper state, or mutable static. Nine binaries migrated; construction/allocation-ordinal/topology/backend subjects keep direct factory calls by the preregistered allowlist. Native Debug, fast-math Release, and CUDA/Release each pass 57/57 (56 baseline + the harness self-test); WSL Clang 18 ASan+UBSan+LSan passes 57/57 with zero findings in 139.14 s. All 43 pre-existing binaries meet or exceed their frozen per-binary assertion floors with no test-case drop (39 identical; Experiment 370→436, ForwardSensitivity 3919→16172, P7G3 6334→7274, Recorder 168→222; CudaSpmBatch unchanged at 433,671/4); floors are re-frozen at these counts. Fifteen registered mutations turn a gate red (hard-coded lanes/nch, flipped discharge sign, omitted A→A/m² conversion, dropped initial sample, dropped final interval, uniform dt, self-comparison, ignored final metric element, RMS/(N−1), naive IEEE-unguarded metric, wrong absolute step time, default lane count, harness inside an allocation-fault window, 50 mV PyBaMM reference perturbation). Adversarial review refuted the by-value-seam, allocation-window, oracle-weakening, and circularity attacks with evidence; its surviving findings were fixed — the structural blacklist now also forbids `EulerLegacy`/`stepper_`/Catch matcher tolerances (each verified RED by injection), `core_ForwardSensitivity_test.cpp` now bitwise-asserts that the accumulated grid re-differences to the exact `dt` the harness re-derives, and the configured-batch no-move invariant (`Recorder`/`CyclerV2` hold a raw `SpmBatch *`) is documented on the seam. **Two preregistered mutations are FALSIFIED and recorded:** the aggregate per-binary floor cannot detect a single removed assertion (harness admission `REQUIRE`s lift counts far above the floors — remedied by re-freezing), and the allocation-window prohibition is enforced by the allocation binaries' ordinals, NOT by the structural gate as preregistered. **M0.8 did not reduce lines** — caller-owned scratch spans cost more per call site than the lambdas they replaced (`core_AsyncRecorder_test.cpp` 1,075→1,118; `core_Experiment_test.cpp` 1,091→1,168; ~+1,190 test lines overall); the MC-1 line debt passes to M0.10, which owns it. TSan was not rerun and no new TSan claim is made; hosted CI is committed, not claimed run. Artifacts: `.claude/reports/p9c4-shared-test-harness-validation-2026-07-13.md` and `.claude/reports/p9c4-pre-refactor-harness-baseline-2026-07-12.md`. |
| 2026-07-13 | M0.9 / 9C-5 public-surface audit | PASSED locally — MC-5 is now a property, not a label. `SpmFactory.hpp` no longer includes `SpmPipeline.hpp`: the parameter blocks and row layouts it needed **by value** moved into `SeiParams`/`SurfaceCrackParams`/`LamParams`/`LithiumPlatingParams`/`SpmBatchLayout`/`AgeingModelMask.hpp`, the factory header names types and `SpmFactory.cpp` names kernels, and `SpmPipelineLayout` is renamed `SpmBatchLayout` (what `SpmBatch::layout()` returns). A public-consumer TU fell from **26 core headers to 21, with all eight ageing/pipeline kernel headers gone** (`clang++ -MM`; registered target 8 → 0 and 20 ± 3). **Corrected overclaim (adversarial review):** the TU is NOT kernel-free — `SpmScalarKernels.hpp` remains reachable via `CompiledCurve`'s inline `linearInterpolate`, because removing it would duplicate physics and break PC-10, which outranks MC-5. Earlier wording said "zero kernels"; that was false and is fixed in CHANGELOG, the report, and the gate's own comment. Every core header carries an `@surface` tag with tier api, support, or internal (16 api / 12 support / rest internal) and `tests/structural/p9c5_public_surface.cmake` enforces R1 classification, R2 the pinned api set, **R3 no api/support header may include an internal one** (which makes the MC-5 property inductive, not merely measured), R4 each api header's pinned declaration count and public types, R5 bindings/docs see api only, R6 the naming lexicon. Naming pass (one verb per concept): `n_lanes()`/`n_rows()` everywhere (`Recorder`, `AsyncRecorder`, `CudaSpmBatch` dropped `nLanes()`/`nRows()`), CUDA shim counters unified to `deviceAllocationCount`/`deviceWideSynchronizationCount`/`deviceArenaBytes`, one `checked_lane_count`/`checked_shape` where three spellings lived; `validate…`→`Status`, `valid…`/`is…`→`bool`, no `get` prefix. **Digit-identical:** all 53 Catch2 binaries report identical assertion and test-case counts to `df2da52` in Debug *and* in a like-for-like fast-math Release baseline build; CTest 57/57 Debug, 57/57 Release, 57/57 CUDA; `CudaSpmBatch` unchanged at 433,671/4. **Twenty registered mutations turn a gate red**, after three rounds of hardening. My own battery found two holes (R4's substring check could not see a renamed public type; the `get`-prefix rule could never fire against whitespace-stripped text). **Adversarial review then broke the gate seven more ways, all now red:** a public method, an enumerator, or a `*Params` field could be added unseen (R4 counted column-0 declarations only, and pinned NO support header — yet the `*Params` structs are the public parameter surface, and a new `SeiParams` field changes the layout of the api type `SpmFactoryInput`); a kernel could be included through a `../` path or an angle bracket (R3 matched quoted, non-`../` includes only); a binding could include `core/detail/…` (R5's regex could not match a nested path). **It also found a real code defect, now fixed:** `SeiParams.hpp` and `SurfaceCrackParams.hpp` used `slide::Status` without including `../types/Status.hpp` and compiled only by luck of include order — a new compiler-verified gate (`structural_test_core_9C5HeaderSelfContained`, CTest 58) now compiles each of the 28 api/support headers standalone. `p9c2` additionally re-pins mask *use* (`for_each_enabled_ageing_model_lane<N>(p.model_mask`), since include-presence alone was weaker than what it defended before. **Recorded, not hidden:** the `*_fast`/`*_ipo` binaries in `build-release/bin` are stale 2026-07-12 artifacts outside the Ninja graph — one fails an assertion and is NOT evidence about current code (the live Release lane is already `-ffast-math` and passes that case); the single Release count that differed from baseline was the `SLIDE_ENABLE_RECORDED_SCALAR_BITS` option, not behaviour — matched like-for-like, the baseline reproduces it exactly, and the two extra assertions are fast-math bit-hash checks that PASS. `SpmScalarKernels.hpp` stays reachable from user code by design (PC-10 beats MC-5 there). The then-recorded CUDA shorthand—MSVC `cl.exe` on `PATH`—was superseded by MQ.1: a fresh lane needs the complete VS x64 developer environment. TSan not rerun, no new TSan claim; hosted CI committed, not claimed run. Artifacts: `.claude/reports/p9c5-public-surface-{preregistration,validation}-2026-07-13.md`, `.claude/designs/m0-9-public-surface.md`. |
| 2026-07-13 | M0.10 / 9C-6 dead-code + line-debt sweep (P9-G5) | PASSED locally — **there was no dead code to remove**, and that is the finding: every namespace-scope entity in `src/core` is referenced (the 17 a naive scan flagged are aggregate members or same-file callees), and `SpectralDiffusionLegacy` stays as the parity oracle. So 9C-6 is a line-debt sweep. The four oversized TUs split along boundaries that proved real: `SpmFactory.cpp` 713 → 494 + **`SpmBatch.cpp` 237** (the batch reaches its implementation only through type-erased pointers, so `SpmBatch.cpp` does not include `SpmPipeline.hpp` at all); `PackSolver.cpp` 910 → 610 + **`PackSolverIterative.cpp` 266** (ladder and relaxation never touch `SolverWorkspace::Impl` — they factorise nothing, which is why they can move); `AsyncRecorder.cpp` 890 → 396 + **`AsyncRecordingCodec.cpp` 452** + **`detail/AsyncRecordingFormat.hpp` 116** (transport vs codec, one shared format); `Recorder.cpp` 721 → 293 + **`RecordingFormat.cpp` 451** (writer and reader stay together — they must agree byte for byte). `SpmPipeline.hpp` 783 → 745 + **`SpmDiffusionRhs.hpp` 64** (every other mechanism's RHS already had its own header). `scripts/status_failure_coverage.py` 1,281 → 135 + the `slide_coverage` package (scanner/session/classify/report/constants) behind the same CLI. **One real duplicate removed:** `checkedAdd`/`checkedMultiply` lived twice with identical semantics → `detail/CheckedArithmetic.hpp`. **Digit-identical:** all 53 binaries identical in Debug and in a like-for-like fast-math Release baseline; CTest **58/58 Debug, Release, and CUDA**; `CudaSpmBatch` unchanged at 433,671/4; the coverage reporter's self-test passes, its census is byte-identical, and `ruff` proves no undefined name crosses the new module seams. **Five mutations, one per new TU, turn a gate red** (three first-attempt mutations were DISCARDED, not counted, because no existing oracle could observe them — a self-consistent CRC, an unexercised overflow, and an `SpmBatch::rhs` path `ExponentialModal` never calls). **Line count GREW 17,363 → 17,595 (+1.3%)** — predicted and registered BEFORE the work, against PLAN's default expectation of a reduction: splitting a TU adds a guard, an include block, and a namespace, and there was no dead code to delete. MC-1's goal is met where it matters: only two `src/core` files exceed 700 lines, both with written justifications — `SpmPipeline.hpp` (745; a template whose composition cannot be split without destroying D-02's zero overhead) and `CyclerV2.cpp` (771; the M0.7 exception). **UNMET, recorded not buried:** the oversized TEST files (`core_ParserAllocation` 1,229, `core_PackSolver` 1,198, `core_Experiment` 1,168, `core_AsyncRecorder` 1,118, `core_ParameterSet` 1,007) are untouched for the second milestone running — M0.8 grew them, M0.10 did not cut them, and they need an explicit owner rather than a third inherited assignment. The Linux llvm-cov lane was NOT re-run; the workflow's path filter now includes `scripts/slide_coverage/**`, but no coverage claim is made. Artifacts: `.claude/reports/p9c6-line-debt-{preregistration,validation}-2026-07-13.md`. |
| 2026-07-14 | PLAN extended to THE superior-stack goal (Fable): Volkan's 16-point directive + integrator-research directive absorbed — goals 15–23; §1.4 validation contract VC-1..VC-5; §3.23 model-family tiers (ECM, semi-empirical/SimSES, lead-acid Schiffer, storage-system), §3.24 `slide.pybamm` drop-in shim, §3.25 integrator & predictive-acceleration research programme; D-43..D-49; ladder extended M13–M19 with v6.0.0 at M19; M1.0 test-file-debt owner box added (closes the M0.10 UNMET carry); M12.0 hardening pass H1 added; Q16 overturned for lead-acid (D-43). Directive items already covered were POINTED at, not duplicated (D-01/D-10/D-16/D-18/D-08/§3.13/§3.8). Independent orthogonal review (Fable critic, same day) returned 12 findings — 1 CRITICAL (this very row had recorded the review as done before it ran), 5 MAJOR (stale header wave numbering; goal-17 GPU promise resting on skippable/missing boxes; VC-2 quantifier unclosed over ECM/Schiffer; M12.3 release-note pointers contradicting D-43/v6; unmarked from-memory VC-1 waiver for Schiffer), 6 MINOR — ALL twelve applied before commit; review artifact: `.claude/reports/plan-revision-orthogonal-review-2026-07-14.md`. | DONE (this revision) |
| 2026-07-14 | Newman instrumentation & design wave added (Fable, on Volkan's "think like John Newman" + "EIS, more optimisation, analytical Jacobians" + "good unit tests" directives): goal 24; §3.26; D-50 (thermodynamic identities as runtime gates), D-51 (EIS by analytic linearisation of the compiled system, optional C_dl rows, never time-domain FFT), D-52 (design mode over ensemble lanes with EXACT forward-sensitivity Jacobians — FD gradients forbidden in shipped optimisers); new milestone **M18** (entropy-production + Bernardi energy-balance gates incl. heat of mixing, voltage-loss decomposition identity, virtual reference electrode + plating-onset event, Newman current-distribution groups → M9 advisor, EIS, Ragone/sizing optimisation with analytically derived benchmark optima, Jacobian service surface, ICA/DVA tool, (opt) PSD/MPM tier, D_s convention doc) with an explicit MC-2/MC-4 unit-test preamble (every box ships oracle-first `tests/unit/core_X_test.cpp`); former M18/M19 renumbered **M19/M20** — v6.0.0 now M20; cross-references updated (goals 20/21 → M19.1/M19.2, goal 23 + D-46 H2 → M20.1, header/§1 → goals 15–24, M13–M20); M7.1 extended to record the BAND(J) lineage + the nondimensionalisation decision. All Newman-school citations marked WEB-VERIFY at their boxes. | DONE (this revision) |
| 2026-07-21 | Implementation & code-quality pass (Volkan-requested, not a ladder box) | DONE (Debug only) — nine-way disjoint fan-out review of `src/core` with an adversarial verifier per finding: **117 raised, 68 survived, 49 refuted**. Landed: the **PC-10 violation** where `Sei`/`Lam`/`LithiumPlating` hand-rolled `(1/T_ref−1/T)/Rg` while nine sites (incl. `CudaSpmRuntime.cu`) already called `spm_scalar::arrheniusFactor` — digit-identical, and provably so for `Dual` because `Dual` declares only `operator/(Dual,Dual)`, so the hand-rolled form was already doing the implicit conversion the cast makes explicit (`SurfaceCrack.hpp` deliberately excluded: different association, would change bits — recorded as remaining PC-10 debt). Also: `SLIDE_ROOT_DIR` was **referenced but never defined anywhere**, so `data/`/`results/` resolved against the CWD and 25 of 53 binaries failed when run directly while CTest passed only because it launches from `<build>/bin` — registered before the fix (25 → 0 failures, counts unchanged) and confirmed; `project_warnings` had **no library consumer**, so `slide_core` and all 53 unit tests compiled with zero `-W` flags (core was already clean bar one `-Wshadow`, so enforcement was free; 22/22 TUs now 0 warnings); Boost joins Eigen/range-v3 as `SYSTEM`; `PathVar` globals `static` → `inline const`. **Two tests that could not fail now can:** the byte-shuffle codec was round-trip-only (green under a self-consistent wrong permutation — demonstrated by reversing plane order in both directions), and `parseValue`'s exponent branch was unexercised though pandas writes `1e-05`. Evidence: 51/53 binaries digit-identical, the two deltas being exactly the added tests, CTest 58/58, every mutation observed red with sources restored. **Recorded honestly:** one 57/58 run and a 6/12 flake loop were **confounded by concurrent builds** and are not evidence; **Release and CUDA were NOT re-run**, so the usual Debug-AND-Release standard is not met by this pass; 60 verified findings remain unapplied and need owners (O(cells²) ladder detection, a line-for-line duplicated topology derivation, `substeps` that does not subdivide `dt`, a thrice-reused relaxation buffer, `crc32` defined twice, a fast-math-unsafe FP gate). Artifacts: `.claude/reports/code-quality-pass-2026-07-21.{md,-survivors.json}`. |
| 2026-07-23 | PLAN extended with the MQ quality wave (Fable, on Volkan's clean/re-derive/hunt directive): milestone **MQ** inserted between M0 and M1 in ladder order — MQ.1 restore the Debug-AND-Release baseline (the 2026-07-21..23 landings were Debug-only), MQ.2 quality-backlog triage owner (all 68 survivors dispositioned: applied/REFUTED/deferred-with-owner; six named landmines may not be silently deferred), MQ.3 repo-hygiene sweep (tracked-tree audit + untracked deletion PROPOSAL, no self-serve deletion), MQ.4 derivation-inventory matrix `docs/derivations/INDEX.md`, MQ.5 independent re-derivation sweep (derive first, then diff against code; discrepancy = failing test first), MQ.6 fresh logic hunt loop-until-dry (gates are code too), MQ.7 structural-only performance hunt, MQ.8 closeout. This closes the previous NEXT row's open question — the quality backlog is pulled forward as MQ.2, not deferred to M12.0/M20.1. Standing operating contract split out to **`AGENTS.md`** at the repo root for autonomous Codex sessions (compression only — PLAN.md remains the single source of truth). | DONE (this revision) |
| 2026-07-23 | MQ.1 three-lane baseline restoration | PASSED locally — fresh tree-local Clang 21.1.8 Debug/ThinLTO, fast-math Release/IPO-off, and nvcc 13.0.48/MSVC 19.50 CUDA Release with host-C++ ThinLTO each discovered and passed 58/58 tests; CUDA executed the frozen 433,671 assertions / 4 cases, and each repo-root direct-CWD SEI witness passed 57 / 3. No source or test changed. The preregistered `cl.exe`-on-`PATH` CUDA recipe was **FALSIFIED** twice (missing manifest tooling, then `LNK1104 kernel32.lib`); a complete x64 `vcvars64.bat` environment is required. The literal before-every-corrective-config source snapshot cadence was not observed, but clean pretest provenance plus a no-op rebuild strongly corroborates that the final tested artifacts match `af898d1`. Known fast-math, assert-boundary, and CUDA test gaps remain MQ.2 inputs; MQ.2 also owns the newly observed deprecated-`-Ofast` and legacy narrowing diagnostics. No timing, sanitizer, coverage, hosted-CI, installed-package, Linux, or macOS claim is made. Artifact: `.claude/reports/mq1-baseline-validation-2026-07-23.md`. | PASSED locally |
| 2026-07-23 | MQ.2 evidence plan registered | The 68-survivor census plus three stable supplemental IDs is frozen before source edits or decisive runs. Seven deeper findings have explicit existing owners (MQ.3/MQ.5/MQ.7); MQ.3's text now expressly owns the repo-wide, content-aware MC-3 contract census rather than an invented implicit deferral. Behavior changes, exact/no-op fixtures, structural counts, mutation restoration, and final three-lane gates are preregistered in `.claude/reports/mq2-quality-backlog-preregistration-2026-07-23.md`. | REGISTERED; MQ.2 remains open |
| 2026-07-23 | MQ.2 factory no-op boundary | `shared-constants-fanout` is APPLIED at `b9c4db1`: one factory owner now propagates the batch constants after each whole-struct mechanism copy, and the NCH=12 padded-arena plus NCH=5 all-ageing hashes remain exact in Debug, fast-math Release, and host-ThinLTO/CUDA builds. The preliminary `homogeneous-init-lane-loop` APPLIED hypothesis is **FALSIFIED after three distinct designs**: live-row fill changed all three Release F0 hashes; first-lane scalar caching changed the negative zero mode by one ULP and propagated into stress history; verbatim lane-zero computation plus row copies preserved NCH=12 but changed the NCH=5 all-ageing hashes from `dc8f92a59dd8d67f / 59085638e06725d3` to `8ce61773119858d1 / 7aaf4294740b7321`. No hash was reblessed; the initialization source is restored and the finding is `DEFERRED — MQ.7` for structural-benefit proof plus an explicit numerical policy. Detailed numbers: `.claude/reports/mq2-quality-backlog-validation-2026-07-23.md`. | FALSIFIED/deferred; MQ.2 remains open |
| 2026-07-23 | MQ.2 direct fused-Euler boundary | `advance-euler-asserts-only` is APPLIED at `e034d23`: exact state/current/output extents, finite context/current/step values, and positive explicit `dt` now reject transactionally before lane-period indexing. The old-source zero-step witness failed 2/8 assertions. The final 21-row invalid matrix plus non-no-op control passes 114/2 in Debug, fast-math Release, and host-ThinLTO/CUDA; the exact factory 1472/8 and all-ageing 384/6 fixtures remain green in all three modes. Oversized-shape, no-op-success, strict-positive-step, and moved-after-`lanePeriod` mutations all turn the gate red or hang before validation; each was explicitly reversed and source hashes restored. Artifact: `.claude/reports/mq2-quality-backlog-validation-2026-07-23.md`. | PASSED locally; MQ.2 remains open |
| 2026-07-23 | MQ.2 O1 observable storage/access | All six O1 findings are APPLIED at `83ab7be`: two scratch owners share one cold checked-extent helper; every observable subspan is pre-guarded and final consumption is exact; transport-cache representation is private without moving scalar physics; modal addressing uses `BatchView::at`; SEI binds its lane index once; duplicate Status includes are removed and structurally rejected; public `TheveninBatchView::lanes()` is now `n_lanes()`. A stale 9C-3 “all top-level headers are public” rule was narrowed only by two exact allowlisted edges, not broadly weakened. Debug/Release/CUDA focused counts are SpmObservables 69/3, SpmElectrical 48/3, SEI 57/3, Ageing 384/6, PackSolver 878/26, and Factory 1472/8; all existing hashes/bands remain exact. The existing recorded fixtures were falsified as a cache-wiring mutation gate (swapped outputs stayed green), so a direct bitwise hit/miss/four-key invalidation oracle was registered; the same mutation then fails 30/48 assertions. Five structural/storage mutations also turn red. Artifact: `.claude/reports/mq2-quality-backlog-validation-2026-07-23.md`. | PASSED locally; MQ.2 remains open |
| 2026-07-23 | MQ.2 A1 ageing scalar ownership | `sei-kinetic-current-duplicated-in-one-function` and supplemental `surfacecrack-arrhenius-association` are APPLIED at `a699da5`/`e34f548`, after oracle-only commit `43c4953`. Five SEI activation expressions and three kinetic expressions now have one scalar owner; a normal function boundary was **FALSIFIED three ways** by the Release/ThinLTO exact gate (same changed ageing/factory hashes each time), while the preregistered caller-expanded fallback retains every existing hash in Debug, fast-math Release, and host-ThinLTO/CUDA. SurfaceCrack model 5 now uses the shared association under an independent nonzero scalar oracle (old/shared relative delta `4.9e-16` in Debug), a temperature-`Dual` centered-FD gate, legacy `1e-12`, and frozen pre-edit bits for models 1–4. Final focused counts are SurfaceCrack 95/4 Debug and 94/4 both Release modes, SEI 57/3, Ageing 384/6, Factory 1472/8; no recorded hash was reblessed. Four ownership/wiring mutations turn red. Artifact: `.claude/reports/mq2-quality-backlog-validation-2026-07-23.md`. | PASSED locally; MQ.2 remains open |
| — | NEXT | Ladder continues at **MQ.2** (disposition all 68 quality survivors plus MQ.1's two registered build-policy findings) → MQ.3 → … → MQ.8, then **M1.0** → M1.1 → … — take the first unticked box, §0.3 protocol. |
