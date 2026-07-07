# SLIDE v4 — SOTA verification of PLAN.md §3.4/§3.4.1/§3.5/§3.8 solver & perf claims

> Adversarial external check, 2026-07-07. Verdict first line per claim, then evidence with links.
> Method: WebSearch/WebFetch only; no builds, no repo edits. Tags: [confirmed] = named source; [inferred] = reasoned, names what would confirm.

---

## C1 — Exponential (matrix-exponential per-mode) propagator for spectral linear diffusion

**Verdict: CONFIRMED (nuance: applies to the *linear* diffusion sub-operator; the full nonlinear cell system still wants a splitting, which the PLAN provides via Strang multirate).**

- For the Chebyshev-discretised solid-phase diffusion the modes are decoupled and **linear**; over a piecewise-constant-flux (CC) substep the per-mode matrix exponential `z_k(t+h)=e^{Dλ_k h}z_k + ((e^{Dλ_k h}−1)/(Dλ_k))B_k j` is the **exact** solution. Exact ⇒ unconditionally stable, and no scheme can beat exact in accuracy-per-cost for that subproblem. Deletes the O(nch⁴) Euler cliff by construction. [confirmed]
- This is **established prior art in battery SPM codes**, not novel: Bizeray, Duncan & Howey evaluate the spectral solid-phase mass-transfer states with "a 1st order exponential integrator approach at each discrete time point" — exactly the PLAN's scheme. Howey's `Spectral_li-ion_SPM` (which the PLAN's P1-G3 already cross-checks) is the MATLAB reference. Bizeray et al., *arXiv:1506.08689* (J. Power Sources 2015); repo github.com/davidhowey/Spectral_li-ion_SPM. [confirmed]
- Is anything "better" (e.g. PyBaMM's adaptive BDF via SUNDIALS IDAKLU)? For the SPM use case (small fixed nch, current steps at *known* times): **no.** Adaptive BDF/IDAKLU targets the full nonlinear coupled DAE and pays nonlinear-Newton + step-adaptivity + error-estimation overhead that a linear, modally-diagonal, event-aligned subproblem does not need. IDAKLU is the correct tool for SPMe/DFN electrolyte coupling (nonlinear, stiff, algebraic) — not for the isolated linear diffusion piece. The PLAN's split (exponential for diffusion, multirate for the slow nonlinear thermal/aging) is the right factorisation. [inferred — would be confirmed by a like-for-like accuracy/cost bench, which PLAN P3-G1 registers]
- Only theoretically "more general" alternative is a Krylov/Leja matrix-free exponential integrator (e.g. *arXiv:2108.13622*) — irrelevant here because for tiny diagonalised modes the closed-form scalar `exp` is optimal; Krylov pays for problems where you cannot diagonalise. No improvement available.

## C2 — Sparse MNA + Newton with per-cell Thevenin linearization (liionpack architecture) for packs

**Verdict: CONFIRMED as the right reference architecture; NUANCED — liionpack is now maintenance-mode, and a 2025 analytical result improves the parallel sub-case (see Better-than-planned).**

- liionpack architecture is exactly as the PLAN states: circuit of current sources / voltage sources (OCV) / resistors, solved by **modified nodal analysis**, staggered against per-cell PyBaMM stepping. JOSS **10.21105/joss.04051** (2022). [confirmed]
- **Status caveat:** liionpack is explicitly "in maintenance mode … no longer actively developed"; last release **v0.4.0, 2025-01-20**; only critical-bug/PyBaMM-compat updates. Its published examples are small (≈16p2s); it does not itself demonstrate 10⁵. So "still SOTA" = it remains the accepted open-source *architecture*, but not an actively-scaling tool. github.com/pybamm-team/liionpack. [confirmed]
- Nothing **demonstrably faster or more scalable for arbitrary topology** was found. Julia ecosystem (JuBat — FEM P2D/SPM/SPMe, agrees with PyBaMM, *SoftwareX 2024*; LiiBRA.jl — reduced-order realisation, *arXiv:2203.17105*; JuliaSim/Dyad Batteries) targets single-cell/ROM, not pack current-distribution at scale. Commercial (Simscape Battery, Ansys) publish no scalability numbers. [confirmed — absence of public claim]
- The PLAN's three-mode design (add Mode B ladder + Mode C relaxation on top of Mode-A MNA) is **more** ambitious than liionpack, which only has the MNA path. Verdict: the mirror is correct and current.

## C3 — Chord/Shamanskii Newton (frozen Jacobian, refresh on ρ>0.5 or k>4)

**Verdict: CONFIRMED — standard; and Broyden+Sherman–Morrison would NOT beat it for ladder MNA (it destroys sparsity).**

- Shamanskii = Newton with Jacobian refreshed every m iterations (m=1 Newton, m=∞ chord). Textbook: Kelley, *Iterative Methods for Linear and Nonlinear Equations*, SIAM 1995, ch. 5. The contraction-monitored refresh rule (linear rate ∝ ‖J_cached−J_true‖) is exactly Kelley's chord-convergence theory. [confirmed]
- **Broyden rank-1 with Sherman–Morrison for a tridiagonal ladder MNA is worse, not better.** SM updates the *dense inverse* → fills in, breaking the O(n) band structure. For a tridiagonal/banded system a direct refactorisation is already O(n) (Thomas), so the chord's *occasional* refactor is cheaper than carrying a densifying secant. Sparsity-preserving quasi-Newton (Schubert / sparse-Broyden, *SIAM J. Sci. Comput.* 0911036) exists but requires a sparse linear solve per iteration anyway — no win over just re-running an O(n) elimination. Broyden's own SM route is documented as not sparsity-preserving (Broyden's method, Wikipedia; Schubert-update literature). [confirmed]
- Corollary: PLAN Mode B (Thomas) rightly skips `fact` entirely; the workspace's Tier-0 (linear ⇒ factor once) is exact-Newton-in-one-step. Design is sound.

## C4 — Waveform relaxation + Baumgarte for 10⁴–10⁵ cells

**Verdict: NUANCED — WR is validated for battery packs and Miekkala–Nevanlinna is the correct theory; the WR+Baumgarte *combination* is a sound but (as far as found) unpublished synthesis; a simpler exact route exists for the parallel sub-case.**

- WR for battery packs is real and cited correctly in spirit: Zhang/… "A computationally efficient implementation of a battery pack electrochemical model using **waveform relaxation**," *J. Energy Storage* 2022 (S2352152X21014304) — WR decomposes the pack DAE into per-cell subproblems solved iteratively, "significant reduction only when a sufficiently large number of cells are connected in parallel." Convergence theory: Miekkala & Nevanlinna, *SIAM J. Sci. Stat. Comput.* 1987, **10.1137/0908046**. Baumgarte constraint stabilisation (ġ+2αg=0) is standard DAE index reduction (Baumgarte 1972). [confirmed]
- **No source found doing WR+Baumgarte *together* for battery packs** — it is a novel-but-reasonable combination of established pieces. Fine as a research contribution; flag it as unproven-in-combo, arbiter (Mode A) already registered in PLAN. [confirmed absence]
- Newer massively-parallel alternatives checked: parareal / parareal-Schwarz-WR (space-time parallel, *ResearchGate 287233530*) — adds complexity and parareal gives limited speedup on dissipative long-horizon problems; not obviously better here. GPU-batched Newton — see C6. None strictly dominates for arbitrary topology. **But** for pure-parallel ladders the analytical-ODE reformulation below removes the algebraic constraint outright (no WR, no Baumgarte). See Better-than-planned.

## C5 — Eigen SparseLU vs SuiteSparse KLU for repeated small–medium MNA solves

**Verdict: CONFIRMED — KLU is the right optional upgrade; note CKTSO/NICSLU as a future parallel step beyond KLU; dense-blocked for n<100 is correct.**

- KLU (Davis & Palamadai Natarajan, ACM TOMS **Algorithm 907**, 2010) is purpose-built for SPICE/MNA matrices with **fixed sparsity and repeated solves**: analyze/order once, refactor-without-pivoting many times — exactly the PLAN's "symbolic factorisation cached at compile(), numeric refactor each Newton iteration." Eigen ships first-class KLU bindings (KLUSupport module), so the optional-dep swap is a one-line solver change. [confirmed]
- Better-than-KLU exists but is heavier and only pays at large n: **NICSLU** (2.08–8.57× over KLU, 1–12 threads, *IEEE TCAD 2013*) and its successor **CKTSO** (2024, *arXiv:2411.14082 / IEEE TCAD 10.1109/TCAD.2024.3506215*) — faster than KLU/NICSLU/MKL-PARDISO and two GPU solvers on 56 circuit matrices up to 5.5M. For the PLAN's regime (n=10–10⁵, nnz≈3n) matrices are tiny/ultra-sparse; KLU is ample and dense-blocked is genuinely faster for n<100 (BLAS-3, no symbolic overhead). Verdict: KLU correct now; CKTSO a documented upgrade path if pack-solve ever dominates. [confirmed]

## C6 — GPU one-cell-per-thread batched ODE stepping (DiffEqGPU claim)

**Verdict: CONFIRMED that the DiffEqGPU 20–100× claim is real and applies to the ~30-state SPM shape; NUANCED — no *published battery-pack* GPU SPM simulation doing this was found (it would be novel).**

- DiffEqGPU / EnsembleGPUKernel (Utkarsh et al., *arXiv:2304.06835*, 2023): vendor-agnostic, one small ODE system per GPU thread, "20–100× faster than the vmap approach in JAX and PyTorch," matching hand-written CUDA-C. A ~30-state SPM cell is exactly the "many small independent stiff/nonstiff ODEs" target — applicable. [confirmed]
- **No published battery-pack GPU precedent** of this exact shape (batched SPM cells, thousands of cells) surfaced; closest hits are unrelated PIC/Monte-Carlo GPU work. So the PLAN's Phase-8 GPU path is well-motivated but would be a **new demonstration**, not a reproduction. The SoA arena (§3.1) is the correct device layout; host-side MNA/relaxation coupling between batch steps (as the PLAN states) is the standard batched-ODE pattern. [confirmed absence of battery precedent — would be confirmed by a targeted lit check at Phase 8]

## C7 — PyBaMM CURRENT API surface (mid-2026)

**Verdict: CONFIRMED, with the diffusivity rename VERIFIED (v25.6.0) — the PLAN's freeze-time caution is valid.**

- **Version:** current stable **26.6.2.0 (2026-06-16)**; scheme is **CalVer YY.M.patch**. github.com/pybamm-team/PyBaMM/releases; pybamm.org/changelog. [confirmed]
- **(a) Experiment grammar:** "(Dis)charge at x A|C|W", "Rest", "Hold at x V"; duration `for <n> seconds|minutes|hours`; termination `until <x> V|A|C-rate`; combine with `or` (whichever first); cycles via tuples/lists of strings. Extended: **custom steps** (v25.4.0), **custom terminations** (v25.10.0), `start_time` (v23.9). Tutorial-5 docs. [confirmed]
- **(b) ParameterValues:** `ParameterValues("Chen2020")`; `["key"]` access; `search()`; `update()` is now **update-insert** — `check_already_exists` **deprecated in v25.12.0** (2026-01); **default constants no longer added on construction** (v25.12.0, mildly breaking). [confirmed]
- **(c) Solution:** dict access `Solution["Terminal voltage [V]"]`, `.entries`, call-interpolation `sol["…"](t)`, `save_data(..., to_format="csv"|"matlab")`; added `.yp` (time-derivatives) and `.observe(symbol)` in v25.12.0. [confirmed]
- **(d) options dict keys:** `"SEI"`, `"SEI porosity change"`, `"SEI on cracks"`, `"lithium plating"`, `"lithium plating porosity change"`, `"particle mechanics"`, `"loss of active material"`. Values e.g. SEI `"solvent-diffusion limited"`; plating `"partially reversible"`; mechanics `("swelling and cracking","swelling only")`; LAM `"stress-driven"|"reaction-driven"`. Per-electrode via **2-tuple** (v23.9). Constraint: `SEI on cracks=true` ⇒ `particle mechanics="swelling and cracking"`. Maps cleanly onto PLAN §3.2 registry. [confirmed]
- **(e) key rename — VERIFIED:** `"Negative electrode diffusivity [m2.s-1]"` → `"Negative particle diffusivity [m2.s-1]"`, **non-breaking with a deprecation warning, v25.6.0 (2025)** (PR #3624; present v23.5→gone by v25.8.0). Other recent renames affecting an absorption table: `"Exchange-current density for lithium plating [A.m-2]"` → `"…for lithium metal electrode [A.m-2]"` (v25.4.0); `"1 + dlnf/dlnc"` → `"Thermodynamic factor"` (v23.3). **Action:** the Chen2020 absorption table must key on the NEW names and carry deprecation-aware aliases for the old ones. [confirmed]

---

## Better-than-planned alternatives found

1. **Analytical ODE reformulation for parallel packs (arXiv:2508.14454, 2025).** Gives *closed-form* current distribution for parallel-connected cells with interconnection resistances, converting the pack **DAE → ODE** — no MNA solve, no WR iteration, **no Baumgarte constraint** at all. ~44% faster than direct DAE at n=135, gains grow with pack size. This is a strictly simpler route than PLAN Mode C for the **pure-parallel** case, and a time-continuous exactification of Mode B. **Not a general replacement** (parallel-only; linear-in-current ohmic drop; needs known R_k; heterogeneity capped at ~1e-4 σ for n>100). Recommendation: cite it in §3.4 and consider adopting it as the Mode-B/parallel fast path; keep WR+Baumgarte (Mode C) for series-parallel/arbitrary topology where this does not apply.
2. **CKTSO (2024)** as a documented post-KLU parallel sparse-LU upgrade (C5) — only if pack-solve becomes the bottleneck at large n; heavier dep, so behind the same optional flag as KLU.

Nothing found overturns the PLAN's architecture (D-05/D-06/D-07/D-18/D-19). The exponential propagator, chord/Shamanskii memory, three-mode pack solver, and PyBaMM-compat surface are all consistent with the state of the art.

---

## Description-layer runtime-cost paragraph (Volkan's speed worry)

**The two new description layers have strong precedent and are the accepted pattern; compile-time/`build()`-time erasure carries no runtime residue *provided the lowering is total*.** The value-type combinator pack tree (§3.4.2) and designated-initializer parameter packs canonicalised at `build()` (§3.11/D-16) are the same "friendly authoring graph, compiled to flat fast code" split used throughout high-performance scientific software: **CasADi** builds a symbolic graph once then compiles/JITs to branch-free functions (the PLAN's stated model); **PyBaMM itself** discretises its symbolic expression tree and lowers to CasADi/JAX/C so the Python tree never runs in the hot loop; **FEniCS/Firedrake** compile UFL forms to C kernels; **Eigen expression templates** collapse a combinator tree to a single fused loop with zero temporaries at compile time; **Halide/TVM** separate algorithm from schedule. In every case the description objects are *gone* before the inner loop. The **only** documented failure mode — where these layers leak runtime cost in comparable libraries — is **incomplete erasure**: retaining a `std::function`/virtual callable or an *interpreted* expression tree in the kernel (opaque indirection, no inlining, no vectorisation). PLAN D-16 pre-empts exactly this by tabulating every injected function to {scalar, SoA row, uniform LUT, separable LUT, inline Arrhenius}, and PC-1/PC-5 assert zero per-step allocation and concrete-types-only kernels. So the pattern is safe *by the invariants already registered*; the one thing to police in review is that no `std::string`/`std::vector`/`std::function` from the description layer survives into `CellBatch` params or `BatchView` — designated-initializer aggregates and combinator trees are pure compile-/build-time and cost nothing at run time. Verdict: **no hidden-cost precedent when erasure is total; this IS the accepted zero-overhead-abstraction pattern.**
