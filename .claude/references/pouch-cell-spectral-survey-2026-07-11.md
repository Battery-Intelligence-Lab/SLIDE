# pouch-cell-spectral — survey digest for the SLIDE M8 port (scout, 2026-07-11)

> Source repo: `C:\D\git\pouch-cell-spectral` (MATLAB, Volkan). Read the named files before designing;
> this digest is a map, not a substitute. Quarantine note: the Howie electro-thermal closure is contested
> upstream (EX780–802) — see PLAN.md D-36.

A MATLAB reduced-order **coupled electro-thermal** spectral solver for a homogenized LFP pouch cell,
benchmarked against COMSOL. Core solver `src/+solver/pouch_cell_phase1.m`; canonical method write-up
`whitepaper/methodology.tex`.

## 1. Mathematical methods

**Chebyshev collocation (through-thickness x).** Multi-domain CGL collocation over 4 stacked 1D subdomains:
neg-CC → neg-electrode → pos-electrode → pos-CC, shared interface nodes, one-sided derivative rows for
conductivity-jump/separator BCs. Files: `src/+spectral/cheby_nodes.m`, `cheby_diffmat.m`,
`cheby_multidomain.m`; `docs/spectral.md`, `methodology.tex` §Step 2. Typical `N_cheby≈14`. Also reused for
the solid-diffusion radial operator (`src/+spectral/solid_diffusion_operator.m`).

**Cosine/Fourier basis (in-plane y,z, electrochemistry).** Neumann cosine modes `cos(nπy/L_y)cos(mπz/L_z)`,
in-plane eigenvalue `μ_nm=(nπ/L_y)²+(mπ/L_z)²`; tab current injection = exact orthogonal projection of tab
indicator functions. File `src/+spectral/cosine_modes.m`; `methodology.tex` §Step 3, §Electrochemical solve.
Governs `∇·i_s=−a_s i_k`, `∇·i_e=+a_s i_k`, `i_k=(i0F/RT)η`, `η=φ_s−φ_e−E`, `∂q/∂t=−a_s i_k/Q_max`
(`methodology.tex` L52-64).

**Robin eigenfunctions (in-plane thermal).** Sturm–Liouville basis `Y_n(ξ)=(Bi/b_n)sin(b_nξ)+cos(b_nξ)`
with roots of `tan(b)=2·Bi·b/(b²−Bi²)`, exactly encoding convective (Robin) edge BCs; thermal decay
`λ_nm=(k_y/ρC_p L_y²)β_n*²+(k_z/ρC_p L_z²)γ_m*²`. Files `src/+spectral/robin_eigenvalues.m`,
`robin_eigenfunctions.m`, `robin_normalization.m`, `robin_tab_projection.m`, `heat2d_eigenfunction.m`.
Bi=0 degenerates to Neumann `nπ`. `docs/spectral.md` L9-15.

**erf / short-time (image) expansions — three distinct uses:**
- *Electrolyte transmission line (diffusion-free Howie)* — `docs/transmission_line_erf.md`. Per-electrode
  SOC/potential reduces exactly to a diffusion equation for local SOC `û''=m(s)²û`,
  `m(s)=√(s/(D_q(1+sτ_ct)))`, effective SOC-diffusivity `D_q=κ(k_U/2)/Q_max≈7.98e-12 m²/s`. Terminal
  impedance `Z_el(s)=√((1+sτ_ct)/(κC_v s))·coth(mL)`. Time domain = image/erfc series (fast small
  θ=D_q t/L²), the **theta-function dual** of the cosine/eigen series (fast large θ). Fix for the
  "reversal-spike + cosine-truncation" residual. Intermediate regime = **Warburg √t** — curvature is SOC
  redistribution, *not* solid diffusion.
- *Tab-heat short-time overlay* — `mittag/theory/tab_heat_erfc.md`. Concentrated tab source ⊗ 2D Gaussian
  heat kernel = closed-form erf-difference response; blended `w(t)θ_erfc+(1−w)θ_eigen` kills Gibbs ringing
  of the truncated Robin series near the tab at pulse onset. Feature-flagged, default off.
- *Spherical-particle diffusion (Jie)* — `docs/spherical_diffusion_erf.md`. Exact Laplace surface
  concentration `ĉ_s=(Ja/Ds)/(qa·coth(qa)−1)`; short-time closed form `ĉ(τ)=e^τ erfc(−√τ)−1`;
  Warburg+curvature series `(2/√π)√τ+τ+…`. Verified EX671 (max rel err 3.2e-6, τ≤0.1).

**Mittag-Leffler folder** (`mittag/`) — learning sandbox, NOT production (`mittag/README.md`,
`mittag/theory/five_objects_table.md`). Separates five conflated objects: ML theorem (pole-sum),
ML function E_{α,β} (fractional — quarantined to `99_aside_fractional/`, "LFP has no published sub-diffusion
in this regime"), erfc-image series, residue grouping (Smith-Rahn-Wang), Padé. Conclusion: production
inversion is **AAA rational approximation** (Nakatsukasa–Sète–Trefethen) — sample `Ĥ(s)`, partial-fraction
`Σ r_j/(s−s_j)`, invert as `Σ r_j e^{s_j t}`. Benchmark `mittag/matlab_1d/benchmark.m` compares six methods
(fd, cheb, eigen_series, erfc_images, aaa, hyperbolic_talbot).

## 2. Models

Coupled electrochemical + thermal. Full heat source `Q=σ|∇φ_s|²+κ|∇φ_e|²+a_s i_k η + a_s i_k T ΔS_rxn`
(`methodology.tex` L68). Phase-1 omits solid-particle diffusion (Howie); SOC is a local 0-D integrator per
(y,z) point. Segregated: T can be frozen while solving φ_s/φ_e/q (`opts.fixed_T`).

- `howie-chu-model/` = reference COMSOL+MATLAB from Chu et al., *"Parameterization of Prismatic LFP Cells
  through a Streamlined Thermal/Electrochemical Model"* (J. Power Sources 2020) — diffusion-free "Howie"
  (20 Ah LFP). Binary COMSOL model, `Main_Opt.m`/`ObjFunc.m`, reference data `SRun4c.mat`/`SRun2c.mat`.
- `jie-lin-model/` = Lin et al., *"Multiscale coupling of surface temperature with solid diffusion in large
  Li-ion pouch cells"* (Nature Comms Eng 2022) — Howie + spherical solid diffusion. `Battery_model.m`,
  `Exp.mat`, `SurfaceT(t=2500s).txt`.

## 3. Maturity / validated results

**Whitepaper headline (canonical claims, `whitepaper/sections/results.tex`):** Iter-23 SOC-robust closure,
canonical COMSOL params, N_y=N_z=10, 2499 s: 6/6 cases PASS, `|T_err|≤0.55 K`. 4C: T_err −0.50 K, V_RMSE
6.9–8.5 mV; 8C: T_err −0.03…−0.14 K, V_RMSE 8.3–11.4 mV (L163-168). The earlier +1.47 K 8C "ceiling" was
**retracted** (V-only pulse-fit artefact, `limitations.tex` L210-216). NULL results retained (all ≤0.015 K):
per-mode κ(T), per-mode γ(T), live-i_rxn α-blend; full pointwise pseudospectral "Option B" (324 BVPs/RHS)
INFEASIBLE under ode15s+numerical-Jacobian. Accepted bar (`ROADMAP.md` L9-10): RMSE < 20 mV / 0.2 K;
"COMSOL-grade" = 22.2 mV / 0.128 K (COMSOL-vs-data residual, EX662).

**Current sprint — CAUTION (quarantine):** Howie voltage/topology REOPENED. EX783: reduced Howie terminal
polarity reversed vs binary; mirror operator (EX787) passes 100 s gates but leaves +0.875 K uniform bias and
57% transfer-NRMS in the CC-interface current field. EX800–802: missing physics = local-T electrochemistry
needing off-diagonal cosine coupling — kinetics γ(T)·η owns transfer saturation (projection 0.9633), κ(T)
owns current/q geometry. Scalar T_avg + diagonal modal-κ structurally incomplete. EX741 T-closure (0.140 K)
is DATA-calibrated R_extra=59.5 µΩ, not first-principles. Production ku/Qmax labelled a "compensating
compatibility fit". **Jie/diffusion side:** erf toolkit verified (EX671) but not yet a validated coupled model.

## 4. Port-worthy algorithmic ideas

- Tensor-product separable bases with per-direction 1D Sturm–Liouville solvers — `robin_*` are
  direction-agnostic; 2D→3D = call a third time (`docs/3d_thermal_design.md` L32-33).
- Per-mode decoupling of the electrochemical BVP: each (n,m) is an independent linear elliptic 2-point BVP
  in x on the 4-domain Chebyshev grid; block Gauss–Seidel (φ_s±↔φ_e) with row scaling (`methodology.tex`
  §Step D). Modes embarrassingly parallel.
- Analytic equilibrium bypass for the singular (0,0)/zero-current mode: η=0 directly (`methodology.tex`
  L240-241).
- Duhamel/exponential thermal stepping: `dc_nm/dt=−λ_nm c_nm+D_nm`, exact
  `c=c0 e^{−λt}+(D/λ)(1−e^{−λt})` (`docs/spectral.md` L53) — congruent with SLIDE's D-07 propagator.
- Heat assembly: reconstruct φ,∇φ on x-nodes × (y,z), form Q_local, trapezoid-integrate through x → Q_2D,
  project onto Robin basis (`src/+physics/heat_source_assembly.m`).
- The OPEN off-diagonal cosine coupling (triad selection `m''∈{|m−m'|,m+m'}`) for κ(T(y,z))/γ(T(y,z)):
  `docs/howie_modal_temperature_coupling.md`, `docs/kappa_per_mode_derivation.md` — the contested frontier.
- Surrogate/bridge layer (`docs/operator_eigenstructure.md`): small state-space `z=[q_bar,x_j,T,xd_l]`,
  poles `−1/τ_j`, `−G/C`, `−Ds μ_l²/Rp²`; DMD/DMDc fitting `src/+solver/dmdc_fit.m`,
  `modal_match_surrogate.m`.
- Finite-pole spherical kernel with DC-preserving renormalisation `w_k=(1/μ_k²)/Σ(1/μ_j²)`
  (`docs/sphere_multipole_diffusion.md`, `src/+spectral/sphere_multipole_kernel.m`).

## 5. Biot regimes / model hierarchy

- **Bi≪1 through-thickness** is the load-bearing simplification: T uniform in x ⇒ thermal PDE is 2D
  in-plane; x-source integrated out before Robin projection (`ROADMAP.md` L124/L171, `methodology.tex` L76;
  full 3D thermal parked — no x-gradient observed).
- In-plane Biot enters as the Robin eigenvalue parameter `Bi_y=hL_y/k_y`, `Bi_z=hL_z/k_z`; Bi→0 Neumann.
- **σ≫κ regime** (σ_el/κ≈4500): φ_s uniform through electrode to O(κ/σ)≈2e-4 — enables the
  transmission-line reduction (`docs/transmission_line_erf.md` §1, A1).
- Hierarchy tiers: (i) full pointwise pseudospectral (infeasible in MATLAB); (ii) diagonal modal Galerkin
  (production, COMSOL-grade); (iii) transmission-line/erf closed forms; (iv) lumped state-space surrogate.
  Layer-reduction identity: per-layer current I/N_layers with single-layer thickness ≡ N_layers²-scaled
  homogenized conductivity under x_h=N_layers·x_s (`methodology.tex` L339-340).
- Π-group regime map is *specified but unstarted* upstream — `ROADMAP.md` item 9, `GOAL.md` Step 4.1,
  intended `docs/nondimensional_analysis.md`. SLIDE M9 realises this intent.

## Pre-port reading list

`docs/architecture.md`, `docs/solver.md`, `docs/physics.md`, `docs/spectral_construction_tutorial.md`,
`docs/model_reduction_techniques.md`, `docs/analytical_rom_design.md`,
`docs/model_equation_audit_2026-07-10.md` (the reopen's evidence base).
