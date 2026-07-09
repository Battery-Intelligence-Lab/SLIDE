# Chebyshev spectral discretisation — mathematical audit (2026-07-09, Fable)

Scope: `src/cells/Cell_SPM/Model_SPM.hpp` (constructor, lines 65–201; `cumsummat`, 215–261) on branch `Claude`.
Purpose: (1) full derivation trail of what the code implements; (2) defects/fragilities; (3) improvements for the
v4 per-batch model build (§3.2 of PLAN.md). Everything tagged [confirmed] cites a line or a derivation below;
[inferred] names what would confirm it.

## 1. The PDE and the two transforms

Solid diffusion in a spherical particle, Fickian, isotropic:

    ∂c/∂t = D · (1/r²) ∂/∂r ( r² ∂c/∂r ),   r ∈ (0, R)                       (1)

Units: [D]=m²/s, [c]=mol/m³ → both sides mol m⁻³ s⁻¹. ✓
BCs: symmetry ∂c/∂r|₀ = 0; surface flux  D ∂c/∂r|_R = −j  with j the molar flux [mol m⁻² s⁻¹]
(sign: j>0 = delithiation for the code's convention; consumed via `molarFlux` in `Cell_SPM_dstate.cpp:47-53`).

**Transform 1 (kill the 1/r² singularity):** u(r,t) = r·c(r,t). Substituting,

    ∂u/∂t = D ∂²u/∂r²,   u(0,t) = 0                                          (2)

— the plain 1-D heat equation. The surface BC becomes, using c_r = (r·u_r − u)/r²:

    u_r(R) − u(R)/R = −jR/D                                                   (3)

**Transform 2 (nondimensional radius):** x = r/R ∈ [0,1]; ∂²/∂r² = (1/R²)∂²/∂x²; BC (3): u_x(1) − u(1) = −jR²/D.
The 1/R² is folded into the discrete operator at `Model_SPM.hpp:128,138` (`A_/(R²)`); the runtime multiplies by
D(T) (`Cell_SPM_dstate.cpp:47-53`), so eigenvalue units are [D·A] = m²/s · 1/m² = 1/s. ✓ (matches PLAN §3.5.)

**Parity:** c is even in r (symmetry), so u = r·c is ODD. This is what the "folding" exploits.

## 2. What the constructor builds, step by step [confirmed by inspection]

1. **Nodes** (line 70): xm(i) = sin((N−i)π/(2N)) ≡ cos(iπ/(2N)), i = 0…2N, N = nch+1 — Chebyshev–Gauss–Lobatto
   points on [−1,1] in Trefethen's symmetric sin form. `xch` = the nch positive interior nodes (excl. surface x=1
   and centre x=0).
2. **First differentiation matrix** (77–94): standard CGL weights (c_i/c_j)(−1)^{i+j}/(x_i−x_j), diagonal by
   negative row sum (exactness on constants — Trefethen `cheb.m`). Note: node values use the sin form but the
   differences DX use cos(iθ)−cos(jθ) (identical analytically; the sin form's symmetric-rounding benefit is
   partially forfeited — cosmetic).
3. **Second derivative** (98–115): Welfert/Weideman–Reddy in-place recursion
   D₂(i,j) = 2·(C_ij·D₁(i,i) − D₁(i,j))/(x_i−x_j), diagonal by row sum — standard `poldif.m` form.
4. **Odd folding** (117–119): for odd u, u(x_{2N−k}) = −u(x_k) and u(x_N)=u(0)=0, so the 2N+1-point operator
   restricts to the N values at x_0..x_{N−1} (surface + interior):
   DN2 = DM2[:, 0:N] − DM2[:, N+1:2N+1] reversed; same for the surface-derivative row DN1. The centre column
   drops (u=0 there). This halves the operator and enforces u(0)=0 EXACTLY — this is the right construction.
5. **BC elimination** (121–126): the Robin BC (3) is solved for the surface value:
   u_x(1) = u(1) + φ  with φ ∝ flux ⇒  Σ_j DN1(0,j)u_j = u_0 + φ ⇒ u_0 = (Σ_{j≥1} DN1(0,j)u_j + φ)/(1 − DN1(0,0)).
   Substituting into the interior rows of DN2 gives the standard eliminated system
   A_ = DN2[1:,1:] + DN2[1:,0]·DN1[0,1:]/temp,  B_ = DN2[1:,0]/temp,  C_ = DN1[0,1:]/temp,  D_ = 1/temp
   with temp = 1 − DN1(0,0). This is a (static condensation) Schur complement of the surface node. [confirmed]
6. **Scaling** (128–146): A → A_/R² per electrode; output maps: c_surf row = C_/R (c = u/r at r=R... note
   c_surf = u_surf/R), interior rows diag(1/r_i); feedthrough D(0) = R·D_ (surface value's direct flux term).
7. **Modal decomposition** (148–177): A_ is NONSYMMETRIC (folding + condensation destroy symmetry), so
   `Eigen::EigenSolver` (general) is used; z = V⁻¹u, A = Re(Λ), B = V⁻¹B_, C = C_·V; V then overwritten by V⁻¹
   (`.eval()` fix, §2.2(b)). Imag parts warn-to-stderr if ratio > 1e-10 (150–154) but execution continues on
   `.real()` truncation.
8. **Zero mode** (179–196): min |λ|/max|λ| per electrode forced to exactly 0 (the mass mode — see §3); indices
   asserted equal across electrodes (assert vanishes in Release — known gap, PLAN §2.2).
9. **Centre node** (198–199): c(0) = lim u/r = u_r(0). The code computes c_centre = cc_coeff·(Cc·c + flux·R/D)
   with Cc = EVEN fold of the SURFACE derivative row and cc_coeff = −1/DM1(N) = −0.5·(−1)^N (the §2.2(a) sign
   fix). **This path is NOT re-derived in this audit** — it was root-caused and fixed against MATLAB in a prior
   session [PLAN §2.2], but the identity (surface-row even fold ↦ centre derivative) is exactly the kind of
   clever step an external oracle must cover: P1-G3 MUST include the centre-node value (the Carslaw & Jaeger
   series gives c(0,t) as well as c(R,t)).
10. **`cumsummat`** (215–261): Clenshaw–Curtis cumulative-integration matrix (values → Chebyshev coefficients →
    integrated coefficients → values; Chebfun `cumsummat` idiom); row 0 zeroed (integral from the left end).
    Consumer: stress integrals (`getDaiStress`). Cold-ish path; not audited further here.

## 3. The free analytic oracle the tests never used: eigenvalues are roots of tan μ = μ

The homogeneous problem behind A_ is (2) with u(0)=0 and the ZERO-FLUX surface BC u_x(1) = u(1).
Try u = sin(μx): u_xx = −μ²u, so the eigenvalue is λ = −μ² (dimensionless; −μ²D/R² physically). The BC:

    μ cos μ = sin μ   ⇔   tan μ = μ                                           (4)

Roots: μ₀ = 0 (u ∝ x ⇒ c = u/x = const — the MASS mode: this is why one zero eigenvalue exists and why forcing
it to exactly 0 is correct), μ₁ ≈ 4.493409457909064, μ₂ ≈ 7.725251836937707, μ₃ ≈ 10.904121659428899, …
(Newton on (4) converges in ~5 iterations from μ_k ≈ (k+½)π).

**Therefore:** the sorted nonzero eigenvalues of the discrete A_·R² must converge SPECTRALLY (exponentially in
nch) to −μ_k². The low modes should match to near machine precision at nch = 8–12. This validates the operator
(steps 1–8) directly, independently of any trajectory, any initial condition, and any output map — different
mathematics from both the parity harness and the Carslaw & Jaeger trajectory oracle. It also kills the
`findZeroEigenvalue` risk class outright (a spurious second near-zero eigenvalue can never hide).

Registered form (goes into P1-G3, band written before the run):
  for nch ∈ {5, 8, 12}: rel err |λ_k·R² + μ_k²|/μ_k² < 1e-10 for k = 1..⌈nch/2⌉; the zero mode exact by
  construction; ALL eigenvalues real to ratio < 1e-12 and strictly negative (except the forced zero).
  [Band ASSUMED from spectral-convergence experience; tighten after first run — record the actual drift.]

Trajectory oracle (P1-G3 as already planned) — constant surface flux, Carslaw & Jaeger §9.5 / Crank §6.3:
  c(x,τ)−c₀ = (jR/D)·[ 3τ + ½(5x²−3)/5 − (2/x)·Σ_k sin(μ_k x)/(μ_k² sin μ_k) · e^{−μ_k²τ} ],  τ = Dt/R²
  (dimensionless form; the series eigenvalues are the SAME μ_k — the two oracles share (4), which is fine
  because they test different objects: (4)-roots test A; the series tests A,B,C,D together, including x→0.)

## 4. Findings (ranked)

| # | Finding | Where | Class |
|---|---------|-------|-------|
| C1 | `EigenSolver::info()` never checked; failed convergence would be consumed silently | 148, 162 | v4-hardening |
| C2 | Imag-part contamination WARNS to stderr and continues on `.real()` truncation; v4 `build()` must return Status failure (cold path can afford it) | 150–154, 164–168 | v4-hardening |
| C3 | `assert(zero_pos == zero_neg)` and the zero-index selection are Release-silent; superseded by the §3 eigenvalue check at build | 189–196 | v4-hardening |
| C4 | Centre-node identity (step 9) unverified by any external oracle; P1-G3 must check c(0,t) | 198–199 | test gap |
| C5 | Conditioning: CGL D² entries and A_ grow O(N⁴); the NONSYMMETRIC eigensolve + explicit V⁻¹ amplify roundoff as nch grows (f64 fine at nch ≤ ~15 [inferred — confirm with the §3 oracle's error-vs-nch curve]); model build must stay f64 even for an f32 GPU state instantiation | 148–177 | numeric |
| C6 | Static singleton `makeModel()` = one geometry for ALL cells (already in PLAN §2.1; v4 per-batch build fixes) | 203–207 | known |
| C7 | Stale comment "326 elements" (struct is templated); sin/cos node-form mixing (cosmetic) | 26, 70/83 | trivial |

## 5. Improvements for v4 `build()` (ranked; slot-compatible per PLAN §3.2)

1. **Adopt the §3 analytic-eigenvalue check as a BUILD-TIME gate** (not only a unit test): after eigensolve,
   verify (4)-roots to a registered tolerance; Status-fail otherwise. Catches C1–C3, C5 in one check, at zero
   hot-path cost. Also fixes cross-electrode zero-index pairing by construction (match modes BY VALUE against
   μ_k², never by solver ordering).
2. **Optional: skip the numeric eigensolve entirely.** Λ is known from (4); eigenfunctions are sin(μ_k x)/norm —
   B and C can be built by collocating the ANALYTIC eigenbasis (Galerkin-in-eigenfunctions). Removes the
   nonsymmetric-solver and V-conditioning concerns for good; the discrete collocation C-map is then only an
   interpolation choice. Candidate for high-nch/DFN-era work, not Phase 1.
3. **Mass-mode consistency check:** dz₀/dt = B₀·j must reproduce d(total Li)/dt = −3j/(R·1) exactly in the
   discrete model (one dot product at build; complements P3-G2 system-level conservation).
4. **Alternative discretisations for the same slot** (user directive 2026-07-09: "Chebyshev and maybe some
   others"): (a) conservative FVM on r (PyBaMM's default — valuable for like-for-like PyBaMM parity in P7-G1);
   (b) 3-parameter polynomial (parabolic) profile — ultra-cheap low-fidelity tier for 10⁵-cell screening;
   (c) eigenfunction/Duhamel method (= improvement 2 taken to its conclusion); (d) Legendre–Galerkin (symmetric
   by construction, kills C5 differently). Each is one `StateSpec`+kernel+registry entry (§3.2 axes).

References: Trefethen, *Spectral Methods in MATLAB* (2000), chs. 6–7; Weideman & Reddy, ACM TOMS 26 (2000)
(`poldif.m`); Bizeray, Zhao, Duncan, Howey, J. Power Sources 296 (2015) 400–412 (the SPM spectral formulation
this implements); Carslaw & Jaeger, *Conduction of Heat in Solids* 2e §9.5; Crank, *Mathematics of Diffusion*
2e §6.3; Chebfun `cumsummat`.
