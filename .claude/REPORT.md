# Chebyshev Discretisation Investigation Report

## Summary

The Chebyshev spectral discretisation in SLIDE (`Model_SPM.hpp`) had **three bugs** preventing it from working for `nch != 5`. The matrix construction (differentiation matrices, state-space A/B/C/D, integration matrix Q) was correct for arbitrary `nch`, but runtime issues in eigenvalue detection, Eigen aliasing, and a parity-dependent sign error in the centre concentration formula broke it.

## Background

SLIDE uses Chebyshev spectral collocation to discretise the solid diffusion PDE in spherical particles. The implementation is a C++ translation of the MATLAB reference (`matlab/get_model_vk_slide.m`), itself inspired by [Spectral_li-ion_SPM](https://github.com/davidhowey/Spectral_li-ion_SPM).

The continuous PDE (in spherical coordinates with variable transformation u = r*c):

```
du/dt = D_s * d^2u/dr^2
```

is discretised into a state-space model:

```
dz/dt = D * A * z + B * j    (nch eigenspace state equations)
c     = C * z + D * j/D_s    (nch+1 output equations for surface + inner nodes)
```

where `z` is the twice-transformed concentration (first r*c, then eigendecomposition), `j` is the molar flux boundary condition, and nch is the number of inner Chebyshev nodes.

## Root Causes

### BUG 1 (Critical): Zero eigenvalue detection — absolute vs relative threshold

**Location:** `src/cells/Cell_SPM/Model_SPM.hpp`, lines 162-164 (original)

```cpp
A[pos].array().abs().minCoeff(&zero);      // finds min |eigenvalue|
A[pos](zero) = A[neg](zero) = 0.0;        // applies SAME index to both
```

**Two problems:**

1. **Absolute vs relative comparison.** The code found the eigenvalue with the smallest absolute value, assuming it was the zero eigenvalue. For `nch=5`, the eigenvalue spectrum happens to be well-separated enough that this works. For other `nch` values, the eigenvalue spectrum is denser, and the "near-zero" eigenvalue may not be the absolute minimum due to numerical noise. The MATLAB reference normalises by `max(abs(d))` first (relative comparison), which is robust across all `nch`.

2. **Shared index across electrodes.** A single `zero` index from `A[pos]` was applied to `A[neg]`. Since `Eigen::EigenSolver` runs independently on each matrix and does not guarantee eigenvalue ordering, the zero eigenvalue could theoretically be at different indices in the two decompositions.

### BUG 2 (Critical): Eigen aliasing on `V.inverse()` for small matrices

**Location:** `src/cells/Cell_SPM/Model_SPM.hpp`, line ~175 (original)

```cpp
V[pos] = V[pos].inverse();   // aliasing: undefined behaviour for matrices size 2-4
V[neg] = V[neg].inverse();
```

For `nch=5`, the matrix is 5×5 — Eigen uses a general path that happens to avoid aliasing. For `nch=3` (3×3 matrix) or `nch=4` (4×4), Eigen uses a specialised small-matrix inverse that evaluates lazily, and writing the result back to the same matrix triggers an aliasing assertion in debug mode and produces incorrect results in release mode.

**Fix:** Add `.eval()` to force eager evaluation before assignment:
```cpp
V[pos] = V[pos].inverse().eval();
V[neg] = V[neg].inverse().eval();
```

### BUG 3 (Critical): Centre concentration sign error for even nch

**Location:** `src/cells/Cell_SPM/Cell_SPM.cpp`, line 196 (original)

```cpp
concentration[nch + 1] = -0.5 * (cpt + molarFlux * R / Dt);
```

The hardcoded `-0.5` coefficient is derived from the entry `DM1(N)` of the first-order Chebyshev differentiation matrix at the centre node, where `N = nch + 1`. The correct formula is:

```
c_centre = -1/DM1(N) * (Cc . c + flux*R/D)
```

Since `DM1(N) = 2 * (-1)^N`:
- **Even N** (odd nch = 3, 5, 7, ...): `DM1(N) = 2`, coefficient = `-0.5` ✓
- **Odd N** (even nch = 4, 6, 8, 10, ...): `DM1(N) = -2`, coefficient = `+0.5` ✗ (code uses -0.5)

**Fix:** Store the actual coefficient `cc_coeff = -1.0 / DM1(N)` in the model and use it instead of the hardcoded `-0.5`.

## Items Confirmed NOT Bugs

| Item | Conclusion |
|------|-----------|
| `Cell_SPM.cpp:180` `//!< Problem here!!!!!!! #TODO` | Loop is dimensionally correct. C matrix has N=nch+1 rows, D has N elements, z has nch elements. Comment updated. |
| "Model.Dn and Dp created as nch elements but require nch+1 elements" | Naming confusion: `D[dom]` (state-space output vector) has N=nch+1 elements (correct). `st.Dp()`/`st.Dn()` are scalar diffusion constants (different thing). |
| Differentiation matrix construction | C++ matches MATLAB line-for-line: Chebyshev nodes, DX formula, boundary coefficient C_vk, iterative higher-order computation. Correct for arbitrary nch. |
| Symmetry exploitation (DN1, DN2) | `leftCols(N) - rightCols(N).rowwise().reverse()` correctly implements the odd-symmetry reduction matching MATLAB's `P = fliplr(eye(N))`. |
| Cc centre-node matrix | Uses `+` sign (even symmetry for concentration), matching MATLAB. The matrix itself is correct; only the scalar coefficient was wrong. |
| cumsummat integration matrix | Template parameter mapping is correct: `cumsummat<M>()` with M=2*(nch+1). |

## Test Coverage

New test file `tests/unit/Chebyshev_test.cpp` with 10 test cases parameterised over nch = {3, 5, 7, 10} (37 test cases total, 397 assertions):

| Test | What it validates |
|------|-------------------|
| Node positions | xch values in (0,1), strictly decreasing |
| Eigenvalues non-positive | All A[dom] eigenvalues <= 0 |
| Zero eigenvalue at stored index | Exactly one zero per electrode, both at `model.zero` |
| Uniform concentration round-trip | setC logic -> getC logic recovers uniform c at all nodes |
| Centre concentration via Cc | Cc with `cc_coeff` produces correct centre value for uniform c |
| Zero time-derivative | dz/dt = 0 for uniform c at zero current |
| Q matrix first row zero | Q(0,:) = 0 (integral boundary condition) |
| Linear profile recovery | State-space model recovers c(r) = r/R |
| nch=5 backward compatibility | Kokam NMC setC -> getC round-trip produces correct concentration |
| Matrix dimensions | All matrices have correct compile-time sizes |

## Files Changed

| File | Change |
|------|--------|
| `src/cells/Cell_SPM/Model_SPM.hpp` | Fixed zero eigenvalue detection, Eigen aliasing, added `cc_coeff` member, imaginary-part warning |
| `src/cells/Cell_SPM/Cell_SPM.cpp` | Use `M->cc_coeff` instead of hardcoded `-0.5`; updated misleading "Problem here" comment |
| `src/settings/settings.hpp` | Replaced "DON'T CHANGE" warning with documentation of valid ranges (3-15) |
| `src/settings/constants.hpp` | Fixed integer overflow in `TIME_INF` (was `250 * 365 * 24 * 3600`, now `250.0 * ...`) |
| `tests/unit/Cell_SPM_test.cpp` | Updated z-state check to be eigenvector-ordering-independent |
| `develop/TODO.md` | Marked two Chebyshev items as resolved |
| `tests/unit/Chebyshev_test.cpp` | **New:** 10 parameterised unit tests across nch = {3, 5, 7, 10} |
| `tests/unit/CMakeLists.txt` | Registered new test target |

## Verification

```
All tests passed (397 assertions in 37 test cases)
```

The fix enables `nch` values in the range 3-15 (and potentially beyond). Existing Cell_SPM tests pass with the same physical outputs (surface concentrations, voltage, time stepping).
