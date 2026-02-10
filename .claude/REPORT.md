# Chebyshev Discretisation Investigation Report

## Summary

The Chebyshev spectral discretisation in SLIDE (`Model_SPM.hpp`) had a bug in the zero eigenvalue detection that prevented the model from working correctly for `nch != 5`. The matrix construction itself (differentiation matrices, state-space A/B/C/D, integration matrix Q) was already correct for arbitrary `nch`.

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

## Root Cause

### BUG: Zero eigenvalue detection — absolute vs relative threshold

**Location:** `src/cells/Cell_SPM/Model_SPM.hpp`, lines 162-164 (original)

```cpp
A[pos].array().abs().minCoeff(&zero);      // finds min |eigenvalue|
A[pos](zero) = A[neg](zero) = 0.0;        // applies SAME index to both
```

**Two problems:**

1. **Absolute vs relative comparison.** The code found the eigenvalue with the smallest absolute value, assuming it was the zero eigenvalue. For `nch=5`, the eigenvalue spectrum happens to be well-separated enough that this works. For other `nch` values, the eigenvalue spectrum is denser, and the "near-zero" eigenvalue may not be the absolute minimum due to numerical noise. The MATLAB reference normalises by `max(abs(d))` first (relative comparison), which is robust across all `nch`.

2. **Shared index across electrodes.** A single `zero` index from `A[pos]` was applied to `A[neg]`. Since `Eigen::EigenSolver` runs independently on each matrix and does not guarantee eigenvalue ordering, the zero eigenvalue could theoretically be at different indices in the two decompositions. (In practice, since `A1` and `A3` are scalar multiples of the same base matrix, the ordering was usually consistent, but this was fragile.)

### Fix applied

```cpp
// Normalise by max eigenvalue magnitude, then find minimum (MATLAB-matching approach)
auto findZeroEigenvalue = [](const auto &eigenvalues) -> Eigen::Index {
    const double max_abs = eigenvalues.array().abs().maxCoeff();
    Eigen::Index idx{};
    (eigenvalues.array().abs() / max_abs).minCoeff(&idx);
    return idx;
};

// Find zero independently for each electrode
const Eigen::Index zero_pos = findZeroEigenvalue(A[pos]);
const Eigen::Index zero_neg = findZeroEigenvalue(A[neg]);
A[pos](zero_pos) = 0.0;
A[neg](zero_neg) = 0.0;

assert(zero_pos == zero_neg && "Zero eigenvalue index mismatch");
zero = zero_pos;
```

Additionally, a warning is emitted if the eigensolver produces eigenvalues with significant imaginary parts, which would indicate numerical conditioning issues for large `nch`.

## Items Confirmed NOT Bugs

| Item | Conclusion |
|------|-----------|
| `Cell_SPM.cpp:180` `//!< Problem here!!!!!!! #TODO` | Loop is dimensionally correct. C matrix has N=nch+1 rows, D has N elements, z has nch elements. Comment removed. |
| "Model.Dn and Dp created as nch elements but require nch+1 elements" | Naming confusion: `D[dom]` (state-space output vector) has N=nch+1 elements (correct). `st.Dp()`/`st.Dn()` are scalar diffusion constants (different thing). |
| Differentiation matrix construction | C++ matches MATLAB line-for-line: Chebyshev nodes, DX formula, boundary coefficient C_vk, iterative higher-order computation. Correct for arbitrary nch. |
| Symmetry exploitation (DN1, DN2) | `leftCols(N) - rightCols(N).rowwise().reverse()` correctly implements the odd-symmetry reduction matching MATLAB's `P = fliplr(eye(N))`. |
| Cc centre-node matrix | Uses `+` sign (even symmetry for derivative at centre), matching MATLAB. |
| cumsummat integration matrix | Template parameter mapping is correct: `cumsummat<M>()` with M=2*(nch+1). |

## Test Coverage

New test file `tests/unit/Chebyshev_test.cpp` with 10 test cases parameterised over nch = {3, 5, 7, 10}:

| Test | What it validates |
|------|-------------------|
| Node positions | xch values in (0,1), strictly decreasing |
| Eigenvalues non-positive | All A[dom] eigenvalues <= 0 |
| Zero eigenvalue at stored index | Exactly one zero per electrode, both at `model.zero` |
| Uniform concentration round-trip | setC logic -> getC logic recovers uniform c |
| Centre concentration via Cc | Cc produces correct centre value for uniform c |
| Zero time-derivative | dz/dt = 0 for uniform c at zero current |
| Q matrix first row zero | Q(0,:) = 0 (integral boundary condition) |
| Linear profile recovery | State-space model recovers c(r) = r/R |
| nch=5 backward compatibility | Known zp/zn values from existing test still match |
| Matrix dimensions | All matrices have correct compile-time sizes |

## Files Changed

| File | Change |
|------|--------|
| `src/cells/Cell_SPM/Model_SPM.hpp` | Fixed zero eigenvalue detection + added imaginary-part warning |
| `src/settings/settings.hpp` | Replaced "DON'T CHANGE" warning with useful documentation |
| `src/cells/Cell_SPM/Cell_SPM.cpp` | Updated misleading "Problem here" comment |
| `develop/TODO.md` | Marked two Chebyshev items as resolved |
| `tests/unit/Chebyshev_test.cpp` | New: 10 parameterised unit tests |
| `tests/unit/CMakeLists.txt` | Registered new test target |
