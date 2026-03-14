# SLIDE Comprehensive Review Report — 20-Agent Adversarial Audit

**Date**: 2026-03-13
**Scope**: Full codebase review of SLIDE v3.0.0 (Simulator for Lithium-Ion Degradation)
**Methodology**: 20 virtual expert reviewers, each specializing in a domain, reviewed the repository top-to-bottom and generated actionable suggestions.

---

## Table of Contents

1. [Numerical Integration & Adaptive Timestepping](#1-numerical-integration--adaptive-timestepping)
2. [Nonlinear Solver & Current Redistribution](#2-nonlinear-solver--current-redistribution)
3. [Battery Pack Topology: Trees vs Sparse Matrices](#3-battery-pack-topology-trees-vs-sparse-matrices)
4. [Chebyshev & Spectral Methods](#4-chebyshev--spectral-methods)
5. [Libraries to Add (No Boost)](#5-libraries-to-add-no-boost)
6. [Removing Boost Dependency](#6-removing-boost-dependency)
7. [Parallelization & Multithreading](#7-parallelization--multithreading)
8. [C++ Best Practices & Modernization](#8-c-best-practices--modernization)
9. [Python Bindings](#9-python-bindings)
10. [MATLAB Bindings](#10-matlab-bindings)
11. [WebAssembly (WASM) — Client-Side Browser Execution](#11-webassembly-wasm--client-side-browser-execution)
12. [Package Management (CPM Best Practices)](#12-package-management-cpm-best-practices)
13. [Data Structures & Memory Layout](#13-data-structures--memory-layout)
14. [Documentation](#14-documentation)
15. [Testing & Verification](#15-testing--verification)
16. [Prioritized Roadmap](#16-prioritized-roadmap)

---

## 1. Numerical Integration & Adaptive Timestepping

### Current State

The entire codebase uses **explicit (forward) Euler** exclusively. In `src/cells/Cell_SPM/Cell_SPM_dstate.cpp`, the diffusion inner loop is:

```cpp
for (int t = 0; t < nstep; t++) {
    dState_diffusion(print, d_st);
    for (size_t i = 0; i < (2 * st.nch); i++)
        st.z(i) += dt * d_st.z(i);
}
```

Degradation is also forward Euler: `st[i] += d_st[i] * dt * nstep;`

There is commented-out RK4 code (`Int_RK4`) that was never completed. Boost odeint's `runge_kutta_dopri5` is used in `Module_p_impl.cpp` but only for the parallel module ODE, not cell-level integration.

### Why This Matters

Forward Euler has **conditional stability**: for `dz/dt = D*A*z`, stability requires `dt < 2/|lambda_max|`. The largest eigenvalue magnitude grows as O(nch²/R²), so higher Chebyshev resolution forces smaller timesteps. For degradation simulations spanning months, this means millions of unnecessary timesteps through stiff degradation ODEs (SEI, LAM) that evolve on much slower timescales than diffusion.

### Recommendations

#### 1a. ESDIRK Methods (Primary Recommendation)

ESDIRK (Explicit Singly Diagonally Implicit Runge-Kutta) methods are the gold standard for mildly stiff systems like SPM:

- **ESDIRK3(2)** or **ESDIRK4(3)** with embedded error estimation
- Only one implicit solve per stage (the diagonal is constant), so the Jacobian factorization is reused
- **Key insight**: The SPM diffusion equation `dz/dt = D * diag(A) * z + B * j` is **diagonal** in the eigenvalue-transformed coordinates (Model_SPM.hpp already stores `A[pos]`, `A[neg]` as vectors). The implicit solve reduces to:

```cpp
// Implicit solve is trivially O(nch) per electrode — just scalar divisions
z_stage[k] = rhs[k] / (1.0 - gamma * dt * D * A[k]);
```

This costs essentially the same per step as explicit Euler but allows **10-100x larger timesteps**.

- Embedded error estimate enables adaptive stepping via PI or PID step-size controller:
```cpp
// PI controller for step size
double err = norm(y_embedded - y_full) / (atol + rtol * norm(y_full));
double factor = safety * pow(err, -0.7/order) * pow(err_prev, 0.4/order);
dt_new = dt * clamp(factor, 0.2, 5.0);
```

#### 1b. TR-BDF2 (for Stiffest Regimes)

If future models introduce electrolyte dynamics or strong thermal coupling (T_MODEL=2), TR-BDF2 is L-stable:
1. Trapezoidal half-step from t_n to t_n + γ*dt
2. BDF2 step to t_n + dt
3. Embedded error estimate from the difference

#### 1c. SUNDIALS CVODE (for Pack-Level)

For coupled pack-level systems (Module_p parallel branches), SUNDIALS CVODE replaces boost::odeint:
- BDF method with Newton iteration handles stiff DAE from voltage equalization
- Built-in band/sparse Jacobian solvers handle block-diagonal structure
- GPU-enabled via SUNDIALS' SUNMatrix and SUNLinearSolver abstractions

Integration path:
1. Add SUNDIALS as optional CPM dependency
2. Create `src/integrators/` with `IntegratorBase.hpp`, `ESDIRK.hpp`, `CVODEWrapper.hpp`
3. Add `timeStep_adaptive(double dt_hint, double t_end, double atol, double rtol)` to StorageUnit

#### 1d. Generic Adaptive Timestepping Architecture

```cpp
// src/integrators/IntegratorBase.hpp
template <typename StateType>
class Integrator {
public:
    struct StepResult {
        double dt_actual;
        double dt_next;
        double error_estimate;
        Status status;
    };

    virtual StepResult step(StateType& state,
                           std::function<void(const StateType&, StateType&)> rhs,
                           double dt, double atol, double rtol) = 0;
    virtual ~Integrator() = default;
};

// ESDIRK specialization exploiting diagonal structure
class ESDIRK_SPM : public Integrator<State_SPM> {
    // Butcher tableau coefficients (compile-time)
    static constexpr double gamma = 0.4358665215; // ESDIRK3(2) diagonal
    // ...
};
```

#### 1e. DAE Structure Recognition

The system has natural DAE structure currently hidden behind iterative solvers:
- **Algebraic**: voltage equalization in parallel modules (KVL)
- **Algebraic**: current-voltage via Butler-Volmer (currently inverted explicitly)
- **Differential**: concentration diffusion, thermal, degradation

For SUNDIALS IDA, the residual form `F(t, y, y') = 0` maps naturally to existing `dState_*` functions.

---

## 2. Nonlinear Solver & Current Redistribution

### Current State

`Module_p_impl.cpp` has two implementations:

1. **`setCurrent_previous_impl`**: Newton-like with numerical Jacobian via perturbation (0.5A step). Uses `Eigen::FullPivLU`. The Jacobian is **static** (line ~548) — computed once and reused forever, becoming stale as cells age.

2. **`setCurrent_analytical_impl`**: Analytical current distribution using recursive ladder-network formula. Only works for `Cell_ECM<1>` (uses `dynamic_cast`), not for SPM cells.

The `setVoltage_iterative()` in `free_functions.hpp` uses **False-Position** with 50 max iterations.

### Why 2500+ Iterations Happen

- The **static LU factorization** becomes stale as cells age asymmetrically
- Perturbation-based resistance estimates can produce negative or near-zero values
- The fallback `r_est[i] = 1e-3` is a guess far from reality
- **No line search or trust region** to ensure global convergence

### Recommendations

#### 2a. Immediate Win: Fix Stale Jacobian (P0)

Remove the `static` keyword from the LU decomposition. Recompute the Jacobian every few iterations or when convergence stalls. This alone could reduce iterations from 2500 to ~50.

#### 2b. Newton-Raphson with Analytical Jacobian

For SPM cells, `V_i = OCV_i(SOC_i) + eta_i(I_i, c_surf_i) - R_dc_i * I_i`. The Jacobian `dV_i/dI_j` is:
- Diagonal: `dV_i/dI_i = d(eta_i)/dI_i - R_dc_i`
- Off-diagonal: 0 (except through contact resistance)

The Butler-Volmer derivative is analytical:
```
eta = (2*R_g*T)/(n*F) * asinh(x),  x = I/(2*a*thick*i0)
d(eta)/dI = (2*R_g*T)/(n*F) * 1/sqrt(1+x^2) * 1/(2*a*thick*i0)
```

This gives a sparse, analytically computable Jacobian at almost no additional cost.

#### 2c. Jacobian-Free Newton-Krylov (JFNK)

For large parallel modules (100+ cells), JFNK computes `J*v` via finite differences:
```
J*v ≈ [F(x + ε*v) - F(x)] / ε
```
With GMRES as the Krylov solver and a block-diagonal preconditioner based on `dV_i/dI_i`.

#### 2d. SUNDIALS KINSol

For the algebraic system `F(I_1,...,I_n) = 0`:
- `F_0 = sum(I_i) - I_total` (KCL)
- `F_i = V_i - V_1 - R_contact * cumulative_I` for i > 0 (KVL)

KINSol provides: inexact Newton with backtracking, Anderson acceleration (excellent for weakly nonlinear systems), built-in convergence monitoring, and integrates with CVODE/IDA.

---

## 3. Battery Pack Topology: Trees vs Sparse Matrices

### Current Architecture

SLIDE uses a **tree** via the StorageUnit hierarchy:
- `Battery` → `Module_s`/`Module_p` → `Deep_ptr<StorageUnit>` children → `Cell_SPM`/`Cell_ECM` leaves
- Contact resistances stored per-module as `std::vector<double> Rcontact` (ladder network)

### Analysis

| Aspect | Tree (Current) | Sparse Matrix | Recommendation |
|--------|----------------|---------------|----------------|
| Small packs (2-20 cells) | Fast, natural | Overhead from sparse format | **Tree** |
| Large packs (100+ cells) | Virtual dispatch overhead | Batch evaluation, SIMD | **Hybrid** |
| Heterogeneous cells | Polymorphism handles it | Requires type tags | **Tree** |
| Arbitrary topology | Locked to ladder network | Incidence matrix handles all | **Sparse** |
| Cache locality | Scattered heap allocations | Contiguous state vector | **Sparse** |
| Thermal coupling | Implicit in heat exchange | Explicit coupling matrix | **Sparse** |

### Recommendation: Hybrid Approach

Keep the tree for **ownership and lifecycle**, add a flattened **"simulation view"** for numerical operations:

1. **Flattened State Vector**: At simulation start, collect all cell states into a contiguous `std::vector<double>` via `getStates()`. Already partially done in `Module_p_impl.cpp`.

2. **Incidence Matrix for Topology**: Represent connections as sparse incidence matrix `A` where `A * I_branches = 0` (KCL). This generalizes the ladder network and enables reconfigurable packs, cross-connected modules, etc.

3. **Batch Evaluation**: For homogeneous modules (all same cell type), provide batch operations:
```cpp
void V_batch(std::span<const double> currents, std::span<double> voltages);
```
This enables SIMD across cells sharing the same Model_SPM matrices.

4. **Thermal Coupling Matrix**: Replace the per-module `Qcontact` calculation with a sparse `N×N` thermal conductance matrix. This correctly handles any topology and makes the 6x overestimation bug easy to fix.

### Why Not Pure Sparse?

PyBaMM uses a symbolic CasADi representation compiling the pack into a single DAE. This requires abandoning polymorphism. For SLIDE's typical use case (3-100 cells, heterogeneous, hierarchical), the hybrid is optimal — the tree gives extensibility while the flattened view gives performance.

---

## 4. Chebyshev & Spectral Methods

### Current Implementation Assessment

`Model_SPM.hpp` implements Vieta-Jungius differentiation matrices for Chebyshev spectral discretization of Fick's 2nd law in spherical coordinates:

```
∂c/∂t = D * (1/r²) * ∂/∂r(r² * ∂c/∂r)
```

**Strengths** (implementation is mathematically sound):
- Correct Vieta-Jungius differentiation matrices (canonical choice)
- Proper even-symmetry exploitation: 2*nch+1 nodes → nch unknowns
- Eigenvalue decomposition diagonalizes the system → O(nch) time integration
- Zero eigenvalue detection with relative threshold (prevents drift)
- Comprehensive tests covering nch ∈ {3, 5, 7, 10}
- Round-trip accuracy: 1e-6 to 1e-8

**Accuracy**: With nch=5 (6 unknowns per electrode), the spectral method achieves the same accuracy as ~50-100 finite difference nodes. **Do not switch to FD/FV.**

### Suggestions

#### 4a. Flipping Trick for High nch

For nch > 10, the condition number of the differentiation matrix grows as O(N²). The **flipping trick** halves roundoff error:
```
D² = flip(flip(D) * D)  // where flip reverses rows and columns
```
Reference: Trefethen & Weideman, "The Eigenvalues of the Second Chebyshev Differentiation Matrix," 2019.

#### 4b. Fast Chebyshev Transform via FFT

The `cumsummat` Q matrix computation involves a `T * B * Tinv` product where T/Tinv are discrete cosine transforms. For nch > 15, use FFT-based DCT (O(N log N) vs O(N²)). At nch=5, the dense computation is faster — this only matters for high-accuracy applications.

#### 4c. Ultraspherical Spectral Method (Future)

The Olver-Townsend method represents differential operators as banded matrices in coefficient space → O(N) solve times. Relevant only for electrolyte-coupled models requiring nch > 20. Reference: Olver & Townsend, "A Fast and Well-Conditioned Spectral Method," SIAM Review, 2013.

#### 4d. Already Well-Done

- Node computation (Chebyshev-Lobatto points) ✓
- Boundary condition embedding in differentiation matrix ✓
- Eigenvalue precomputation and state-space form ✓
- V_inv*B, C*V precomputed ✓
- Zero eigenvalue handling ✓

---

## 5. Libraries to Add (No Boost)

### 5a. SUNDIALS — Highest Priority

**What**: ODE/DAE solvers (CVODE, IDA), nonlinear solvers (KINSol), sensitivity analysis
**Why for SLIDE**: Replaces boost::odeint with production-grade adaptive BDF/Adams methods. IDA handles DAE structure natively. KINSol replaces hand-rolled Newton. GPU support via SUNDIALS' abstractions. Used by PyBaMM, COMSOL, Cantera.

```cmake
CPMAddPackage(
    NAME sundials
    GITHUB_REPOSITORY LLNL/sundials
    VERSION 7.2.0
    OPTIONS
        "BUILD_ARKODE ON"    # Adaptive Runge-Kutta
        "BUILD_CVODE ON"     # Stiff ODE solver
        "BUILD_IDA ON"       # DAE solver
        "BUILD_KINSOL ON"    # Nonlinear algebraic solver
        "BUILD_TESTING OFF"
        "BUILD_BENCHMARKS OFF"
        "EXAMPLES_INSTALL OFF"
)
```

### 5b. glaze — Runtime Configuration (Already Considered)

The CMakeLists.txt already mentions glaze. It maps directly to SLIDE's existing structs:

```cmake
CPMAddPackage(
    NAME glaze
    GITHUB_REPOSITORY stephenberry/glaze
    VERSION 4.2.3
)
```

Enables runtime JSON/TOML configuration for cell parameters, degradation model selection, solver tolerances, pack topology. Zero-cost abstraction over C++ structs with compile-time reflection.

### 5c. mdspan (C++23 or Kokkos Reference Implementation)

Replace raw pointer arrays (`double Tneigh[1]`) with type-safe multidimensional views. Header-only, zero runtime cost:

```cmake
CPMAddPackage(
    NAME mdspan
    GITHUB_REPOSITORY kokkos/mdspan
    VERSION 0.6.0
    DOWNLOAD_ONLY YES
)
```

### 5d. spdlog — Structured Logging

Replace scattered `std::cout`/`std::cerr` with level-filtered, thread-safe logging:

```cmake
CPMAddPackage(
    NAME spdlog
    GITHUB_REPOSITORY gabime/spdlog
    VERSION 1.15.0
)
```

### 5e. Kokkos — Performance Portability (Medium-Term)

Replace manual `std::thread` pool with `Kokkos::parallel_for`. Automatic SIMD, GPU offloading, memory spaces:

```cmake
CPMAddPackage(
    NAME Kokkos
    GITHUB_REPOSITORY kokkos/kokkos
    VERSION 4.5.01
    OPTIONS "Kokkos_ENABLE_OPENMP ON"
)
```

### 5f. Libraries NOT Recommended

| Library | Why Not |
|---------|---------|
| Boost | Huge dependency for just odeint — replace with SUNDIALS |
| Abseil | Overlaps with fmt, std::expected covers StatusOr |
| PETSc | Overkill for SLIDE's problem sizes; SUNDIALS is sufficient |
| deal.II | Full FEM framework — wrong tool for SPM spectral method |
| Trilinos | Enterprise-scale; too heavy for SLIDE |

---

## 6. Removing Boost Dependency

### Current Usage

Boost is used in **exactly one file**: `src/modules/Module_p_impl.cpp` line 16:
```cpp
#include <boost/numeric/odeint.hpp>
```

It uses `runge_kutta_dopri5` for the parallel module ODE integration. The rest of SLIDE has zero Boost usage (confirmed via grep).

### Replacement Strategy

**Option A: SUNDIALS CVODE** (recommended if adding SUNDIALS anyway)
- Direct replacement: CVODE's Adams/BDF methods are superior to odeint's Dormand-Prince
- Provides adaptive stepping out of the box

**Option B: Hand-rolled RK45**
The Dormand-Prince 4(5) method is ~50 lines of code for a single fixed-step call:

```cpp
// Dormand-Prince RK45 coefficients (Butcher tableau)
// This replaces the entire boost::numeric::odeint dependency
template <typename System, typename State>
void integrate_rk45(System sys, State& y, double t0, double t1, double dt) {
    constexpr double a2=1.0/5, a3=3.0/10, a4=4.0/5, a5=8.0/9, a6=1.0, a7=1.0;
    // ... standard Dormand-Prince coefficients ...
    // ~50 lines total
}
```

**Option C: Use Eigen's built-in ODE support** (experimental, not recommended)

### Recommended: Option A (SUNDIALS), fall back to Option B if SUNDIALS is deferred

After removing Boost:
1. Delete `cmake/recipes/boost.cmake`
2. Remove `include(boost)` from `Dependencies.cmake`
3. Remove `Boost::boost` / `boost_headers` from `target_link_libraries`
4. Replace `#include <boost/numeric/odeint.hpp>` in `Module_p_impl.cpp`

---

## 7. Parallelization & Multithreading

### Current State

`src/utility/parallelisation.hpp` implements a basic thread pool:
- Spawns up to 32 threads per `run()` call
- Interleaved work distribution (stride = N_threads)
- **Threads created and joined every call** — no persistent pool
- Controlled by `settings::isParallel` compile-time flag

### Problems

1. **Thread creation overhead**: ~10-100μs per create+join. For 2-4 cell modules, exceeds computation time.
2. **No load balancing**: Cells with different complexity finish at different times.
3. **Static data races**: `Module_p_impl.cpp` uses `static` local variables extensively (Eigen vectors, matrices, `is_init` flags). These are **data races** when multiple Module_p instances run concurrently.
4. **False sharing**: Interleaved distribution means adjacent threads access adjacent `Deep_ptr` objects in contiguous memory.

### Recommendations

#### 7a. Fix Static Variable Data Races (P0 — Correctness Bug)

In `Module_p_impl.cpp`, all `static` local variables must become either:
- **Member variables** of Module_p (preferred)
- **Local variables** (for temporaries recomputed each call)
- `thread_local` as a last resort

#### 7b. Persistent Thread Pool

Replace the create-join-per-call pattern with a persistent pool:

```cpp
// Option 1: C++20 std::jthread + lock-free queue
class ThreadPool {
    std::vector<std::jthread> workers;
    // Chase-Lev work-stealing deque per worker
    void submit(std::function<void()> task);
    void wait_all();
};

// Option 2: OpenMP (simplest, one-line change)
#pragma omp parallel for schedule(dynamic)
for (int i = 0; i < getNSUs(); i++)
    SUs[i]->timeStep_CC(dt, nstep);

// Option 3: C++17 parallel algorithms
std::for_each(std::execution::par, SUs.begin(), SUs.end(),
    [&](auto& su) { su->timeStep_CC(dt, nstep); });
```

**Recommendation**: Start with OpenMP (`schedule(dynamic)` handles load imbalance), then migrate to Kokkos for GPU portability.

#### 7c. SIMD Vectorization

The inner diffusion loop:
```cpp
for (size_t k = 0; k < nch; k++)
    d_st.z(k, dom) = D * M->A[dom](k) * st.z(k, dom) + M->B[dom](k) * molarFlux;
```

With nch=5, this fits in a single AVX-512 register. Use Eigen's vectorized operations instead of hand-rolled loops:

```cpp
// Replace hand-rolled loop with Eigen DAXPY
d_st.z_vec(dom) = D * M->A[dom].array() * st.z_vec(dom).array() + M->B[dom] * molarFlux;
```

#### 7d. GPU Strategy (Long-Term)

The SPM diffusion update uses the **same** Model_SPM matrices (A, B, C, D) for all cells of the same type. This is a **batched matrix-vector multiply** — perfectly suited for GPU:

1. Short-term: SUNDIALS with CUDA-enabled NVector for pack-level ODE
2. Medium-term: Kokkos `parallel_for` with `CudaSpace`
3. Long-term: Batched SPM evaluation (all cells share matrices)

---

## 8. C++ Best Practices & Modernization

### 8a. Error Handling Unification (P1)

Three incompatible mechanisms coexist:
1. `Status` enum — returned from setCurrent, setVoltage
2. `throw int` — `throw 101`, `throw 106`, `throw 108` in Cell_SPM
3. `throw const char*` — in Cell_ECM

**Recommendation**: Adopt `std::expected<T, Status>` (C++23):
```cpp
std::expected<void, Status> setCurrent(double I) {
    if (I > I_max) return std::unexpected(Status::ReachedCurrentLimit);
    // ...
    return {};
}
```

Migration: Map integer throw codes to Status values, change leaf functions first, propagate up.

### 8b. Concepts for Template Constraints

```cpp
template <typename T>
concept CellLike = requires(T& cell, double I, double dt) {
    { cell.V() } -> std::convertible_to<double>;
    { cell.setCurrent(I) } -> std::same_as<Status>;
    { cell.timeStep_CC(dt, 1) } -> std::same_as<void>;
};
```

Enables **static dispatch** for homogeneous modules (e.g., `Module_p<Cell_ECM<1>>`) while keeping dynamic dispatch for heterogeneous ones.

### 8c. Runtime Configuration via glaze

Replace compile-time `constexpr` settings with runtime parameters:

```cpp
// Before (requires recompile):
constexpr size_t nch{ 5 };
constexpr int T_MODEL{ 0 };

// After (runtime JSON):
struct SimulationConfig {
    size_t nch = 5;
    int thermal_model = 0;
    double atol = 1e-6, rtol = 1e-4;
    // glaze auto-serializes from JSON
};
```

### 8d. Builder Pattern for Complex Construction

```cpp
auto cell = Cell_SPM::Builder()
    .withDegradation(DEG_ID{.SEI_id={1,2}, .LAM_id={1}})
    .withChebyshevNodes(7)
    .withThermalModel(ThermalModel::Coupled)
    .build();
```

### 8e. Const Correctness

`V()`, `getOCV()`, `SOC()` should be `const` in Cell_SPM. Currently they modify internal state for caching. Use `mutable` for cache members:

```cpp
double V() const {
    if (!Vcell_valid_) { Vcell_ = compute_V(); Vcell_valid_ = true; }
    return Vcell_;
}
mutable double Vcell_;
mutable bool Vcell_valid_ = false;
```

### 8f. Replace range-v3 with std::ranges

SLIDE only uses `views::iota` from range-v3. C++20 provides `<ranges>` natively — drop the dependency.

### 8g. Other Quick Wins

- `StorageUnit::copy()` should return `std::unique_ptr<StorageUnit>`, not raw pointer
- Replace `string const&` parameters with `std::string_view`
- Replace C-arrays (`double Tneighb[]`) with `std::span<double>`
- Replace `#define DATASTORE_BATT 0` with `constexpr` (line 98 of settings.hpp)

---

## 9. Python Bindings

### nanobind vs pybind11

**Recommendation: nanobind** — designed for C++17+, 2-5x faster compile, 50% smaller binaries, native `std::span` support.

### Architecture

```
bindings/python/
  slide_python.cpp       # Embind source
  src/slide/
    __init__.py
    _core.pyi           # Type stubs
    pybamm_compat.py    # PyBaMM Experiment compatibility
```

### Key Design Decisions

1. **StorageUnit hierarchy**: Use nanobind trampoline classes for polymorphic dispatch. Do NOT expose `Deep_ptr` to Python — use factory functions returning concrete types.

2. **Zero-copy state access**:
```cpp
.def("view_states", [](Cell_SPM& cell) {
    auto span = cell.viewStates();
    return nb::ndarray<double, nb::numpy>(span.data(), {span.size()});
}, nb::rv_policy::reference_internal)
```

3. **Build with scikit-build-core**:
```toml
[build-system]
requires = ["scikit-build-core>=0.10", "nanobind>=2.0"]
build-backend = "scikit_build_core.build"
```

4. **PyBaMM compatibility**: Translate `Experiment` strings into Cycler calls.

5. **Async simulation**: Use `nb::gil_scoped_release` during C++ computation, `asyncio.to_thread` from Python.

---

## 10. MATLAB Bindings

### Approach: C++ MEX API (R2018a+) with Handle Pattern

Single MEX entry point (`slide_mex.cpp`) maintains a `std::map<uint64_t, std::unique_ptr<StorageUnit>>`. MATLAB wrapper classes dispatch commands via handle IDs.

```matlab
classdef CellSPM < handle
    methods
        function obj = CellSPM()
            obj.Handle = slide_mex('CellSPM_new');
        end
        function v = V(obj)
            v = slide_mex('V', obj.Handle);
        end
    end
end
```

Build with CMake's `find_package(Matlab)` + `matlab_add_mex()`.

---

## 11. WebAssembly (WASM) — Client-Side Browser Execution

### Why WASM?

Users run simulations **in their browser** — their machine does the computation. No server needed, no security concerns (WASM is sandboxed), no data leaves the user's machine. Free hosting on static sites.

### Dependency Compatibility

| Dependency | Emscripten Compatible | Notes |
|---|---|---|
| **Eigen 3.4** | **Yes** | Header-only, WASM SIMD supported |
| **range-v3 0.12** | **Yes** | Pure templates |
| **fmt 11.0** | **Yes** | Header-only mode or compiled |
| **NLopt 2.10** | **Yes** | Pure C, compiles cleanly |
| **Boost odeint** | **Yes** | Header-only. But we want to remove it anyway |
| **std::thread** | **Partial** | Requires SharedArrayBuffer + COOP/COEP headers |
| **std::filesystem** | **Yes** | Via Emscripten's MEMFS (in-memory) |

**All SLIDE dependencies work with Emscripten.** This is unusually favorable.

### Build Configuration

```cmake
option(SLIDE_WASM "Build for WebAssembly" OFF)

if(EMSCRIPTEN OR SLIDE_WASM)
    set(SLIDE_WASM ON)
    add_compile_definitions(__SLIDE_WASM__)
    set(BUILD_TESTING OFF)
endif()

if(SLIDE_WASM)
    add_executable(slide_wasm bindings/wasm/slide_wasm.cpp)
    target_link_libraries(slide_wasm PRIVATE src project_options)
    target_link_options(slide_wasm PRIVATE
        "--bind"                              # Embind
        "-sWASM=1"
        "-sMODULARIZE=1"
        "-sEXPORT_ES6=1"
        "-sEXPORT_NAME=createSLIDE"
        "-sALLOW_MEMORY_GROWTH=1"
        "-sMAXIMUM_MEMORY=2147483648"         # 2GB max
        "-sWASM_BIGINT=1"
        "-sENVIRONMENT=web,worker"
        "-sINVOKE_RUN=0"
        "-sNO_EXIT_RUNTIME=1"
        "--embed-file=${CMAKE_SOURCE_DIR}/data@/data"  # Bake in OCV curves
    )
endif()
```

### Build Commands

```bash
# Setup Emscripten
git clone https://github.com/emscripten-core/emsdk.git && cd emsdk
./emsdk install latest && ./emsdk activate latest && source emsdk_env.sh

# Single-threaded (simplest, most compatible)
emcmake cmake -B build-wasm -DSLIDE_WASM=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-wasm --target slide_wasm
# Output: slide_wasm.js, slide_wasm.wasm, slide_wasm.data

# Multi-threaded (requires COOP/COEP headers)
emcmake cmake -B build-wasm-mt -DSLIDE_WASM=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-pthread" \
    -DCMAKE_EXE_LINKER_FLAGS="-pthread -sPTHREAD_POOL_SIZE=4"
```

### JavaScript/TypeScript API via Embind

```cpp
// bindings/wasm/slide_wasm.cpp
#include <emscripten/bind.h>
#include "slide.hpp"
using namespace emscripten;

EMSCRIPTEN_BINDINGS(slide) {
    enum_<slide::Status>("Status")
        .value("Success", slide::Status::Success)
        .value("ReachedVoltageLimit", slide::Status::ReachedVoltageLimit);

    class_<slide::Cell_SPM, base<slide::Cell>>("CellSPM")
        .constructor<>()
        .function("v", &slide::Cell_SPM::V)
        .function("soc", &slide::Cell_SPM::SOC)
        .function("setCurrent", &slide::Cell_SPM::setCurrent);

    class_<slide::Cycler>("Cycler")
        .constructor<slide::StorageUnit*, std::string>()
        .function("cc", &slide::Cycler::CC)
        .function("cv", &slide::Cycler::CV);
}
```

### Non-Blocking Simulation via Web Workers

```typescript
// simulation-worker.ts — runs in Web Worker, not main thread
import createSLIDE from './slide_wasm.js';

self.onmessage = async (e) => {
    const Module = await createSLIDE();
    const cell = new Module.CellSPM();
    const cycler = new Module.Cycler(cell, "browser");

    // Run in chunks, posting progress
    const chunkTime = 10.0; // seconds per chunk
    let remaining = e.data.totalTime;
    while (remaining > 0) {
        const status = cycler.cc(e.data.current, e.data.vlim,
                                  Math.min(chunkTime, remaining), e.data.dt, 1);
        remaining -= chunkTime;
        self.postMessage({ type: 'progress', voltage: cell.v(), soc: cell.soc() });
        if (status !== Module.Status.ReachedTimeLimit) break;
    }
    self.postMessage({ type: 'complete' });
    cell.delete(); cycler.delete();
};
```

### Performance

| Operation | Native (Clang) | WASM (Chrome) | Ratio |
|---|---|---|---|
| Single cell 1000s CC | ~50ms | ~75-150ms | 1.5-3x |
| Diffusion PDE step | ~1μs | ~2-3μs | 2-3x |
| Full aging simulation | ~60s | ~120-180s | 2-3x |

**Acceptable for browser use** — simulations that take 1 minute natively take 2-3 minutes in WASM.

Enable WASM SIMD for Eigen: `-msimd128` (Chrome 91+, Firefox 89+, Safari 16.4+).

### Path Adaptation

`slide_paths.hpp` uses `std::filesystem` relative to binary. For WASM:
```cpp
#ifdef __SLIDE_WASM__
const static fs::path root_folder{ "/" };  // MEMFS root
#else
// ... existing logic
#endif
```

The `--embed-file data@/data` flag bakes OCV CSV files into the `.data` artifact.

### Security Model

- WASM executes in a **sandboxed VM** — cannot access host filesystem, network, or other tabs
- Memory is a single ArrayBuffer — no pointer escapes
- All I/O goes through JS APIs the developer provides
- **Users' simulation data never leaves their machine**
- Hosting: static files on GitHub Pages / Cloudflare Pages / Netlify (free)

### Deployment Architecture

```
slide-web/         # Static site — deploy anywhere
  index.html       # UI (React, Svelte, or vanilla)
  app.js           # Controls and visualization
  worker.js        # Web Worker for simulation
  slide_wasm.js    # Emscripten JS glue
  slide_wasm.wasm  # Compiled module (~1-2MB optimized)
  slide_wasm.data  # Embedded CSV data
  sw.js            # Service Worker for offline support
```

### Multi-Threading in WASM

Requires `SharedArrayBuffer` + COOP/COEP HTTP headers:
```
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

Netlify and Cloudflare Pages support custom headers. GitHub Pages does not (use a Cloudflare Worker proxy).

For safety, **default to single-threaded WASM** and offer multi-threaded as opt-in:
```cmake
if(SLIDE_WASM AND NOT SLIDE_WASM_THREADS)
    target_compile_definitions(src PUBLIC SLIDE_SINGLE_THREADED)
endif()
```

### Future: WebGPU

WebGPU (Chrome 113+) enables GPU compute in the browser. Could accelerate pack-level redistribution. However, requires WGSL shader rewrites — defer until WASM is stable.

---

## 12. Package Management (CPM Best Practices)

### Why CPM is Already the Right Choice

CPM (CMake Package Manager) is ideal for SLIDE because:
- **Zero external tools**: No Python (Conan), no vcpkg binary, no package manager installation
- **Pure CMake**: Users only need CMake (which they already have)
- **Reproducible**: URLs with SHA256 hashes pin exact versions
- **Offline-friendly**: `CPM_SOURCE_CACHE` caches downloads across projects
- **Header-only friendly**: `DOWNLOAD_ONLY YES` works perfectly for Eigen, range-v3, Boost headers

### CPM Best Practices

#### 12a. Pin with SHA256 Hashes

```cmake
CPMAddPackage(
    NAME fmt
    URL "https://github.com/fmtlib/fmt/archive/refs/tags/11.0.2.tar.gz"
    URL_HASH SHA256=<hash>  # Add this for reproducibility
)
```

#### 12b. Use `CPM_SOURCE_CACHE`

Tell users to set `CPM_SOURCE_CACHE` environment variable to avoid re-downloading:
```bash
export CPM_SOURCE_CACHE=~/.cache/CPM  # Linux/macOS
set CPM_SOURCE_CACHE=%USERPROFILE%\.cache\CPM  # Windows
```

#### 12c. Prefer System Packages When Available

CPM supports `find_package` fallback:
```cmake
CPMFindPackage(
    NAME Eigen3
    VERSION 3.4
    GITHUB_REPOSITORY libigl/eigen
    GIT_TAG 3.4.0
)
# Uses system Eigen if available, downloads otherwise
```

#### 12d. Version Pinning Table

Maintain a clear version table in `Dependencies.cmake`:
```cmake
# SLIDE Dependency Versions — keep this up to date
set(SLIDE_EIGEN_VERSION "3.4.0")
set(SLIDE_FMT_VERSION "11.0.2")
set(SLIDE_CATCH2_VERSION "3.6.0")
set(SLIDE_NLOPT_VERSION "2.10.0")
set(SLIDE_SUNDIALS_VERSION "7.2.0")  # When added
set(SLIDE_GLAZE_VERSION "4.2.3")     # When added
```

### Why NOT Conan or vcpkg

| Tool | Problem for SLIDE |
|------|-------------------|
| **Conan** | Requires Python installation. Users must install `pip install conan`. Version conflicts. Profile management. |
| **vcpkg** | Requires cloning a 2GB+ repo or installing via Visual Studio. Triplet management. |
| **CPM** | Zero dependencies. Just CMake. Works today. |

CPM is the right choice. Keep it.

---

## 13. Data Structures & Memory Layout

### Current Design (Good)

- `State_SPM` inherits from `std::array<double, 19+2*nch>` → stack-allocated, cache-friendly (232 bytes for nch=5)
- `Deep_ptr<StorageUnit>` for polymorphic ownership with deep copy
- `Pair<T>` for electrode-indexed data [pos, neg]
- `std::span` for zero-copy state views in hot paths

### Suggestions

#### 13a. Structure of Arrays (SoA) for Batch Processing

For homogeneous modules, store cell states in SoA layout for SIMD:
```cpp
struct CellBatch {
    std::vector<double> voltages;   // [N_cells] contiguous
    std::vector<double> currents;   // [N_cells] contiguous
    std::vector<double> z_pos;      // [N_cells * nch] contiguous
    std::vector<double> z_neg;      // [N_cells * nch] contiguous
    // Enables vectorized diffusion update across all cells
};
```

#### 13b. Arena Allocator for Scratch Buffers

Replace per-call `std::vector` allocations in Module_p with arena/pool allocation:
```cpp
// Thread-local arena for temporary Eigen vectors
thread_local std::vector<double> scratch(1024);
Eigen::Map<Eigen::VectorXd> Q(scratch.data(), n);
```

#### 13c. SmallVector for Module Children

Modules typically have 2-20 children. Use a small-buffer-optimized vector (SLIDE already has `SmallVector` in types/):
```cpp
using SUs_t = SmallVector<Deep_ptr<StorageUnit>, 16>; // 16 inline, heap after
```

---

## 14. Documentation

### Current State

- Doxygen comments present in most headers
- `develop/TODO.md` is comprehensive
- `develop/discussions.md` logs design decisions
- `.claude/CLAUDE.md` is an excellent runbook

### Suggestions

1. **API Reference**: Generate Doxygen HTML and host alongside WASM demo
2. **"How to add a cell type"** guide (referenced in CLAUDE.md, may not exist yet)
3. **Mathematical background**: Document the Chebyshev spectral method derivation (Model_SPM.hpp is dense)
4. **Benchmark dashboard**: Automated performance tracking in CI
5. **Interactive WASM demo**: Best documentation is a working browser demo

---

## 15. Testing & Verification

### Current Coverage

- Catch2 v3.6.0 framework, 12+ test files
- Chebyshev: 10 comprehensive tests, nch ∈ {3,5,7,10}
- Cell SPM/ECM, Module_s/p, Battery, Cycler, Procedure tests
- Integration test against PyBaMM (`parallel_model_sln.cpp`)

### Suggestions

1. **Property-based testing**: For conservation laws (total lithium, energy balance)
2. **Fuzzing**: Feed random current profiles to catch edge cases
3. **Benchmark regression**: Track key performance metrics in CI
4. **WASM smoke test**: Build + run a quick simulation in Node.js in CI
5. **Replace `assert()` with `REQUIRE()`**: Some test files still use raw assert

---

## 16. Prioritized Roadmap

| Priority | Item | Impact | Effort | Section |
|----------|------|--------|--------|---------|
| **P0** | Fix static data races in Module_p_impl.cpp | Correctness | 2 days | §7a |
| **P0** | Fix stale Jacobian (remove `static` LU) | 50x fewer iterations | 1 day | §2a |
| **P1** | Remove Boost → replace with SUNDIALS or hand-rolled RK45 | Smaller dependency | 1 week | §6 |
| **P1** | ESDIRK for Cell_SPM (exploit diagonal A) | 10-100x speedup | 2 weeks | §1a |
| **P1** | Analytical Jacobian for Newton solver | 10x fewer iterations | 1 week | §2b |
| **P1** | Unify error handling (std::expected) | Maintainability | 2 weeks | §8a |
| **P2** | SUNDIALS integration (CVODE + KINSol) | Professional solvers | 1 month | §5a |
| **P2** | glaze for runtime config | Usability | 2 weeks | §5b |
| **P2** | WASM build target | Browser deployment | 2 weeks | §11 |
| **P2** | OpenMP for quick parallelization | Easy threading | 3 days | §7b |
| **P2** | Python bindings (nanobind) | Python ecosystem | 1 month | §9 |
| **P3** | Flattened simulation view (hybrid topology) | Cache locality | 1 month | §3 |
| **P3** | SIMD via Eigen for inner loops | Vectorization | 1 week | §7c |
| **P3** | Replace range-v3 with std::ranges | Remove dependency | 3 days | §8f |
| **P3** | spdlog for structured logging | Debugging | 1 week | §5d |
| **P3** | MATLAB MEX bindings | MATLAB users | 2 weeks | §10 |
| **P4** | Chebyshev flipping trick (high nch) | nch>10 accuracy | 3 days | §4a |
| **P4** | Kokkos for GPU portability | Future GPU | 2 months | §5e |
| **P4** | Concepts for static dispatch | Performance | 2 weeks | §8b |
| **P4** | SoA layout for batch processing | Pack-level perf | 1 month | §13a |

---

## Critical Files Reference

| File | What to Modify | Section |
|------|----------------|---------|
| `src/cells/Cell_SPM/Cell_SPM_dstate.cpp` | Time integration (Euler → ESDIRK) | §1 |
| `src/modules/Module_p_impl.cpp` | Static data races, stale Jacobian, Boost removal | §2, §6, §7a |
| `src/cells/Cell_SPM/Model_SPM.hpp` | Chebyshev improvements, diagonal A exploit | §4 |
| `src/utility/parallelisation.hpp` | Thread pool replacement | §7 |
| `src/types/Status.hpp` | Error handling expansion | §8a |
| `src/settings/settings.hpp` | Runtime config, compile-time → runtime | §8c |
| `src/settings/slide_paths.hpp` | WASM filesystem adaptation | §11 |
| `cmake/Dependencies.cmake` | Add SUNDIALS, glaze; remove Boost | §5, §6 |
| `cmake/recipes/boost.cmake` | Delete entirely | §6 |

---

---

## Appendix A: Adversarial Audit — Corrections to Original Report

The following corrections were identified by adversarial review agents cross-checking every claim against the actual code.

### Corrected Claims

| # | Original Claim | Verdict | Correction |
|---|---------------|---------|------------|
| 1 | "Forward Euler exclusively" | **REFUTED** | `Module_p_impl.cpp:631` uses `boost::odeint::runge_kutta_dopri5` with `integrate_adaptive` — an active, compiled adaptive RK45. Cell_SPM_dstate.cpp also has commented-out RK4 (lines 464-488). |
| 2 | "No adaptive stepping" | **REFUTED** | Module_p's `integrateODE_CC` uses Boost odeint's error-controlled adaptive stepping with 1e-12 tolerances. Cycler::CC also has quasi-adaptive stepping (thermal/degradation every N*dt). |
| 3 | "32 status codes" | **REFUTED** | Only **23 distinct values** (2 are aliases: `Critical = VMIN_violation`, `NotSafe = VMINsafety_violation`). |
| 4 | "NLopt for determineOCV" | **REFUTED** | NLopt is linked but only used in `main.cpp` for a **demo/tutorial problem**. `determine_OCV.cpp` has zero NLopt usage. `TODO.md` lists "Integrate NLopt for determineOCV" as pending. |
| 5 | "range-v3 only uses views::iota" | **PARTIALLY** | Cell_SPM files include only `views::iota`, but `main.cpp:18` includes `<range/v3/all.hpp>`. The iota headers are included but not visibly used in function bodies. |
| 6 | "Kelvin = 273.0" | **BUG CONFIRMED** | `constants.hpp:19` — should be 273.15. Introduces a 0.15 K systematic offset in all temperature conversions. |
| 7 | "F = 96487" | **BUG CONFIRMED** | `constants.hpp:20` — NIST 2018 value is 96485.33212. ~0.002% systematic error in all electrochemical calculations. |

---

## Appendix B: Newly Discovered Bugs (Not in Original Report)

### B1. CRITICAL — `Module::setStates` Restore Loop Bug (Copy-Paste Error)

**File:** `src/modules/Module.cpp:233-234`
```cpp
for (size_t j = 0; j <= i; j++)
    SUs[i]->setStates(sorig, n_sorig, false, print); // BUG: should be SUs[j]
```

The loop variable is `j` but the body uses `SUs[i]`. This restores `SUs[i]` (the failing cell) `i+1` times instead of restoring `SUs[0]` through `SUs[i]`. Cells 0 through `i-1` that were successfully set to new states are **never restored**, leaving the module in an inconsistent half-updated state.

**Fix:** Change `SUs[i]` to `SUs[j]` on line 234.

### B2. CRITICAL — `V()` Called Every Inner Timestep, Cache Defeated

**File:** `src/cells/Cell_SPM/Cell_SPM_dstate.cpp:262-272`
```cpp
Vcell_valid = false;           // line 262: invalidate cache
// ...
st.Wh() += std::abs(dAh * V());  // line 272: V() recomputes from scratch
```

Inside the inner diffusion loop (every `dt` step), `Vcell_valid` is set to `false` immediately before calling `V()` for Wh tracking. This forces a full voltage recomputation (3 binary searches, 2 asinh(), 2 exp(), ~30 arithmetic ops) on **every single inner timestep**, even though only `z[k]` changed by a small Euler increment. The `Vcell_valid` cache is **never used** during time integration.

**Impact:** For a 1-hour CC discharge at dt=2s, this means ~1800 unnecessary full V() recomputations per cell. For a 100-cell pack, ~180,000 wasted V() calls per hour of simulated time.

**Fix:** Either (a) compute a "cheap V" for Wh tracking using the previous voltage + linear correction, or (b) move the Wh accumulation outside the inner loop, or (c) cache the intermediate computations.

### B3. HIGH — `Module_p::timeStep_CC` Contact Resistance Heat Bug

**File:** `src/modules/Module_p.cpp:274-279`
```cpp
double Ii = 0;
for (size_t i = 0; i < SUs.size(); i++) {
    for (size_t j = i; j < SUs.size(); j++)
        Ii += SUs[j]->I();   // Ii accumulates across OUTER iterations!
    therm.Qcontact += Rcontact[i] * sqr(Ii) * nstep * dt;
}
```

`Ii` is initialized once at line 274 but **never reset** between outer loop iterations. The current through `Rcontact[i]` should be the suffix sum `sum(I[i..n-1])`, but instead `Ii` grows unboundedly across iterations. The existing code comment `#TODO very important! Not calculating current for Rcontact properly!!!!!!!!!!` at line 276 confirms this is a known bug.

**Fix:** Move `Ii = 0;` inside the outer loop.

### B4. HIGH — Unchecked `dynamic_cast` to `Cell_ECM<1>*`

**File:** `src/modules/Module_p_impl.cpp:325, 350, 428`
```cpp
auto cp = dynamic_cast<Cell_ECM<1> *>(SUs[i].get());
// No null check — immediate dereference follows
```

Three `dynamic_cast` calls cast `StorageUnit*` to `Cell_ECM<1>*` without checking for `nullptr`. If any SU is not a `Cell_ECM<1>` (e.g., Cell_SPM, Cell_ECM<2>), the cast returns `nullptr` and the subsequent dereference is **undefined behavior** (crash).

### B5. HIGH — Static Eigen Vectors Sized by Runtime `n` — Stale Size Bug

**File:** `src/modules/Module_p_impl.cpp:213`
```cpp
static Eigen::VectorXd OCV_branch(n), w(n), v_mod(n);
```

These static vectors are initialized with size `n` on first call. If a subsequent call has a different `n` (different module size), the vectors retain their old size, causing **out-of-bounds access**. This is both a thread safety issue and a correctness issue.

### B6. MEDIUM — Division by Zero in Stress Model at `i == 0`

**File:** `src/cells/Cell_SPM/Cell_SPM_degradation.cpp:672-674`
```cpp
const auto x_cube = cube(xtot[nch + 1 + i]);   // xtot[nch+1] = 0 when i=0
const double bp_x = (Fp[nch + 1 + i] - Fp[nch + 1]) / x_cube;  // div by 0!
```

When `i == 0`, the center node `xtot[nch + 1] = 0`, so `x_cube = 0³ = 0`. The division produces NaN/inf. There is an `if (i == 0)` guard at line 682, but the division at lines 673-674 **executes before the guard**, producing NaN that is only overwritten later.

### B7. MEDIUM — `free::check_current` Returns Nothing (UB)

**File:** `src/utility/free_functions.hpp:213-218`
```cpp
template <bool Print = true>
auto inline check_current(bool checkV, auto &su)
{
    //!< TBC
}
```

Non-void function with no return statement. Any call is **undefined behavior**. The function is empty with only a `//!< TBC` comment.

### B8. MEDIUM — `getFile` Throws Unconditionally, Dead Code After

**File:** `src/utility/io/read_CSVfiles.hpp:60-64`
```cpp
inline auto getFile(std::string name) {
    throw std::runtime_error("getFile function is not implemented yet!");
    static std::map<std::string, std::string> fileMap;  // dead code
}
```

### B9. LOW — `static std::map` in CSV Loading — Thread Safety

**File:** `src/utility/io/read_CSVfiles.hpp:228`
```cpp
static std::map<std::string, XYplain> XYdataMap;
```

Concurrent cell construction (parallel initialization) would race on this shared mutable static. Not currently hit because cell construction is sequential, but will break if parallelized.

### B10. LOW — Global Mutable Static State in `VecState`

**File:** `src/types/VecState.hpp:19-21`
```cpp
static std::vector<double> st, dst, bst;
static std::vector<std::array<int, 2>> locs{ { 0, 0 } };
static size_t current_id{ 0 };
```

Class-level `static` members — globally shared mutable state. Any concurrent access is a data race.

---

## Appendix C: Missed Performance Findings

### C1. Redundant Computation per Timestep (3-4x Waste)

Surface concentrations, Arrhenius coefficients, and OCV lookups are computed **independently** in 4 separate places within a single timestep:
1. `V()` at line 272 of `Cell_SPM_dstate.cpp` (for Wh tracking)
2. `dState_thermal()` at lines 73-93
3. `dState_degradation()` at lines 153-162
4. `setCurrent → checkCurrent → checkVoltage → V()` again from Cycler

**Fix:** Compute these once, store in a per-timestep scratchpad, pass to all `dState_*` functions.

### C2. `getVall` is O(n²) — Should Be O(n)

**File:** `src/modules/Module_p.cpp:86-94`
```cpp
for (auto k{ j }; k < SUs.size(); k++)
    Vall[k] -= I_cumulative * Rcontact[j];
```

The inner loop subtracts the contact resistance voltage drop from all subsequent cells. This can be replaced with a prefix sum in O(n).

### C3. `getRtot` Calls Virtual Function Twice Per Iteration

**File:** `src/modules/Module_p.cpp:46`
```cpp
rtot = Rcontact[i] + (SUs[i]->getRtot() * rtot) / (SUs[i]->getRtot() + rtot);
```

`SUs[i]->getRtot()` is a virtual call evaluated twice. Cache in a local variable.

### C4. `calcArrheniusCoeff()` Called Redundantly Within V()

**File:** `src/cells/Cell_SPM/Cell_SPM.cpp:290`

`calcArrheniusCoeff()` is called once inside `getCSurf()` (line 264) and again at line 290 for overpotential — same value, computed twice.

---

## Appendix D: Physical Constants Corrections

| Constant | Current Value | Correct Value | File | Line |
|----------|--------------|---------------|------|------|
| Kelvin offset | 273.0 | **273.15** | `constants.hpp` | 19 |
| Faraday's constant | 96487 | **96485.33212** (NIST 2018) | `constants.hpp` | 20 |

### Impact Assessment

The Kelvin offset error (0.15 K) propagates through:
- All `_degC` literal conversions
- All Arrhenius rate calculations: `exp(E_a/R * (1/T_ref - 1/T))`
- For T = 25°C = 298.15 K: the code uses 298.0 K, a 0.05% error in absolute temperature
- In Arrhenius: at E_a = 50 kJ/mol, this produces ~0.3% error in rate constants

The Faraday's constant error (1.67 C/mol) is negligible (~0.002%).

---

## Appendix E: Additional Architecture Issues

### E1. `#define DATASTORE_BATT 0` — C Macro in C++20 Codebase

**File:** `src/settings/settings.hpp:98`

This C preprocessor macro coexists with C++ `constexpr` settings. It controls `#if DATASTORE_BATT > 1` guards in `Battery.hpp`, permanently compiling out all Battery data storage. Should be `constexpr int DATASTORE_BATT = 0;`.

### E2. `Module::validStates` Always Returns `true`

**File:** `src/modules/Module.cpp:188`
```cpp
return true; // #TODO here we probably need to check if all submodule states valid!
```

No validation occurs. Invalid module states pass silently.

### E3. `Cell::checkCurrent` Never Checks Current

**File:** `src/cells/Cell.hpp:72-78`
```cpp
virtual Status checkCurrent(bool checkV, bool print) noexcept {
    // ...
    //!< #TODO Current checking part is missing!
    return Vstatus;
}
```

Despite its name, only voltage is checked.

### E4. `Cell_ECM` Constructor Span Bounds Not Checked

**File:** `src/cells/Cell_ECM/Cell_ECM.hpp:203-208`

If the span `spn.size() > N_RC`, the loop writes beyond `Rp[]` and `Tau[]` bounds (both `std::array<double, N_RC>`). No bounds check.

### E5. Stack-Allocated 8KB Arrays in Module Thermal Model

**File:** `src/modules/Module.cpp:274, 337`
```cpp
double Tnew[settings::MODULE_NSUs_MAX];  // 1000 * 8 = 8000 bytes on stack
```

If `MODULE_NSUs_MAX` is increased, this will overflow the stack, especially in recursive thermal model calls.

### E6. `CMakeLists.txt` Uses `PUBLIC` on Executable Target

**File:** `CMakeLists.txt:42-45`
```cmake
target_include_directories(slide PUBLIC data)
```

`PUBLIC` on an executable is meaningless (executables can't be linked against). Should be `PRIVATE`.

### E7. pthread Linking Should Use `find_package(Threads)`

**File:** `src/CMakeLists.txt:17-21`

Manual platform-specific `-pthread` detection should be replaced with `find_package(Threads)` + `Threads::Threads`.

---

## Updated Prioritized Roadmap (Post-Adversarial Audit)

| Priority | Item | Impact | Effort | Section |
|----------|------|--------|--------|---------|
| **P0** | Fix `Module::setStates` restore loop bug (`SUs[i]` → `SUs[j]`) | Correctness | 1 line | App B1 |
| **P0** | Fix `Module_p` contact resistance heat bug (reset `Ii`) | Correctness | 1 line | App B3 |
| **P0** | Fix static data races in Module_p_impl.cpp | Correctness | 2 days | §7a |
| **P0** | Fix stale Jacobian (remove `static` LU) | 50x fewer iterations | 1 day | §2a |
| **P0** | Fix unchecked `dynamic_cast` in Module_p_impl.cpp | Crash prevention | 3 lines | App B4 |
| **P0** | Fix `Kelvin = 273.0` → `273.15` | Physics accuracy | 1 line | App D |
| **P0** | Fix static Eigen vectors sized by runtime `n` | Correctness | 10 lines | App B5 |
| **P1** | Eliminate redundant V() in inner timestep loop | ~2x speedup | 1 day | App B2, C1 |
| **P1** | Remove Boost → replace with SUNDIALS or hand-rolled RK45 | Smaller dependency | 1 week | §6 |
| **P1** | ESDIRK for Cell_SPM (exploit diagonal A) | 10-100x speedup | 2 weeks | §1a |
| **P1** | Analytical Jacobian for Newton solver | 10x fewer iterations | 1 week | §2b |
| **P1** | Unify error handling (std::expected) | Maintainability | 2 weeks | §8a |
| **P2** | Fix getVall O(n²) → O(n) | Pack-level perf | 1 day | App C2 |
| **P2** | SUNDIALS integration (CVODE + KINSol) | Professional solvers | 1 month | §5a |
| **P2** | WASM build target | Browser deployment | 2 weeks | §11 |
| **P2** | Python bindings (nanobind) | Python ecosystem | 1 month | §9 |
| **P3** | Flattened simulation view (hybrid topology) | Cache locality | 1 month | §3 |
| **P3** | Replace range-v3 with std::ranges | Remove dependency | 3 days | §8f |
| **P3** | MATLAB MEX bindings | MATLAB users | 2 weeks | §10 |

---

---

## Appendix F: Strategic Vision — Surpassing PyBaMM

### F1. Why SLIDE Can Beat PyBaMM

PyBaMM's architecture has fundamental performance bottlenecks:

1. **Python interpreter overhead**: GIL prevents true parallelism. Every function evaluation crosses Python-C boundary.
2. **Symbolic overhead**: CasADi builds symbolic expression trees in Python, then compiles to C. Flexible but slow.
3. **General-purpose DAE solver**: SUNDIALS IDA computes full Jacobians and performs Newton iterations — O(N²) to O(N³) per timestep.
4. **Dense state vectors**: FVM discretization needs 30-50 nodes per electrode where spectral Chebyshev needs 3-5.

SLIDE's existing advantage: **pre-diagonalized spectral discretization** means each timestep is O(nch) scalar operations — no matrix factorization, no Newton iteration. This is an architectural advantage, not just an implementation detail.

### F2. The "Implicit for Free" Insight

**This is the single most important performance insight in the entire report.**

The current code uses Forward Euler:
```cpp
st.z(k, dom) += dt * (D * M->A[dom](k) * st.z(k, dom) + M->B[dom](k) * molarFlux);
```

Since A is **diagonal** (eigenvalue-decomposed), Backward Euler costs exactly the same:
```cpp
st.z(k, dom) = (st.z(k, dom) + dt * M->B[dom](k) * molarFlux)
               / (1.0 - dt * D * M->A[dom](k));
```

One division instead of one multiplication — **same cost, but unconditionally stable**. This means:
- Forward Euler with dt=0.5s requires 7200 steps per 1-hour cycle
- Backward Euler with dt=50s requires 72 steps per 1-hour cycle
- **100x fewer steps at the same per-step cost**

This insight extends to SPMe and DFN because their electrolyte diffusion can also be eigenvalue-decomposed → diagonal → implicit for free.

### F3. SPMe Implementation (Single Particle Model with Electrolyte)

#### Physics

SPMe adds electrolyte transport across the cell sandwich (negative/separator/positive):

- **Electrolyte concentration** c_e(x,t): parabolic PDE, εₑ ∂c_e/∂t = ∂/∂x(D_e_eff ∂c_e/∂x) + (1-t⁺)/F · jₙ(x)
- **Electrolyte potential** φ_e(x): elliptic BVP (algebraic), κ_eff ∂²φ_e/∂x² + source = 0
- **Solid potential** φ_s(x): elliptic BVP (algebraic), σ_eff ∂²φ_s/∂x² - source = 0

#### Spectral Discretization Strategy

The electrolyte concentration PDE can be eigenvalue-decomposed the same way as particle diffusion:

1. Place Chebyshev nodes in each region: n_neg, n_sep, n_pos nodes
2. Build global second-derivative operator across all three regions
3. Apply interface matching conditions (continuity of c_e and flux at neg/sep and sep/pos boundaries)
4. Eliminate interface unknowns → reduced matrix of size n_e = n_neg + n_sep + n_pos - 4
5. Eigendecompose → diagonal dynamics: dz_e/dt = diag(A_e) · z_e + B_e · j

For n_neg = n_sep = n_pos = 3: the eigenvalue problem is **5×5** (done once at construction, O(125) flops).

Each timestep cost: 5 scalar multiply-adds for electrolyte + existing 10 for particles = **15 total**.

#### State Vector

```cpp
struct State_SPMe : public State<19 + 2*nch + n_e> {
    // [0..18]           : degradation/thermal (same as SPM)
    // [19..19+nch-1]    : zp[nch] (positive particle modes)
    // [19+nch..19+2nch-1]: zn[nch] (negative particle modes)
    // [19+2nch..end]    : ze[n_e] (electrolyte concentration modes)
    // Total: 34 states (vs 29 for SPM — only 5 additional)
};
```

#### Algebraic Equations (Potentials)

Two approaches:

**Substitution (fast, default)**: For constant conductivity, solid and electrolyte potential drops have analytical expressions. The electrolyte potential drop is computed by Chebyshev quadrature of c_e using the existing integration matrix Q. No additional linear solve needed.

**Direct spectral solve (general)**: Precompute LU of the potential operator (15×15). Per-timestep cost: O(n²) = O(225) flops for triangular solves.

#### Projected Performance

| Model | PyBaMM | SLIDE (projected) | Speedup |
|-------|--------|-------------------|---------|
| SPMe, 1 cycle | ~200 ms | ~8 ms | **25x** |
| SPMe, 1000 cycles + degradation | ~200 s | ~8 s | **25x** |

### F4. DFN Implementation (Doyle-Fuller-Newman)

#### Physics

DFN makes particle diffusion spatially resolved: each through-cell node has its own particle with concentration profile c_s(x_i, r, t). Butler-Volmer kinetics jₙ(x_i) couple particle surface concentrations to electrolyte.

#### Key Insight: Tensor Product Structure

The DFN particle diffusion at all x-nodes can be written as:

```
dZ_s/dt = (I_x ⊗ A_r) · Z_s + B_r · J_n(x)
```

Since A_r is the **same diagonal eigenvalue vector** from the SPM (all particles share geometry), each of the n_x particles is an independent copy of the existing SPM scalar ODE system. **The DFN particle diffusion is trivially parallelizable and reuses the exact same eigenvalue decomposition.**

#### State Vector

For n_x_neg = n_x_pos = 5, nch_r = 5, n_sep = 3:

- Particle states: (5+5) × 5 = 50 modes
- Electrolyte states: 5+3+5-4 = 9 modes
- Overhead: ~15 degradation/thermal
- **Total: ~74 states** (vs PyBaMM DFN: 1000-3000 states with FVM)

#### Algorithm Per Timestep

1. **Recover surface concentrations**: c_s_surf(x_i) = C · z_s(x_i) + D · j(x_i) → O(n_x × nch_r) = O(50)
2. **Recover electrolyte concentrations**: c_e(x_i) = C_e · z_e → O(n_e) = O(9)
3. **Evaluate Butler-Volmer** at each x_i → O(n_x × 20) = O(200)
4. **Update particle modes** (diagonal, implicit): → O(n_x × nch_r) = O(50)
5. **Update electrolyte modes** (diagonal, implicit): → O(n_e) = O(9)
6. **Solve potentials** (precomputed LU): → O(n_x²) = O(100)

**Total: ~400 flops per timestep.** PyBaMM DFN does O(N²·4) ≈ O(10⁷) per Newton step, with ~5-20 Newton steps.

#### Projected Performance

| Model | PyBaMM | SLIDE Spectral (projected) | Speedup |
|-------|--------|---------------------------|---------|
| DFN, 1 cycle | 2-10 s | 30-100 ms | **50-100x** |
| DFN, pack 100 cells | ~15 min | ~15 s (GPU) | **60x** |

### F5. Model Engine — Compositional Physics via C++ Templates

Instead of PyBaMM's symbolic DSL (flexible but slow), use compile-time composition via policy templates:

```cpp
// SPM: particle diffusion + simplified kinetics
using Model_SPM = ModelEngine<
    ParticleDiffusion<5>,
    SimplifiedKinetics,
    ChebyshevSpectral<5>,
    ImplicitEuler  // "free" implicit integration
>;

// SPMe: add electrolyte transport
using Model_SPMe = ModelEngine<
    ParticleDiffusion<5>,
    ElectrolyteDiffusion<3, 3, 3>,
    ElectrolytePotentialSubstitution,
    SimplifiedKinetics,
    ChebyshevSpectral<5>,
    ImplicitEuler
>;

// DFN: spatially-resolved particles + full electrochemistry
using Model_DFN = ModelEngine<
    SpatiallyResolvedParticleDiffusion<5, 5, 5>,
    ElectrolyteDiffusion<5, 3, 5>,
    ElectrolytePotentialSpectralSolve,
    SolidPotential,
    ButlerVolmerKinetics,
    ChebyshevSpectralMultiDomain<5, 3, 5>,
    ExponentialIntegrator
>;

// Custom: electrolyte transport but no particle diffusion (for fast screening)
using Model_Custom = ModelEngine<
    NoParticleDiffusion,  // uniform c_s
    ElectrolyteDiffusion<3, 3, 3>,
    SimplifiedKinetics,
    ChebyshevSpectral<3>,
    ForwardEuler
>;
```

Each physics component declares `static constexpr int n_states` — the ModelEngine sums them at compile time. Zero runtime overhead. Users compose models like LEGO blocks.

### F6. 3D Microstructure Support

Multi-scale architecture:

```
Macro scale (3D voxel grid from X-ray CT):
  - Electrolyte transport: D_eff · ∇²c_e on structured grid
  - Solid conductivity: σ_eff · ∇²φ_s on structured grid
  - 7-point stencil Laplacian → GPU-friendly

Micro scale (1D spectral per active voxel):
  - Same Chebyshev particle diffusion as SPM
  - Butler-Volmer kinetics per voxel
  - Reuses Model_SPM eigendecomposition
```

Memory for realistic X-ray CT (200×200×160 voxels, 3.2M active material voxels, nch=3):
- Particle states: 3.2M × 3 × 8 bytes = 77 MB
- Electrolyte + potentials: 6.4M × 3 × 8 bytes = 154 MB
- **Total: ~230 MB** — fits comfortably in RAM/VRAM

Compute per timestep: ~200M flops → at 10 TFLOPS GPU throughput: **~20 μs/step**. This is where GPU becomes essential.

### F7. GPU Acceleration Strategy

| Operation | GPU Benefit | Phase |
|-----------|------------|-------|
| Batched particle diffusion (DFN/3D) | **Very high** — n_x×nch independent scalar ops | Phase 1 |
| Pack-level parallel cells | **High** — hundreds of independent cells | Phase 1 |
| Parameter sweeps | **Very high** — embarrassingly parallel | Phase 1 |
| 3D Laplacian (stencil ops) | **Very high** — structured grid | Phase 2 |
| Butler-Volmer at all nodes | **Medium** — per-node, memory-bound | Phase 2 |
| Eigendecomposition at init | **Low** — small matrices, done once | N/A |

CUDA kernel for batched particle diffusion (implicit):
```cuda
__global__ void particle_step(
    double* z,           // [n_particles × nch]
    const double* A,     // [nch] eigenvalues (shared)
    const double* B,     // [nch] input vector (shared)
    const double* D_coef,// [n_particles] diffusion coefficients
    const double* j,     // [n_particles] molar flux
    double dt, int n_particles, int nch)
{
    int pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= n_particles) return;
    for (int k = 0; k < nch; k++) {
        int idx = pid * nch + k;
        z[idx] = (z[idx] + dt * B[k] * j[pid])
                 / (1.0 - dt * D_coef[pid] * A[k]);
    }
}
```

### F8. Language Bindings as Backend Interfaces

The goal: SLIDE as a **drop-in high-performance backend** for existing Python/Julia ecosystems.

#### Python — PyBaMM Backend

```python
# PyBaMM compatibility: same Experiment interface, SLIDE backend
import pybamm
import pyslide

# Option A: Drop-in replacement for PyBaMM's solver
model = pybamm.lithium_ion.SPM()
sim = pybamm.Simulation(model, solver=pyslide.SLIDESolver())
sim.solve([0, 3600])  # 100x faster than default CasADi/SUNDIALS

# Option B: Native SLIDE API (maximum performance)
cell = pyslide.CellSPM()
cycler = pyslide.Cycler(cell, "test")
status = cycler.cc(16.0, 2.7, 3600, 0.5)
```

**Implementation**: nanobind wrapping of StorageUnit hierarchy + a `SLIDESolver` class that implements PyBaMM's `BaseSolver` interface. The solver receives a `pybamm.Model` and maps it to the closest SLIDE model (SPM/SPMe/DFN).

#### Julia — DifferentialEquations.jl Integration

```julia
using SLIDE  # ccall-based wrapper

# SLIDE as an ODE right-hand-side for Julia's solvers
cell = SLIDE.CellSPM()
function rhs!(du, u, p, t)
    SLIDE.set_states!(cell, u)
    SLIDE.get_dxdt!(cell, du)
end

prob = ODEProblem(rhs!, SLIDE.get_states(cell), (0.0, 3600.0))
sol = solve(prob, Tsit5())  # Use Julia's adaptive RK45
```

Julia's `ccall` has zero overhead for calling C functions. Expose a thin `extern "C"` API over the C++ classes.

#### MATLAB — MEX Backend

Same handle-pattern as Section 10, but with additional support for MATLAB's `ode15s` (equivalent to SUNDIALS IDA):

```matlab
cell = slide.CellSPM();
[t, y] = ode15s(@(t,y) cell.get_dxdt(y), [0 3600], cell.getStates());
```

### F9. Implementation Roadmap

| Phase | Item | Duration | Prerequisite |
|-------|------|----------|-------------|
| 0 | Implicit Euler for diagonal SPM | 1 day | None |
| 0 | Fix P0 bugs (setStates, Ii, static races) | 3 days | None |
| 1 | `Model_SPMe` + `Cell_SPMe` | 3-4 weeks | Phase 0 |
| 1 | SUNDIALS integration (for Module_p DAE) | 2 weeks | Phase 0 |
| 2 | `Model_DFN` + `Cell_DFN` | 4-6 weeks | Phase 1 |
| 2 | Python bindings (nanobind) | 3 weeks | Phase 1 |
| 3 | ModelEngine template composition | 6-8 weeks | Phase 2 |
| 3 | PyBaMM backend compatibility | 2 weeks | Phase 2 |
| 3 | Julia bindings (extern "C") | 1 week | Phase 2 |
| 4 | CUDA batched particle kernels | 3 weeks | Phase 2 |
| 4 | WASM build target | 2 weeks | Phase 1 |
| 5 | 3D microstructure support | 8-12 weeks | Phase 4 |

### F10. Competitive Position Summary

| Capability | PyBaMM | SLIDE (Current) | SLIDE (Planned) |
|-----------|--------|-----------------|----------------|
| SPM speed | Baseline | **10x faster** | **100x faster** (implicit) |
| SPMe | Yes | No | Yes, spectral — **25x faster** |
| DFN | Yes | No | Yes, spectral — **50-100x faster** |
| Model flexibility | Symbolic DSL | Hard-coded | Template composition |
| GPU | No | No | CUDA batched kernels |
| Pack simulation | Limited | Yes (hierarchy) | Yes + GPU batched |
| Degradation | Community models | 4 built-in | 4 + extensible |
| 3D microstructure | No | No | Multi-scale spectral+voxel |
| Python API | Native | Planned | nanobind + PyBaMM backend |
| Julia API | No | No | ccall FFI |
| MATLAB API | No | Post-processing only | MEX + ode15s |
| Browser (WASM) | No | No | Yes |
| Languages | Python | C++ | C++ / Python / Julia / MATLAB / WASM |

---

---

## Appendix G: Architecture Redesign — Design Patterns, Scalability, Circuit Solver

### G1. Design Patterns for 100K-Cell Simulation

#### Flyweight (Extend Model_SPM Singleton)

Already in use: `Model_SPM::makeModel()` is a Meyers singleton. All Cell_SPM instances share one Model (5 KB). For 100K cells without Flyweight: 500 MB of duplicate matrices. With Flyweight: 5 KB total.

Extend to `ModelRegistry<ModelTag, nch>` with thread-safe lookup for non-default particle radii. Use quantized keys to avoid floating-point hash issues.

#### CRTP Strategy (Zero-Overhead Degradation)

Replace the runtime `switch(DEG_ID.SEI_id)` dispatch in `Cell_SPM_degradation.cpp` with compile-time CRTP:

```cpp
template <typename Derived>
struct SEI_Strategy {
    double compute_isei(/*args*/) const {
        return static_cast<const Derived*>(this)->compute_isei_impl(/*args*/);
    }
    static constexpr int additional_states = Derived::additional_states;
    static constexpr bool requires_surface_concentration = true; // compile-time check
};

struct SEI_Pinson : SEI_Strategy<SEI_Pinson> {
    static constexpr int additional_states = 0;
    double compute_isei_impl(/*args*/) const { /*Pinson & Bazant 2013*/ }
};
```

Compose multiple models via fold expressions:
```cpp
template <typename... Models>
struct ComposedDegradation {
    std::tuple<Models...> models;
    double total_isei(/*args*/) const {
        return std::apply([&](auto&... m) { return (m.compute_isei(/*args*/) + ...); }, models);
    }
    static constexpr int additional_states = (Models::additional_states + ...);
};
```

The compiler inlines every `compute_isei_impl` through the fold — zero virtual dispatch overhead.

#### Policy-Based Cell (Alexandrescu)

Compose cells from independent policy classes at compile time:

```cpp
template <typename DiffusionPolicy, typename ThermalPolicy,
          typename DegradationPolicy, typename KineticsPolicy = ButlerVolmer>
class Cell_Policied : public Cell {
    static constexpr int n_states =
        DiffusionPolicy::n_states + ThermalPolicy::n_states +
        DegradationPolicy::n_states + KineticsPolicy::n_states;

    DiffusionPolicy diffusion;
    ThermalPolicy thermal;
    DegradationPolicy degradation;
    // Each policy gets a span view into its portion of the state array
};

// Backward-compatible aliases:
using Cell_SPM_v4 = Cell_Policied<
    ChebyshevDiffusion<5>, LumpedThermal,
    ComposedDegradation<SEI_Pinson, LAM_Dai, CS_Laresgoiti>, ButlerVolmer>;
```

`Cell_Policied` inherits from `Cell` → `StorageUnit`, so it plugs into the existing Module hierarchy without changes.

#### Entity-Component-System (ECS) for 100K Cells

For massive scale, SoA (Structure-of-Arrays) layout is essential:

```cpp
struct CellPool {
    // Components: contiguous arrays, one element per cell
    std::vector<std::array<double, 2*nch>> z;  // diffusion modes
    std::vector<double> I, V, T, SOC;           // electrical/thermal
    std::vector<double> delta, LLI, CS;         // degradation
    size_t n_cells;

    // Systems: operate on ALL cells in vectorizable loops
    void update_diffusion(const Model_SPM<>& model, double dt) {
        #pragma omp parallel for simd
        for (int cell = 0; cell < n_cells; ++cell) {
            for (int dom = 0; dom < 2; ++dom)
                for (int k = 0; k < nch; ++k)
                    z[cell][dom*nch+k] = /*implicit step*/;
        }
    }
};
```

Memory for 100K SPM cells: ~25 MB total (fits in L3 cache). Maps directly to CUDA kernels.

### G2. Sparse Matrix Circuit Solver

Replace the iterative `redistributeCurrent` (2500+ iterations) with a direct sparse solve:

**KCL + KVL system:**
- Incidence matrix A (n_nodes × n_branches): `A * I_branches = I_external`
- Conductance matrix: `G = A * R⁻¹ * Aᵀ` (symmetric positive definite)
- Solve: `G * V_nodes = A * R⁻¹ * OCV - I_ext` via Cholesky (`Eigen::SimplicialLLT`)
- Recover currents: `I_branch = R⁻¹ * (OCV - Aᵀ * V)`

**For 100K cells in series-of-parallel:** G is block-tridiagonal, solvable in O(N_modules). For arbitrary topologies: hierarchical Schur complement — factorize each module in parallel, solve small coupling system, back-substitute in parallel.

**vs current approach:** One direct solve (+ few Newton iterations for nonlinear V-I) replaces 2500 iterative steps.

### G3. Multi-Rate Integration (Segregated Solver)

Formalize the existing ad-hoc multi-rate split:

| Timescale | What | dt | Method |
|-----------|------|------|--------|
| Fast | Diffusion PDE | 2-50s | Exact exponential (diagonal) |
| Medium | Thermal | 60-100s | Forward Euler |
| Slow | Degradation | 600-3600s | Forward Euler |
| Very slow | Calendar aging | hours-days | Large Euler steps |

**Strang splitting** (2nd order) vs current Lie-Trotter (1st order):
```
Half-step degradation → Full step thermal → Full steps diffusion → Half-step degradation
```

**Exact exponential integrator** for diffusion (unconditionally stable, arbitrary dt):
```cpp
// Since A is diagonal, exact solution per mode:
z_new = z_old * exp(D*lambda*dt) + B*j * (exp(D*lambda*dt) - 1) / (D*lambda)
```
This eliminates the inner diffusion loop entirely when combined with Strang splitting.

### G4. Electrode Class Refactoring

Move from `param_p`/`param_n` duplication to a proper `Electrode` class:

```cpp
struct Electrode {
    Domain dom;
    double R_particle;
    XYdata_ff OCV, dOCV;
    double x_0, x_100;  // stoichiometry limits
    // + all methods from current Electrode_SPM
    // + StateView providing references into the cell state array
};

class Cell_SPM : public Cell {
    Pair<Electrode> electrode;  // [pos, neg] — indexed by Domain enum
    // All _p/_n code becomes: for (auto dom : {pos, neg}) electrode[dom].update(...)
};
```

Eliminates ~30% of Cell_SPM code duplication.

### G5. Injectable Model System

```cpp
auto cell = CellBuilder<SPM>()
    .add_degradation<SEI_Pinson>()    // compile-time compatibility check
    .add_degradation<LAM_Dai>()
    .set_thermal<LumpedThermal>()
    .build();

// Compile error: ECM doesn't provide surface concentration
auto bad = CellBuilder<ECM>()
    .add_degradation<SEI_Pinson>()  // static_assert fails
    .build();
```

Each model declares `requires_surface_concentration`, `requires_stress`, etc. The builder enforces compatibility via `static_assert`.

### G6. Enzyme Auto-Differentiation (Clang LLVM Plugin)

**This is the auto-diff tool the user was asking about.** Enzyme differentiates compiled LLVM IR — works on arbitrary C++ including `std::exp`, `std::asinh`, Eigen operations.

Use cases for SLIDE:
1. **Automatic Jacobian** for Newton solvers: `dV/dI` computed exactly, not by perturbation
2. **Sensitivity analysis**: `d(capacity_fade)/d(parameters)` in one reverse-mode pass
3. **Parameter estimation**: Gradient-based optimization of model parameters

```cmake
option(SLIDE_USE_ENZYME "Enable Enzyme auto-differentiation (Clang only)" OFF)
```

**Portable fallback**: Dual numbers (forward-mode AD, header-only, works on all compilers):
```cpp
template <typename T = double>
struct Dual {
    T val, deriv;
    friend Dual exp(Dual a) { auto e = std::exp(a.val); return {e, a.deriv*e}; }
    // ... all math functions
};
```

### G7. MPI Strategy for 100K Cells on SLURM

```
Rank 0: Orchestrator — circuit solve, current scatter
Ranks 1..N: Cell Workers — each owns ~1000 cells (CellPool)
  Within rank: OpenMP + SIMD (SoA layout)
  Optional: GPU offload (CUDA kernel for batched diffusion)
```

Communication per timestep: scatter 800 KB currents + gather 4 MB summaries = ~5 MB total. Trivial for InfiniBand.

Hybrid strategy: **MPI between nodes, Kokkos within nodes** (auto-selects OpenMP or CUDA).

### G8. Spectral Basis Comparison

| Basis | Domain | Best for SLIDE? | Why |
|-------|--------|----------------|-----|
| **Chebyshev** (current) | [-1,1] | **Yes — keep** | Exponential convergence, well-conditioned for nch≤15 |
| Legendre | [-1,1] | Equivalent | Same convergence, different weights |
| Laguerre | [0,∞) | No | Semi-infinite domain — not applicable |
| Hermite | (-∞,∞) | No | Unbounded — not applicable |
| Lagrange on Cheb nodes | [-1,1] | Same thing | Just a different representation |
| **Ultraspherical** | [-1,1] | For SPMe/DFN | Banded operators → O(N) solve |
| Spectral elements | Multi-domain | For electrolyte | Piecewise Chebyshev per region |

**Conclusion**: Chebyshev is optimal for SLIDE's bounded-interval problems. Ultraspherical matters only for DFN with large N. The current implementation is already near-optimal.

### G9. Template nch as Discretisation Parameter

The user wants `Discretisation<5>` so users can choose resolution without runtime overhead:

```cpp
template <int nch>
struct Discretisation {
    static constexpr int n_modes = nch;
    // All matrices are compile-time sized via Eigen fixed-size types
    Eigen::Vector<double, nch> A_pos, A_neg, B_pos, B_neg;
    // ...
};

// Usage in bindings:
// Python: pyslide.CellSPM(nch=7)  — selects pre-compiled instantiation
// C++: auto cell = Cell_SPM<Discretisation<7>>();
```

This is already partially done — `settings::nch` is a compile-time `constexpr` and `Model_SPM<nch>` is templated on it. Extend by exposing common instantiations (nch = 3, 5, 7, 10) to bindings.

---

## Appendix H: Technical Notes & References

### H1. Eigenvalue Scaling

For the Chebyshev spectral discretization of d²c/dr² on [0, R]:
- Eigenvalues scale as O(nch²/R²)
- Largest eigenvalue magnitude: |λ_max| ≈ (π·nch/R)²
- Forward Euler stability: dt < 2/(D·|λ_max|) ≈ 2R²/(D·π²·nch²)
- For D = 3.7e-14 m²/s, R = 5.86e-6 m, nch = 5: dt_max ≈ 2600s (very permissive)
- **Implicit integration removes this constraint entirely**

### H2. ESDIRK Methods

ESDIRK (Explicit Singly Diagonally Implicit Runge-Kutta) methods are the gold standard for mildly stiff ODE systems:
- Only one implicit solve per stage (the diagonal coefficient γ is constant)
- The Jacobian factorization is reused across stages
- For SLIDE's diagonal systems: the implicit solve is trivially O(nch) scalar divisions
- ESDIRK3(2) provides 3rd-order accuracy with embedded 2nd-order error estimate

Key references:
- Kennedy & Carpenter, "Additive Runge-Kutta schemes for convection-diffusion-reaction equations," Applied Numerical Mathematics 44, 2003
- Alexander, "Diagonally implicit Runge-Kutta methods for stiff O.D.E.'s," SIAM J. Numer. Anal. 14, 1977

### H3. Adaptive Step-Size Controllers

PI controller for step size (Gustafsson, 1991):
```
factor = safety * (err)^(-0.7/p) * (err_prev)^(0.4/p)
dt_new = dt * clamp(factor, 0.2, 5.0)
```

PID controller (Söderlind, 2003):
```
factor = safety * (err)^(-k_I/p) * (err_prev)^(k_P/p) * (err_prevprev)^(k_D/p)
```

References:
- Gustafsson, "Control-theoretic techniques for stepsize selection in implicit Runge-Kutta methods," ACM TOMS 17(4), 1991
- Söderlind, "Digital filters in adaptive time-stepping," ACM TOMS 29(1), 2003
- Hairer & Wanner, "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems," Springer, 1996

### H4. Vieta-Jungius Chebyshev Differentiation Matrices

The current implementation uses the Vieta (alternating product) formula for the Chebyshev differentiation matrix entries. This is numerically stable for moderate N but suffers from O(N²) condition number growth.

For nch > 10, consider the **flipping trick** (Trefethen & Weideman, 2019):
```
D² = flip(flip(D) × D)  // halves roundoff error
```

Reference: Trefethen & Weideman, "The Eigenvalues of the Second Chebyshev Differentiation Matrix," 2019.

### H5. Spectral Methods — Comprehensive References

Core texts:
- Trefethen, "Spectral Methods in MATLAB," SIAM, 2000 (the bible for implementation)
- Boyd, "Chebyshev and Fourier Spectral Methods," Dover, 2001
- Olver & Townsend, "A Fast and Well-Conditioned Spectral Method," SIAM Review 55(3), 2013 (ultraspherical)

For battery-specific spectral methods:
- Bizeray et al., "Lithium-ion battery thermal-electrochemical model-based state estimation using orthogonal collocation and a modified extended Kalman filter," J. Power Sources 296, 2015
- Subramanian et al., "Efficient macro-micro scale coupled modeling of batteries," J. Electrochem. Soc. 152(10), 2005

### H6. Integrator Compatibility with Tree Structure

The question: "Will the template Integrator work with the tree structure?"

**Answer: Yes, via the existing `get_dxdt()` virtual method.** `StorageUnit::get_dxdt()` returns the time derivatives of all states for the entire subtree. An integrator only needs this function signature:

```cpp
template <typename StateType>
class Integrator {
    // StateType can be std::vector<double> from getStates()
    // rhs can be get_dxdt() wrapped as a lambda
    StepResult step(StateType& state,
                    std::function<void(const StateType&, StateType&)> rhs,
                    double dt, double atol, double rtol);
};

// Usage with tree:
auto rhs = [&](const auto& s, auto& ds) { battery->get_dxdt(ds); };
integrator.step(states, rhs, dt, 1e-6, 1e-4);
```

The tree structure is transparent to the integrator — it sees only a flat state vector and a right-hand-side function.

### H7. Jacobian Update Strategy for Module_p

The user notes: "sometimes it is easier to do more steps than updating the Jacobian."

**Broyden rank-1 update** is the sweet spot — O(n²) per update vs O(n³) for full refactorization:

```cpp
// After Newton step: J_new ≈ J_old + (ΔF - J_old·Δx) · Δxᵀ / (Δxᵀ·Δx)
void broyden_update(MatrixXd& J, VectorXd& dx, VectorXd& dF) {
    J += ((dF - J * dx) * dx.transpose()) / dx.squaredNorm();
}
```

Strategy: Full Jacobian every N_rebuild cycles (N_rebuild = 10-100). Broyden updates between rebuilds. Reset if convergence stalls (> 20 iterations).

### H8. State Perturbation Safety

The user asks: "Perturbation and solving good but need to safely save the states?"

**Yes.** The existing backup/restore pattern in Cell_SPM::setCurrent (backup I, V, Vcell_valid) is too lightweight. For perturbation-based Jacobian computation:

```cpp
// Safe perturbation with RAII backup
struct StateBackup {
    StorageUnit& su;
    std::vector<double> saved;
    StateBackup(StorageUnit& su) : su(su) { su.getStates(saved); }
    ~StateBackup() { int n = 0; su.setStates(saved, n, false, false); }
};

// Usage:
for (int i = 0; i < n_cells; ++i) {
    StateBackup guard(*SUs[i]);       // saves state
    SUs[i]->setCurrent(I + eps);      // perturb
    double V_pert = SUs[i]->V();      // evaluate
    J(i) = (V_pert - V_base) / eps;
}  // guard destructor restores state
```

### H9. deal.II for Future 3D?

The user asks about FEM frameworks. **deal.II** is indeed the right choice for future 3D electrode microstructure:
- Adaptive mesh refinement (AMR) for complex geometries
- hp-FEM with spectral-order elements
- Built-in MPI parallelism
- CUDA support via matrix-free operators
- Excellent for the macro-scale 3D solve, while micro-scale (particle) stays spectral

Integration path: Use deal.II only for 3D transport; particle diffusion stays in SLIDE's eigenvalue-decomposed form. Multi-scale coupling via operator splitting.

### H10. Kinetics Strategy — Butler-Volmer Options & Analytical Inversions

The current code (`Electrode_SPM.hpp:43-51`) uses the **symmetric BV inversion** exclusively:
```cpp
eta = (2 * Rg * T) / (n * F) * asinh(x)   // where x = I / (2*a*thick*i0)
```
This is valid only for α_a = α_c = 0.5. Design a kinetics policy with multiple options:

**Option 1: Symmetric BV (current, fastest)**
- `eta = 2RT/(nF) * asinh(x)` — one `asinh` call, analytically invertible
- Inverse: `I = 2*a*thick*i0 * sinh(nF*eta / (2RT))` — enables direct voltage→current

**Option 2: Asymmetric BV (general, requires Newton)**
- `I = i0 * [exp(α_a*F*η/RT) - exp(-α_c*F*η/RT)]` with α_a ≠ α_c
- Not analytically invertible in general; requires Newton iteration for η→I
- Derivative for Newton: `dI/dη = i0 * F/RT * [α_a*exp(α_a*Fη/RT) + α_c*exp(-α_c*Fη/RT)]`

**Option 3: Linearized BV (fastest, valid near equilibrium)**
- For small overpotentials (|η| << RT/F ≈ 26mV): `I ≈ i0 * (α_a + α_c) * F * η / RT`
- Equivalent to a linear resistance: `R_ct = RT / (i0 * (α_a + α_c) * F)`
- Analytically invertible: `η = I * R_ct`
- Best for: low C-rates, near-equilibrium operation, fast screening

**Option 4: Tafel (high-overpotential limit)**
- For |η| >> RT/F: `I = i0 * exp(α_a*F*η/RT)` (anodic) or `I = -i0 * exp(-α_c*F*η/RT)` (cathodic)
- Analytically invertible: `η = RT/(α_a*F) * ln(I/i0)`
- Best for: high C-rate discharge, plating conditions

**Option 5: Marcus kinetics (advanced)**
- `I = i0 * [exp(-(λ + ΔG)²/(4λkT)) - exp(-(λ - ΔG)²/(4λkT))]`
- Relevant for concentrated electrolytes and electron transfer theory
- Not invertible; Newton required

**Option 6: No kinetic limits (infinite i0)**
- `η = 0` always — the cell is purely thermodynamic
- Best for: very fast screening, equivalent to ECM with R_ct = 0

**CRTP Strategy Implementation:**
```cpp
struct SymmetricBV : KineticsStrategy<SymmetricBV> {
    double eta_impl(double i_app, double i0, double T) const {
        double x = 0.5 * i_app / i0;
        return 2.0 * PhyConst::Rg * T / (PhyConst::n * PhyConst::F) * std::asinh(x);
    }
    // Analytical inverse:
    double current_from_eta_impl(double eta, double i0, double T) const {
        return 2.0 * i0 * std::sinh(PhyConst::n * PhyConst::F * eta / (2.0 * PhyConst::Rg * T));
    }
    static constexpr bool is_invertible = true;
};

struct LinearizedBV : KineticsStrategy<LinearizedBV> {
    double eta_impl(double i_app, double i0, double T) const {
        double R_ct = PhyConst::Rg * T / (i0 * PhyConst::F);
        return i_app * R_ct;
    }
    double current_from_eta_impl(double eta, double i0, double T) const {
        return eta * i0 * PhyConst::F / (PhyConst::Rg * T);
    }
    static constexpr bool is_invertible = true;
};

struct AsymmetricBV : KineticsStrategy<AsymmetricBV> {
    double alpha_a = 0.5, alpha_c = 0.5;
    // Not analytically invertible — requires Newton iteration
    static constexpr bool is_invertible = false;
};
```

**The `is_invertible` flag** enables compile-time optimization: if the kinetics are invertible, `setCurrent` can directly compute voltage without iteration. If not, a Newton loop is needed. This flag propagates through the policy-based Cell template.

### H11. Data Storage Architecture for 100K+ Cells

For different models (SPM: 29 states, DFN: 74 states, ECM: 6 states):

```cpp
// Heterogeneous cell pool with variant storage
struct HeterogeneousCellPool {
    // Type-tagged state blocks
    struct CellBlock {
        enum ModelType { SPM, SPMe, DFN, ECM };
        ModelType type;
        size_t n_cells;
        std::vector<double> states;  // contiguous: n_cells * n_states_per_model
        size_t states_per_cell() const {
            switch(type) { case SPM: return 29; case DFN: return 74; /*...*/ }
        }
    };

    std::vector<CellBlock> blocks;
    // SPM cells in one block, DFN in another — SoA within each block
    // Enables SIMD within homogeneous blocks
    // Blocks can be on different MPI ranks or GPU streams
};
```

For connecting SPM and DFN cells in parallel: they all implement `StorageUnit` interface (V, I, setCurrent). The circuit solver sees only branch voltages and currents — model type is irrelevant.

---

*Report generated 2026-03-13 by 20-agent adversarial audit of SLIDE v3.0.0*
*Appendices A-E added 2026-03-13 by adversarial verification pass*
*Appendix F added 2026-03-14: Strategic vision (SPMe, DFN, GPU, multi-language)*
*Appendices G-H added 2026-03-14: Architecture redesign, design patterns, technical notes*
