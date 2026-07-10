# PyBaMM Technical Analysis Report

**Date:** 2026-03-16
**Purpose:** Detailed technical analysis of PyBaMM (Python Battery Mathematical Modelling) to inform SLIDE's Python bindings and API design.

> **Note on sources:** WebSearch/WebFetch tools were unavailable during generation. This report is based on extensive knowledge of PyBaMM up to v24.x (mid-2025). Key references are listed at the end. The user should verify any details against the latest docs at https://docs.pybamm.org and https://github.com/pybamm-team/PyBaMM.

---

## 1. Architecture & Design Patterns

### High-Level Structure

PyBaMM is organized around a **symbolic computation pipeline**:

```
Model Definition (symbolic) --> Discretization --> Solver --> Solution
```

The key insight is that models are defined as **systems of symbolic equations** (using PyBaMM's own expression tree), then discretized onto meshes, then handed to numerical solvers. This separation of concerns is the central design pattern.

### Package Structure

```
pybamm/
  __init__.py
  models/               # Battery model definitions
    full_battery_models/ # Complete models (SPM, SPMe, DFN, etc.)
      lithium_ion/       # Li-ion specific models
      lead_acid/         # Lead-acid models
    submodels/           # Composable physics submodels
      particle/          # Particle diffusion
      electrolyte_diffusion/
      electrolyte_conductivity/
      electrode/
      interface/         # Reaction kinetics (Butler-Volmer, etc.)
      thermal/           # Thermal models
      sei/               # SEI growth models
      active_material/   # LAM models
      ...
    base_model.py        # BaseModel class
  parameters/            # Parameter handling
    parameter_sets.py    # Registry of parameter sets
    parameter_values.py  # ParameterValues class
    process_parameter_data.py
  geometry/              # Domain geometries
  meshes/                # Mesh generation
  spatial_methods/       # Finite volume, spectral, etc.
  discretisations/       # Discretisation class
  solvers/               # ODE/DAE solvers
    casadi_solver.py     # CasADi-based solver (primary)
    idaklu_solver.py     # SUNDIALS IDA/KLU (fast C++ solver)
    scipy_solver.py      # SciPy wrapper
    base_solver.py
  experiment/            # Experiment protocol definition
    experiment.py        # Experiment class
    step/                # Step definitions (charge, discharge, rest, etc.)
  expression_tree/       # Symbolic math engine
    symbol.py            # Base Symbol class
    binary_operators.py
    unary_operators.py
    concatenations.py
    state_vector.py
    parameter.py
    scalar.py
    variable.py
    ...
  plotting/              # Visualization
  input/                 # Input parameter data files
  simulation.py          # High-level Simulation class
  batch_study.py         # Parameter sweeps
```

### Key Design Patterns

1. **Expression Tree / Symbolic DAG**: All model equations are built as a directed acyclic graph of `Symbol` nodes. This enables automatic differentiation, code generation (to CasADi), and optimization passes before numerical solve.

2. **Submodel Composition**: Models are assembled from interchangeable submodels. Each submodel defines its own variables, equations (algebraic and ODE), boundary conditions, and initial conditions. The parent model collects and merges these.

3. **Options Dictionary Pattern**: Models accept an `options` dict that controls which submodels are selected:
   ```python
   model = pybamm.lithium_ion.SPM(
       options={
           "thermal": "lumped",
           "SEI": "solvent-diffusion limited",
           "loss of active material": "stress-driven",
           "particle": "Fickian diffusion",
       }
   )
   ```

4. **Separation of Symbolic and Numeric**: The model definition phase is purely symbolic. `ParameterValues.process_model()` substitutes numeric values. `Discretisation.process_model()` converts spatial operators to discrete matrices. Only then does a solver touch numerics.

5. **Domain System**: Variables and equations are tagged with domains (`"negative electrode"`, `"separator"`, `"positive electrode"`, `"negative particle"`, `"positive particle"`, etc.) which control mesh assignment and discretization.

### Key Classes

| Class | Module | Role |
|-------|--------|------|
| `Symbol` | `expression_tree.symbol` | Base node in symbolic expression tree |
| `Variable` | `expression_tree.variable` | Unknown variable to solve for |
| `Parameter` | `expression_tree.parameter` | Named parameter (substituted later) |
| `InputParameter` | `expression_tree.input_parameter` | Runtime-changeable parameter |
| `BaseModel` | `models.base_model` | Abstract base for all models |
| `BaseBatteryModel` | `models.full_battery_models.base_battery_model` | Battery-specific base |
| `BaseSubModel` | `models.submodels.base_submodel` | Base for composable submodels |
| `ParameterValues` | `parameters.parameter_values` | Parameter substitution engine |
| `Discretisation` | `discretisations.discretisation` | Spatial discretization |
| `BaseSolver` | `solvers.base_solver` | Abstract solver interface |
| `CasadiSolver` | `solvers.casadi_solver` | Primary production solver |
| `IDAKLUSolver` | `solvers.idaklu_solver` | Fast C++ SUNDIALS solver |
| `Simulation` | `simulation` | High-level orchestrator |
| `Experiment` | `experiment.experiment` | Cycling protocol definition |
| `Solution` | `solvers.solution` | Result container with lazy evaluation |

---

## 2. Model Framework

### Model Hierarchy

```
BaseModel
  BaseBatteryModel
    lithium_ion.BaseModel
      lithium_ion.SPM          # Single Particle Model
      lithium_ion.SPMe         # SPM with electrolyte dynamics
      lithium_ion.DFN          # Doyle-Fuller-Newman (full P2D)
      lithium_ion.MPM          # Many Particle Model
      lithium_ion.MSMR         # Multi-Species Multi-Reaction
    lead_acid.BaseModel
      lead_acid.LOQS           # Leading Order Quasi-Static
      lead_acid.Full
```

### How Model Composition Works

Each `BaseBatteryModel` subclass defines a `set_submodels()` method that assembles physics:

```python
class SPM(BaseModel):
    def set_submodels(self):
        self.submodels["external circuit"] = pybamm.external_circuit.ExplicitCurrentControl(...)
        self.submodels["negative electrode potential"] = pybamm.electrode.ohm.LeadingOrder(...)
        self.submodels["positive electrode potential"] = pybamm.electrode.ohm.LeadingOrder(...)
        self.submodels["negative particle"] = pybamm.particle.FickianDiffusion(...)
        self.submodels["positive particle"] = pybamm.particle.FickianDiffusion(...)
        self.submodels["negative interface"] = pybamm.kinetics.ButlerVolmer(...)
        self.submodels["positive interface"] = pybamm.kinetics.ButlerVolmer(...)
        self.submodels["electrolyte conductivity"] = pybamm.electrolyte_conductivity.LeadingOrder(...)
        # etc...
```

Each submodel independently provides:
- `get_fundamental_variables()` - variables this submodel introduces
- `get_coupled_variables()` - derived variables that depend on other submodels
- `set_rhs()` - ODE right-hand sides (dy/dt = f(y))
- `set_algebraic()` - algebraic constraints (0 = g(y))
- `set_boundary_conditions()` - BCs for spatial operators
- `set_initial_conditions()` - initial values
- `set_events()` - termination events (voltage cutoffs, etc.)

The parent model calls `build_model()` which iterates over all submodels, collects their contributions, and assembles a complete system.

### Symbolic Math (Expression Tree)

PyBaMM has its **own symbolic math engine** (not SymPy). The expression tree:

- Every mathematical expression is a tree of `Symbol` nodes
- Operators: `Addition`, `Multiplication`, `Division`, `Power`, `MatrixMultiplication`
- Functions: `Exponential`, `Log`, `Tanh`, `Arcsinh`, etc.
- Spatial operators: `Gradient`, `Divergence`, `Laplacian`, `BoundaryValue`, `BoundaryGradient`
- Special: `StateVector` (indexes into the solution vector), `Variable`, `Parameter`

Example of how a PDE gets defined symbolically:
```python
c = pybamm.Variable("c", domain="negative particle")
N = -D * pybamm.grad(c)          # Flux
dcdt = -pybamm.div(N)            # Conservation: dc/dt = -div(N)
model.rhs = {c: dcdt}
```

### CasADi Integration

CasADi is used as a **backend for efficient numerical evaluation and automatic differentiation**:

1. PyBaMM's expression tree is converted to CasADi's `SX` or `MX` symbolic types
2. CasADi generates optimized C code for function evaluation and Jacobians
3. The `CasadiSolver` wraps CasADi's integrator interface (which calls SUNDIALS CVODES/IDAS internally)

The conversion happens in `pybamm.CasadiConverter` which walks the PyBaMM expression tree and maps each node to its CasADi equivalent.

### Model Differences (SPM vs SPMe vs DFN)

| Model | Particle | Electrolyte | Electrode Potential | Typical States |
|-------|----------|-------------|---------------------|----------------|
| SPM | Fickian (r-direction) | Uniform | Leading-order | ~20-40 |
| SPMe | Fickian (r-direction) | 1D diffusion + migration | Composite | ~100-200 |
| DFN | Fickian (r-direction, x-coupled) | Full 1D | Full Ohm's law | ~500-2000 |

SPM uses a single representative particle per electrode (like SLIDE). SPMe adds electrolyte concentration and potential variation. DFN is the full pseudo-2D model with x-coupled particles.

---

## 3. Solver Strategy

### Available Solvers

| Solver | Backend | Type | When to Use |
|--------|---------|------|-------------|
| `CasadiSolver` | CasADi + SUNDIALS | ODE/DAE | Default, reliable, good for most models |
| `IDAKLUSolver` | SUNDIALS IDA with KLU sparse solver (C++) | DAE | Fastest for large models (DFN) |
| `ScipySolver` | `scipy.integrate.solve_ivp` | ODE | Simple, good for debugging |
| `CasadiAlgebraicSolver` | CasADi rootfinding | Algebraic | Steady-state / algebraic-only systems |
| `JaxSolver` | JAX | ODE | GPU acceleration, differentiable simulation |

### Stiff System Handling

- **CasadiSolver** uses SUNDIALS CVODES (for ODEs, BDF method) or IDAS (for DAEs, BDF method). BDF is ideal for stiff systems. Adaptive timestepping with error control.
- **IDAKLUSolver** is a custom C++ wrapper around SUNDIALS IDA with the KLU sparse direct linear solver. It operates on the sparse Jacobian structure directly, making it extremely efficient for large sparse systems (DFN models with ~1000+ states).
- Jacobians are computed via **CasADi automatic differentiation** -- no finite differences needed. This is critical for stiff systems where accurate Jacobians dramatically improve convergence.

### Spatial Discretization

PyBaMM uses the **Finite Volume Method** as its primary spatial discretization:

```python
# In discretisations/discretisation.py
disc = pybamm.Discretisation(mesh, spatial_methods)
disc.process_model(model)  # Replaces spatial operators with discrete versions
```

**Spatial methods available:**
- `pybamm.FiniteVolume` -- primary method, used for particle and electrolyte domains
- `pybamm.SpectralVolume` -- higher-order spectral volume method
- Chebyshev spectral methods for particle domain (option)

**Mesh types:**
- `pybamm.Uniform1DSubMesh` -- uniform spacing
- `pybamm.Exponential1DSubMesh` -- exponential clustering near boundaries
- `pybamm.MeshGenerator` -- configurable

The discretization replaces symbolic `grad`, `div` etc. with matrix operations:
- `grad(c)` becomes a difference matrix `G @ c_vec`
- `div(N)` becomes a divergence matrix `D @ N_vec`
- Boundary conditions modify these matrices appropriately

### Event Handling

Solvers support event detection (root-finding) for:
- Voltage cutoffs (min/max)
- Capacity limits
- Temperature limits
- Custom user-defined events

Events are defined as `pybamm.Event(name, expression)` where the expression crosses zero at the event.

---

## 4. Experiment API

### The `Experiment` Class

This is one of PyBaMM's most user-friendly features. It allows defining complex cycling protocols with a readable syntax:

```python
experiment = pybamm.Experiment(
    [
        (
            "Discharge at C/10 for 10 hours or until 3.3 V",
            "Rest for 1 hour",
            "Charge at 1 A until 4.1 V",
            "Hold at 4.1 V until 50 mA",
            "Rest for 1 hour",
        )
    ] * 3,  # Repeat 3 cycles
    period="1 minute",  # Data recording interval
)
```

### Step Types (as of v24.x)

The experiment parser understands natural-language-like strings:

**Current control:**
- `"Discharge at 1C for 1 hour"`
- `"Charge at 2 A until 4.2 V"`
- `"Discharge at C/5 for 30 minutes or until 2.5 V"`

**Voltage control:**
- `"Hold at 4.2 V until C/50"`
- `"Hold at 4.2 V for 1 hour"`
- `"Hold at 4.2 V until 10 mA or for 2 hours"` (compound termination)

**Power control:**
- `"Discharge at 1 W for 2 hours"`
- `"Charge at 0.5 W until 4.2 V"`

**Rest:**
- `"Rest for 30 minutes"`

**Drive cycle:**
- `"Run US06"` (loads predefined drive cycle data)
- `pybamm.step.current(current_data, duration=600)`

### Programmatic Step API (newer, more flexible)

```python
experiment = pybamm.Experiment(
    [
        pybamm.step.c_rate(-1, duration="1 hour"),          # Charge at 1C
        pybamm.step.voltage(4.2, duration="1 hour"),         # Hold at 4.2V
        pybamm.step.rest(duration="30 minutes"),             # Rest
        pybamm.step.c_rate(0.5, termination="2.5 V"),       # Discharge at C/2
    ]
)
```

### How Experiments Execute

The `Simulation` class processes experiments step-by-step:

1. Each step is a separate solver call
2. The final state of one step becomes the initial condition of the next
3. Each step can have different operating conditions (current/voltage/power control changes which equations are active)
4. `InputParameter` values are updated between steps
5. Solutions from all steps are concatenated into a single `Solution` object

### Termination Conditions

- Time-based: `"for 1 hour"`, `"for 3600 seconds"`
- Voltage: `"until 4.2 V"`, `"until 2.5 V"`
- Current: `"until 50 mA"`, `"until C/50"`
- Compound: `"until 50 mA or for 1 hour"` (first condition met wins)

---

## 5. Parameter Management

### ParameterValues Class

```python
param = pybamm.ParameterValues("Chen2020")
# or
param = pybamm.ParameterValues({
    "Nominal cell capacity [A.h]": 5.0,
    "Current function [A]": 5.0,
    "Negative electrode thickness [m]": 85.2e-6,
    # ... hundreds of parameters
})
```

`ParameterValues` is essentially a dictionary that maps parameter name strings to values (scalars, functions, or data interpolants). It processes a model by walking the expression tree and replacing `Parameter` nodes with their numeric values.

### Built-in Parameter Sets

PyBaMM ships parameter sets as Python modules in `pybamm/input/parameters/`:

| Parameter Set | Chemistry | Source |
|--------------|-----------|--------|
| `Chen2020` | NMC/Graphite (LGM50) | Chen et al. 2020 |
| `Marquis2019` | LCO/Graphite | Marquis et al. 2019 |
| `Ecker2015` | NMC/Graphite | Ecker et al. 2015 |
| `Mohtat2020` | NMC/Graphite | Mohtat et al. 2020 |
| `OKane2022` | NMC/Graphite (with degradation) | O'Kane et al. 2022 |
| `Ai2020` | LCO/Graphite | Ai et al. 2020 |
| `Ramadass2004` | LCO/Graphite | Ramadass et al. 2004 |
| `Xu2019` | NMC811/Si-Graphite | Xu et al. 2019 |

As of v24.x, parameter sets moved to a separate package: **`pybamm-parameters`** or are accessible via BPX format.

### Parameter Organization

Each parameter set defines:
- **Cell geometry**: electrode thicknesses, particle radii, separator thickness, etc.
- **Transport properties**: diffusivities D(c,T), conductivities, transference number
- **Thermodynamics**: OCV curves U(SOC) as functions or interpolated data
- **Kinetics**: exchange current density functions j0(c,T), activation energies
- **Thermal**: thermal conductivities, heat capacities, heat transfer coefficients
- **Degradation**: SEI parameters, LAM parameters, lithium plating parameters

Parameters can be:
- **Scalars**: `"Negative electrode thickness [m]": 85.2e-6`
- **Functions**: `"Negative electrode diffusivity [m2.s-1]": D_n_function`  (Python callable taking `(sto, T)`)
- **Data interpolants**: `"Negative electrode OCP [V]": pybamm.Interpolant(sto_data, ocp_data, ...)`

### BPX (Battery Parameter Exchange) Format

PyBaMM supports the BPX JSON standard for parameter exchange:
```json
{
  "Header": {"BPX": 0.3, "Model": "SPM"},
  "Parameterisation": {
    "Cell": {"Nominal cell capacity [A.h]": 5.0, ...},
    "Negative electrode": {"Particle radius [m]": 5.86e-6, ...},
    "Positive electrode": {...},
    "Electrolyte": {...}
  }
}
```

### Parameter Validation

- Parameters are validated at `process_model()` time
- Missing parameters raise `KeyError` with helpful messages
- Function parameters are checked for correct signature (number of arguments)
- Units are embedded in parameter names by convention: `"[m]"`, `"[A.h]"`, `"[m2.s-1]"`
- The BPX schema provides JSON Schema validation

---

## 6. Output / Post-processing

### Solution Object

```python
sim = pybamm.Simulation(model, experiment=experiment, parameter_values=param)
sol = sim.solve()
```

The `Solution` object is the primary output container:

```python
# Access variables by name (lazy evaluation)
voltage = sol["Terminal voltage [V]"]        # Returns ProcessedVariable
time = sol["Time [s]"]
soc = sol["Discharge capacity [A.h]"]
c_n = sol["Negative particle concentration [mol.m-3]"]

# Get numpy arrays
t = sol["Time [s]"].entries                  # 1D array
V = sol["Terminal voltage [V]"].entries       # 1D array
c = sol["Negative particle concentration [mol.m-3]"].entries  # 2D+ array (space x time)

# Evaluate at specific times
V_at_100s = sol["Terminal voltage [V]"](t=100)

# For spatial variables, evaluate at specific positions
c_at_surface = sol["Negative particle concentration [mol.m-3]"](r=1, t=100)
```

### Lazy Evaluation / ProcessedVariable

Variables are **not computed until accessed**. The `Solution` stores the state vector trajectory `y(t)` and the CasADi functions to compute any output variable. When you access `sol["Terminal voltage [V]"]`, it:

1. Looks up the symbolic expression for that variable
2. Evaluates the CasADi function with the stored `y(t)` values
3. Caches the result

This is memory-efficient -- only requested variables are computed.

### Plotting

```python
# Quick plot (built-in)
sim.plot()

# Custom plot
sol.plot(["Terminal voltage [V]", "Current [A]", "Negative electrode SOC"])

# Or use the QuickPlot class
plot = pybamm.QuickPlot(sol, output_variables=["Terminal voltage [V]"])
plot.dynamic_plot()  # Interactive slider for time

# Comparison plots (multiple solutions)
pybamm.QuickPlot([sol1, sol2], labels=["Model A", "Model B"])
```

### Data Export

```python
# To pandas DataFrame
df = sol.get_data_dataframe(
    ["Time [s]", "Terminal voltage [V]", "Current [A]"]
)

# Save to CSV
sol.save_data("output.csv", ["Time [s]", "Terminal voltage [V]"])

# Save full solution (pickle)
sol.save("solution.pkl")
sol = pybamm.load("solution.pkl")
```

### Cycle Summary Variables

For experiments with multiple cycles:
```python
sol.summary_variables["Capacity [A.h]"]          # Per-cycle capacity
sol.summary_variables["Loss of lithium inventory [%]"]
sol.summary_variables["Change in loss of active material [%]"]
```

---

## 7. Extensibility

### Adding a New Submodel

1. Create a class inheriting from `pybamm.BaseSubModel`
2. Implement required methods:

```python
class MyNewSEI(pybamm.BaseSubModel):
    def get_fundamental_variables(self):
        L_sei = pybamm.Variable("SEI thickness [m]", domain="negative electrode")
        return {"SEI thickness [m]": L_sei}

    def get_coupled_variables(self, variables):
        # Compute derived quantities using variables from other submodels
        j_sei = ...
        variables.update({"SEI interfacial current density [A.m-2]": j_sei})
        return variables

    def set_rhs(self, variables):
        L_sei = variables["SEI thickness [m]"]
        dLdt = ...  # Your growth equation
        self.rhs = {L_sei: dLdt}

    def set_initial_conditions(self, variables):
        L_sei = variables["SEI thickness [m]"]
        self.initial_conditions = {L_sei: pybamm.Scalar(5e-9)}
```

3. Register it in the model's options system
4. Select it via `options={"SEI": "my-new-model"}`

### Adding a New Full Model

Inherit from `BaseBatteryModel`, define `set_submodels()`, and select which submodels to use.

### Adding New Parameters

1. Create a Python module or BPX JSON file
2. Define all required parameters
3. Load with `pybamm.ParameterValues("path/to/params.py")` or `pybamm.ParameterValues(bpx="params.json")`

### Custom Experiments

Users can define custom `step` types or use functional forms:
```python
def custom_current(t):
    return 2.0 * np.sin(2 * np.pi * t / 3600)

pybamm.step.current(custom_current, duration="1 hour")
```

### Extensibility Assessment

**Strengths:**
- Extremely modular submodel system -- swap individual physics easily
- Expression tree means new equations "just work" with existing solvers
- Options dictionary pattern is user-friendly
- Adding degradation models is well-supported

**Weaknesses:**
- Adding fundamentally new model structures (not submodel-swappable) requires deep knowledge of the framework
- The symbolic layer has a learning curve for contributors
- Performance-critical additions may need C++ (IDAKLUSolver extensions)

---

## 8. Performance

### Known Bottlenecks

1. **Model creation/processing**: Building the expression tree, discretizing, and converting to CasADi can take 5-30 seconds for complex models. This is a one-time cost per model configuration.

2. **DFN solve time**: The full DFN model with many nodes can be slow (seconds to minutes for a single cycle). The IDAKLUSolver (C++) is 5-10x faster than CasadiSolver for large DAE systems.

3. **Experiment step transitions**: Each step in an experiment requires re-initialization of the solver. For experiments with many short steps, this overhead accumulates.

4. **Memory for long simulations**: Storing full state trajectories for degradation simulations (1000+ cycles) can use significant memory.

### Performance Strategies

1. **IDAKLUSolver**: Written in C++ with pybind11, wraps SUNDIALS IDA with KLU sparse direct solver. Exploits Jacobian sparsity. 5-10x faster than pure Python/CasADi path for DFN.

2. **JaxSolver**: Enables GPU acceleration and JIT compilation. Useful for parameter fitting (differentiable simulation via JAX autodiff).

3. **Model simplification**: SPM is ~100x faster than DFN. SPMe is a good middle ground. Use the simplest model that captures needed physics.

4. **Adaptive timestepping**: All solvers use adaptive step size control. BDF methods with error estimation.

5. **Batch studies**: `pybamm.BatchStudy` runs multiple parameter combinations. Parallelism via multiprocessing or concurrent.futures.

### Parallelism

- **Within a single simulation**: SUNDIALS solvers can use BLAS-level parallelism (OpenBLAS/MKL) for linear algebra. No explicit thread-level parallelism in PyBaMM's Python layer.
- **Across simulations**: `BatchStudy` and user-level parallelism via multiprocessing. No shared-memory concerns because each simulation is independent.
- **liionpack** (separate package): Provides pack-level simulation by running many cell models in parallel using MPI or multiprocessing, communicating via electrical network constraints.
- **GPU**: JaxSolver can offload to GPU. Experimental but promising for parameter sweeps.

### Typical Performance Numbers (approximate)

| Model | 1 cycle solve | 1000 cycles | Notes |
|-------|--------------|-------------|-------|
| SPM (CasADi) | ~0.1-0.5s | ~2-8 min | Fast, good for degradation |
| SPMe (CasADi) | ~0.5-2s | ~10-30 min | |
| DFN (CasADi) | ~2-10s | hours | Use IDAKLUSolver |
| DFN (IDAKLU) | ~0.5-2s | ~15-30 min | Recommended for DFN |
| SPM (Jax, GPU) | ~0.05s | ~1 min | After JIT warmup |

---

## 9. Testing & CI

### Test Framework

- **pytest** as the test runner
- Tests in `tests/unit/` and `tests/integration/`
- Extensive use of `pytest.mark.parametrize` for testing multiple model configurations

### Test Structure

```
tests/
  unit/
    test_models/
      test_full_battery_models/
        test_lithium_ion/
          test_spm.py
          test_spme.py
          test_dfn.py
      test_submodels/
    test_solvers/
    test_parameters/
    test_discretisations/
    test_expression_tree/
    test_experiment/
    test_simulation.py
  integration/
    test_models/      # Full solve + compare to known solutions
```

### CI Pipeline

- **GitHub Actions** for CI
- Runs on: Ubuntu, macOS, Windows
- Python versions: 3.9-3.12+
- Test matrix includes:
  - With/without optional dependencies (CasADi, JAX, IDAKLU)
  - Different solver backends
  - Linting (ruff)
  - Type checking (mypy, partially)
  - Documentation builds (Sphinx)
  - Coverage reporting (Codecov)
- **Nox** as the task runner (replacement for tox)
- Pre-commit hooks for formatting

### Coverage

- Aim for high coverage (>90% on core modules)
- Integration tests solve full models and compare to reference solutions
- Regression tests catch numerical drift

---

## 10. Community & Ecosystem

### Core Team & Governance

- Developed at University of Oxford (initially), now multi-institutional
- Part of the **NumFOCUS** affiliated projects
- Lead developers: Valentin Sulzer (original author), Robert Timms, Scott Marquis, Martin Robinson, Ferran Brosa Planella, and others
- Open governance model with regular contributor meetings

### Key Publications

1. **Sulzer et al. (2021)**: "Python Battery Mathematical Modelling (PyBaMM)" -- JORS paper, the primary citation
2. **Sulzer et al. (2019)**: "Faster Lead-Acid Battery Simulations from Porous-Electrode Theory" -- JES
3. **Marquis et al. (2019)**: "An Asymptotic Derivation of a Single Particle Model with Electrolyte" -- JES
4. **O'Kane et al. (2022)**: Degradation modeling paper

### Ecosystem Packages

| Package | Purpose | Relationship |
|---------|---------|-------------|
| **liionpack** | Pack-level simulation | Uses PyBaMM cells in a network model |
| **PyBOP** | Bayesian Parameter Optimization | Parameter fitting using PyBaMM models |
| **pybamm-cookiecutter** | Template for creating PyBaMM submodel packages | Standardized extension pattern |
| **BPX** | Battery Parameter Exchange standard | JSON schema for parameter sharing |
| **tec-reduced-order-models** | Thermal-electrochemical ROMs | Built on PyBaMM |
| **pybamm-eis** | Electrochemical Impedance Spectroscopy | EIS simulation using PyBaMM |

### Integration Points

- **CasADi**: Core dependency for symbolic-to-numeric conversion and AD
- **SUNDIALS (CVODES/IDAS)**: ODE/DAE solvers (via CasADi or IDAKLU)
- **JAX**: Optional, for GPU acceleration and differentiable simulation
- **Matplotlib**: Plotting
- **pandas/NumPy**: Data handling
- **SciPy**: Fallback solver, interpolation
- **pybind11**: C++ extensions (IDAKLUSolver)

### liionpack Details (Pack-Level Simulation)

Since SLIDE also does pack-level simulation, this is particularly relevant:

```python
import liionpack as lp

# Define a network (e.g., 4s3p = 12 cells)
netlist = lp.setup_circuit(Np=3, Ns=4, Rb=1e-3, Rc=1e-3, Ri=1e-3)

# Use PyBaMM model for each cell
parameter_values = pybamm.ParameterValues("Chen2020")

# Can vary parameters per cell
parameter_values_list = [param.copy() for _ in range(12)]
for i, pv in enumerate(parameter_values_list):
    pv["Current function [A]"] = ...  # Varies per cell

# Solve
output = lp.solve(
    netlist=netlist,
    parameter_values=parameter_values_list,
    experiment=experiment,
    output_variables=["Terminal voltage [V]", "Current [A]"],
    initial_soc=0.5,
)
```

liionpack solves the electrical network (Kirchhoff's laws) at each timestep, distributing currents to individual PyBaMM cell models. It supports:
- Series/parallel configurations
- Cell-to-cell variability (different parameters per cell)
- Thermal coupling (basic)
- MPI parallelism for many cells

---

## 11. Key Takeaways for SLIDE

### What SLIDE Can Learn from PyBaMM

1. **Experiment API is excellent** -- the natural-language string parsing and programmatic `step` API are very user-friendly. SLIDE's Python bindings should aim for similar ergonomics.

2. **Options dictionary pattern** -- selecting model variants via a simple dict is much better than compile-time configuration. SLIDE could adopt this for its Python interface.

3. **Solution object with lazy evaluation** -- storing the state trajectory and computing output variables on demand is memory-efficient and flexible.

4. **Submodel composition** -- PyBaMM's approach of assembling models from swappable submodels is elegant but comes at the cost of complexity. SLIDE's compile-time approach is faster but less flexible at runtime.

5. **Parameter management via BPX** -- adopting the BPX standard would make SLIDE interoperable with PyBaMM's parameter ecosystem.

6. **IDAKLUSolver shows C++ is needed for performance** -- PyBaMM's fastest solver is a C++ extension. This validates SLIDE's approach of having the core in C++ with Python bindings.

### Where SLIDE Has Advantages

1. **Raw performance**: C++ core is inherently faster than PyBaMM's Python + CasADi codegen path. SLIDE should be 10-100x faster for equivalent models.

2. **Pack-level simulation**: SLIDE has native Module/Battery hierarchy. PyBaMM needs the separate liionpack package with network solvers.

3. **Memory efficiency**: C++ arrays vs Python objects. Critical for 1000+ cycle degradation studies.

4. **No symbolic overhead**: SLIDE compiles equations directly. No expression tree build + discretize + codegen pipeline.

### Recommended Python Binding Strategy

Based on this analysis, SLIDE's Python bindings should:

1. **Mimic the Experiment API**: Support both string-based and programmatic step definitions
2. **Expose a Solution-like object**: With named variable access, numpy integration, and lazy evaluation where possible
3. **Support BPX parameter format**: For interoperability
4. **Provide a Simulation class**: High-level orchestrator similar to `pybamm.Simulation`
5. **Use pybind11 or nanobind**: As planned -- this is the same approach PyBaMM uses for IDAKLUSolver

---

## Key References & URLs

- **GitHub Repository**: https://github.com/pybamm-team/PyBaMM
- **Documentation**: https://docs.pybamm.org/en/latest/
- **API Reference**: https://docs.pybamm.org/en/latest/source/api/index.html
- **Examples/Tutorials**: https://docs.pybamm.org/en/latest/source/examples/index.html
- **JORS Paper**: Sulzer et al., "Python Battery Mathematical Modelling (PyBaMM)", Journal of Open Research Software, 9(1), 14, 2021. DOI: 10.5334/jors.309
- **PyPI**: https://pypi.org/project/pybamm/
- **BPX Standard**: https://github.com/FaradayInstitution/BPX
- **liionpack**: https://github.com/pybamm-team/liionpack
- **PyBOP**: https://github.com/pybamm-team/PyBOP
- **CasADi**: https://web.casadi.org/
- **SUNDIALS**: https://computing.llnl.gov/projects/sundials

> **Disclaimer**: This analysis is based on knowledge up to mid-2025. PyBaMM is under active development and details may have changed. Verify against current documentation before making design decisions.
