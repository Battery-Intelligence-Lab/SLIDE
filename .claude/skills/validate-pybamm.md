# /validate-pybamm - Cross-Validate Against PyBaMM

Compares SLIDE simulation results against PyBaMM for validation.

## Usage

```
/validate-pybamm [--model=SPM|ECM] [--test=TESTNAME]
```

## Background

PyBaMM is a widely-used Python battery modeling framework. Cross-validation ensures SLIDE's physics implementation matches established references.

## Comparison Scripts Location

```
benchmark/python/
├── slide_defaults.py          # Default parameter generation
├── test_single_cell_SPM.py    # SPM single cell validation
├── test_single_cell_ECM.py    # ECM validation (if exists)
└── test_multi_cell.py         # Module validation (if exists)
```

## Steps

### 1. Run SLIDE Simulation

Export results to CSV for comparison:

```cpp
// In test code or main
auto cell = make_Cell_SPM();
Cycler cycler(cell.get(), "pybamm_compare");

// Run protocol
cycler.CC(1.0, 2.7, 3600, 1.0, 10);  // 1C discharge

// Data is saved to results folder
```

### 2. Run PyBaMM Equivalent

```python
# benchmark/python/test_single_cell_SPM.py
import pybamm
import numpy as np

# Create model
model = pybamm.lithium_ion.SPM()
param = pybamm.ParameterValues("Kokam_Marquis2019")

# Run simulation
sim = pybamm.Simulation(model, parameter_values=param)
sim.solve([0, 3600])

# Extract results
t = sim.solution["Time [s]"].entries
V = sim.solution["Voltage [V]"].entries
```

### 3. Compare Results

```python
# Load SLIDE results
slide_data = np.loadtxt("results/pybamm_compare.csv", delimiter=",")
slide_t = slide_data[:, 0]
slide_V = slide_data[:, 1]

# Interpolate to common time points
from scipy.interpolate import interp1d
pybamm_interp = interp1d(t, V)
pybamm_V_at_slide_t = pybamm_interp(slide_t)

# Calculate error
error = np.abs(slide_V - pybamm_V_at_slide_t)
max_error = np.max(error)
mean_error = np.mean(error)

print(f"Max voltage error: {max_error*1000:.2f} mV")
print(f"Mean voltage error: {mean_error*1000:.2f} mV")

# Assert tolerance
assert max_error < 0.01, f"Voltage mismatch: {max_error} V"
```

## Validation Criteria

| Quantity | Tolerance | Notes |
|----------|-----------|-------|
| Voltage | < 10 mV | Over full discharge |
| Capacity | < 1% | Measured at same cutoff |
| Temperature | < 1 K | If thermal model enabled |

## Common Discrepancies

1. **OCV curve differences**: Ensure same OCV data source
2. **Discretization effects**: SLIDE uses Chebyshev (nch=5), PyBaMM configurable
3. **Time stepping**: Different ODE solvers may give slightly different results
4. **Parameter values**: Verify all parameters match (diffusivity, kinetics, geometry)

## Parameter Mapping

Key parameters to align:

| SLIDE | PyBaMM | Description |
|-------|--------|-------------|
| `Dp` | `D_p` | Positive diffusivity |
| `Dn` | `D_n` | Negative diffusivity |
| `kp` | `k_p` | Positive rate constant |
| `kn` | `k_n` | Negative rate constant |
| `Rp` | `R_p` | Positive particle radius |
| `Rn` | `R_n` | Negative particle radius |

## Reporting

Document validation results in `benchmark/README.md`:

```markdown
## PyBaMM Validation (2026-01-29)

| Test | SLIDE Version | PyBaMM Version | Max Error | Status |
|------|---------------|----------------|-----------|--------|
| SPM 1C discharge | v3.0.0 | v24.1 | 5.2 mV | PASS |
| SPM CC-CV | v3.0.0 | v24.1 | 8.1 mV | PASS |
```

## Troubleshooting

- **Large discrepancies**: Check parameter values, especially OCV curves
- **Different capacity**: Verify active material amounts and electrode areas
- **Temperature mismatch**: Check thermal parameters (Cp, convection coefficient)
