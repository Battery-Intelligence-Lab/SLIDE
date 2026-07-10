---
layout: default
title: Python quickstart
nav_order: 3
---

# Python quickstart

The Python API deliberately follows PyBaMM's setup shape while running SLIDE's compiled spectral core.

<!-- doc-test:python -->
```python
import numpy as np

import slide

parameters = slide.ParameterValues("Chen2020")
experiment = slide.Experiment(
    "Discharge at 1 C for 60 seconds",
    period="10 seconds",
)
simulation = slide.Simulation(
    slide.lithium_ion.SPM({"nch": 12}),
    experiment=experiment,
    parameter_values=parameters,
)
solution = simulation.solve(initial_soc=0.8)
voltage = solution["Terminal voltage [V]"].entries

assert solution.t.shape == (7,)
assert voltage.shape == (7,)
assert np.isfinite(voltage).all()
print(f"samples={solution.t.size} final_voltage={voltage[-1]:.6f}")
```

`Solution` variables expose `.entries` and interpolation. Use `solution.plot(...)` or `slide.plot(...)` only after installing the optional plotting dependency. Read the [compatibility table](compatibility.html#pybamm-shaped-api-known-gaps) before passing arbitrary PyBaMM models, meshes, solvers, or callbacks.
