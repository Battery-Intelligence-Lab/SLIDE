---
layout: default
title: MATLAB quickstart
nav_order: 4
---

# MATLAB quickstart

Build `slide_mex` first, then run this from the repository root. If `SLIDE_ROOT` is set, the script uses that checkout instead of `pwd`.

<!-- doc-test:matlab -->
```matlab
slideRoot = getenv("SLIDE_ROOT");
if strlength(slideRoot) == 0
    slideRoot = pwd;
end
addpath(fullfile(slideRoot, "matlab"));

parameters = slide.ParameterValues("Chen2020");
experiment = slide.Experiment( ...
    "Discharge at 1 C for 60 seconds", period="10 seconds");
simulation = slide.Simulation( ...
    slide.lithium_ion.SPM(struct("nch", 12)), ...
    experiment=experiment, parameter_values=parameters);
solution = solve(simulation, initial_soc=0.8);
voltage = variable(solution, "Terminal voltage [V]").entries;

assert(isequal(size(solution.t), [7, 1]));
assert(isequal(size(voltage), [7, 1]));
assert(all(isfinite(voltage)));
fprintf("samples=%d final_voltage=%.6f\n", numel(solution.t), voltage(end));
```

The `+slide` package mirrors the Python names where MATLAB syntax permits: `Experiment`, `ParameterValues`, `lithium_ion.SPM`, `Simulation`, `Solution`, processed variables, `varied`, and device discovery through `slide.availableDevices()`.
