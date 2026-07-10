from __future__ import annotations

import numpy as np
import pytest

import slide

pybop = pytest.importorskip("pybop")
from slide.pybop import fitting_problem


D_NAME = "Negative particle diffusivity [m2.s-1]"
R_NAME = "Contact resistance [Ohm]"


def test_p7g2_pybop_gradient_fit_recovers_synthetic_truth():
    truth = {D_NAME: 4.0e-14, R_NAME: 1.5e-3}
    experiment = slide.Experiment(
        "Discharge at 1 C for 1800 seconds", period="10 seconds"
    )
    truth_parameters = slide.ParameterValues("Chen2020")
    truth_parameters.update(truth)
    truth_solution = slide.Simulation(
        slide.SPM({"nch": 12}),
        experiment=experiment,
        parameter_values=truth_parameters,
    ).solve(initial_soc=0.8)

    fit_parameters = pybop.Parameters(
        {
            D_NAME: pybop.Parameter(
                initial_value=2.5e-14,
                bounds=[1.0e-14, 8.0e-14],
                transformation=pybop.ScaledTransformation(1e15),
            ),
            R_NAME: pybop.Parameter(
                initial_value=2.5e-3,
                bounds=[1.0e-5, 5.0e-3],
                transformation=pybop.ScaledTransformation(1e3),
            ),
        }
    )
    fit_values = slide.ParameterValues("Chen2020")
    fit_values["Initial state-of-charge"] = 0.8
    simulation = slide.Simulation(
        slide.SPM({"nch": 12}), experiment=experiment, parameter_values=fit_values
    )
    problem = fitting_problem(
        simulation,
        fit_parameters,
        truth_solution.t,
        truth_solution["Voltage [V]"].entries,
    )
    optimiser = pybop.SciPyMinimize(
        problem,
        options=pybop.SciPyMinimizeOptions(
            method="L-BFGS-B",
            jac=True,
            maxiter=500,
            tol=1e-15,
            solver_options={"ftol": 1e-15, "gtol": 1e-12, "maxls": 100},
        ),
    )
    result = optimiser.run()
    fitted = result.best_inputs
    assert abs(fitted[D_NAME] / truth[D_NAME] - 1.0) <= 0.02
    assert abs(fitted[R_NAME] / truth[R_NAME] - 1.0) <= 0.005
