"""PyBOP 25.11 adapter for SLIDE's exact forward-sensitivity surface."""

from __future__ import annotations

from collections.abc import Sequence
from copy import copy
from typing import Any

import numpy as np
import pybop

from . import Simulation

# PyBOP 25.11 still calls the NumPy <2 alias in BaseCost.stack_sensitivities.
# Keep the compatibility seam local to the optional adapter until upstream
# replaces it with vstack (the operations are identical for this use).
if not hasattr(np, "row_stack"):
    np.row_stack = np.vstack


class Simulator(pybop.BaseSimulator):
    """Wrap a fixed-CC :class:`slide.Simulation` as a PyBOP simulator."""

    def __init__(
        self,
        simulation: Simulation,
        parameters: pybop.Parameters | dict[str, pybop.Parameter],
    ) -> None:
        super().__init__(
            parameters if isinstance(parameters, pybop.Parameters) else pybop.Parameters(parameters)
        )
        self.simulation = simulation
        self._output_variables = ["Voltage [V]"]

    @property
    def has_sensitivities(self) -> bool:
        return True

    def set_output_variables(self, target: list[str]) -> None:
        unsupported = set(target) - {"Voltage [V]", "Terminal voltage [V]"}
        if unsupported:
            raise ValueError(f"SLIDE's PyBOP adapter supports voltage, not {sorted(unsupported)}")
        self._output_variables = list(target)

    def batch_solve(
        self,
        inputs: list[dict[str, float]],
        calculate_sensitivities: bool = False,
    ) -> list[pybop.Solution]:
        results: list[pybop.Solution] = []
        names = self.parameters.names
        for values in inputs:
            values = dict(values or {})
            if calculate_sensitivities:
                slide_solution = self.simulation.simulateS1(names, inputs=values)
                sensitivities = slide_solution.sensitivities
            else:
                slide_solution = self.simulation.solve(inputs=values)
                sensitivities = None
            result = pybop.Solution(inputs=values)
            result.set_solution_variable("Time [s]", slide_solution.t)
            voltage = slide_solution["Voltage [V]"].entries
            for variable in self._output_variables:
                result.set_solution_variable(
                    variable,
                    np.asarray(voltage),
                    sensitivities=sensitivities,
                )
            results.append(result)
        return results

    def copy(self) -> "Simulator":
        result = copy(self)
        result.parameters = self.parameters
        return result


def fitting_problem(
    simulation: Simulation,
    parameters: pybop.Parameters | dict[str, pybop.Parameter],
    time: Sequence[float],
    voltage: Sequence[float],
) -> pybop.Problem:
    """Build a voltage SSE problem using the native `simulateS1` gradients."""

    dataset = pybop.Dataset(
        {
            "Time [s]": np.asarray(time, dtype=float),
            "Voltage [V]": np.asarray(voltage, dtype=float),
        },
        variables=["Time [s]", "Voltage [V]"],
    )
    cost = pybop.SumSquaredError(dataset, target="Voltage [V]")
    return pybop.Problem(Simulator(simulation, parameters), cost)


__all__ = ["Simulator", "fitting_problem"]
