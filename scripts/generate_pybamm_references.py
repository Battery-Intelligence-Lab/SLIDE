"""Regenerate the committed Phase-7 PyBaMM 26.6.2.0 SPM references.

This script is intentionally not part of C++ CI. It pins and asserts the
upstream version so a deliberate regeneration exposes upstream drift instead
of silently moving the parity oracle.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pybamm


PYBAMM_VERSION = "26.6.2.0"
TUTORIAL_STEPS = [
    "Discharge at C/10 for 10 hours or until 3.3 V",
    "Rest for 1 hour",
    "Charge at 1 A until 4.1 V",
    "Hold at 4.1 V until 50 mA",
    "Rest for 1 hour",
]


def write_trace(path: Path, solution: pybamm.Solution) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("segment", "time_s", "local_time_s", "current_A", "voltage_V"))
        for segment, cycle in enumerate(solution.cycles):
            time = np.asarray(cycle["Time [s]"].entries, dtype=float)
            voltage = np.asarray(cycle["Terminal voltage [V]"].entries, dtype=float)
            current = np.asarray(cycle["Current [A]"].entries, dtype=float)
            if (
                time.ndim != 1
                or voltage.shape != time.shape
                or current.shape != time.shape
                or not np.all(np.diff(time) > 0)
            ):
                raise RuntimeError("PyBaMM returned an invalid reference segment")
            local_time = time - time[0]
            writer.writerows(
                (
                    segment,
                    format(t, ".17g"),
                    format(local_t, ".17g"),
                    format(i, ".17g"),
                    format(v, ".17g"),
                )
                for t, local_t, i, v in zip(time, local_time, current, voltage, strict=True)
            )


def solve(steps: list[str], period_seconds: float) -> pybamm.Solution:
    model = pybamm.lithium_ion.SPM()
    experiment = pybamm.Experiment(steps, period=f"{period_seconds:g} seconds")
    parameters = pybamm.ParameterValues("Chen2020")
    solver = pybamm.IDAKLUSolver(rtol=1e-10, atol=1e-12)
    simulation = pybamm.Simulation(
        model,
        experiment=experiment,
        parameter_values=parameters,
        solver=solver,
    )
    return simulation.solve(initial_soc=1.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tests" / "reference",
    )
    args = parser.parse_args()
    if pybamm.__version__ != PYBAMM_VERSION:
        raise RuntimeError(
            f"reference generation requires PyBaMM {PYBAMM_VERSION}, got {pybamm.__version__}"
        )
    args.output.mkdir(parents=True, exist_ok=True)
    scenarios = {
        "pybamm_chen2020_c050.csv": (["Discharge at C/50 for 36000 seconds"], 60.0),
        "pybamm_chen2020_c1.csv": (["Discharge at 1 C for 3000 seconds"], 10.0),
        "pybamm_chen2020_tutorial5.csv": (TUTORIAL_STEPS, 60.0),
    }
    for filename, (steps, period) in scenarios.items():
        write_trace(args.output / filename, solve(steps, period))
        print(f"wrote {args.output / filename}")


if __name__ == "__main__":
    main()
