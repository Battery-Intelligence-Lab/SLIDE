"""Generate the exact Python-wheel trace used by the MATLAB parity gate."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import slide


SLIDE_VERSION = "4.0.0.dev0"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            Path(__file__).resolve().parents[1]
            / "tests"
            / "reference"
            / "slide_python_chen2020_1c_600.csv"
        ),
    )
    args = parser.parse_args()
    if slide.__version__ != SLIDE_VERSION:
        raise RuntimeError(
            f"reference generation requires SLIDE {SLIDE_VERSION}, got {slide.__version__}"
        )

    simulation = slide.Simulation(
        slide.SPM({"nch": 12}),
        experiment=slide.Experiment(
            "Discharge at 1 C for 600 seconds", period="10 seconds"
        ),
    )
    solution = simulation.solve(initial_soc=0.8)
    time = np.asarray(solution.t, dtype=float)
    current = np.asarray(solution["Current [A]"].entries, dtype=float)
    voltage = np.asarray(solution["Voltage [V]"].entries, dtype=float)
    if not (
        time.ndim == 1
        and current.shape == time.shape
        and voltage.shape == time.shape
        and np.all(np.diff(time) > 0)
        and np.all(np.isfinite(voltage))
    ):
        raise RuntimeError("SLIDE returned an invalid reference trace")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("time_s", "current_A", "voltage_V"))
        writer.writerows(
            (format(t, ".17g"), format(i, ".17g"), format(v, ".17g"))
            for t, i, v in zip(time, current, voltage, strict=True)
        )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
