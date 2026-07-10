"""Compare regenerated PyBaMM references without platform-sensitive byte diffs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("generated", type=Path)
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tests" / "reference",
    )
    args = parser.parse_args()

    names = (
        "pybamm_chen2020_c050.csv",
        "pybamm_chen2020_c1.csv",
        "pybamm_chen2020_tutorial5.csv",
    )
    limits = np.asarray([0.0, 1e-6, 1e-6, 1e-10, 1e-9])
    failed = False
    for name in names:
        reference = np.loadtxt(args.reference / name, delimiter=",", skiprows=1)
        generated = np.loadtxt(args.generated / name, delimiter=",", skiprows=1)
        if reference.shape != generated.shape:
            print(f"{name}: shape drift {reference.shape} -> {generated.shape}")
            failed = True
            continue
        errors = np.max(np.abs(reference - generated), axis=0)
        print(
            f"{name}: max segment/time/local/current/voltage drift "
            + ", ".join(f"{value:.3g}" for value in errors)
        )
        if np.any(errors > limits):
            failed = True
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
