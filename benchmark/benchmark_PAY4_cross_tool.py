"""Reproducible PAY-4 cross-tool positioning benchmark.

The default driver compares a Release SLIDE build with PyBaMM 26.x/IDAKLU and
launches this file under a separate liionpack environment for the pack cases.
All subprocesses emit one JSON object so that results can be archived without
scraping progress bars or human-readable logs.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def _timed(function):
    begin = time.perf_counter()
    value = function()
    return time.perf_counter() - begin, value


def _json_process(command: list[str]) -> dict[str, Any]:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    try:
        value = json.loads(completed.stdout)
    except json.JSONDecodeError:
        value = None
    if isinstance(value, dict):
        return value
    for line in reversed(completed.stdout.splitlines()):
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            return value
    raise RuntimeError(
        f"no JSON object from {command!r}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )


def _pybamm_case(repetitions: int) -> dict[str, Any]:
    import pybamm

    if not pybamm.__version__.startswith("26."):
        raise RuntimeError(f"PAY-4 requires PyBaMM 26.x, found {pybamm.__version__}")
    instructions = [
        "Discharge at 1 C until 2.5 V",
        "Charge at 1 C until 4.2 V",
        "Hold at 4.2 V until C/20",
    ]

    def build():
        model = pybamm.lithium_ion.SPM()
        parameters = pybamm.ParameterValues("Chen2020")
        experiment = pybamm.Experiment(instructions, period="10 seconds")
        solver = pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8)
        simulation = pybamm.Simulation(
            model,
            parameter_values=parameters,
            experiment=experiment,
            solver=solver,
        )
        simulation.build_for_experiment(initial_soc=1.0)
        return simulation

    setup_s, simulation = _timed(build)
    first_s, solution = _timed(
        lambda: simulation.solve(initial_soc=1.0, showprogress=False)
    )
    warm_times: list[float] = []
    for _ in range(repetitions):
        elapsed, solution = _timed(
            lambda: simulation.solve(initial_soc=1.0, showprogress=False)
        )
        warm_times.append(elapsed)
    return {
        "tool": "pybamm",
        "version": pybamm.__version__,
        "solver": "IDAKLUSolver",
        "rtol": 1e-6,
        "atol": 1e-8,
        "setup_s": setup_s,
        "first_solve_s": first_s,
        "warm_solve_median_s": statistics.median(warm_times),
        "warm_solve_min_s": min(warm_times),
        "warm_solve_max_s": max(warm_times),
        "samples": len(solution.t),
        "duration_s": float(solution.t[-1]),
        "terminal_voltage_V": float(solution["Terminal voltage [V]"].entries[-1]),
    }


def _liionpack_worker(series: int, parallel: int, steps: int) -> dict[str, Any]:
    # Imported only in the isolated worker: current liionpack pins an older PyBaMM.
    import liionpack as lp
    import numpy as np
    import pandas
    import pybamm

    holder: dict[str, Any] = {}

    def setup():
        holder["netlist"] = lp.setup_circuit(
            Np=parallel,
            Ns=series,
            Ri=1e-4,
            Rc=1e-10,
            Rb=1e-10,
            Rt=1e-10,
            I=5.0 * parallel,
            V=4.0,
        )
        holder["experiment"] = pybamm.Experiment(
            [f"Discharge at {5.0 * parallel:g} A for {steps * 10} seconds"],
            period="10 seconds",
        )
        holder["parameters"] = pybamm.ParameterValues("Chen2020")
        holder["manager"] = lp.CasadiManager()
        holder["kwargs"] = {
            "netlist": holder["netlist"],
            "sim_func": None,
            "parameter_values": holder["parameters"],
            "experiment": holder["experiment"],
            "inputs": None,
            "output_variables": None,
            "initial_soc": 0.8,
            "nproc": 1,
            "node_termination_func": None,
        }
        holder["manager"].solve(**holder["kwargs"], setup_only=True)

    def solve():
        manager = holder["manager"]
        manager.global_step = 0
        for index, protocol in enumerate(manager.protocol_steps):
            termination = manager.terminations[index]
            if termination == []:
                termination = 0.0
            manager._step_solve_step(
                protocol,
                termination,
                manager.step_types[index],
                None,
            )
        return manager.step_output()

    # liionpack writes progress bars to both streams even with its logger muted.
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        setup_s, _ = _timed(setup)
        solve_s, output = _timed(solve)
    voltage = np.asarray(output["Pack terminal voltage [V]"])
    return {
        "tool": "liionpack",
        "version": getattr(lp, "__version__", "unknown"),
        "pybamm_version": pybamm.__version__,
        "pandas_version": pandas.__version__,
        "manager": "CasadiManager",
        "series": series,
        "parallel": parallel,
        "cells": series * parallel,
        "steps_requested": steps,
        "samples": int(voltage.size),
        "setup_s": setup_s,
        "solve_s": solve_s,
        "terminal_voltage_V": float(voltage[-1]),
        "connector_resistance_ohm": 1e-10,
        "cell_resistance_ohm": 1e-4,
    }


def _default_executable(name: str) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    config = "Release" if os.name == "nt" else ""
    return ROOT / "bin" / config / f"{name}{suffix}"


def _driver(args: argparse.Namespace) -> dict[str, Any]:
    slide_single = _json_process(
        [str(args.slide_single), "--repetitions", str(args.repetitions)]
    )
    pybamm_result = _pybamm_case(args.repetitions)
    single = {
        "slide": slide_single,
        "pybamm": pybamm_result,
        "per_solve_speedup": (
            pybamm_result["warm_solve_median_s"] / slide_single["solve_median_s"]
        ),
        "cold_setup_plus_solve_speedup": (
            (pybamm_result["setup_s"] + pybamm_result["first_solve_s"])
            / (slide_single["setup_s"] + slide_single["solve_median_s"])
        ),
    }

    packs = []
    for series, parallel in ((16, 4), (1, 100)):
        slide_pack = _json_process(
            [
                str(args.slide_pack),
                "--series",
                str(series),
                "--parallel",
                str(parallel),
                "--steps",
                str(args.pack_steps),
                "--repetitions",
                str(args.repetitions),
            ]
        )
        liionpack_result = _json_process(
            [
                str(args.liionpack_python),
                str(Path(__file__).resolve()),
                "--liionpack-worker",
                "--series",
                str(series),
                "--parallel",
                str(parallel),
                "--pack-steps",
                str(args.pack_steps),
            ]
        )
        packs.append(
            {
                "slide": slide_pack,
                "liionpack": liionpack_result,
                "solve_speedup": (
                    liionpack_result["solve_s"] / slide_pack["solve_median_s"]
                ),
                "conservative_solve_speedup": (
                    liionpack_result["solve_s"] / slide_pack["solve_max_s"]
                ),
                "setup_speedup": (
                    liionpack_result["setup_s"] / slide_pack["setup_s"]
                ),
                "terminal_voltage_difference_V": abs(
                    liionpack_result["terminal_voltage_V"]
                    - slide_pack["terminal_voltage_V"]
                ),
            }
        )

    return {
        "benchmark": "PAY-4",
        "timing_qualification": (
            "Targets are positioning measurements, not gates; rerun on a quiet machine. "
            "liionpack records all outputs during solve while the SLIDE pack harness advances state only."
        ),
        "environment": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "protocol": {
            "single": [
                "Discharge at 1 C until 2.5 V",
                "Charge at 1 C until 4.2 V",
                "Hold at 4.2 V until C/20",
            ],
            "pack": f"600 s equivalent: {args.pack_steps} x 10 s CC at 5 A/cell",
            "parameters": "Chen2020; nch=12 in SLIDE; initial SOC 1.0 single, 0.8 pack",
        },
        "single": single,
        "packs": packs,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--liionpack-worker", action="store_true")
    parser.add_argument("--series", type=int, default=16)
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--pack-steps", type=int, default=60)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument(
        "--slide-single",
        type=Path,
        default=_default_executable("benchmark_PAY4_slide_single"),
    )
    parser.add_argument(
        "--slide-pack",
        type=Path,
        default=_default_executable("benchmark_PAY4_slide_pack"),
    )
    parser.add_argument("--liionpack-python", type=Path)
    parser.add_argument("--output", type=Path)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.liionpack_worker:
        result = _liionpack_worker(args.series, args.parallel, args.pack_steps)
    else:
        if args.liionpack_python is None:
            raise SystemExit("--liionpack-python is required by the PAY-4 driver")
        result = _driver(args)
    encoded = json.dumps(result, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.write_text(encoded + "\n", encoding="utf-8")
    print(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
