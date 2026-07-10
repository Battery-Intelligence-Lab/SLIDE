"""PyBaMM-shaped Python interface to the compiled SLIDE v4 core."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, MutableMapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any
import csv
import re

import numpy as np

from . import _slide_core

__all__ = [
    "Experiment",
    "ParameterValues",
    "ProcessedVariable",
    "Simulation",
    "Solution",
    "SPM",
    "available_devices",
    "lithium_ion",
    "plot",
    "varied",
]
__version__ = "4.0.0.dev0"


_TIME_UNITS = {
    "s": 1.0,
    "sec": 1.0,
    "second": 1.0,
    "seconds": 1.0,
    "min": 60.0,
    "minute": 60.0,
    "minutes": 60.0,
    "h": 3600.0,
    "hr": 3600.0,
    "hour": 3600.0,
    "hours": 3600.0,
}

_PARAMETER_ALIASES = {
    "Negative electrode diffusivity [m2.s-1]": "Negative particle diffusivity [m2.s-1]",
    "Positive electrode diffusivity [m2.s-1]": "Positive particle diffusivity [m2.s-1]",
    "Exchange-current density for lithium plating [A.m-2]": (
        "Exchange-current density for lithium metal electrode [A.m-2]"
    ),
    "1 + dlnf/dlnc": "Thermodynamic factor",
}


def _parameter_name(name: str) -> str:
    return _PARAMETER_ALIASES.get(name, name)


def _seconds(value: float | str) -> float:
    if isinstance(value, (int, float)):
        result = float(value)
    else:
        match = re.fullmatch(
            r"\s*([+]?(?:\d+(?:\.\d*)?|\.\d+))\s*([A-Za-z]+)\s*", value
        )
        if match is None or match.group(2).lower() not in _TIME_UNITS:
            raise ValueError(f"invalid time interval: {value!r}")
        result = float(match.group(1)) * _TIME_UNITS[match.group(2).lower()]
    if not np.isfinite(result) or result <= 0:
        raise ValueError("time interval must be finite and positive")
    return result


@dataclass(frozen=True)
class Varied:
    """A cold per-lane parameter vector used by an ensemble simulation."""

    values: np.ndarray

    def __post_init__(self) -> None:
        values = np.asarray(self.values, dtype=float)
        if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
            raise ValueError("varied values must be a finite, non-empty 1-D array")
        object.__setattr__(self, "values", values.copy())


def varied(values: Iterable[float]) -> Varied:
    """Mark explicit values as lanes of one SoA parameter-study batch."""

    return Varied(np.asarray(list(values), dtype=float))


class ParameterValues(MutableMapping[str, Any]):
    """PyBaMM-named SI parameter values backed by Chen2020 or a BPX file."""

    def __init__(self, values: str | Path | MutableMapping[str, Any] = "Chen2020"):
        if isinstance(values, MutableMapping):
            self.source = "Chen2020"
            self._base = dict(_slide_core.parameter_values(self.source))
            self._overrides: dict[str, Any] = {}
            self.update(values)
        else:
            self.source = str(values)
            self._base = dict(_slide_core.parameter_values(self.source))
            self._overrides = {}

    @classmethod
    def create_from_bpx(cls, filename: str | Path) -> "ParameterValues":
        return cls(filename)

    def __getitem__(self, key: str) -> Any:
        key = _parameter_name(key)
        return self._overrides[key] if key in self._overrides else self._base[key]

    def __setitem__(self, key: str, value: Any) -> None:
        key = _parameter_name(key)
        if isinstance(value, Varied):
            self._overrides[key] = value
            return
        if isinstance(value, dict) and set(value) == {"x", "y"}:
            x = np.asarray(value["x"], dtype=float)
            y = np.asarray(value["y"], dtype=float)
            if x.ndim != 1 or x.size < 2 or x.shape != y.shape:
                raise ValueError(f"invalid curve parameter {key!r}")
            self._overrides[key] = {"x": x.tolist(), "y": y.tolist()}
            return
        number = float(value)
        if not np.isfinite(number):
            raise ValueError(f"parameter {key!r} must be finite")
        self._overrides[key] = number

    def __delitem__(self, key: str) -> None:
        key = _parameter_name(key)
        if key in self._overrides:
            del self._overrides[key]
        elif key in self._base:
            raise KeyError("base parameter entries cannot be removed; override them instead")
        else:
            raise KeyError(key)

    def __iter__(self) -> Iterator[str]:
        return iter(dict.fromkeys((*self._base, *self._overrides)))

    def __len__(self) -> int:
        return len(set(self._base) | set(self._overrides))

    def copy(self) -> "ParameterValues":
        result = object.__new__(ParameterValues)
        result.source = self.source
        result._base = self._base.copy()
        result._overrides = self._overrides.copy()
        return result

    def scalar_overrides(self) -> dict[str, float]:
        result: dict[str, float] = {}
        for name, value in self._overrides.items():
            if isinstance(value, Varied):
                continue
            if isinstance(value, dict):
                raise NotImplementedError("Python curve overrides are not yet accepted by the native boundary")
            result[name] = float(value)
        return result

    def varied_overrides(self) -> dict[str, np.ndarray]:
        return {
            name: value.values
            for name, value in self._overrides.items()
            if isinstance(value, Varied)
        }


class Experiment:
    """A sequence of PyBaMM experiment instruction strings."""

    def __init__(
        self,
        operating_conditions: str | Sequence[str | Sequence[str]],
        period: float | str = "1 minute",
        *,
        temperature: float | str | None = None,
        termination: str | Sequence[str] | None = None,
    ) -> None:
        del temperature, termination  # accepted for source compatibility; model options own them
        if isinstance(operating_conditions, str):
            operating_conditions = [operating_conditions]
        self.steps = self._flatten(operating_conditions)
        if not self.steps:
            raise ValueError("an experiment needs at least one operating condition")
        self.period = _seconds(period)

    @classmethod
    def _flatten(cls, values: Sequence[str | Sequence[str]]) -> list[str]:
        result: list[str] = []
        for value in values:
            if isinstance(value, str):
                result.append(value)
            else:
                result.extend(cls._flatten(value))
        return result


class SPM:
    """Precompiled single-particle model selection."""

    def __init__(self, options: MutableMapping[str, Any] | None = None, name: str = "SPM"):
        self.options = dict(options or {})
        self.name = name


lithium_ion = SimpleNamespace(SPM=SPM)


def _mask(value: Any, choices: dict[str, int], name: str) -> int:
    if value in (None, False, "none", "false"):
        return 0
    if isinstance(value, bool):
        return 1 if value else 0
    key = str(value).strip().lower()
    if key not in choices:
        raise ValueError(f"unsupported {name} option {value!r}; choose from {sorted(choices)}")
    return choices[key]


def _compile_options(options: MutableMapping[str, Any]) -> dict[str, int]:
    recognised = {
        "nch",
        "particle mesh points",
        "thermal",
        "sei",
        "sei porosity change",
        "particle mechanics",
        "loss of active material",
        "lithium plating",
    }
    unknown = {str(name).lower() for name in options} - recognised
    if unknown:
        raise ValueError(f"unsupported model option(s): {', '.join(sorted(unknown))}")
    normal = {str(name).lower(): value for name, value in options.items()}
    nch = int(normal.get("nch", normal.get("particle mesh points", 8)))
    if nch not in (5, 8, 12):
        raise ValueError("nch must select a compiled registry entry: 5, 8, or 12")
    thermal = str(normal.get("thermal", "isothermal")).lower()
    if thermal not in ("isothermal", "lumped"):
        raise ValueError("thermal must be 'isothermal' or 'lumped'")
    return {
        "nch": nch,
        "thermal": int(thermal == "lumped"),
        "sei_mask": _mask(
            normal.get("sei"),
            {
                "reaction limited": 1,
                "constant": 1,
                "solvent-diffusion limited": 2,
                "electron-migration limited": 4,
                "interstitial-diffusion limited": 8,
            },
            "SEI",
        ),
        "sei_porosity": int(
            str(normal.get("sei porosity change", "false")).lower() == "true"
        ),
        "crack_mask": _mask(
            normal.get("particle mechanics"),
            {"swelling and cracking": 1, "swelling only": 1},
            "particle mechanics",
        ),
        "lam_mask": _mask(
            normal.get("loss of active material"),
            {"stress-driven": 1, "reaction-driven": 2},
            "loss of active material",
        ),
        "plating": int(
            str(normal.get("lithium plating", "none")).lower()
            not in ("none", "false")
        ),
    }


def available_devices() -> tuple[str, ...]:
    """Return compiled execution backends. CUDA is added by the Phase-8 build when available."""

    return ("cpu",)


class ProcessedVariable:
    """A sampled solution field with PyBaMM-style entries and interpolation."""

    def __init__(self, name: str, time: np.ndarray, entries: np.ndarray):
        self.name = name
        self.t = np.asarray(time, dtype=float)
        self.entries = np.asarray(entries, dtype=float)
        self.data = self.entries

    def __call__(self, t: float | Sequence[float], **_: Any) -> np.ndarray | float:
        query = np.asarray(t, dtype=float)
        if self.entries.ndim == 1:
            result = np.interp(query, self.t, self.entries)
        else:
            result = np.stack(
                [np.interp(query, self.t, self.entries[:, lane]) for lane in range(self.entries.shape[1])],
                axis=-1,
            )
        return float(result) if result.ndim == 0 else result

    def plot(self, **kwargs: Any) -> Any:
        import matplotlib.pyplot as plt

        return plt.plot(self.t, self.entries, **kwargs)


class Solution:
    """Eager time series with lazy PyBaMM-compatible processed-variable views."""

    def __init__(self, data: MutableMapping[str, Any]):
        self.t = np.asarray(data["time"], dtype=float)
        self.termination = str(data["termination"])
        self.segment = int(data["segment"])
        self._fields = {
            "Time [s]": np.asarray(data["time"], dtype=float),
            "Terminal voltage [V]": np.asarray(data["voltage"], dtype=float),
            "Voltage [V]": np.asarray(data["voltage"], dtype=float),
            "Current [A]": np.asarray(data["current"], dtype=float),
        }
        self.all_models: list[Any] = []
        self.sensitivities: dict[str, np.ndarray] = {}
        self.cycles = [self]

    def __getitem__(self, name: str) -> ProcessedVariable:
        if name not in self._fields:
            raise KeyError(f"SLIDE did not record {name!r}; available: {sorted(self._fields)}")
        return ProcessedVariable(name, self.t, self._fields[name])

    def plot(self, output_variables: str | Sequence[str] = "Terminal voltage [V]", **kwargs: Any) -> Any:
        names = [output_variables] if isinstance(output_variables, str) else list(output_variables)
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(len(names), 1, squeeze=False)
        for axis, name in zip(axes[:, 0], names):
            axis.plot(self.t, self._fields[name])
            axis.set_xlabel("Time [s]")
            axis.set_ylabel(name)
        if kwargs.get("show_plot", False):
            plt.show()
        return figure

    def save_data(
        self,
        filename: str | Path,
        variables: Sequence[str] | None = None,
        *,
        to_format: str = "csv",
        short_names: MutableMapping[str, str] | None = None,
    ) -> None:
        names = list(variables or ("Time [s]", "Terminal voltage [V]", "Current [A]"))
        labels = [dict(short_names or {}).get(name, name) for name in names]
        arrays = [self._fields[name] for name in names]
        path = Path(filename)
        format_name = to_format.lower()
        if format_name == "csv":
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(labels)
                writer.writerows(zip(*arrays, strict=True))
        elif format_name in ("matlab", "mat"):
            from scipy.io import savemat

            savemat(path, {label: values for label, values in zip(labels, arrays, strict=True)})
        else:
            raise ValueError("to_format must be 'csv' or 'matlab'")


class Simulation:
    """Build and solve one precompiled SLIDE model."""

    def __init__(
        self,
        model: SPM | None = None,
        *,
        experiment: Experiment | None = None,
        parameter_values: ParameterValues | None = None,
        solver: Any = None,
        geometry: Any = None,
        submesh_types: Any = None,
        var_pts: Any = None,
        spatial_methods: Any = None,
        device: str = "cpu",
    ) -> None:
        if any(value is not None for value in (solver, geometry, submesh_types, var_pts, spatial_methods)):
            raise NotImplementedError(
                "SLIDE selects precompiled numerical methods; custom solver/mesh objects are out of scope"
            )
        self.model = model or SPM()
        self.experiment = experiment
        self.parameter_values = parameter_values or ParameterValues("Chen2020")
        self.device = str(device).lower()
        if self.device not in available_devices():
            raise RuntimeError(
                f"device={device!r} is unavailable in this build; compiled devices: {available_devices()}"
            )
        self.solution: Solution | None = None

    def solve(
        self,
        t_eval: Sequence[float] | None = None,
        *,
        initial_soc: float | None = None,
        inputs: MutableMapping[str, float] | None = None,
        **_: Any,
    ) -> Solution:
        parameters = self.parameter_values.copy()
        if initial_soc is not None:
            parameters["Initial state-of-charge"] = initial_soc
        if inputs:
            parameters.update(inputs)
        if parameters.varied_overrides():
            varied_values = parameters.varied_overrides()
        if self.experiment is None:
            if t_eval is None or len(t_eval) < 2:
                raise ValueError("supply an Experiment or at least two t_eval samples")
            duration = float(t_eval[-1]) - float(t_eval[0])
            experiment = Experiment(f"Rest for {duration} seconds", period=np.min(np.diff(t_eval)))
        else:
            experiment = self.experiment
        sample_step = experiment.period
        if t_eval is not None and len(t_eval) >= 2:
            differences = np.diff(np.asarray(t_eval, dtype=float))
            if np.any(differences <= 0):
                raise ValueError("t_eval must be strictly increasing")
            sample_step = float(np.min(differences))
        if parameters.varied_overrides():
            native = _slide_core.solve_ensemble(
                parameters.source,
                parameters.scalar_overrides(),
                {name: values.tolist() for name, values in varied_values.items()},
                _compile_options(self.model.options),
                experiment.steps,
                sample_step,
            )
            lanes = int(native["n_lanes"])
            native["voltage"] = np.asarray(native["voltage"], dtype=float).reshape(-1, lanes)
            native["current"] = np.asarray(native["current"], dtype=float).reshape(-1, lanes)
        else:
            native = _slide_core.solve_experiment(
                parameters.source,
                parameters.scalar_overrides(),
                _compile_options(self.model.options),
                experiment.steps,
                sample_step,
            )
        self.solution = Solution(native)
        return self.solution

    def plot(self, output_variables: str | Sequence[str] = "Terminal voltage [V]", **kwargs: Any) -> Any:
        if self.solution is None:
            raise RuntimeError("call solve() before plot()")
        return self.solution.plot(output_variables, **kwargs)


def plot(solution: Solution, output_variables: str | Sequence[str] = "Terminal voltage [V]", **kwargs: Any) -> Any:
    return solution.plot(output_variables, **kwargs)
