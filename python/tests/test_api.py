from __future__ import annotations

import csv

import numpy as np
import pytest

import slide


def test_parameter_values_are_dict_like_and_aliases_round_trip():
    values = slide.ParameterValues("Chen2020")
    assert values["Nominal cell capacity [A.h]"] == 5.0
    assert values["Negative electrode diffusivity [m2.s-1]"] == 3.3e-14
    values["Initial state-of-charge"] = 0.7
    assert values["Initial state-of-charge"] == 0.7
    copy = values.copy()
    copy["Initial state-of-charge"] = 0.6
    assert values["Initial state-of-charge"] == 0.7


def test_tutorial_five_style_experiment_and_solution_surface(tmp_path):
    model = slide.lithium_ion.SPM()
    experiment = slide.Experiment(
        [
            "Discharge at C/10 for 10 hours or until 3.3 V",
            "Rest for 1 hour",
            "Charge at 1 A until 4.1 V",
            "Hold at 4.1 V until 50 mA",
            "Rest for 1 hour",
        ],
        period="5 minutes",
    )
    simulation = slide.Simulation(model, experiment=experiment)
    solution = simulation.solve()
    voltage = solution["Terminal voltage [V]"]
    assert solution.t.size > 10
    assert voltage.entries.shape == solution.t.shape
    assert np.all(np.isfinite(voltage.entries))
    assert voltage(solution.t[3]) == pytest.approx(voltage.entries[3])
    assert solution.termination in ("event", "final time")

    filename = tmp_path / "solution.csv"
    solution.save_data(filename, to_format="csv")
    with filename.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    assert rows[0] == ["Time [s]", "Terminal voltage [V]", "Current [A]"]
    assert len(rows) == solution.t.size + 1


def test_device_preflight_and_input_validation():
    assert slide.available_devices() == ("cpu",)
    with pytest.raises(RuntimeError, match="unavailable"):
        slide.Simulation(slide.SPM(), device="gpu")
    with pytest.raises(ValueError, match="registry"):
        slide.Simulation(slide.SPM({"nch": 7}), experiment=slide.Experiment("Rest for 1 second")).solve()


def test_varied_marker_is_defensive():
    marker = slide.varied([0.4, 0.5, 0.6])
    values = slide.ParameterValues("Chen2020")
    values["Initial state-of-charge"] = marker
    np.testing.assert_array_equal(
        values.varied_overrides()["Initial state-of-charge"], [0.4, 0.5, 0.6]
    )
    with pytest.raises(ValueError):
        slide.varied([])


def test_varied_lanes_are_one_ensemble_batch():
    values = slide.ParameterValues("Chen2020")
    values["Initial state-of-charge"] = slide.varied([0.45, 0.55, 0.65])
    experiment = slide.Experiment("Discharge at 1 C for 60 seconds", period="10 seconds")
    solution = slide.Simulation(
        slide.SPM(), experiment=experiment, parameter_values=values
    ).solve()
    voltage = solution["Terminal voltage [V]"].entries
    assert voltage.shape == (7, 3)
    assert np.all(np.diff(voltage[0]) > 0)
    assert np.all(voltage[-1] < voltage[0])
