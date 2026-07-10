from __future__ import annotations

import csv
from datetime import datetime, timedelta
from pathlib import Path

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


def test_p7g1_tutorial_five_segment_local_parity():
    reference = np.loadtxt(
        Path(__file__).parents[2] / "tests" / "reference" / "pybamm_chen2020_tutorial5.csv",
        delimiter=",",
        skiprows=1,
    )
    steps = [
        "Discharge at C/10 for 10 hours or until 3.3 V",
        "Rest for 1 hour",
        "Charge at 1 A until 4.1 V",
        "Hold at 4.1 V until 50 mA",
        "Rest for 1 hour",
    ]
    solution = slide.Simulation(
        slide.SPM({"nch": 12}), experiment=slide.Experiment(steps, period="60 seconds")
    ).solve(initial_soc=1.0)
    voltage = solution["Terminal voltage [V]"].entries
    errors = []
    for segment in range(len(steps)):
        reference_mask = reference[:, 0] == segment
        slide_mask = solution.sample_segment == segment
        start = 0.0 if segment == 0 else solution.t[solution.sample_segment == segment - 1][-1]
        local_time = solution.t[slide_mask] - start
        comparable = (
            reference_mask
            & (reference[:, 2] >= local_time[0] - 1e-9)
            & (reference[:, 2] <= local_time[-1] + 1e-9)
        )
        errors.extend(
            np.interp(reference[comparable, 2], local_time, voltage[slide_mask])
            - reference[comparable, 4]
        )
    errors = np.asarray(errors)
    assert np.max(np.abs(errors)) <= 15e-3
    assert np.sqrt(np.mean(errors**2)) <= 8e-3


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


def test_simulate_s1_exposes_true_forward_sensitivities():
    names = slide.sensitivity_parameters()
    assert len(names) == 10
    selected = ["Negative particle diffusivity [m2.s-1]", "Contact resistance [Ohm]"]
    simulation = slide.Simulation(
        slide.SPM({"nch": 12}),
        experiment=slide.Experiment("Discharge at 1 C for 60 seconds", period="10 seconds"),
    )
    solution = simulation.simulateS1(selected, initial_soc=0.8)
    assert set(solution.sensitivities) == set(selected)
    assert all(values.shape == solution.t.shape for values in solution.sensitivities.values())
    np.testing.assert_allclose(
        solution.sensitivities["Contact resistance [Ohm]"], -5.0, rtol=0, atol=1e-14
    )
    assert solution["Voltage [V]"].sensitivities.keys() == solution.sensitivities.keys()


def test_custom_steps_and_terminations_run_through_the_core():
    explicit = slide.step.CustomStepExplicit(
        lambda variables: 2.5,
        duration=20,
        period=10,
        direction="discharge",
    )
    custom = slide.Simulation(
        experiment=slide.Experiment(explicit)
    ).solve(initial_soc=0.55)
    standard = slide.Simulation(
        experiment=slide.Experiment("Discharge at 2.5 A for 20 seconds", period=10)
    ).solve(initial_soc=0.55)
    np.testing.assert_array_equal(custom.t, standard.t)
    np.testing.assert_array_equal(
        custom["Current [A]"].entries, standard["Current [A]"].entries
    )
    np.testing.assert_array_equal(
        custom["Voltage [V]"].entries, standard["Voltage [V]"].entries
    )

    termination = slide.step.CustomTermination(
        "five-second event",
        lambda variables: 5.0 - variables["Local time [s]"],
    )
    event = slide.Simulation(
        experiment=slide.Experiment(
            slide.step.c_rate(
                0.2,
                duration=20,
                period=10,
                termination=termination,
            )
        )
    ).solve(initial_soc=0.55)
    assert event.termination == "event"
    assert event.termination_detail == "five-second event"
    assert event.t[-1] == pytest.approx(5.0, abs=1e-11)

    implicit = slide.Simulation(
        experiment=slide.Experiment(
            slide.step.CustomStepImplicit(
                lambda variables: variables["Voltage [V]"] - 3.8,
                duration=1,
                period=1,
                direction="charge",
            )
        )
    ).solve(initial_soc=0.55)
    np.testing.assert_allclose(implicit["Voltage [V]"].entries, 3.8, atol=2e-10)
    assert np.all(implicit["Current [A]"].entries < 0)

    differential = slide.Simulation(
        experiment=slide.Experiment(
            slide.step.CustomStepImplicit(
                lambda variables: 1.0,
                duration=2,
                period=1,
                direction="discharge",
                control="differential",
            )
        )
    ).solve(initial_soc=0.55)
    np.testing.assert_array_equal(differential["Current [A]"].entries, [1.0, 1.0, 2.0])


def test_start_time_cuts_steps_and_inserts_rest():
    origin = datetime(2023, 1, 1, 8, 0, 0)
    experiment = slide.Experiment(
        [
            slide.step.string("Rest for 1 hour", start_time=origin),
            slide.step.string(
                "Rest for 10 minutes", start_time=origin + timedelta(minutes=30)
            ),
            slide.step.string(
                "Rest for 30 minutes", start_time=origin + timedelta(hours=1)
            ),
            slide.step.rest("1 hour"),
        ],
        period="10 minutes",
    )
    solution = slide.Simulation(experiment=experiment).solve(initial_soc=0.55)
    assert solution.t[-1] == 9000.0
    assert 3000.0 in solution.t
    assert 3600.0 in solution.t
    np.testing.assert_array_equal(solution["Current [A]"].entries, 0.0)

    with pytest.raises(ValueError, match="first experiment step"):
        slide.Experiment(
            [
                "Rest for 1 hour",
                slide.step.string("Rest for 1 hour", start_time=origin),
            ]
        )
