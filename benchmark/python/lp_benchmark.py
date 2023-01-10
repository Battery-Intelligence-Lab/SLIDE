# This file is directly taken from liionpack benchmarks folder
# Please see https://github.com/pybamm-team/liionpack/blob/develop/benchmarks/benchmarks.py
# for further details.

import liionpack as lp
import pybamm


class BasicBenchmark:
    def setup(self):
        self.sim = lp.basic_simulation()

    def time_solve_model(self):
        self.sim.solve([0, 1800])


class SmallPack:

    timeout = 60

    def setup(self):
        self.netlist = lp.setup_circuit(Np=2, Ns=1, Rb=1e-4, Rc=1e-2)
        self.parameter_values = pybamm.ParameterValues("Chen2020")
        self.experiment = pybamm.Experiment(
            [
                "Discharge at 2 A for 5 minutes",
            ],
            period="10 seconds",
        )

    def time_discharge_1cpu(self):
        _ = lp.solve(
            netlist=self.netlist.copy(),
            parameter_values=self.parameter_values.copy(),
            experiment=self.experiment,
            initial_soc=0.5,
            nproc=1,
        )

    def time_discharge_2cpu(self):
        _ = lp.solve(
            netlist=self.netlist.copy(),
            parameter_values=self.parameter_values.copy(),
            experiment=self.experiment,
            initial_soc=0.5,
            nproc=2,
        )


class MediumPack:

    timeout = 120

    def setup(self):
        self.netlist = lp.setup_circuit(Np=32, Ns=10, Rb=1e-4, Rc=1e-2)
        self.parameter_values = pybamm.ParameterValues("Chen2020")
        self.experiment = pybamm.Experiment(
            [
                "Discharge at 32 A for 5 minutes",
            ],
            period="10 seconds",
        )

    def time_discharge_1cpu(self):
        _ = lp.solve(
            netlist=self.netlist.copy(),
            parameter_values=self.parameter_values.copy(),
            experiment=self.experiment,
            initial_soc=0.5,
            nproc=1,
        )

    def time_discharge_2cpu(self):
        _ = lp.solve(
            netlist=self.netlist.copy(),
            parameter_values=self.parameter_values.copy(),
            experiment=self.experiment,
            initial_soc=0.5,
            nproc=2,
        )


class LargePack:

    timeout = 600

    def setup(self):
        self.netlist = lp.setup_circuit(Np=64, Ns=10, Rb=1e-4, Rc=1e-2)
        self.parameter_values = pybamm.ParameterValues("Chen2020")
        self.experiment = pybamm.Experiment(
            [
                "Discharge at 64 A for 5 minutes",
            ],
            period="10 seconds",
        )
        I_app = 64
        self.long_experiment = pybamm.Experiment(
            [
                f"Charge at {I_app} A for 20 minutes",
                "Rest for 15 minutes",
                f"Discharge at {I_app} A for 20 minutes",
                "Rest for 30 minutes",
            ]
            * 3,
            period="10 seconds",
        )

    def time_discharge_1cpu(self):
        _ = lp.solve(
            netlist=self.netlist.copy(),
            parameter_values=self.parameter_values.copy(),
            experiment=self.experiment,
            initial_soc=0.5,
            nproc=1,
        )

    def time_discharge_2cpu(self):
        _ = lp.solve(
            netlist=self.netlist.copy(),
            parameter_values=self.parameter_values.copy(),
            experiment=self.experiment,
            initial_soc=0.5,
            nproc=2,
        )

    def time_long_cycle_2cpu(self):
        _ = lp.solve(
            netlist=self.netlist.copy(),
            parameter_values=self.parameter_values.copy(),
            experiment=self.long_experiment,
            initial_soc=0.5,
            nproc=2,
        )
