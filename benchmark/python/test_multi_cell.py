import liionpack as lp
import numpy as np
import pybamm

import importlib
import slide_defaults
import os
importlib.reload(slide_defaults)
os.chdir(os.path.dirname(os.path.abspath(__file__)))

model = pybamm.lithium_ion.SPM()

params = slide_defaults.get_default_param()

# Generate the netlist
netlist = lp.setup_circuit(
    Np=2, Ns=1, Rb=1e-4, Rc=5e-2, Ri=9e-2, V=3.2, I=80.0)

output_variables = [
    'X-averaged total heating [W.m-3]',
    'Volume-averaged cell temperature [K]',
    'X-averaged negative particle surface concentration [mol.m-3]',
    'X-averaged positive particle surface concentration [mol.m-3]',
]

# Heat transfer coefficients

# Cycling experiment, using PyBaMM
experiment = pybamm.Experiment([
    "Discharge at 16 A for 30 minutes"],
    period="10 seconds")


# Solve pack
output = lp.solve(netlist=netlist,
                  parameter_values=params,
                  experiment=experiment,
                  output_variables=output_variables)


# Draw the circuit at final state
lp.draw_circuit(netlist, cpt_size=1.0, dpi=150, node_spacing=2.5)
lp.plot_output(output)
