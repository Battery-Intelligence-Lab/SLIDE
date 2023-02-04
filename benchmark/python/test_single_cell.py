import os
import pybamm
import matplotlib.pyplot as plt
import numpy as np
import time

import importlib
import slide_defaults

importlib.reload(slide_defaults)

os.chdir(os.path.dirname(os.path.abspath(__file__)))

model = pybamm.lithium_ion.SPM()

params = slide_defaults.get_default_param()
slide_OCV = np.genfromtxt(
    "../../results/PyBAMM_1_CC_Crate_Cell_SPM_cellData.csv", delimiter=",", skip_header=3)


experiment = pybamm.Experiment([
    "Discharge at 16 A for 30 minutes"])

# experiment = pybamm.Experiment([
#     "Charge at 10 A for 15 minutes",
#     "Rest for 15 minutes",
#     "Discharge at 10 A for 15 minutes",
#     "Rest for 30 minutes"])


sim = pybamm.Simulation(model, parameter_values=params, experiment=experiment)

start_time = time.time()
solution = sim.solve([0, 3600])
end_time = time.time()

print(f"Elapsed time: {(end_time - start_time):.3f} seconds")


t = solution["Time [s]"]
V = solution["Terminal voltage [V]"]
Tc = solution['Cell temperature [C]']

plt.plot(t.data, V.data)
plt.plot(slide_OCV[1:, 4], slide_OCV[1:, 1], '--')

plt.figure()

plt.plot(t.data, Tc.data[0, :])
plt.plot(slide_OCV[1:, 4], slide_OCV[1:, 3]-273.15, '--')
