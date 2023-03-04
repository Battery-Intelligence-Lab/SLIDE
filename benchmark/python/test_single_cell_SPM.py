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
    "../../results/PyBAMM_2_CC_Crate_Cell_SPM_cellData.csv", delimiter=",", skip_header=3)


experiment = pybamm.Experiment([
    "Discharge at 16 A for 40 seconds"])

# experiment = pybamm.Experiment([
#     "Charge at 10 A for 15 minutes",
#     "Rest for 15 minutes",
#     "Discharge at 10 A for 15 minutes",
#     "Rest for 30 minutes"])

params['Current function [A]'] = 1

# experiment=experiment
sim = pybamm.Simulation(model, parameter_values=params)

start_time = time.time()
solution = sim.solve(t_eval=np.linspace(0, 3600, 4000))
end_time = time.time()

print(f"Elapsed time: {(end_time - start_time):.3f} seconds")

t = solution["Time [s]"]
V = solution["Terminal voltage [V]"]
overpotential = solution["X-averaged battery reaction overpotential [V]"]
ocv = solution["X-averaged battery open circuit voltage [V]"]
Tc = solution['Cell temperature [C]']

#fig, ax = pybamm.plot_voltage_components(solution)
plt.figure()
plt.plot(t.data, V.data)
plt.plot(slide_OCV[:, 4], slide_OCV[:, 1], '--')
plt.gca().legend(('PyBAMM-V', 'SLIDE-V'))

plt.figure()
plt.plot(t.data, ocv.data)
plt.plot(slide_OCV[:, 4], slide_OCV[:, 7], '--')
plt.gca().legend(('PyBAMM-OCV', 'SLIDE-OCV'))


plt.figure()
plt.plot(t.data, overpotential.data)
plt.gca().legend(('PyBAMM-overpotential', 'SLIDE-overpotential'))


surf_p = solution['Positive particle surface concentration']

plt.figure()
plt.plot(t.data, overpotential.data)
plt.gca().legend(('PyBAMM-overpotential', 'SLIDE-overpotential'))


plt.show()


# plt.figure()

# plt.plot(t.data, Tc.data[0, :])
# plt.plot(slide_OCV[1:, 4], slide_OCV[1:, 3]-273.15, '--')
# plt.gca().legend(('PyBAMM','SLIDE'))
# plt.show()
