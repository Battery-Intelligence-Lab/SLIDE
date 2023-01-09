import pybamm
import matplotlib.pyplot as plt
import numpy as np

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

model = pybamm.lithium_ion.SPM()
params = pybamm.ParameterValues('Chen2020')


def Positive_exchange_current(c_e, c_s_surf, c_s_max, T):
    n = 1
    kp = 5e-11  # rate constant of main reaction at positive electrode at reference temperature
    Tref = 25 + 273.15
    kp_T = 58000
    ArrheniusCoeff = (1 / Tref - 1 / T) / pybamm.constants.R
    kpt = kp * pybamm.exp(kp_T * ArrheniusCoeff)
    i0p = kpt * n * pybamm.constants.F * \
        (c_e*c_s_surf*(c_s_max - c_s_surf)) ** 0.5

    return i0p


def Negative_exchange_current(c_e, c_s_surf, c_s_max, T):
    n = 1
    # rate constant of main reaction at positive electrode at reference temperature
    kn = 1.7640e-11
    Tref = 25 + 273.15
    kn_T = 20000
    ArrheniusCoeff = (1 / Tref - 1 / T) / pybamm.constants.R
    knt = kn * pybamm.exp(kn_T * ArrheniusCoeff)
    i0n = knt * n * pybamm.constants.F * \
        (c_e*c_s_surf*(c_s_max - c_s_surf)) ** 0.5

    return i0n


Tref = 25 + 273.15
Dp_T = 29000
Dn_T = 35000.0 / 5.0
Rg = 8.314


Dp_ref = 8e-14
Dn_ref = 7e-14

ArrheniusCoeff = (1.0 / Tref - 1.0 / 288.0) / Rg

Dpt = Dp_ref * np.exp(Dp_T * ArrheniusCoeff)
Dnt = Dn_ref * np.exp(Dn_T * ArrheniusCoeff)

params['Positive electrode exchange-current density [A.m-2]'] = Positive_exchange_current
params['Negative electrode exchange-current density [A.m-2]'] = Negative_exchange_current

params['Maximum concentration in positive electrode [mol.m-3]'] = 51385
params["Maximum concentration in negative electrode [mol.m-3]"] = 30555

params['Positive particle radius [m]'] = 8.5e-6
params['Negative particle radius [m]'] = 1.25e-5

params['Positive electrode diffusivity [m2.s-1]'] = Dpt
params['Negative electrode diffusivity [m2.s-1]'] = Dnt

params['Positive electrode thickness [m]'] = 70e-6
params['Negative electrode thickness [m]'] = 73.5e-6
params['Separator thickness [m]'] = (
    1.685e-4 - 73.5e-6 - 70e-6)  # Total minus electrodes

params['Positive electrode active material volume fraction'] = 0.5
params['Negative electrode active material volume fraction'] = 0.5

params['Initial concentration in positive electrode [mol.m-3]'] = 0.689332*51385
params['Initial concentration in negative electrode [mol.m-3]'] = 0.479283*30555

params['Positive electrode conductivity [S.m-1]'] = 1/2.8e-3
params['Negative electrode conductivity [S.m-1]'] = 1/2.8e-3

params['Initial inner SEI thickness [m]'] = 1e-9
params['Initial outer SEI thickness [m]'] = 1e-9

params['Nominal cell capacity [A.h]'] = 16
params['Electrode height [m]'] = 0.1
params['Electrode width [m]'] = 31*0.1*2  # 31 layers double layer


params['Reference temperature [K]'] = 25 + 273.15
params['Initial temperature [K]'] = 288
params['Ambient temperature [K]'] = 288

params['Typical current [A]'] = 16
params['Upper voltage cut-off [V]'] = 4.2
params['Lower voltage cut-off [V]'] = 2.7

params['SEI resistivity [Ohm.m]'] = 2037.4 * 50

params["Current function [A]"] = 16.0  # 16.0/50

slide_Un = np.genfromtxt("../../data/Kokam_OCV_C.csv", delimiter=",")
slide_Up = np.genfromtxt("../../data/Kokam_OCV_NMC.csv", delimiter=",")
slide_Un[0, 0] = 0
slide_Up[0, 0] = 0


slide_OCV = np.genfromtxt(
    "../../results/PyBAMM_1_CC_Cell_SPM_cellData.csv", delimiter=",", skip_header=3)


def Up_fun(x):
    return pybamm.Interpolant(slide_Up[:, 0], slide_Up[:, 1], x, 'slide_Up')


def Un_fun(x):
   # return pybamm.Interpolant(slide_Up[:, 0], slide_Up[:, 1],x)
    return pybamm.Interpolant(slide_Un[:, 0], slide_Un[:, 1], x, 'slide_Un')


params.update({"Negative electrode OCP [V]": Un_fun,
               "Positive electrode OCP [V]": Up_fun})

sim = pybamm.Simulation(model, parameter_values=params)

solution = sim.solve([0, 3600])

t = solution["Time [s]"]
V = solution["Terminal voltage [V]"]
Tc = solution['Cell temperature [C]']

plt.plot(t.data, V.data)
plt.plot(slide_OCV[1:, 4], slide_OCV[1:, 1], '--')

plt.figure()

plt.plot(t.data, Tc.data[0, :])
plt.plot(slide_OCV[1:, 4], slide_OCV[1:, 3]-273.15, '--')


# x = pybamm.linspace(0, 1, 1000)  # sto

# U_n = params["Negative electrode OCP [V]"]
# U_p = params["Positive electrode OCP [V]"]


# fig, ax = plt.subplots(1, 3, figsize=(12, 4))
# #ax[0].plot(x.entries, U_n(x))
# ax[0].plot(slide_Un[:, 0], slide_Un[:, 1])
# ax[0].set(xlabel="sto [-]", ylabel="$U_n$ [V]")
# ax[1].plot(x.entries, U_p(x).entries)
# ax[1].set(xlabel="sto [-]", ylabel="$U_p$ [V]")
# ax[2].plot(x.entries, U_p(x).entries - U_n(1-x).entries)
# ax[2].set(xlabel="sto [-]", ylabel="$U$ [V]")
# plt.tight_layout()
