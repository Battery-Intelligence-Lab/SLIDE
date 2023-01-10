import pybamm
import matplotlib.pyplot as plt
import numpy as np

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

slide_Un = np.genfromtxt("../../data/Kokam_OCV_C.csv", delimiter=",")
slide_Up = np.genfromtxt("../../data/Kokam_OCV_NMC.csv", delimiter=",")
slide_Un[0, 0] = 0
slide_Up[0, 0] = 0


def Up_fun(x):
    return pybamm.Interpolant(slide_Up[:, 0], slide_Up[:, 1], x, 'slide_Up')


def Un_fun(x):
    return pybamm.Interpolant(slide_Un[:, 0], slide_Un[:, 1], x, 'slide_Un')


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


def Positive_diffusivity(sto, T):
    D_ref = 8e-14
    Dp_T = 29000
    T_ref = 25 + 273.15
    arrhenius = pybamm.exp(Dp_T / pybamm.constants.R * (1 / T_ref - 1 / T))
    return D_ref * arrhenius


def Negative_diffusivity(sto, T):
    D_ref = 7e-14
    Dp_T = 35000.0 / 5.0
    T_ref = 25 + 273.15
    arrhenius = pybamm.exp(Dp_T / pybamm.constants.R * (1 / T_ref - 1 / T))
    return D_ref * arrhenius


def get_default_param():
    params = pybamm.ParameterValues('Chen2020')

    params['Positive electrode exchange-current density [A.m-2]'] = Positive_exchange_current
    params['Negative electrode exchange-current density [A.m-2]'] = Negative_exchange_current

    params['Maximum concentration in positive electrode [mol.m-3]'] = 51385
    params["Maximum concentration in negative electrode [mol.m-3]"] = 30555

    params['Positive particle radius [m]'] = 8.5e-6
    params['Negative particle radius [m]'] = 1.25e-5

    params['Positive electrode diffusivity [m2.s-1]'] = Positive_diffusivity
    params['Negative electrode diffusivity [m2.s-1]'] = Negative_diffusivity

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

    params.update({"Negative electrode OCP [V]": Un_fun,
                   "Positive electrode OCP [V]": Up_fun})

    return params
