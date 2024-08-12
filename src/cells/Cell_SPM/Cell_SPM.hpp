/**
 * @file Cell_SPM.hpp
 * @brief Header file for Cell_SPM class
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 2021
 */

#pragma once

#include "State_SPM.hpp" //!< class that represents the state of a cell, the state is the collection of all time-varying conditions of the battery
#include "Electrode_SPM.hpp"
#include "Model_SPM.hpp" //!< defines a struct with the values for the matrices used in the spatial discretisation of the diffusion PDE
#include "param/param_SPM.hpp"
#include "../Cell.hpp"
#include "../../utility/utility.hpp"   // Do not remove they are required in cpp files.
#include "settings.hpp" // Do not remove they are required in cpp files.
#include "OCVcurves.hpp"
#include "Pair.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <memory>

namespace slide {

//!< State related functions
void validState(State_SPM &s, State_SPM &s_ini);

class Cell_SPM : public Cell
{
public:
  DEG_ID deg_id; //!< structure with the identification of which degradation model(s) to use #TODO may be protected.
  using sigma_type = std::array<double, settings::nch + 2>;

protected:                 //!< protected such that child classes can access the class variables
  State_SPM st{}, s_ini{}; //!< the battery current/initial state, grouping all parameter which change over the battery's lifetime (see State_SPM.hpp)

  Pair<Electrode_SPM> electrode{
    Electrode_SPM{
      .dom = pos, //!< Positive electrode:
      //!< parameters of the main li-insertion reaction
      .k = 5e-11,   //!< rate constant of main reaction at reference temperature
      .k_T = 58000, //!< activation energy for the Arrhenius relation of k
      .D_T = 29000, //!< activation energy for the Arrhenius relation of D
      .Cmax = 51385 //!< maximum lithium concentration in the cathode [mol m-3]  value for NMC
    },

    Electrode_SPM{
      .dom = neg, //!< Negative electrode
      //!< parameters of the main li-insertion reaction
      .k = 1.764e-11,       //!< rate constant of main reaction at reference temperature
      .k_T = 20000,         //!< activation energy for the Arrhenius relation of k
      .D_T = 35000.0 / 5.0, //!< activation energy for the Arrhenius relation of D
      .Cmax = 30555,        //!< maximum lithium concentration in the anode [mol m-3] value for C
    },

  };

  //!< Thermal model parameters
  double Therm_Qgen{};             //!< total heat generation since the last update [J]
  double Therm_Qgentot{};          //!< variable for unit testing, total heat generation since the beginning of this cell's life [J]
  double Therm_time{};             //!< time since the last update of the thermal model
  double T_env{ settings::T_ENV }; //!< environment temperature [K]
  double T_ref{ 25.0_degC };       //!< reference temperature [K]

  //!< Qch: 90 gives very good cooling, as if there is a fan pointed at the cell. values of 30-60 are representative for a cell on a shelf without forced cooling
  //!< 40 is representative for a cell on a shelf without forced cooling
  double Qch{ 45 };   //!< convective heat transfer coefficient per volume [W K-1 m-3]
  double rho{ 1626 }; //!< density of the battery
  double Cp{ 750 };   //!< thermal capacity of the battery #TODO = units missing.

  double xp_0{ 0.983999588653496 }, xp_100{ 0.400145394039564 }, xn_0{ 0.029397569380507 }, xn_100{ 0.932469496648387 }; // Lithium fractions at 0 to 100 percent.

  //!< Geometric parameters
  param::Geometry_SPM geo{};
  //!< other geometric parameters are part of State because they can change over the battery's lifetime

  param::StressParam sparam{ param::def::StressParam_Kokam }; //!< Stress parameters.

  //!< Constants and parameters for the SEI growth model //!< #TODO if we can make these static and speed gain.
  double rsei{ 2037.4 * 50 }; //!< specific resistance times real surface area of the SEI film [Ohm m] ? #TODO if unit is correct. aging fit.
  double nsei{ 1 };           //!< number of electrons involved in the SEI reaction [-]
  double alphasei{ 1 };       //!< charge transfer coefficient of the SEI reaction [-]
  double OCVsei{ 0.4 };       //!< equilibrium potential of the SEI side reaction [V]
  double rhosei{ 100e3 };     //!< partial molar volume of the SEI layer [m3 mol-1]
  double c_elec0{ 4.541e-3 }; //!< bulk concentration of the electrolyte molecule participating in the SEI growth (e.g. EC) [mol m-3]
  double Vmain{ 13 };         //!< partial molar volume of the main reaction, see Ashwin et al, 2016
  double Vsei{ 64.39 };       //!< partial molar volume of the SEI side reaction, see Ashwin et al., 2016

  param::SEIparam sei_p{ param::def::SEIparam_Kokam }; //!< structure with the fitting parameters of the different SEI growth models

  //!< surface crack parameters & constants
  param::CSparam csparam{}; //!< structure with the fitting parameters of the different crack growth models

  //!< LAM parameters & constants
  double OCVnmc{ 4.1 };                                //!< equilibrium potential of the NMC dissolution side reaction [V]
  param::LAMparam lam_p{ param::def::LAMparam_Kokam }; //!< structure with the fitting parameters of the different LAM models

  //!< Li-plating parameters & constants //!< #TODO if we can make these static and speed gain.
  double npl{ 1 };      //!< number of electrons involved in the plating reaction [-]
  double alphapl{ 1 };  //!< charge transfer constant for the plating reaction [-]
  double OCVpl{ 0 };    //!< OCV of the plating reaction [V]
  double rhopl{ 10e6 }; //!< density of the plated lithium layer
  param::PLparam pl_p;  //!< structure with the fitting parameters of the different plating models

  //!< Matrices for spatial discretisation of the solid diffusion model
  Model_SPM<settings::nch> *M{ Model_SPM<settings::nch>::makeModel() };

  //!< OCV curves
  OCVcurves OCV_curves;

  bool Vcell_valid{ false };

  //!< Functions
  double calcSurfaceConcentration(double molarFlux, double Dt, Domain dom);

  inline double calcArrheniusCoeff() { return (1 / T_ref - 1 / st.T()) / PhyConst::Rg; } //!< Calculates Arrhenius coefficient.
  std::pair<double, DPair> calcMolarFlux();                                              //!< Calculate molar flux

  //!< void setStates(State_SPM &&states);											  //!< set the cell's states to the states in the array

  //!< degradation models
  void SEI(double OCVnt, double etan, double *isei, double *den);                                                                             //!< calculate the effect of SEI growth
  void CS(double OCVnt, double etan, double *isei_multiplyer, double *dCS, double *dDn);                                                      //!< calculate the effect of surface crack growth
  void LAM(bool critical, double zp_surf, double etap, double *dthickp, double *dthickn, double *dap, double *dan, double *dep, double *den); //!< calculate the effect of LAM
  double LiPlating(double OCVnt, double etan);

  //!< state space model
  void dState_diffusion(bool print, State_SPM &d_state);   //!< just diffusion PDE
  void dState_thermal(bool print, double &dQgen);          //!< calculate the heat generation
  void dState_degradation(bool print, State_SPM &d_state); //!< calculate the effect of lithium plating
  void dState_all(bool print, State_SPM &d_state);         //!< individual functions are combined in one function to gain speed.

  //!< thermal model
  double thermalModel_cell();
  double thermalModel_coupled(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim);

  //!< cell to cell variations // #TODO why do we store variations?
  double var_cap{ 1 };    //!< relative factor increasing the capacity of the cell
  double var_R{ 1 };      //!< relative factor increasing the DC resistance
  double var_degSEI{ 1 }; //!< relative factor to speed up or slow down the rate of SEI growth
  double var_degLAM{ 1 }; //!< relative factor to speed up or slow down the rate of LAM

public:
  //!< Constructor
  Cell_SPM(OCVcurves OCV_curves_) : OCV_curves(OCV_curves_) {}
  Cell_SPM(std::string IDi, const DEG_ID &degid, double capf, double resf, double degfsei, double degflam);

  Cell_SPM(); //!< Default constructor.

  Cell_SPM(Model_SPM<settings::nch> *M_ptr) : M(M_ptr) {}

  //!< getters
  double T() noexcept override { return st.T(); }   //!< returns the uniform battery temperature in [K]
  double getTenv() const noexcept { return T_env; } //!< get the environmental temperature [K]

  Status setCurrent(double Inew, bool checkV = true, bool print = true) override;
  Status setSOC(double SOCnew, bool checkV = true, bool print = true) override;

  auto &getStateObj() { return st; }
  auto setStateObj(const State_SPM &st_new)
  {
    st = st_new;
    Vcell_valid = false;
  }

  std::array<double, 4> getVariations() const noexcept override { return { var_cap, var_R, var_degSEI, var_degLAM }; } // #TODO : deprecated will be deleted.

  void getTemperatures(double *Tenv, double *Tref) noexcept //!< get the environmental and reference temperature
  {
    /*
     * Function to get the environmental and reference temperatures
     * OUT
     * Tenv 	environmental temperature [K]
     * Tref 	reference temperature [K] at which cell parameters are measured
     */
    *Tenv = T_env;
    *Tref = T_ref;
  }

  //!< thermal model
  double getThermalSurface() override;
  double thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim) override;
  double thermal_getTotalHeat(); //!< function for unit testing

  //!< double getR() noexcept //!< get the total cell DC resistance
  //!< {
  //!< 	/*
  //!< 	 * Return the total cell DC resistance [Ohm]
  //!< 	 * The total cell resistance is the sum of the resistance of the electrodes and the SEI layer
  //!< 	 * the resistance of the electrodes increases due to loss of active material (LAM)
  //!< 	 * the resistance of the SEI layer increases as the layer becomes thicker
  //!< 	 * it is assumed the plated lithium does not add resistance
  //!< 	 */
  //!< 	return (rsei * st.delta() + st.getR(elec_surf));
  //!< }

  double getRdc() noexcept; //!< calculate the resistance from every component
  double getRtot() override { return getRdc(); }

  double getAnodeSurface() noexcept { return st.an() * st.thickn() * geo.elec_surf; } //!< get the anode pure surface area (without cracks) product of the effective surface area (an) with the electrode volume

  double I() const override { return st.I(); } //!< get the cell current [A]  positive for discharging
  double V() override;

  bool getCSurf(DPair &c_surf, bool print);                                                                                      //!< get the surface concentrations
  std::array<double, settings::nch + 2> getC(Domain dom) noexcept;                                                               //!< get the concentrations at all nodes
  int getVoltage(bool print, double *V, double *OCVp, double *OCVn, double *etap, double *etan, double *Rdrop, double *Temp);    //!< get the cell's voltage
  int getVoltage_ne(bool print, double *V, double *OCVp, double *OCVn, double *etap, double *etan, double *Rdrop, double *Temp); //!< get the cell's voltage noexcept

  void getDaiStress(double *sigma_p, double *sigma_n, sigma_type &sigma_r_p, sigma_type &sigma_r_n, sigma_type &sigma_t_p, sigma_type &sigma_t_n, sigma_type &sigma_h_p, sigma_type &sigma_h_n) noexcept; //!< get the stresses at all nodes according to Dai's stress model
  void updateDaiStress() noexcept;                                                                                                                                                                        //!< updated the stored stress values for Dai's stress model
  void getLaresgoitiStress(bool print, double *sigma_n);                                                                                                                                                  //!< get the stresses at all nodes according to Laresgoiti's stress model
  void updateLaresgoitiStress(bool print);                                                                                                                                                                //!< update the stored stress values for Laresgoiti's stress model

  //!< setters
  //!< void setVlimits(double VMAX, double VMIN); //!< set the voltage limits of the cell
  void setT(double T) override; //!< set the cell's temperature
  void setTenv(double Tenv);    //!< set the environmental temperature
  //!< void setStates(const State_SPM &si, double I);		  //!< set the cell's states to the states in the State object and the cell current to the given value
  void setC(DPair lifrac); //!< set the concentrations to the given (uniform) concentration
  //!< void setCurrent(bool critical, bool check, double I); //!< set the cell's current to the specified value -> From old slide.
  void peekVoltage(double I); //!< Peeks voltage for state for given I

  //!< State related functions

  //!< void validState()
  //!< {
  //!< 	slide::validState(st, s_ini);
  //!< }
  ThroughputData getThroughputs() override { return { st.time(), st.Ah(), st.Wh() }; }

  void overwriteCharacterisationStates(DPair D_, double ri)
  {
    //!< Overwrite both current and initial states.
    st.overwriteCharacterisationStates(D_, ri);
    s_ini.overwriteCharacterisationStates(D_, ri);
  }

  void overwriteGeometricStates(DPair thick_, DPair e_, DPair a_)
  {
    //!< Overwrite both current and initial states.
    st.overwriteGeometricStates(thick_, e_, a_);
    s_ini.overwriteGeometricStates(thick_, e_, a_);
  }

  double calculateIntegral();
  //!< time integration
  //!< void integratorStep(bool print, double dti, bool blockDegradation); //!< step forward in time using specified numerical integrator

  //!< void ETI_electr(bool print, double I, double dti, bool blockDegradation, bool pos_); //!< step forward with only one electrode using forward Euler time integration

  //!< Base integrators:

  //!< State_SPM::states_type Int_FWEuler(bool print, double dti, bool blockDegradation); //!< Forward euler integrator.
  //!< State_SPM::states_type Int_RK4(bool print, double dti, bool blockDegradation);	   //!< Runge-Kutta 4 integrator.

  //!< Calculate the time derivatives of the states at the actual cell current (state-space model)
  State_SPM::states_type dState(bool print, State_SPM &d_state);

  //!< Utility
  void checkModelparam(); //!< check if the inputs to the MATLAB code are the same as the ones here in the C++ code

  //!< --- slide-pack functions --- //
  void getStates(getStates_t s) override { s.insert(s.end(), st.begin(), st.end()); } //!< returns the states of the cell collectively.
  std::span<double> viewStates() override { return std::span<double>(st.begin(), st.end()); }
  double getOCV() override;
  Status setStates(setStates_t sSpan, bool checkV, bool print) override;
  bool validStates(bool print = true) override;
  inline double SOC() override { return calculateIntegral(); }
  void timeStep_CC(double dt, int steps = 1) override;

  Cell_SPM *copy() override { return new Cell_SPM(*this); }
  //!< Obsolete functions (do not use):
  void ETI(bool print, double dti, bool blockDegradation); //!< step forward in time using forward Eurler time integration

  //!< Functions for fitting:
  void setOCVcurve(const Pair<std::string> &name);                                      //!< sets the OCV curve of the cell to the given value
  void setInitialConcentration(DPair cmax, DPair lifrac);                               //!< sets the initial concentration
  void setGeometricParameters(double capnom, double elec_surf, DPair e_, DPair thick_); //!< sets the geometric parameters related to the amount of active material
  //!< void setRamping(double Istep, double tstep);																	  //!< sets the ramping parameters

  void setCharacterisationParam(DPair D_, DPair k_, double Rdc); //!< sets the parameters related to the characterisation of the cell
};
} // namespace slide