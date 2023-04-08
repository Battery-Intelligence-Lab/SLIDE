/*
 * Cell_SPM_dstate.cpp
 *
 * Implements the functions for the parent class of the Cells
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Cell_SPM.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <utility>

namespace slide {
void Cell_SPM::dState_diffusion(bool print, State_SPM &d_st)
{
  /*
   * Calculate just the diffusion PDE
   * IN
   * print 	boolean indicating if we want to print error messages or not
   * 				if true, error messages are printed
   * 				if false no error messages are printed (but the error will still be thrown)
   * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
   * blockDegAndTherm 	if true, battery degradation is ignored (i.e. the time derivatives of those states are 0)
   *
   * OUT
   * dstates	change in the states
   * 		dzp			time derivative of the transformed concentration at the positive inner nodes of the positive electrode (dzp/dt)
   * 		dzn			time derivative of the transformed concentration at the positive inner nodes of the negative electrode (dzn/dt)
   *
   */

  const auto nch = st.nch;

  auto [Dpt, Dnt] = calcDiffConstant();
  auto [i_app, jp, jn] = calcMolarFlux(); //!< current density, molar flux on the pos/neg particle

  //!< Calculate the effect of the main li-reaction on the (transformed) concentration
  for (size_t j = 0; j < nch; j++)
    d_st.zp(j) = (Dpt * M->Ap[j] * st.zp(j) + M->Bp[j] * jp); //!< dz/dt = D * A * z + B * j

  //!< loop for each row of the matrix-vector product A * z
  for (size_t j = 0; j < nch; j++)                            //!< A is diagonal, so the array M->A has only the diagonal elements
    d_st.zn(j) = (Dnt * M->An[j] * st.zn(j) + M->Bn[j] * jn); //!< dz/dt = D * A * z + B * j

  d_st.SOC() += -I() / (Cap() * 3600); //!< dSOC state of charge
}

void Cell_SPM::dState_thermal(bool print, double &dQgen)
{
  /*
   * Calculate the internal heat generation
   * IN
   * print 	boolean indicating if we want to print error messages or not
   * 				if true, error messages are printed
   * 				if false no error messages are printed (but the error will still be thrown)
   * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
   * blockDegAndTherm 	if true, battery degradation is ignored (i.e. the time derivatives of those states are 0)
   *
   * OUT
   * dQgen 	change in generated heat
   */

  //!< Calculcate the lithium fractions at the surface of the particles
  auto [Dpt, Dnt] = calcDiffConstant();
  auto [i_app, jp, jn] = calcMolarFlux(); //!< current density, molar flux on the pos/neg particle
  auto [cps, cns] = calcSurfaceConcentration(jp, jn, Dpt, Dnt);

  //!< check if the surface concentration is within the allowed range
  //!< 	0 < cp < Cmaxpos
  //!< 	0 < cn < Cmaxneg
  if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg) //!< Do not delete if you cannot ensure zp/zn between 0-1
  {
    if (print) {
      std::cerr << "ERROR in Cell_SPM::dState: concentration out of bounds. the positive lithium fraction is "
                << cps / Cmaxpos << " and the negative lithium fraction is " << cns / Cmaxneg;
      std::cerr << "they should both be between 0 and 1.\n";
    }
    throw 101;
  }

  const double zp_surf = (cps / Cmaxpos); //!< lithium fraction (0 to 1)
  //!< const double zn_surf = (cns / Cmaxneg);

  //!< Calculate the overpotentials if needed
  auto [etap, etan] = calcOverPotential(cps, cns, i_app);

  constexpr int electr = 0;    //!< #TODO why don't we use in the new one?
  if constexpr (electr == 2) { //!< only consider negative electrode, ignore the positive electrode
    jp = 0;
    Dpt = 0;
    etap = 0;
  } else if constexpr (electr == 1) { //!< only consider positive electrode, ignore the negative electrode
    jn = 0;
    Dnt = 0;
    etan = 0;
  }

  //!< Calculate the entropic coefficient
  const bool bound = true;                                               //!< in linear interpolation, throw an error if you are outside of the allowed range of the data
  const double dOCV = OCV_curves.dOCV_tot.interp(zp_surf, print, bound); //!< entropic coefficient of the entire cell OCV [V K-1]

  //!< temperature model
  //!< Calculate the thermal sources/sinks/transfers per unit of volume of the battery
  //!< The battery volume is given by the product of the cell thickness and the electrode surface
  const double Qrev = -I() * T() * dOCV;    //!< reversible heat due to entropy changes [W]
  const double Qrea = I() * (etan - etap);  //!< reaction heat due to the kinetics [W]
  const double Qohm = I() * I() * getRdc(); //!< Ohmic heat due to electrode resistance [W]
  //!< const double Qc = -Qch * SAV * (st.T() - T_env); 				//!< cooling with the environment [W m-2]

  //!< total heat generation and cooling in W
  dQgen = (Qrev + Qrea + Qohm);
  //!< dQcool = (Qc) * (L*elec_surf);

  //!< dT = 1/(rho*Cp)*(Qrev+Qrea+Qohm+Qc);					//!< dT		cell temperature
}

void Cell_SPM::dState_degradation(bool print, State_SPM &d_state)
{
  /*
   * calculate the effects of degradation
   *
   * IN
   * print 	boolean indicating if we want to print error messages or not
   * 				if true, error messages are printed
   * 				if false no error messages are printed (but the error will still be thrown)
   * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
   *
   * OUT
   * dstates	change in the states
   * 		dzp			time derivative of the transformed concentration at the positive inner nodes of the positive electrode (dzp/dt)
   * 		dzn			time derivative of the transformed concentration at the positive inner nodes of the negative electrode (dzn/dt)
   * 		ddelta 		time derivative of the SEI thickness [m s-1] (ddelta/dt)
   * 		dLLI 		time derivative of the lost lithium inventory [C s-1] (dLLI/dt)
   * 		dthickp 	time derivative of the thickness of the positive electrode [m s-1] (dthickp/dt), <0 (dthickp/dt)
   * 		dthickn		time derivative of the thickness of the negative electrode [m s-1] (dthickn/dt), <0 (dthickn/dt)
   * 		dep			time derivative of the volume fraction in the positive electrode [s-1] (dep/dt)
   * 		den			time derivative of the volume fraction in the negative electrode [s-1] (den/dt)
   * 		dap			time derivative of the effective surface area of the positive electrode [m2 m-3 s-1] (dap/dt)
   * 		dan			time derivative of the effective surface area of the negative electrode [m2 m-3 s-1] (dan/dt)
   * 		dCS 		time derivative of the crack surface [m2 s-1], dCS/st > 0 (dCS/dt)
   * 		dDp 		time derivative of the diffusion constant at reference temperature of the positive electrode [m s-1 s-1] (dDp/dt)
   * 		dDn			time derivative of the diffusion constant at reference temperature of the negative electrode [m s-1 s-1] (dDn/dt)
   * 		drdc_p		time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
   * 		drdc_n		time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
   * 		drdc_cc		time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
   * 		ddelta_pl 	time derivative of the thickness of the plated lithium layer [m s-1] (ddelta_pl/dt)
   * 		dSoC
   * 		dT			time derivative of the battery temperature [K s-1] (dT/dt)
   * 		dI
   */
  using namespace PhyConst;
  using settings::nch;

  //!< Calculcate the lithium fractions at the surface of the particles
  auto [Dpt, Dnt] = calcDiffConstant();
  auto [i_app, jp, jn] = calcMolarFlux(); //!< current density, molar flux on the pos/neg particle
  auto [cps, cns] = calcSurfaceConcentration(jp, jn, Dpt, Dnt);

  //!< check if the surface concentration is within the allowed range
  //!< 	0 < cp < Cmaxpos
  //!< 	0 < cn < Cmaxneg
  if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg) //!< Do not delete if you cannot ensure zp/zn between 0-1
  {
    if (print) {
      std::cerr << "ERROR in Cell_SPM::dState: concentration out of bounds. the positive lithium fraction is "
                << cps / Cmaxpos << " and the negative lithium fraction is " << cns / Cmaxneg
                << " they should both be between 0 and 1.\n";
    }
    throw 101;
  }

  const double zp_surf = (cps / Cmaxpos); //!< lithium fraction (0 to 1)
  const double zn_surf = (cns / Cmaxneg);

  //!< Calculate the overpotentials if needed
  auto [etap, etan] = calcOverPotential(cps, cns, i_app);

  const bool bound = true;
  //!< calculate the anode potential (needed for various degradation models)
  const double dOCVn = OCV_curves.dOCV_neg.interp(zp_surf, print, bound); //!< entropic coefficient of the anode potential [V K-1]
  const double OCV_n = OCV_curves.OCV_neg.interp(zn_surf, print, bound);  //!< anode potential [V]
  const double OCVnt = OCV_n + (st.T() - T_ref) * dOCVn;                  //!< anode potential at the cell's temperature [V]

  //!< SEI growth
  double isei;        //!< current density of the SEI growth side reaction [A m-2]
  double den_sei;     //!< decrease in volume fraction due to SEI growth [s-1]
  double dznsei[nch]; //!< additional diffusion in the anode due to isei

  SEI(OCVnt, etan, &isei, &den_sei); //!< Throws but not wrapped in try-catch since only appears here.

  //!< Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF)
  for (size_t j = 0; j < nch; j++)
    dznsei[j] = (M->Bn[j] * isei / (nsei * F));

  //!< crack growth leading to additional exposed surface area
  double isei_multiplyer; //!< relative increase in isei due to additional SEI growth on the extra exposed surface area [-]
  double dCS;             //!< increase in crack surface area [m2 s-1]
  double dDn;             //!< change in negative diffusion constant [m s-1 s-1]
  double dznsei_CS[nch];  //!< additional diffusion in the anode due to extra SEI growth on the crack surface

  CS(OCVnt, etan, &isei_multiplyer, &dCS, &dDn); //!< Throws but not wrapped in try-catch since only appears here.

  //!< crack surface leads to extra SEI growth because the exposed surface area increases.
  //!< (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF + isei_CS/nF)
  const double isei_CS = isei * isei_multiplyer; //!< extra SEI side reaction current density due to the crack surface area [A m-2]
  for (int j = 0; j < nch; j++)
    dznsei_CS[j] = (M->Bn[j] * isei_CS / (nsei * F));

  //!< loss of active material LAM
  double dthickp, dthickn, dap, dan, dep, den; //!< change in geometric parameters describing the amount of active material

  LAM(print, zp_surf, etap, &dthickp, &dthickn, &dap, &dan, &dep, &den); //!< Throws but not wrapped in try-catch since only appears here.

  //!< lithium plating
  const double ipl = LiPlating(OCVnt, etan); //!< current density of the plating side reaction [A m-2]
  double dzn_pl[nch];                        //!< additional diffusion in the anode due to ipl

  //!< Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + ipl/nF)
  for (size_t j = 0; j < nch; j++)
    dzn_pl[j] = (M->Bn[j] * ipl / (npl * F)); //!< #TODO (npl * F) division is unnecessary it is already multiple in the function.

  //!< output
  for (size_t j = 0; j < nch; j++)                           //!< dzp += 0 //!< dzp should be added from diffusion function
    d_state.zn(j) += (dznsei[j] + dznsei_CS[j] + dzn_pl[j]); //!< dzn		jtot = jn + isei/nF + isei_CS/nF + ipl/nF

  d_state.delta() = isei / (nsei * F * rhosei);                                   //!< ddelta	SEI thickness
  d_state.LLI() = (isei + isei_CS + ipl) * geo.elec_surf * st.thickn() * st.an(); //!< dLLI		lost lithium
  d_state.thickp() = dthickp;                                                     //!< dthickp 	electrode thickness
  d_state.thickn() = dthickn;                                                     //!< dthickn
  d_state.ep() = dep;                                                             //!< dep		volume fraction of active material
  d_state.en() = den + den_sei;                                                   //!< den
  d_state.ap() = dap + 3 / geo.Rp * d_state.ep();                                 //!< dap		effective surface are, a = 3 e/R
  d_state.an() = dan + 3 / geo.Rn * d_state.en();                                 //!< dan
  d_state.CS() = dCS;                                                             //!< dCS		surface area of the cracks
  d_state.Dp() = 0;                                                               //!< dDp 		diffusion constant
  d_state.Dn() = dDn;                                                             //!< dDn
  d_state.rDCp() = 0;                                                             //!< drdc_p 	cathode resistance
  d_state.rDCn() = 0;                                                             //!< drdc_n 	anode resistance
  d_state.rDCcc() = 0;                                                            //!< drdc_cc 	current collector resistance
  d_state.delta_pl() = ipl / (npl * F * rhopl);                                   //!< ddelta_pl thickness of the plated lithium
}

void Cell_SPM::timeStep_CC(double dt, int nstep)
{
  /*
   * take a number of time steps with a constant current
   *
   * The diffusion model is resolved to every time step of dt
   * The thermal and degradation models are resolved once, for dt*nstep seconds
   *
   */

  //!< check the time step is positive
  if (dt < 0) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Cell_ECM::timeStep_CC, the time step dt must be 0 or positive, but has value " << dt << '\n';
    throw 10;
  }

  const bool print = true;

  //!< Update the stress values stored in the attributes with the stress of the previous time step
  sparam.s_dai_p_prev = sparam.s_dai_p;     //!< Dai's stress in the positive particle in the previous time step
  sparam.s_dai_n_prev = sparam.s_dai_n;     //!< Dai's stress in the negative particle in the previous time step
  sparam.s_lares_n_prev = sparam.s_lares_n; //!< Laresgoiti's stress in the negative particle in the previous time step
  sparam.s_dt = nstep * dt;

  //!< *********************************************  Resolve the diffusion model for every dt time step ****************************************************************************************
  const auto dth = dt / 3600.0;
  for (int t = 0; t < nstep; t++) {
    slide::State_SPM d_st{};
    //!< Calculate the time derivatives
    dState_diffusion(print, d_st);

    //!< forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
    for (size_t i = 0; i < (2 * st.nch); i++)
      st.z(i) += dt * d_st.z(i);

    st.SOC() += dt * d_st.SOC();

    Vcell_valid = false;
    sparam.s_dai_update = false;
    sparam.s_lares_update = false;


    const auto dAh = st.I() * dth;
    //!< increase the cumulative variables of this cell
    if constexpr (settings::data::storeCumulativeData) {
      st.time() += dt;
      st.Ah() += std::abs(dAh);
      st.Wh() += std::abs(dAh * V());
    }
  }

  //!< **************************************************** Calculate the thermal model once for the nstep * dt time period *****************************************************************
  if (!blockDegAndTherm) {

    //!< Calculate the internal heat generation of the cell
    double dQgen;
    try {
      dState_thermal(print, dQgen);
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in SPM cell " << getFullID() << ", error " << e
                  << " when calculating the time derivatives of the thermal model, throwing it on.\n";
      throw e;
    }

    Therm_Qgen += dQgen * dt * nstep;
    Therm_time += dt * nstep;

    //!< If this cell has a parent module, this parent will call the thermal model with the correct parameters
    //!< which will include heat exchange with the cell's neighbours and cooling from the cooling system of the module.

    //!< if there is no parent, this cell is a stand-alone cell.
    //!< We assume convective cooling with an environment
    //!< If there is no parent, assume we update T every nstep*dt. So update the temperature now
    if (!parent) //!< else it is the responsibility of the parent to call the thermal model function with the correct settings
    {
      double Tneigh[1] = { T_env };
      double Kneigh[1] = { Qch };                 //!< so cooling Qc = Qch * SAV * dT / (rho*cp) = Qch * A * dT / (rho*cp)
      double Atherm[1] = { getThermalSurface() }; //!< calculate the surface of this cell

      const auto new_T = thermalModel(1, Tneigh, Kneigh, Atherm, Therm_time);
      setT(new_T); //!< #TODO in slide-pack it does not check the temperature.
    }

    //!< ******************************************************* Calculate the degradation model once for the nstep * dt time period **************************************************************

    //!< Calculate the stress values stored in the attributes for the stress in this time step
    if (sparam.s_dai) //!< only a few degradation models actually need the stress according to Dai, so only calculate it if needed
      updateDaiStress();
    if (sparam.s_lares) //!< only a few degradation models need the stress according to Laresgoiti
      updateLaresgoitiStress(true);

    slide::State_SPM d_st{};
    //!< Calculate the time derivatives
    try {
      dState_degradation(print, d_st);
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in SPM cell " << getFullID()
                  << " when calculating the time derivatives of the degradation: "
                  << e << ", throwing it on.\n";
      throw e;
    }

    //!< forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
    //!< degradation accumulation
    for (size_t i = 0; i < st.size(); i++)
      st[i] += d_st[i] * dt * nstep;
  }
}

//!< void Cell_SPM::ETI(bool print, double dti, bool blockDegradation)
//!< {
//!< 	/*
//!< 	 * Performs forward Euler time integration over one time step of dti seconds
//!< 	 * s(t+1) = s(t) + ds/dt * dti
//!< 	 *
//!< 	 * IN
//!< 	 * print 			boolean indicating if we want to print error messages or not
//!< 	 * 					if true, error messages are printed
//!< 	 * 					if false no error messages are printed (but the error will still be thrown)
//!< 	 * 					we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
//!< 	 * dti 				time step over which time integradation should be done [s]
//!< 	 * 					it should be small enough to ensure numerical stability and accuracy
//!< 	 * 					1-5 seconds in usually ok (depending on the magnitude of the current, the temperature, etc.)
//!< 	 * blockDegradation if true, degradation is not accounted for in this time step
//!< 	 */

//!< 	//!< Calculate the time derivatives
//!< 	auto states = st;

//!< 	//!< calculate time derivatives, electr = 0 to account for both electrodes (i.e. cycle the full cell)
//!< 	const slide::states_type dstates = dState(print, blockDegradation, 0); //!< arrays with the state and dstate/dt

//!< 	//!< forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
//!< 	for (int i = 0; i < settings::ns; i++) //!< loop for all states
//!< 		states[i] += dti * dstates[i];

//!< 	setStates(std::move(states)); //!< store new states, checks if they are illegal (throws an error in that case)

//!< }

//!< void Cell_SPM::ETI_electr(bool print, double I, double dti, bool blockDegradation, bool pos)
//!< {
//!< 	/*
//!< 	 * Performs forward Euler time integration over one time step of dti seconds
//!< 	 * Considers only one electrode, i.e. you can do half-cell cycling with only one electrode
//!< 	 * This allows to calculate the half-cell OCV curves (versus a reference electrode of metallic lithium)
//!< 	 *
//!< 	 * USE THIS FUNCTION WITH CARE:
//!< 	 * It doesn't check the input current
//!< 	 * It might get the cell to an illegal situation
//!< 	 * You can't undo the effects of this function, unless you had stored the states before calling this function
//!< 	 *
//!< 	 * IN
//!< 	 * print 	boolean indicating if we want to print error messages or not
//!< 	 * 				if true, error messages are printed
//!< 	 * 				if false no error messages are printed (but the error will still be thrown)
//!< 	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
//!< 	 * I				current applied during this time step [A]
//!< 	 * 						> 0 for discharge
//!< 	 * 						< 0 for charge
//!< 	 * dti				time step [s]
//!< 	 * blockDegradation if true, degradation is not accounted for in this time step. Must be true
//!< 	 * pos 				if true, the positive electrode is considered (half-cell cycling with the positive electrode)
//!< 	 * 					if false, the negative electrode is considered (half-cell cycling with the negative electrode)
//!< 	 *
//!< 	 * THROWS
//!< 	 * 109 				illegal input parameters
//!< 	 */

//!< 	if (!blockDegradation)
//!< 	{
//!< 		std::cerr << "ERROR in Cell_SPM::ETI_electr, you are cycling only one electrode but want to account for degradation. This is not allowed.\n";
//!< 		//!< half-cell cycling is only allowed if blockDegradation is true (i.e. you ignore degradation)
//!< 		//!< this is because some of the degradation mechanisms could give NaN or Inf due to a divide by 0
//!< 		//!< So you can only call this function with 'true' for blockDegradation
//!< 		throw 109;
//!< 	}

//!< 	//!< Set the specified current
//!< 	Icell = I; //!< don't call setCurrent because that function ramps the current on both electrodes
//!< 			   //!< so here 'cheat it' and directly set the current
//!< 			   //!< this means you avoid the checks done in setCurrent, so you don't know if the current is feasible or not

//!< 	//!< Calculate the time derivatives
//!< 	auto states = st;

//!< 	//!< calculate time derivatives of the positive or negative electrode
//!< 	const slide::states_type dstates = pos ? dState(print, blockDegradation, 1) : dState(print, blockDegradation, 2); //!< arrays with the state and dstate/dt

//!< 	//!< forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
//!< 	for (int i = 0; i < settings::ns; i++) //!< loop for all states
//!< 		states[i] += dti * dstates[i];

//!< 	setStates(std::move(states)); //!< store new states

//!< }

//!< void Cell_SPM::integratorStep(bool print, double dti, bool blockDegradation)
//!< {
//!< 	/*
//!< 	 * Performs integration with specified solver over one time step of dti seconds
//!< 	 *
//!< 	 * IN
//!< 	 * print 			boolean indicating if we want to print error messages or not
//!< 	 * 					if true, error messages are printed
//!< 	 * 					if false no error messages are printed (but the error will still be thrown)
//!< 	 * 					we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
//!< 	 * dti 				time step over which time integradation should be done [s]
//!< 	 * 					it should be small enough to ensure numerical stability and accuracy
//!< 	 * 					1-5 seconds in usually ok (depending on the magnitude of the current, the temperature, etc.)
//!< 	 * blockDegradation if true, degradation is not accounted for in this time step
//!< 	 */

//!< 	State_SPM new_states{Int_FWEuler(print, dti, blockDegradation)};

//!< 	setStates(std::move(new_states)); //!< store new states, checks if they are illegal (throws an error in that case)

//!< }

//!< Base integrators:
//!< Be careful since they change the state values!

//!< slide::states_type Cell_SPM::Int_FWEuler(bool print, double dti, bool blockDegradation)
//!< {
//!<     //!< Get current states:
//!<     auto states = st;

//!<     //!< calculate time derivatives, electr = 0 to account for both electrodes (i.e. cycle the full cell)
//!<     const auto dstates = dState(print, blockDegradation, 0); //!< arrays with the state and dstate/dt

//!<     //!< forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
//!<     for (int i = 0; i < settings::ns; i++) //!< loop for all states
//!<         states[i] += dti * dstates[i];

//!<     return states;
//!< }

//!< slide::states_type Cell_SPM::Int_RK4(bool print, double dti, bool blockDegradation)
//!< {
//!<     const auto Iprev = I();

//!<     //!< Runge kutta 4
//!<     auto y1 = st;
//!<     const auto k1 = dState(settings::printBool::printCrit, blockDegradation, 0);

//!<     const auto y2 = arrSum(y1, k1, 1.0, 0.5 * dti);
//!<     setStates(State_SPM{y2}, Iprev);
//!<     const auto k2 = dState(settings::printBool::printCrit, blockDegradation, 0);

//!<     const auto y3 = arrSum(y1, k2, 1.0, 0.5 * dti);
//!<     setStates(State_SPM{y3}, Iprev);
//!<     const auto k3 = dState(settings::printBool::printCrit, blockDegradation, 0);

//!<     const auto y4 = arrSum(y1, k3, 1.0, dti);
//!<     setStates(State_SPM{y4}, Iprev); //!< #TODO: eliminate temporary objects State_SPM{}
//!<     const auto k4 = dState(settings::printBool::printCrit, blockDegradation, 0);

//!<     for (size_t i = 0; i < y1.size(); i++)
//!<         y1[i] += (1.0 / 6.0) * dti * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);

//!<     return y1;
//!< }

} // namespace slide
