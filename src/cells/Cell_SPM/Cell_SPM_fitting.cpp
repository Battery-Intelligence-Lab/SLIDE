/*
 * CellFit.cpp
 *
 * This cell is used to fit the cycling parameters of a cell to data.
 * It is done by the functions in determineCharacterisation.cpp
 * This cell should never be used by other functions (e.g. it should never be used for degradation simulations)
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

namespace slide {
//!< Cell_Fit::Cell_Fit(const slide::Model_SPM &MM, int verbosei) #TODO important for ageing to set no ageing.
//!< 	: Cell(settings::path::Kokam::namepos, settings::path::Kokam::nameneg, settings::path::Kokam::nameentropicC, settings::path::Kokam::nameentropicCell)
//!< {
//!< 	/*
//!< 	 * Standard constructor to initialise the battery parameters.
//!< 	 * Nothing should be changed in this constructor unless the user changes the radius of the particles.
//!< 	 * The values used for fitting can be changed using the functions defined in this class
//!< 	 *
//!< 	 * IN
//!< 	 * M 			structure of the type Model, with the matrices of spatial discretisation for solid diffusion
//!< 	 * verbosei		integer indicating how verbose the simulation has to be.
//!< 	 * 				The higher the number, the more output there is.
//!< 	 * 				Recommended value is 1, only use higher values for debugging
//!< 	 * 				From 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
//!< 	 * 					0 	almost no messages are printed, only in case of critical errors related to illegal parameters
//!< 	 * 					1 	error messages are printed in case of critical errors which might crash the simulation
//!< 	 * 					2 	all error messages are printed, whether the simulation can recover from the errors or not
//!< 	 * 					3 	on top of the output from 2, a message is printed every time a function in the Cycler and BasicCycler is started and terminated
//!< 	 * 					4 	on top of the output from 3, the high-level flow of the program in the Cycler is printed (e.g. 'we are going to discharge the cell')
//!< 	 * 					5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
//!< 	 * 					6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
//!< 	 * 					7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated
//!< 	 *
//!< 	 * THROWS
//!< 	 * 110 	the matrices for the solid diffusion discretisation, produced by MATLAB, are wrong
//!< 	 * 111 	the OCV curves are too long
//!< 	 */

//!< 	verbose = verbosei;

//!< 	//!< maximum concentrations
//!< 	Cmaxpos = 51385;
//!< 	Cmaxneg = 30555;
//!< 	C_elec = 1000; //!< standard concentration of 1 molar

//!< 	//!< constants
//!< 	n = 1;

//!< 	//!< Cell parameters
//!< 	nomCapacity = 2.7; //!< 2.7 Ah capacity
//!< 	Vmax = 4.2;		   //!< value for an NMC/C cell
//!< 	Vmin = 2.7;		   //!< value for an NMC/C cell
//!< 	Icell = 0;		   //!< initial cell current is 0A
//!< 	dIcell = 0.1;	   //!< ramp at 0.1A
//!< 	dt_I = 1e-3;	   //!< ramp at 1ms so changing the current goes fast
//!< 					   //!< now changing the current takes 0.01 second per A

//!< 	//!< thermal parameters
//!< 	T_ref = PhyConst::Kelvin + 25;
//!< 	T_env = PhyConst::Kelvin + 25;
//!< 	Qch = 90; //!< 90 gives very good cooling, as if there is a fan pointed at the cell. values of 30-60 are representative for a cell on a shelf without forced cooling
//!< 	rho = 1626;
//!< 	Cp = 750;

//!< 	//!< geometry
//!< 	L = 1.6850e-4;
//!< 	Rp = 8.5e-6;  //!< do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied MATLAB scripts. See the word document '2 overview of the code', section 'MATLAB setup before running the C++ code'
//!< 	Rn = 1.25e-5; //!< do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied MATLAB scripts. See the word document '2 overview of the code', section 'MATLAB setup before running the C++ code'
//!< 	SAV = 252.9915;
//!< 	elec_surf = 0.0982; //!< OCV fitting parameter

//!< 	//!< Stress parameters
//!< 	sparam = param::def::StressParam_Fit;

//!< 	//!< main Li reaction
//!< 	kp = 5e-11; //!< characterisation fitting parameter (at Tref)
//!< 	kp_T = 58000;
//!< 	kn = 1.7640e-11; //!< characterisation fitting parameter (at Tref)
//!< 	kn_T = 20000;
//!< 	//!< The diffusion coefficients at reference temperature are part of 'State'.
//!< 	//!< The values are set in the block of code below ('Initialise state variables')
//!< 	Dp_T = 29000;
//!< 	Dn_T = 35000;

//!< 	//!< spatial discretisation of the solid diffusion PDE
//!< 	M = MM;
//!< 	checkModelparam(); //!< check if the inputs to the MATLAB code are the same as the ones here in the C++ code

//!< 	//!< Initialise state variables
//!< 	slide::State_SPM::z_type up, un;
//!< 	double fp, fn, T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, R, delta_pl;
//!< 	double Rdc = 0.0102;																			  //!< DC resistance of the total cell in Ohm, characterisation fitting parameter
//!< 	fp = 0.689332;																					  //!< lithium fraction in the cathode at 50% soc [-], OCV fitting parameter
//!< 	fn = 0.479283;																					  //!< lithium fraction in the anode at 50% soc [-], OCV fitting parameter
//!< 	T = PhyConst::Kelvin + 25;																		  //!< cell temperature
//!< 	delta = 1e-9;																					  //!< SEI thickness. Start with a fresh cell, which has undergone some formation cycles so it has an initial SEI layer.
//!< 																									  //!< never start with a value of 0, because some equations have a term 1/delta, which would give nan or inf
//!< 																									  //!< so this will give errors in the code
//!< 	LLI = 0;																						  //!< lost lithium. Start with 0 so we can keep track of how much li we lose while cycling the cell
//!< 	thickp = 70e-6;																					  //!< thickness of the positive electrode, OCV fitting parameter
//!< 	thickn = 73.5e-6;																				  //!< thickness of the negative electrode, OCV fitting parameter
//!< 	ep = 0.5;																						  //!< volume fraction of active material in the cathode, OCV fitting parameter
//!< 	en = 0.5;																						  //!< volume fraction of active material in the anode, OCV fitting parameter
//!< 	ap = 3 * ep / Rp;																				  //!< effective surface area of the cathode, the 'real' surface area is the product of the effective surface area (a) with the electrode volume (elec_surf * thick)
//!< 	an = 3 * en / Rn;																				  //!< effective surface area of the anode
//!< 	CS = 0.01 * an * elec_surf * thickn;															  //!< initial crack surface. Start with 1% of the real surface area
//!< 	Dp = 8e-14;																						  //!< diffusion constant of the cathode at reference temperature, characterisation fitting parameter
//!< 	Dn = 7e-14;																						  //!< diffusion constant of the anode at reference temperature, characterisation fitting parameter
//!< 	R = Rdc * ((thickp * ap * elec_surf + thickn * an * elec_surf) / 2);							  //!< specific resistance of the combined electrodes, see State::iniStates
//!< 	delta_pl = 0;																					  //!< thickness of the plated lithium layer. You can start with 0 here
//!< 	s_ini.initialise(up, un, T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, R, delta_pl); //!< set the states, with a random value for the concentration
//!< 	s = s_ini;																						  //!< set the states, with a random value for the concentration

//!< 	//!< Check if this was a valid state
//!< 	try
//!< 	{
//!< 		validState();
//!< 	}
//!< 	catch (int e)
//!< 	{
//!< 		std::cout << "Error in State::initialise, one of the states has an illegal value, throwing an error\n";
//
//!< 		throw 12;
//!< 	}

//!< 	setC(fp, fn); //!< set the lithium concentration

//!< 	//!< SEI parameters
//!< 	nsei = 1;
//!< 	alphasei = 1;
//!< 	OCVsei = 0.4;
//!< 	rhosei = 100e3;
//!< 	rsei = 2037.4;
//!< 	Vmain = 13.0;
//!< 	Vsei = 64.39;
//!< 	c_elec0 = 4.541 / 1000;
//!< 	//!< fitting parameters of the models
//!< 	sei_p.sei1k = 0;
//!< 	sei_p.sei1k_T = 0;
//!< 	sei_p.sei2k = 0;
//!< 	sei_p.sei2k_T = 0;
//!< 	sei_p.sei2D = 0;
//!< 	sei_p.sei2D_T = 0;
//!< 	sei_p.sei3k = 0;
//!< 	sei_p.sei3k_T = 0;
//!< 	sei_p.sei3D = 0;
//!< 	sei_p.sei3D_T = 0;
//!< 	sei_p.sei_porosity = 0;

//!< 	//!< surface cracking
//!< 	//!< fitting parameters of the models
//!< 	csparam.CS1alpha = 0;
//!< 	csparam.CS2alpha = 0;
//!< 	csparam.CS3alpha = 0;
//!< 	csparam.CS4alpha = 0;
//!< 	csparam.CS4Amax = 0;
//!< 	csparam.CS5k = 0;
//!< 	csparam.CS5k_T = 0;
//!< 	csparam.CS_diffusion = 0;

//!< 	//!< LAM
//!< 	OCVnmc = 4.1;
//!< 	//!< fitting parameters
//!< 	lam_p.lam1p = 0;
//!< 	lam_p.lam1n = 0;
//!< 	lam_p.lam2ap = 0;
//!< 	lam_p.lam2bp = 0;
//!< 	lam_p.lam2an = 0;
//!< 	lam_p.lam2bn = 0;
//!< 	lam_p.lam2t = 0;
//!< 	lam_p.lam3k = 0;
//!< 	lam_p.lam3k_T = 0;
//!< 	lam_p.lam4p = 0;
//!< 	lam_p.lam4n = 0;

//!< 	//!< li-plating parameters
//!< 	npl = 1;
//!< 	alphapl = 1;
//!< 	OCVpl = 0;
//!< 	rhopl = 10000e3;
//!< 	//!< fitting parameters
//!< 	pl_p.pl1k = 0;
//!< 	pl_p.pl1k_T = 0;

//!< 	//!< degradation identifiers: no degradation
//!< deg.SEI_id.add_model(0); //!< no SEI growth, there is 1 SEI model (namely '0')
//!< deg.SEI_porosity = 0;	 //!< don't decrease the volume fraction

//!< deg.CS_id.add_model(0); //!< no surface cracks //!< there is 1 model (that there are no cracks)
//!< deg.CS_diffusion = 0;	//!< don't decrease the diffusion coefficient

//!< deg.LAM_id.add_model(0); //!< no LAM //!< there is 1 LAM model
//!< deg.pl_id = 0;			 //!< no lithium plating
//!< 	//!< Check if we will have to calculate the stress
//!< 	sparam.s_dai = false;	//!< according to Dai's stress model
//!< 	sparam.s_lares = false; //!< according to Laresgoiti's stress model
//!< }

void Cell_SPM::setOCVcurve(const std::string &namepos, const std::string &nameneg)
{
  /*
   * Function to load the specified OCV curves
   *
   * IN
   * namepos 	name of the CSV file in which the cathode OCV curve is defined
   * nameneg 	name of the CSV file in which the anode OCV curve is defined
   *
   */

  //!< Store the OCV curves
  try {
    OCV_curves.OCV_neg.setCurve(PathVar::data / nameneg); //!< the OCV curve of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
    OCV_curves.OCV_pos.setCurve(PathVar::data / namepos);
  } catch (int e) {
    //!< std::cout << "Throw test: " << 32 << '\n';
    std::cout << "ERROR in Cell_SPM::setOCVcurve when loading the OCV curves from the CSV files: "
              << e << ".\n";

    throw e;
  }
}

void Cell_SPM::setInitialConcentration(double cmaxp, double cmaxn, double lifracp, double lifracn)
{
  /*
   * Function to specify the lithiun fraction from which we have to start
   *
   * IN
   * cmaxp 	maximum li-concentration in the cathode [mol m-3]
   * cmaxn 	maximum li-concentration in the anode [mol m-3]
   * lifracp 	lithium fraction in the cathode at 50% SOC
   * lifracn 	lithium fraction in the anode at 50% SOC
   */

  //!< Store the maximum concentrations
  Cmaxpos = cmaxp;
  Cmaxneg = cmaxn;

  //!< set the concentration
  try {
    setC(lifracp, lifracn);
  } catch (int e) {
    //!< std::cout << "Throw test: " << 33 << '\n';
    std::cout << "Error in Cell_SPM::setInitialConcentration when setting the concentration: "
              << e << ".\n";
    throw e;
  }
}

void Cell_SPM::setGeometricParameters(double capnom, double surf, double ep, double en, double thickp, double thickn)
{
  /*
   * Function to set the parameters relating to the total amount of active material present in the cell.
   *
   * IN
   * capnom 		nominal capacity of the cell [Ah]
   * surf 		surface area of the electrodes [m2]
   * ep 			volume fraction of active material in the cathode [-]
   * en 			volume fraction of active material in the anode [-]
   * thickp 		thickness of the cathode [m]
   * thickn 		thickness of the anode [m]
   */

  setCapacity(capnom);
  geo.elec_surf = surf;

  //!< If you change the volume fraction, you also have to change the effective surface area
  double ap = 3 * ep / geo.Rp;
  double an = 3 * en / geo.Rn;

  //!< overwrite the original geometric parameters, use cell version
  overwriteGeometricStates(thickp, thickn, ep, en, ap, an);

  Vcell_valid = false;
}

void Cell_SPM::setCharacterisationParam(double Dp, double Dn, double kpi, double kni, double Rdc)
{
  /*
   * Function to set the parameters related to the characterisation of the cell
   *
   * IN
   * Dp	diffusion constant of the cathode at rate temperature [m s-1]
   * Dn 	diffusion constant of the anode at rate temperature [m s-1]
   * kp 	rate constant of the main lithium insertion reaction at the cathode at reference temperature
   * kn 	rate constant of the main lithiuim insertion reaction at the anode at reference temperature
   * Rdc 	DC resistance of the cell [Ohm]
   */

  //!< store the rate constants in the cell
  kp = kpi;
  kn = kni;

  //!< calculate the specific electrode resistance from the total DC resistance
  double r = Rdc * ((st.thickp() * st.ap() * geo.elec_surf + st.thickn() * st.an() * geo.elec_surf) / 2);

  //!< overwrite the cycling related parameters
  overwriteCharacterisationStates(Dp, Dn, r);
}

//!< void Cell_SPM::setRamping(double Istep, double tstep)
//!< {
//!< 	/*
//!< 	 * Function to set the ramping parameters of the cell.
//!< 	 * This allows to decrease them if you are trying weird things (e.g. very small diffusion constants) with the cell
//!< 	 *
//!< 	 * IN
//!< 	 * Istep 	maximum change in current per tstep [A], > 0
//!< 	 * tstep 	time step by which you change the current, [s], > 0
//!< 	 */

//!< 	dIcell = Istep;
//!< 	dt_I = tstep;
//!< }
} // namespace slide