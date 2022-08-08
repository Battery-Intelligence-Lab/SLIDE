/*
 * cell_user.cpp
 *
 * One of the child classes that implements a cell where the user can set this/her own parameters.
 * The initial values are the ones for the high-power Kokam cell
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "cell_user.hpp"

#include <cmath>
#include <iostream>

#include "read_CSVfiles.hpp"
#include "state.hpp"
#include "param/param_default.hpp"

namespace slide
{

	Cell_user::Cell_user(const slide::Model_SPM &MM, int verbosei)
		: Cell(settings::path::Kokam::namepos, settings::path::Kokam::nameneg, settings::path::Kokam::nameentropicC, settings::path::Kokam::nameentropicCell)

	{
		/*
		 * Standard constructor to initialise the battery parameters
		 *
		 * IN
		 * M 			structure of the type Model, with the matrices of spatial discretisation for solid diffusion
		 * verbosei		integer indicating how verbose the simulation has to be.
		 * 				The higher the number, the more output there is.
		 * 				Recommended value is 1, only use higher values for debugging
		 * 				From 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
		 * 					0 	almost no messages are printed, only in case of critical errors related to illegal parameters
		 * 					1 	error messages are printed in case of critical errors which might crash the simulation
		 * 					2 	all error messages are printed, whether the simulation can recover from the errors or not
		 * 					3 	on top of the output from 2, a message is printed every time a function in the Cycler and BasicCycler is started and terminated
		 * 					4 	on top of the output from 3, the high-level flow of the program in the Cycler is printed (e.g. 'we are going to discharge the cell')
		 * 					5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
		 * 					6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
		 * 					7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated
		 *
		 * THROWS
		 * 110 	the matrices for the solid diffusion discretisation, produced by MATLAB, are wrong
		 * 111 	the OCV curves are too long
		 */

		verbose = verbosei;

		// maximum concentrations
		Cmaxpos = 51385; // value for NMC
		Cmaxneg = 30555; // value for C
		C_elec = 1000;	 // standard concentration of 1 molar

		// constants
		n = 1;

		// Cell parameters
		nomCapacity = 2.7;
		Vmax = 4.2;	  // value for an NMC/C cell
		Vmin = 2.7;	  // value for an NMC/C cell
		Icell = 0;	  // initial cell current is 0A
		dIcell = 1.0; // ramp at 1A
		dt_I = 1e-2;  // ramp at 10ms so changing the current goes fast
					  // now changing the current takes 0.01 second per A

		// thermal parameters
		T_ref = 273 + 25;
		T_env = 273 + 25;
		Qch = 90; // 90 gives very good cooling, as if there is a fan pointed at the cell. values of 30-60 are representative for a cell on a shelf without forced cooling
		rho = 1626;
		Cp = 750;

		// geometry
		L = 1.6850e-4;
		Rp = 8.5e-6;  // do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied MATLAB scripts. See the word document '2 overview of the code', section 'MATLAB setup before running the C++ code'
		Rn = 1.25e-5; // do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied MATLAB scripts. See the word document '2 overview of the code', section 'MATLAB setup before running the C++ code'
		SAV = 252.9915;
		elec_surf = 0.0982; // OCV fitting parameter

		// Stress parameters
		sparam = param::def::StressParam_User;

		// main Li reaction
		kp = 5e-11; // characterisation fitting parameter (at Tref)
		kp_T = 58000;
		kn = 1.7640e-11; // characterisation fitting parameter (at Tref)
		kn_T = 20000;
		// The diffusion coefficients at reference temperature are part of 'State'.
		// The values are set in the block of code below ('Initialise state variables')
		Dp_T = 29000;
		Dn_T = 35000;

		// spatial discretisation of the solid diffusion PDE
		M = MM;
		checkModelparam(); // check if the inputs to the MATLAB code are the same as the ones here in the C++ code
		// Initialise state variables
		State_SPM::z_type up, un;
		double fp, fn, T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, R, delta_pl;
		double Rdc = 0.0102;																			  // DC resistance of the total cell in Ohm, characterisation fitting parameter
		fp = 0.689332;																					  // lithium fraction in the cathode at 50% soc [-], OCV fitting parameter
		fn = 0.479283;																					  // lithium fraction in the anode at 50% soc [-], OCV fitting parameter
		T = PhyConst::Kelvin + 25;																		  // cell temperature
		delta = 1e-9;																					  // SEI thickness. Start with a fresh cell, which has undergone some formation cycles so it has an initial SEI layer.
																										  // never start with a value of 0, because some equations have a term 1/delta, which would give nan or inf
																										  // so this will give errors in the code
		LLI = 0;																						  // lost lithium. Start with 0 so we can keep track of how much li we lose while cycling the cell
		thickp = 70e-6;																					  // thickness of the positive electrode, OCV fitting parameter
		thickn = 73.5e-6;																				  // thickness of the negative electrode, OCV fitting parameter
		ep = 0.5;																						  // volume fraction of active material in the cathode, OCV fitting parameter
		en = 0.5;																						  // volume fraction of active material in the anode, OCV fitting parameter
		ap = 3 * ep / Rp;																				  // effective surface area of the cathode, the 'real' surface area is the product of the effective surface area (a) with the electrode volume (elec_surf * thick)
		an = 3 * en / Rn;																				  // effective surface area of the anode
		CS = 0.01 * an * elec_surf * thickn;															  // initial crack surface. Start with 1% of the real surface area
		Dp = 8e-14;																						  // diffusion constant of the cathode at reference temperature, characterisation fitting parameter
		Dn = 7e-14;																						  // diffusion constant of the anode at reference temperature, characterisation fitting parameter
		R = Rdc * ((thickp * ap * elec_surf + thickn * an * elec_surf) / 2);							  // specific resistance of the combined electrodes, see State::iniStates
		delta_pl = 0;																					  // thickness of the plated lithium layer. You can start with 0 here
		s_ini.initialise(up, un, T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, R, delta_pl); // set the states, with a random value for the concentration
		s = s_ini;																						  // set the states, with a random value for the concentration

		// Check if this was a valid state
		try
		{
			validState();
		}
		catch (int e)
		{
			std::cout << "Error in State::initialise, one of the states has an illegal value, throwing an error.\n";
			std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
			throw 12;
		}
		setC(fp, fn); // set the lithium concentration

		// SEI parameters
		nsei = 1;
		alphasei = 1;
		OCVsei = 0.4;
		rhosei = 100e3;
		rsei = 2037.4;
		Vmain = 13.0;
		Vsei = 64.39;
		c_elec0 = 4.541 / 1000;
		// fitting parameters of the models
		sei_p = param::def::SEIparam_User;

		// surface cracking
		// fitting parameters of the models
		csparam.CS1alpha = 4.25e-5;
		csparam.CS2alpha = 6.3e-7;
		csparam.CS3alpha = 2.31e-16;
		csparam.CS4alpha = 4.3306e-8;
		csparam.CS4Amax = 5 * getAnodeSurface(); // assume the maximum crack surface is 5 times the initial anode surface
		csparam.CS5k = 1e-18;
		csparam.CS5k_T = -127040;
		csparam.CS_diffusion = 2;

		// LAM
		OCVnmc = 4.1;
		// fitting parameters
		lam_p = param::def::LAMparam_User;

		// li-plating parameters
		npl = 1;
		alphapl = 1;
		OCVpl = 0;
		rhopl = 10000e3;
		// fitting parameters
		plparam.pl1k = 4.5e-10;
		plparam.pl1k_T = -2.014008e5;

		// degradation identifiers: no degradation
		deg.SEI_id.add_model(0); // no SEI growth, there is 1 SEI model (namely '0')
		deg.SEI_porosity = 0;	 // don't decrease the volume fraction

		deg.CS_id.add_model(0); // no surface cracks // there is 1 model (that there are no cracks)
		deg.CS_diffusion = 0;	// don't decrease the diffusion coefficient

		deg.LAM_id.add_model(0); // no LAM // there is 1 LAM model
		deg.pl_id = 0;			 // no lithium plating

		// Check if we will have to calculate the stress according to Dai's stress model
		for (auto cs_id : deg_id.CS_id)
			sparam.s_dai = sparam.s_dai || cs_id == 2;

		for (auto lam_id : deg_id.LAM_id)
			sparam.s_dai = sparam.s_dai || lam_id == 1;

		// check if we need to calculate the stress according to Laresgoiti's stress model
		for (auto cs_id : deg_id.CS_id)
			sparam.s_lares = sparam.s_lares || cs_id == 1;
	}

	Cell_user::Cell_user(const slide::Model_SPM &M, const DEG_ID &deg_id, int verbosei) : Cell_user(M, verbosei)
	{
		/*
		 * constructor to initialise the degradation parameters
		 *
		 * IN
		 * M 		structure of the type Model, with the matrices of spatial discretisation for solid diffusion
		 * degid 	structure of the type DEG_ID, with the identifications of which degradation model(s) to use
		 */

		// Check if we will have to calculate the stress according to Dai's stress model
		for (auto cs_id : deg_id.CS_id)
			sparam.s_dai = sparam.s_dai || cs_id == 2;

		for (auto lam_id : deg_id.LAM_id)
			sparam.s_dai = sparam.s_dai || lam_id == 1;

		// check if we need to calculate the stress according to Laresgoiti's stress model
		for (auto cs_id : deg_id.CS_id)
			sparam.s_lares = sparam.s_lares || cs_id == 1;
	}

} // namespace slide
