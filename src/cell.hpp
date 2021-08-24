/*
 * Cell.hpp
 *
 * Header file for the parent Class for all Cells.
 *
 * It also defines a number of structs:
 * 		DEG_ID 		settings about which degradation models should be used for the simulations
 * 		SEIparam 	fitting parameters of the various models for SEI growth
 * 		CSparam 	fitting parameters of the various models for crach growth on the surface of the electrodes
 * 		LAMparam 	fitting parameters of the various models for loss of active material
 * 		PLparam 	fitting parameters of the various models for lithium plating
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <vector>
#include <array>
#include <iostream>

#include "state.hpp" // class that represents the state of a cell, the state is the collection of all time-varying conditions of the battery
#include "model.h"	 // defines a struct with the values for the matrices used in the spatial discretisation of the diffusion PDE
#include "read_CSVfiles.h"
#include "constants.hpp"
#include "interpolation.h"
#include "param/cell_param.hpp"

//#include <string>

using sigma_type = std::array<double, settings::nch + 2>;

// Free functions:

// State related functions
void validState(slide::State &s, slide::State &s_ini);

class Cell
{

protected: // protected such that child classes can access the class variables
	// battery states
	slide::State s;		  // the battery state, grouping all parameter which change over the battery's lifetime (see State.hpp)
	slide::State s_ini;	  // battery initial state, moved to here.
	double Icell;		  // current which is currently running through the cell [A]
						  //		> 0 for discharge, < 0 for charge
	struct DEG_ID deg_id; // structure with the identification of which degradation model(s) to use
	int verbose;		  // integer deciding how verbose the simulation should be
						  // The higher the number, the more output there is.
						  // Recommended level is 1, only use higher levels for debugging
						  // From level 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
						  // 	0 	almost no messages are printed, only in case of critical errors related to illegal parameters
						  // 	1 	error messages are printed in case of critical errors which might crash the simulation
						  // 	2 	all error messages are printed, whether the simulation can recover from the errors or not
						  // 	3 	on top of the output from 2, a message is printed every time a function in the Cycler and BasicCycler is started and terminated
						  // 	4 	on top of the output from 3, the high-level flow of the program in the Cycler is printed (e.g. 'we are going to discharge the cell')
						  // 	5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
						  // 	6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
						  // 	7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated

	// Battery model constants
	double Cmaxpos; // maximum lithium concentration in the cathode [mol m-3]
	double Cmaxneg; // maximum lithium concentration in the anode [mol m-3]
	double C_elec;	// Li- concentration in electrolyte [mol m-3]
	double n;		// number of electrons involved in the main reaction [-]
	double T_ref;	// reference temperature [K]

	// parameters of the main li-insertion reaction
	double kp;	 // rate constant of main reaction at positive electrode at reference temperature
	double kp_T; // activation energy for the Arrhenius relation of kp
	double kn;	 // rate constant of main reaction at negative electrode at reference temperature
	double kn_T; // activation energy for the Arrhenius relation of kn
	// The diffusion constants at reference temperature are part of State because they can change over the battery's lifetime
	double Dp_T; // activation energy for the Arrhenius relation of Dp
	double Dn_T; // activation energy for the Arrhenius relation of Dn

	// Thermal model parameters
	double T_env; // environment temperature [K]
	double Qch;	  // convective heat transfer coefficient per volume [W K-1 m-3]
	double rho;	  // density of the battery
	double Cp;	  // thermal capacity of the battery

	// Cell usage parameters
	double nomCapacity; // nominal battery capacity [Ah]
	double Vmax;		// maximum cell voltage [V]
	double Vmin;		// minimal cell voltage [V]
	double dIcell;		// maximum (absolute value of the) change in current per time step dt_I [A s-1],  used for ramping the current up or down
	double dt_I;		// time step used to ramp the current up or down [s]
						// e.g. if dIcell = 0.1 and dt_I = e-3, then the current is ramped up or down at 0.1A per milisecond

	// Geometric parameters
	double L;		  // thickness of the cell [m]
	double elec_surf; // geometric surface area of the electrodes (electrode height * electrode width) [m2]
	double SAV;		  // surface area to volume-ratio of the cell [m2/m3]
	double Rp;		  // radius of the positive sphere of the Single Particle model [m]
	double Rn;		  // radius of the negative sphere of the Single Particle model [m]
	// other geometric parameters are part of State because they can change over the battery's lifetime

	struct StressParam sparam; // Stress parameters.

	// Constants and parameters for the SEI growth model
	double nsei;			  // number of electrons involved in the SEI reaction [-]
	double alphasei;		  // charge transfer coefficient of the SEI reaction [-]
	double OCVsei;			  // equilibrium potential of the SEI side reaction [V]
	double rhosei;			  // partial molar volume of the SEI layer [m3 mol-1]
	double Rsei;			  // resistance of the SEI layer [Ohm/m]
	double c_elec0;			  // bulk concentration of the electrolyte molecule participating in the SEI growth (e.g. EC) [mol m-3]
	double Vmain;			  // partial molar volume of the main reaction, see Ashwin et al, 2016
	double Vsei;			  // partial molar volume of the SEI side reaction, see Ashwin et al., 2016
	struct SEIparam seiparam; // structure with the fitting parameters of the different SEI growth models

	// surface crack parameters & constants
	struct CSparam csparam; // structure with the fitting parameters of the different crack growth models

	// LAM parameters & constants
	double OCVnmc;			  // equilibrium potential of the NMC dissolution side reaction [V]
	struct LAMparam lamparam; // structure with the fitting parameters of the different LAM models

	// Li-plating parameters & constants
	double npl;				// number of electrons involved in the plating reaction [-]
	double alphapl;			// charge transfer constant for the plating reaction [-]
	double OCVpl;			// OCV of the plating reaction [V]
	double rhopl;			// density of the plated lithium layer
	struct PLparam plparam; // structure with the fitting parameters of the different plating models

	// Matrices for spatial discretisation of the solid diffusion model
	struct slide::Model M;

	// OCV curves
	struct OCVcurves OCV_curves;
	// Functions

	// Get and set the cell's state using an array-representation
	// This is for internal use only, the states can be accessed externally using the State-representation
	const auto &getStates() // get an array with the cell's states
	{
		/*
	 	* Returns all states which describe the status of the cell in the array.
	 	* OUT
	 	* states array of length ns (defined on top of Constants.hpp)
	 	*/
		return s.getStates_arr();
	}

	void setStates(slide::states_type &&states); // set the cell's states to the states in the array

	// degradation models
	void SEI(double OCVnt, double etan, double *isei, double *den);																				// calculate the effect of SEI growth
	void CS(double OCVnt, double etan, double *isei_multiplyer, double *dCS, double *dDn);														// calculate the effect of surface crack growth
	void LAM(bool critical, double zp_surf, double etap, double *dthickp, double *dthickn, double *dap, double *dan, double *dep, double *den); // calculate the effect of LAM
	void LiPlating(double OCVnt, double etan, double *isei);																					// calculate the effect of lithium plating

	// Calculate the time derivatives of the states at the actual cell current (state-space model)
	slide::states_type dState(bool critical, bool blockDegradation, int electr);

public:
	// Constructor
	Cell(const std::string &_namepos, int _OCV_pos_n,
		 const std::string &_nameneg, int _OCV_neg_n,
		 const std::string &_nameentropicC, int _dOCV_neg_n,
		 const std::string &_nameentropicCell, int _dOCV_tot_n)

		: OCV_curves(_namepos, _OCV_pos_n,
					 _nameneg, _OCV_neg_n,
					 _nameentropicC, _dOCV_neg_n,
					 _nameentropicCell, _dOCV_tot_n)
	{
	}

	Cell() = default;											  // Default constructor.
																  // getters
	double getNominalCap() const noexcept { return nomCapacity; } // returns nominal capacity of the cell [Ah]	(e.g. to convert from Crate to Amperes)
	double getVmax() const noexcept { return Vmax; }			  // returns the maximum voltage of this cell [V] above which operation is not allowed
	double getVmin() const noexcept { return Vmin; }			  // returns the minimum voltage of this cell [V] below which operation is not allowed
	double getT() noexcept { return s.get_T(); }				  // returns the uniform battery temperature in [K]
	double getTenv() const noexcept { return T_env; }			  // get the environmental temperature [K]

	void getTemperatures(double *Tenv, double *Tref) noexcept // get the environmental and reference temperature
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

	double getR() noexcept // get the total cell DC resistance
	{
		/*
	 	* Return the total cell DC resistance [Ohm]
	 	* The total cell resistance is the sum of the resistance of the electrodes and the SEI layer
	 	* the resistance of the electrodes increases due to loss of active material (LAM)
	 	* the resistance of the SEI layer increases as the layer becomes thicker
	 	* it is assumed the plated lithium does not add resistance
	 	*/
		return (Rsei * s.get_delta() + s.getR(elec_surf));
	}

	double getAnodeSurface() noexcept { return s.get_an() * s.get_thickn() * elec_surf; } // get the anode pure surface area (without cracks) product of the effective surface area (an) with the electrode volume

	double getI() { return Icell; } // get the cell current [A]  positive for discharging

	void getStates(slide::State &si, double *I);

	void getCSurf(double *cps, double *cns);																					 // get the surface concentrations
	void getC(double cp[], double cn[]);																						 // get the concentrations at all nodes
	bool getVoltage(bool print, double *V, double *OCVp, double *OCVn, double *etap, double *etan, double *Rdrop, double *Temp); // get the cell's voltage

	void getDaiStress(double *sigma_p, double *sigma_n, sigma_type &sigma_r_p, sigma_type &sigma_r_n, sigma_type &sigma_t_p, sigma_type &sigma_t_n,
					  sigma_type &sigma_h_p, sigma_type &sigma_h_n) noexcept; // get the stresses at all nodes according to Dai's stress model
	void updateDaiStress() noexcept;										  // updated the stored stress values for Dai's stress model
	void getLaresgoitiStress(bool print, double *sigma_n);					  // get the stresses at all nodes according to Laresgoiti's stress model
	void updateLaresgoitiStress(bool print);								  // update the stored stress values for Laresgoiti's stress model

	// setters
	void setVlimits(double VMAX, double VMIN);		  // set the voltage limits of the cell
	void setT(double T);							  // set the cell's temperature
	void setTenv(double Tenv);						  // set the environmental temperature
	void setStates(const slide::State &si, double I); // set the cell's states to the states in the State object and the cell current to the given value
	void setC(double cp0, double cn0);				  // set the concentrations to the given (uniform) concentration
	void setI(bool critical, bool check, double I);	  // set the cell's current to the specified value

	// State related functions
	void validState() { ::validState(s, s_ini); }
	void overwriteCharacterisationStates(double Dpi, double Dni, double ri)
	{
		// Overwrite both current and initial states.
		s.overwriteCharacterisationStates(Dpi, Dni, ri);
		s_ini.overwriteCharacterisationStates(Dpi, Dni, ri);
	}

	void overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani)
	{
		// Overwrite both current and initial states.
		s.overwriteGeometricStates(thickpi, thickni, epi, eni, api, ani);
		s_ini.overwriteGeometricStates(thickpi, thickni, epi, eni, api, ani);
	}

	// time integration
	void ETI(bool print, double dti, bool blockDegradation);							// step forward in time using forward Eurler time integration
	void ETI_electr(bool print, double I, double dti, bool blockDegradation, bool pos); // step forward with only one electrode using forward Euler time integration

	// Utility
	void checkModelparam(); // check if the inputs to the Matlab code are the same as the ones here in the C++ code
};