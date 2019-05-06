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

#ifndef SRC_CELL_HPP_
#define SRC_CELL_HPP_

#include <string>
#include <cstring>
#include "State.hpp" 			// class that represents the state of a cell, the state is the collection of all time-varying conditions of the battery
#include "Model.h"				// defines a struct with the values for the matrices used in the spatial discretisation of the diffusion PDE

#define printCellFunctions 7	// threshold of verbose of when to print messages at the start and end of functions of the Cell

using namespace std;

// Define a structure with the identifications of which degradation model(s) to use
struct DEG_ID{
	// identifiers for degradation models
	// Each array is made with a length of 'len', which is the maximum number of models per mechanisms.
	// If the user wants to use more models, you only have to increase the value of 'len' to the number you want to use
	// and change the '10' in the definition of the arrays
	int len = 10;		// length of the arrays with identifications of which models to use
	int SEI_id[10];		// array with identifications to decide which SEI models to use. Max length 10
							/* 				0	no SEI growth
							 * 				1 	kinetic model only (Tafel kinetics)
							 * 				2 	Pinson&Bazant model: linear diffusion + Tafel kinetics
							 * 				3	Christensen and Newman model
							 */
	int SEI_n;			// SEI_N 	number of SEI models to use (length of SEI_ID)
	int SEI_porosity;	// integer deciding whether we reduce the active volume fraction due to SEI growth
							/* 				0	don't reduce it
							 * 				1	use correlation from Ashwin et al. 2016
							 */
	int CS_id[10];		// array with identifications for which model to use for surface cracking. Max length 10
							/* 				0 	no surface cracking
							 * 				1 	Laresgoiti's stress + crack growth model
							 * 				2 	Dai stress model + Laresgoiti crack growth
							 * 				3 	model based on Deshpande and Bernardi, 2017
							 * 				4 	model from Barai et al
							 * 				5 	model from Ekstrom et al
							 */
	int CS_n;			// number of surface crack growth models to use (length of CS_ID)
	int CS_diffusion;	// integer deciding whether we reduce the negative diffusion constant due to surface cracks
							/* 				0 	don't decrease diffusion
							 * 				1	decrease according to Barai et al. 2015
							 */
	int LAM_id[10];		// array with the integers deciding which models is to be used for loss of active material. Max length 10
							/* 				0 	no LAM
							 * 				1	Dai's stress model and Laresgoiti's correlation to get LAM
							 * 				2	delacourt's	correlation between abs(j) and porosity
							 * 				3 	Kindermann's model for cathode dissolution: tafel kinetics for increased porosity
							 * 				4 	Narayanrao's correlation which decreases the effective surface area proportionally to itself and j
							 */
	int LAM_n;			// number of LAM models to be used (length of LAM_id)
	int pl_id;			// integer deciding which model is to be used for li-plating
							/* 				0 	no plating
							 * 				1	Yang et al thermodynamic plating (Tafel kinetics)
							 */
};

// Define a structure with the fitting parameters of the SEI growth models (SEI)
struct SEIparam{

	double sei1k;		// rate parameter of the SEI side reaction in the 1st SEI model
	double sei1k_T;		// activation energy of sei1k

	double sei2k; 		// rate parameter of the SEI side reaction in the 2nd SEI model
	double sei2k_T; 	// activation energy of sei2k
	double sei2D; 		// diffusion constant of the SEI layer in the 2nd SEI model
	double sei2D_T; 	// activation energy of sei2D

	double sei3k; 		// rate parameter of the SEI side reaction in the 3rd SEI model
	double sei3k_T; 	// activation energy of sei3k
	double sei3D; 		// diffusion constant of the SEI layer in the 3rd SEI model
	double sei3D_T; 	// activation energy of sei3D

	double sei_porosity;// proportionality constant between the SEI growth and the decrease in volume fraction of active material
							// (because the SEI layer blocks the pores at the surface)
};

// Define a structure with the fitting parameters of the surface crack growth models (CS)
struct CSparam{

	double CS1alpha;	// fitting parameter of the 1st surface crack growth model

	double CS2alpha;	// fitting parameter of the 2nd surface crack growth model

	double CS3alpha;	// fitting parameter of the 3rd surface crack growth model

	double CS4Amax; 	// maximum crack growth surface for the 4th surface crack growth model
	double CS4alpha;	// fitting parameter of the 4th surface crack growth model

	double CS5k;		// rate parameter of the 5th surface crack growth model at reference temperature
	double CS5k_T;		// activation energy of CS5k

	double CS_diffusion; // fitting parameter to decrease the diffusion constant due to surface cracks
};

// Define a structure with the fitting parameters of the models for loss of active material (LAM)
struct LAMparam{
	double lam1p; 		// fitting parameter for the positive electrode for the 1st LAM model
	double lam1n; 		// fitting parameter for the negative electrode for the 1st LAM model

	double lam2ap; 		// fitting parameter 1 at reference temperature for the positive electrode for the 2nd LAM model
	double lam2bp; 		// fitting parameter 2 at reference temperature for the positive electrode for the 2nd LAM model
	double lam2an; 		// fitting parameter 1 at reference temperature for the negative electrode for the 2nd LAM model
	double lam2bn; 		// fitting parameter 2 at reference temperature for the negative electrode for the 2nd LAM model
	double lam2t;		// activation energy for all the parameters of the 2nd LAM model

	double lam3k; 		// rate constant at reference temperature for the cathode dissolution side reaction
	double lam3k_T;		// activation energy for lam3k

	double lam4p; 		// fitting parameter for the positive electrode for the 4th LAM model
	double lam4n; 		// fitting parameter for the negative electrode for the 4th LAM model
};

// Define a structure with the fitting parameters of the li-plating models (PL)
struct PLparam{
	double pl1k;		// rate constant of the li-plating side reaction at reference temperature in the 1st model
	double pl1k_T;		// activation energy of pl1k
};

class Cell {

protected: 						// protected such that child classes can access the class variables
	// battery states
	State s;					// the battery state, grouping all parameter which change over the battery's lifetime (see State.hpp)
	double Icell;				// current which is currently running through the cell [A]
								//		> 0 for discharge, < 0 for charge
	struct DEG_ID deg_id; 		// structure with the identification of which degradation model(s) to use
	int verbose;				// integer deciding how verbose the simulation should be
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
		double Cmaxpos; 		// maximum lithium concentration in the cathode [mol m-3]
		double Cmaxneg; 		// maximum lithium concentration in the anode [mol m-3]
		double C_elec; 			// Li- concentration in electrolyte [mol m-3]
		double n;				// number of electrons involved in the main reaction [-]
		double F; 				// Faraday's constant
		double Rg; 				// ideal gas constant
		double T_ref; 			// reference temperature [K]

	// parameters of the main li-insertion reaction
		double kp; 				// rate constant of main reaction at positive electrode at reference temperature
		double kp_T; 			// activation energy for the Arrhenius relation of kp
		double kn; 				// rate constant of main reaction at negative electrode at reference temperature
		double kn_T; 			// activation energy for the Arrhenius relation of kn
		// The diffusion constants at reference temperature are part of State because they can change over the battery's lifetime
		double Dp_T; 			// activation energy for the Arrhenius relation of Dp
		double Dn_T; 			// activation energy for the Arrhenius relation of Dn

	// Thermal model parameters
		double T_env; 			// environment temperature [K]
		double Qch; 			// convective heat transfer coefficient per volume [W K-1 m-3]
		double rho; 			// density of the battery
		double Cp; 				// thermal capacity of the battery

	// Cell usage parameters
		double nomCapacity; 	// nominal battery capacity [Ah]
		double Vmax; 			// maximum cell voltage [V]
		double Vmin; 			// minimal cell voltage [V]
		double dIcell;			// maximum (absolute value of the) change in current per time step dt_I [A s-1],  used for ramping the current up or down
		double dt_I;			// time step used to ramp the current up or down [s]
								// e.g. if dIcell = 0.1 and dt_I = e-3, then the current is ramped up or down at 0.1A per milisecond

	// Geometric parameters
		double L; 				// thickness of the cell [m]
		double elec_surf; 		// geometric surface area of the electrodes (electrode height * electrode width) [m2]
		double SAV; 			// surface area to volume-ratio of the cell [m2/m3]
		double Rp; 				// radius of the positive sphere of the Single Particle model [m]
		double Rn; 				// radius of the negative sphere of the Single Particle model [m]
		// other geometric parameters are part of State because they can change over the battery's lifetime

	// Constants for the stress model
		double omegap;			// partial molar volume of positive electrode [m3 mol-1]
		double omegan;			// partial molar volume of negative electrode [m3 mol-1]
		double Ep;				// Young's modulus of positive electrode [GPa]
		double En;				// Young's modulus of negative electrode [GPa]
		double nup;				// Poisson's ratio of positive electrode [-]
		double nun;				// Poisson's ratio of negative electrode [-]
		// values of the stress are often needed. Because it takes very long to calculate them, we calculate them once and store them so we don't need to repeat the same calculation twice
		bool s_dai;				// do we need to calculate the stress according to Dai's model?
		bool s_lares;			// do we need to calculate the stress according to Laresgoiti's model?
		bool s_dai_update;		// boolean to indicate if Dai's stress are up to date with the battery state at this time step
		bool s_lares_update;	// boolean to indicate if Dai's stress are up to date with the battery state at this time step
		double s_dai_p;			// maximum hydrostatic stress in the positive particle according to Dai's stress model
		double s_dai_n;			// maximum hydrostatic stress in the negative particle according to Dai's stress model
		double s_lares_n;		// stress in the negative particle according to Laresgoiti's stress model
		double s_dai_p_prev;	// maximum hydrostatic stress in the previous time step in the positive particle according to Dai's stress model
		double s_dai_n_prev;	// maximum hydrostatic stress in the previous time step in the negative particle according to Dai's stress model
		double s_lares_n_prev;	// stress in the previous time step in the negative particle according to Laresgoiti's stress model

	// Constants and parameters for the SEI growth model
		double nsei; 			// number of electrons involved in the SEI reaction [-]
		double alphasei; 		// charge transfer coefficient of the SEI reaction [-]
		double OCVsei; 			// equilibrium potential of the SEI side reaction [V]
		double rhosei; 			// partial molar volume of the SEI layer [m3 mol-1]
		double Rsei; 			// resistance of the SEI layer [Ohm/m]
		double c_elec0;			// bulk concentration of the electrolyte molecule participating in the SEI growth (e.g. EC) [mol m-3]
		double Vmain;			// partial molar volume of the main reaction, see Ashwin et al, 2016
		double Vsei; 			// partial molar volume of the SEI side reaction, see Ashwin et al., 2016
		struct SEIparam seiparam;// structure with the fitting parameters of the different SEI growth models

	// surface crack parameters & constants
		struct CSparam csparam; // structure with the fitting parameters of the different crack growth models

	// LAM parameters & constants
		double OCVnmc;			// equilibrium potential of the NMC dissolution side reaction [V]
		struct LAMparam lamparam; // structure with the fitting parameters of the different LAM models

	// Li-plating parameters & constants
		double npl;				// number of electrons involved in the plating reaction [-]
		double alphapl;			// charge transfer constant for the plating reaction [-]
		double OCVpl;			// OCV of the plating reaction [V]
		double rhopl;			// density of the plated lithium layer
		struct PLparam plparam; // structure with the fitting parameters of the different plating models

	// Matrices for spatial discretisation of the solid diffusion model
		struct Model M;

	// OCV curves
		// The OCV curves have to be stored in arrays. The maximum length of the arrays is 1000.
		// If you exceed the length of the arrays, the model will crash
		int OCV_maxLength = 1000; // maximum length of the OCV curves (length of the arrays in which the curves will be stored)
		int OCV_pos_n; 			// number of data points in the OCV curve for the cathode
		double OCV_pos_x[1000];	// lithium fractions of the points of the cathode OCV curve
		double OCV_pos_y[1000];	// voltage vs li/li+ of the points of the cathode OCV curve [V]

		int OCV_neg_n;			// number of data points in the OCV curve for the anode
		double OCV_neg_x[1000];	// lithium fractions of the points of the anode OCV curve
		double OCV_neg_y[1000];	// voltage vs li/li+ of the points of the anode OCV curve [V]

		int dOCV_neg_n;			// number of data points in the entropic coefficient for the anode
		double dOCV_neg_x[1000];// lithium fractions of the points of the anode entropic coefficient curve
		double dOCV_neg_y[1000];// entropic coefficient curve [V K-1]

		int dOCV_tot_n;			// number of data points in the entropic coefficient for the entire cell OCV curve
		double dOCV_tot_x[1000];// cathodic lithium fractions of the points of the entire cell's entropic coefficient
		double dOCV_tot_y[1000];// the entire cell's entropic coefficient [V K-1]

	// Functions

		// Get and set the cell's state using an array-representation
		// This is for internal use only, the states can be accessed externally using the State-representation
		virtual void getStates(int n, double states[]);			// get an array with the cell's states
		virtual void setStates(int n, double states[]);			// set the cell's states to the states in the array

		// degradation models
		virtual void SEI(double OCVnt, double etan, double* isei, double* den); // calculate the effect of SEI growth
		virtual void CS(double OCVnt, double etan, double* isei_multiplyer, double* dCS, double* dDn); // calculate the effect of surface crack growth
		virtual void LAM(bool critical, double zp_surf, double etap, double* dthickp, double* dthickn, double* dap, double* dan, double* dep, double* den); // calculate the effect of LAM
		virtual void LiPlating(double OCVnt, double etan, double* isei); // calculate the effect of lithium plating

		// Calculate the time derivatives of the states at the actual cell current (state-space model)
		virtual void dState(bool critical, bool blockDegradation, int electr, int length, double dstates[]);

public:

	// getters
	virtual double getNominalCap(); 							// get the nominal cell capacity (e.g. to convert from Crate to Amperes)
	virtual double getVmax();									// get the maximum cell voltage
	virtual double getVmin();									// get the minimum cell voltage
	virtual double getT();										// get the cell temperature
	virtual double getTenv();									// get the environmental temperature
	virtual void getTemperatures(double* Tenv, double* Tref); 	// get the environmental and reference temperature
	virtual double getR();										// get the total cell DC resistance
	virtual double getAnodeSurface();							// get the anode surface area (without cracks)
	virtual double getI();										// get the cell's current
	virtual void getStates(std::State& si, double* I);			// get a State object with the cell's states and the cell current
	virtual void getCSurf(double* cps, double* cns); 			// get the surface concentrations
	virtual void getC(int n, double cp[], double cn[]); 		// get the concentrations at all nodes
	virtual bool getVoltage(bool print, double* V, double* OCVp, double* OCVn, double* etap, double* etan, double* Rdrop, double* Temp);// get the cell's voltage
	virtual void getDaiStress(int n, double* sigma_p, double* sigma_n, double sigma_r_p[], double sigma_r_n[],double sigma_t_p[], double sigma_t_n[], double sigma_h_p[], double sigma_h_n[]);// get the stresses at all nodes according to Dai's stress model
	virtual void updateDaiStress();								// updated the stored stress values for Dai's stress model
	virtual void getLaresgoitiStress(bool print,  double* sigma_n); // get the stresses at all nodes according to Laresgoiti's stress model
	virtual void updateLaresgoitiStress(bool print);			// update the stored stress values for Laresgoiti's stress model

	// setters
	virtual void setVlimits(double VMAX, double VMIN);			// set the voltage limits of the cell
	virtual void setT(double T);								// set the cell's temperature
	virtual void setTenv(double Tenv);							// set the environmental temperature
	virtual void setStates(std::State si, double I);			// set the cell's states to the states in the State object and the cell current to the given value
	virtual void setC(double cp0, double cn0);					// set the concentrations to the given (uniform) concentration
	virtual void setI(bool critical, bool check, double I);		// set the cell's current to the specified value

	// time integration
	virtual void ETI(bool print, double dti, bool blockDegradation); // step forward in time using forward Eurler time integration
	virtual void ETI_electr(bool print, double I, double dti, bool blockDegradation, bool pos); // step forward with only one electrode using forward Euler time integration
};

#endif /* SRC_CELL_HPP_ */
