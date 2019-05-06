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

#include "Cell_Fit.hpp"
#include "ReadCSVfiles.h"

#include <cassert>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <math.h>

using namespace std;

Cell_Fit::Cell_Fit(const Model& MM, int verbosei) {
	/*
	 * Standard constructor to initialise the battery parameters.
	 * Nothing should be changed in this constructor unless the user changes the radius of the particles.
	 * The values used for fitting can be changed using the functions defined in this class
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
	 * 110 	the matrices for the solid diffusion discretisation, produced by Matlab, are wrong
	 * 111 	the OCV curves are too long
	 */

	verbose = verbosei;

	// maximum concentrations
	Cmaxpos = 51385;
	Cmaxneg = 30555;
	C_elec = 1000;						// standard concentration of 1 molar

	// constants
	n = 1;
	F = 96487;
	Rg = 8.314;

	// Cell parameters
	nomCapacity = 2.7;
	Vmax = 4.2;							// value for an NMC/C cell
	Vmin = 2.7;							// value for an NMC/C cell
	Icell = 0;							// initial cell current is 0A
	dIcell = 0.1;						// ramp at 0.1A
	dt_I = pow(10,-3);					// ramp at 1ms so changing the current goes fast
										// now changing the current takes 0.01 second per A

	// thermal parameters
	T_ref = 273+25;
	T_env = 273+25;
	Qch = 90;							// 90 gives very good cooling, as if there is a fan pointed at the cell. values of 30-60 are representative for a cell on a shelf without forced cooling
	rho = 1626;
	Cp = 750;

	// geometry
	L = 1.6850*pow(10,-4);
	Rp = 8.5*pow(10,-6);				// do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied Matlab scripts. See the word document '2 overview of the code', section 'Matlab setup before running the C++ code'
	Rn = 1.25*pow(10,-5);				// do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied Matlab scripts. See the word document '2 overview of the code', section 'Matlab setup before running the C++ code'
	SAV = 252.9915;
	elec_surf = 0.0982;					// OCV fitting parameter

	// Stress parameters
	omegap = 2.1 * pow(10,-6);			// value from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
	omegan = 3.17 * pow(10,-6); 		// value from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
	Ep = 138.73;						// value from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
	En = 10;							// value from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
	nup = 0.3;			 				// value from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
	nun = 0.3;							// value from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
	s_dai = false;
	s_lares = false;
	s_dai_update = false;
	s_lares_update = false;
	s_dai_p = 0;
	s_dai_n = 0;
	s_lares_n = 0;
	s_dai_p_prev = 0;
	s_dai_n_prev = 0;
	s_lares_n_prev = 0;

	// main Li reaction
	kp = 5*pow(10,-11);					// characterisation fitting parameter (at Tref)
	kp_T = 58000;
	kn = 1.7640*pow(10,-11);			// characterisation fitting parameter (at Tref)
	kn_T = 20000;
	// The diffusion coefficients at reference temperature are part of 'State'.
	// The values are set in the block of code below ('Initialise state variables')
	Dp_T = 29000;
	Dn_T = 35000;

	// OCV curves
	string namepos = "Kokam_OCV_NMC.csv"; 						// name of the csv file with the cathode OCV curve
	string nameneg = "Kokam_OCV_C.csv";							// name of the csv file with the anode OCV curve
	string nameentropicC = "Kokam_entropic_C.csv"; 				// name of the csv file with the entropic coefficient curve for the anode OCV curve
	string nameentropicCell = "Kokam_entropic_cell.csv";		// name of the csv file with the entropic coefficient curve for the total cell OCV curve
	OCV_maxLength = 1000;
	OCV_pos_n = 49;
	OCV_neg_n = 63;
	dOCV_neg_n = 11;
	dOCV_tot_n = 11;
	// check that the OCV curves are shorter than the allowed range
	bool leng = false;
	if(OCV_pos_n >= OCV_maxLength){
		cerr<<"ERROR in the constructor of Cell_Fit. The cathode OCV curve has "<<OCV_pos_n<<" points while the maximum number is "<<OCV_maxLength<<". Throwing an error"<<endl<<flush;
		leng = true;
	}
	if(OCV_neg_n >= OCV_maxLength){
		cerr<<"ERROR in the constructor of Cell_Fit. The anode OCV curve has "<<OCV_neg_n<<" points while the maximum number is "<<OCV_maxLength<<". Throwing an error"<<endl<<flush;
		leng = true;
	}
	if(dOCV_neg_n >= OCV_maxLength){
		cerr<<"ERROR in the constructor of Cell_Fit. The anode entropic coefficient curve has "<<dOCV_neg_n<<" points while the maximum number is "<<OCV_maxLength<<". Throwing an error"<<endl<<flush;
		leng = true;
	}
	if(dOCV_tot_n >= OCV_maxLength){
		cerr<<"ERROR in the constructor of Cell_Fit. The cell entropic coefficient curve has "<<dOCV_tot_n<<" points while the maximum number is "<<OCV_maxLength<<". Throwing an error"<<endl<<flush;
		leng = true;
	}
	if(leng)
		throw 111;
		// if you get this error and you do want to use OCV curves which are longer than 1000 points, increase the size of the OCV arrays in Cell.hpp and the value of OCV_maxLength
	// Read the csv files
	try{
		loadCSV_2col(nameneg,OCV_neg_n,OCV_neg_x, OCV_neg_y); 					// the OCV curve of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(namepos,OCV_pos_n,OCV_pos_x, OCV_pos_y);					// the OCV curve of the cathode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(nameentropicC, dOCV_neg_n, dOCV_neg_x, dOCV_neg_y); 		// the entropic coefficient of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
		loadCSV_2col(nameentropicCell, dOCV_tot_n, dOCV_tot_x, dOCV_tot_y); 	// the entropic coefficient of the entire cell, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
	}
	catch(int e){
		cout<<"ERROR in the constructor of Cell_Fit when loading the OCV curves from the CSV files: "<<e<<". Throwing it on."<<endl<<flush;
		throw e;
	}


	// spatial discretisation of the solid diffusion PDE
	M = MM;
	// check if the inputs to the Matlab code are the same as the ones here in the C++ code
	// input:
	// 		M.Input[0] has to be the same as nch (defined in State.hpp)
	// 		M.Input[1] has to be the same as Rp (defined earlier in this constructor)
	// 		M.Input[2] has to be the same as Rn (defined earlier in this constructor)
	// 		M.input[3] has to give the location of the 0 eigenvalue
	bool Mnch = (M.Input[0]-nch)/M.Input[0] > pow(10,-10); 		// allow a relative difference of e-10 due to numerical errors
	if (Mnch)
		cerr<<"ERROR in the constructor of Cell_Fit: the value of nch used in the Matlab script "<<M.Input[0] <<" is not the same as the value of nch used in the c++ code "<<nch <<endl<<flush;
	bool Mrp = (M.Input[1]-Rp)/M.Input[1] > pow(10,-10);		// allow a relative difference of e-10 due to numerical errors
	if (Mrp)
		cerr<<"ERROR in the constructor of Cell_Fit: the value of Rp used in the Matlab script "<<M.Input[1] <<" is not the same as the value of Rp used in the c++ code "<<Rp <<endl<<flush;
	bool Mrn = (M.Input[2]-Rn)/M.Input[2] > pow(10,-10);		// allow a relative difference of e-10 due to numerical errors
	if (Mrn)
		cerr<<"ERROR in the constructor of Cell_Fit: the value of Rn used in the Matlab script "<<M.Input[2] <<" is not the same as the value of Rn used in the c++ code "<<Rn <<endl<<flush;
	int a = M.Input[3];
	bool Meig = abs(M.An[a])> pow(10,-10) || abs(M.Ap[a])> pow(10,-10); // allow a relative difference of e-10 due to numerical errors
	if (Meig)
		cerr<<"ERROR in the constructor of Cell_Fit: the row of the 0-eigenvalue is "<<M.Input[3] <<" but that row has a positive eigenvalue of "<<M.Ap[a]<<" and negative eigenvalue of "<<M.An[a]<<". They are not 0." <<endl;
	if (Mnch || Mrp || Mrn  || Meig){
		cout<<"The Matlab script modelSetup.m produces matrices used by the C++ code for the discretisation of the solid diffusion equation."
				" Matlab needs some input parameters, such as the number of nodes and the radius of the particles. These parameters are specified on top of the Matlab scripts."
				" These values are also defined in the C++ code. Of course, both values have to be the same."
				" It turned out this was not the case, so you either have to change the values in the Matlab script or the ones in the C++ code."
				" We are throwing an error."<<endl<<flush;
		throw 110;
	}

	// Initialise state variables
	double up[nch], un[nch];
	double fp, fn, T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, R, delta_pl;
	double Rdc = 0.0102;			// DC resistance of the total cell in Ohm, characterisation fitting parameter
	fp 		= 0.689332;				// lithium fraction in the cathode at 50% soc [-], OCV fitting parameter
	fn 		= 0.479283;				// lithium fraction in the anode at 50% soc [-], OCV fitting parameter
	T		= 273+25;				// cell temperature
	delta 	= pow(10,-9);			// SEI thickness. Start with a fresh cell, which has undergone some formation cycles so it has an initial SEI layer.
										// never start with a value of 0, because some equations have a term 1/delta, which would give nan or inf
										// so this will give errors in the code
	LLI		= 0;					// lost lithium. Start with 0 so we can keep track of how much li we lose while cycling the cell
	thickp	= 70*pow(10,-6);		// thickness of the positive electrode, OCV fitting parameter
	thickn 	= 73.5*pow(10,-6);		// thickness of the negative electrode, OCV fitting parameter
	ep 		= 0.5;					// volume fraction of active material in the cathode, OCV fitting parameter
	en 		= 0.5;					// volume fraction of active material in the anode, OCV fitting parameter
	ap		= 3*ep/Rp;				// effective surface area of the cathode, the 'real' surface area is the product of the effective surface area (a) with the electrode volume (elec_surf * thick)
	an 		= 3*en/Rn;				// effective surface area of the anode
	CS		= 0.01*an*elec_surf*thickn;	// initial crack surface. Start with 1% of the real surface area
	Dp 		= 8*pow(10,-14);		// diffusion constant of the cathode at reference temperature, characterisation fitting parameter
	Dn 		= 7*pow(10,-14);		// diffusion constant of the anode at reference temperature, characterisation fitting parameter
	R		= Rdc * ( (thickp*ap*elec_surf + thickn*an*elec_surf)/2 );// specific resistance of the combined electrodes, see State::iniStates
	delta_pl = 0;					// thickness of the plated lithium layer. You can start with 0 here
	s.initialise(nch, up, un, T, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, R, delta_pl); // set the states, with a random value for the concentration
	setC(fp, fn);					// set the lithium concentration

	// SEI parameters
		nsei = 1;
		alphasei = 1;
		OCVsei = 0.4;
		rhosei = 100*pow(10,3);
		Rsei = 2037.4;
		Vmain = 13.0;
		Vsei = 64.39;
		c_elec0 = 4.541/1000;
		// fitting parameters of the models
		seiparam.sei1k = 0;
		seiparam.sei1k_T = 0;
		seiparam.sei2k = 0;
		seiparam.sei2k_T = 0;
		seiparam.sei2D = 0;
		seiparam.sei2D_T = 0 ;
		seiparam.sei3k = 0;
		seiparam.sei3k_T = 0;
		seiparam.sei3D = 0;
		seiparam.sei3D_T = 0;
		seiparam.sei_porosity = 0;

	// surface cracking
		// fitting parameters of the models
		csparam.CS1alpha = 0;
		csparam.CS2alpha = 0;
		csparam.CS3alpha = 0;
		csparam.CS4alpha = 0;
		csparam.CS4Amax = 0;
		csparam.CS5k = 0;
		csparam.CS5k_T = 0;
		csparam.CS_diffusion = 0;

	// LAM
		OCVnmc = 4.1;
		// fitting parameters
		lamparam.lam1p = 0;
		lamparam.lam1n = 0;
		lamparam.lam2ap = 0;
		lamparam.lam2bp = 0;
		lamparam.lam2an = 0;
		lamparam.lam2bn = 0;
		lamparam.lam2t = 0;
		lamparam.lam3k = 0;
		lamparam.lam3k_T = 0;
		lamparam.lam4p = 0;
		lamparam.lam4n = 0;

	// li-plating parameters
		npl = 1;
		alphapl = 1;
		OCVpl = 0;
		rhopl = 10000*pow(10,3);
		// fitting parameters
		plparam.pl1k = 0;
		plparam.pl1k_T = 0;

	// degradation identifiers: no degradation
		deg_id.SEI_id[0] = 0;			// no SEI growth
		deg_id.SEI_n = 1;				// there is 1 SEI model (namely '0')
		deg_id.SEI_porosity = 0;		// don't decrease the volume fraction
		deg_id.CS_id[0] = 0;			// no surface cracks
		deg_id.CS_n = 1;				// there is 1 model (that there are no cracks)
		deg_id.CS_diffusion = 0;		// don't decrease the diffusion coefficient
		deg_id.LAM_id[0] = 0;			// no LAM
		deg_id.LAM_n = 1;				// there is 1 LAM model
		deg_id.pl_id = 0;				// no lithium plating

	// Check if we will have to calculate the stress
	s_dai = false;						// according to Dai's stress model
	s_lares = false;					// according to Laresgoiti's stress model
}

Cell_Fit::~Cell_Fit() {

}

void Cell_Fit::setOCVcurve(string namepos, string nameneg, int np, int nn){
	/*
	 * Function to load the specified OCV curves
	 *
	 * IN
	 * namepos 	name of the CSV file in which the cathode OCV curve is defined
	 * nameneg 	name of the CSV file in which the anode OCV curve is defined
	 * np 		number of points in the cathode OCV curve
	 * nn 		number of points in the anode OCV curve
	 *
	 * THROWS
	 * 111 		the OCV curves are too long
	 */

	// Update the number of points in each curve
	OCV_pos_n = np;
	OCV_neg_n = nn;

	// check that the OCV curves are shorter than the allowed range
	bool leng = false;
	if(OCV_pos_n >= OCV_maxLength){
		cerr<<"ERROR in Cell_Fit::setOCVcurve. The cathode OCV curve has "<<OCV_pos_n<<" points while the maximum number is "<<OCV_maxLength<<". Throwing an error"<<endl<<flush;
		leng = true;
	}
	if(OCV_neg_n >= OCV_maxLength){
		cerr<<"ERROR in Cell_Fit::setOCVcurve. The anode OCV curve has "<<OCV_neg_n<<" points while the maximum number is "<<OCV_maxLength<<". Throwing an error"<<endl<<flush;
		leng = true;
	}
	if(leng)
		throw 111;
	// if you get this error and you do want to use OCV curves which are longer than 1000 points, increase the size of the OCV arrays in Cell.hpp and the value of OCV_maxLength


	// Store the OCV curves
	try{
		loadCSV_2col(nameneg,OCV_neg_n,OCV_neg_x, OCV_neg_y); 				// the OCV curve of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(namepos,OCV_pos_n,OCV_pos_x, OCV_pos_y);				// the OCV curve of the cathode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
	}
	catch(int e){
		cout<<"ERROR in Cell_Fit::setOCVcurve when loading the OCV curves from the CSV files: "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}
}

void Cell_Fit::setInitialConcentration(double cmaxp, double cmaxn, double lifracp, double lifracn){
	/*
	 * Function to specify the lithiun fraction from which we have to start
	 *
	 * IN
	 * cmaxp 	maximum li-concentration in the cathode [mol m-3]
	 * cmaxn 	maximum li-concentration in the anode [mol m-3]
	 * lifracp 	lithium fraction in the cathode at 50% SoC
	 * lifracn 	lithium fraction in the anode at 50% SoC
	 */

	// Store the maximum concentrations
	Cmaxpos = cmaxp;
	Cmaxneg = cmaxn;

	// set the concentration
	try{
		setC(lifracp, lifracn);
	}
	catch(int e){
		cout<<"Error in Cell_Fit::setInitialConcentration when setting the concentration: "<<e<<". Throwing it on"<<endl;
		throw e;
	}
}

void Cell_Fit::setGeometricParameters(double capnom, double surf, double ep, double en, double thickp, double thickn){
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

	nomCapacity = capnom;
	elec_surf = surf;

	// If you change the volume fraction, you also have to change the effective surface area
	double ap = 3*ep/Rp;
	double an = 3*en/Rn;

	// overwrite the original geometric parameters
	s.overwriteGeometricStates(thickp, thickn, ep, en, ap, an);
}

void Cell_Fit::setCharacterisationParam(double Dp, double Dn, double kpi, double kni, double Rdc){
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

	// store the rate constants in the cell
	kp = kpi;
	kn = kni;

	// calculate the specific electrode resistance from the total DC resistance
	double r = Rdc * ( (s.getThickp()*s.getAp()*elec_surf + s.getThickn()*s.getAn()*elec_surf)/2 );

	// overwrite the cycling related parameters
	s.overwriteCharacterisationStates(Dp, Dn, r);
}

void Cell_Fit::setRamping(double Istep, double tstep){
	/*
	 * Function to set the ramping parameters of the cell.
	 * This allows to decrease them if you are trying weird things (e.g. very small diffusion constants) with the cell
	 *
	 * IN
	 * Istep 	maximum change in current per tstep [A], > 0
	 * tstep 	time step by which you change the current, [s], > 0
	 */

	dIcell = Istep;
	dt_I = tstep;
}

