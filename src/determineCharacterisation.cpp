/*
 * determineRateParameters.cpp
 *
 * The functions below can be used to fit the 'rate parameters' or 'characterisation parameters' of a cell at one temperature (currently 25 degrees):
 * 		diffusion constants (Dp and Dn)
 * 		rate constants of the main li-insertion reaction (kp and kn)
 * 		DC resistance of the cell (r)
 *
 * The user has to supply voltage measurements for a number of CC CV cycles (i.e. charge and discharge, both with a CC and CV phase).
 * These have to come in csv files with two columns
 * 		The first column has to have the charge throughput in Ah, starting at 0 and increasing as the cell is charged or discharged.
 * 		The second column has to have the cell voltage at the corresponding point.
 * In the top-level function ('fitCharacterisationAtReferenceT), the user has to define which cycles this are exactly (i.e. what C rates for the CC phase and which current threshold for the CV)
 *
 * The script then uses a hierarchical search algorithm to find the parameters which minimise the error between the simulated and measured voltage curves.
 * The results are written in a csv file. If the user wants to use these parameters, they have to be copied to the constructor of the respective cell-class.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Cell_Fit.hpp"
#include "Cycler.hpp"
#include "determineCharacterisation.h"
#include "ReadCSVfiles.h"
#include "Interpolation.h"
#include "determineOCV.h"

void CCCV(double Crate, double Ccut, double Tref, double Dp, double Dn, double kp, double kn, double R, struct OCVparam ocvfit,const struct Model& M,
		int nin, double Vsim[][2], double Tsim[][2], int* nsim){
	/*
	 * Function which simulates a full CC CV (dis)charge with the given rate parameters at reference temperature.
	 * The reference temperature is defined in Cell (Tref), and is 25 degrees.
	 * If you have data at a different reference temperature, you have to change the value of Tref in the constructor of Cell_Fit.
	 *
	 * IN
	 * Crate 	the C rate to be used for the CC phase
	 * 				> 0 for discharge
	 * 				< 0 for charge
	 * Ccut 	Crate of the cutoff current below which the CV phase should be stopped, > 0, [-]
	 * Tref 	temperature at which the characterisation is done [K]
	 * Dp		diffusion constant of the positive electrode at reference temperature
	 * Dn		diffusion constant of the negative electrode at reference temperature
	 * kp		rate constant of the lithium insertion at the positive electrode at reference temperature
	 * kn		rate constant of the lithium insertion at the negative electrode at reference temperature
	 * R		DC resistance of the cell
	 * ocvfit 	structure with the parameters determined by the functions in determineOCV.cpp (the struct is defined in determineCharacterisation.h)
	 * M 		structure with the matrices of the spatial discretisation of the solid diffusion PDE
	 * 				defined in Cell.hpp
	 * nin 		length of the arrays provided for feedback
	 *
	 * OUT
	 * Vsim		Matrix with the cell voltage [V] in the second column and the charge throughput [Ah] in the first one
	 * Tsim		Matrix with the cell temperature [K] in the second column and the charge throughput [Ah] in the first one
	 * nsim 	number of data points in the output arrays (length(Vsim))
	 *
	 * THROWS
	 * 10002 	an error occurred while cycling the cell, and reducing the time step didn't help
	 * 10003	one of the parameters is negative, which is not allowed;
	 */

	// *********************************************************** 1 variables ***********************************************************************
	int verbose = 0;																// integer deciding how verbose the simulation should be
																					// The higher the number, the more output there is.
																					// Recommended value is 1, only use higher value for debugging
																					// From 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
																					// 	0 	almost no messages are printed, only in case of critical errors related to illegal parameters
																					// 	1 	error messages are printed in case of critical errors which might crash the simulation
																					// 	2 	all error messages are printed, whether the simulation can recover from the errors or not
																					// 	3 	on top of the output from 2, a message is printed every time a function in the Cycler and BasicCycler is started and terminated
																					// 	4 	on top of the output from 3, the high-level flow of the program in the Cycler is printed (e.g. 'we are going to discharge the cell')
																					// 	5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
																					// 	6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
																					// 	7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated

	// Make a cell of the sub-class Cell_Fit. This sub-class defines some extra functions to change its parameters
	Cell_Fit c1 = Cell_Fit(M, verbose);

	// Check all parameters are positive
	if (Dp <= 0 || Dn <= 0 || kp <= 0 || kn <= 0 || R < 0)
		throw 10003;

	// Set the characterisation parameters of the cell to the ones given as input to this function
	c1.setOCVcurve(ocvfit.namepos, ocvfit.nameneg, ocvfit.np, ocvfit.nn);
	c1.setInitialConcentration(ocvfit.cmaxp, ocvfit.cmaxn, ocvfit.lifracpini, ocvfit.lifracnini);
	c1.setGeometricParameters(ocvfit.cap, ocvfit.elec_surf, ocvfit.ep, ocvfit.en, ocvfit.thickp, ocvfit.thickn);
	c1.setCharacterisationParam(Dp, Dn, kp, kn, R);
	c1.setVlimits(ocvfit.Vmax, ocvfit.Vmin);
	c1.setT(Tref);																	// set the temperature of the cell to the given value
	c1.setTenv(Tref);																// set the environmental temperature to the given value

	// time steps
	double dt = 2.0;																// time step for cycling [s]
	double Istep = 0.1;																// current step for ramping, indicating how fast the current can change per 'ramp time step', [A s-1]
	double tstep = 0.001;															// time step for ramping [s]
																					// the current can change at Istep/tstep, so currently 0.1A per 1 ms.

	// variables
	State s;																		// initial state of the cell, used to recover after an error
	double Iini;																	// initial current of the cell, used to recover after an error
	double Ccut2 = 0.05;															// Crate for the cutoff current when bringing the cell to the initial state (i.e. charge the cell first before you simulate the CCCV discharge)
	double Vset;																	// voltage at which the CCCV cycle should end (i.e. the minimum voltage if you are simulating a discharge)
	bool blockDegradation = true;													// don't account for degradation while doing the cycles
	string ID = "CharacterisationFit";												// identification string for the Cycler
	int timeCycleData = -1; 														// time interval at which cycling data has to be recorded [s]
																						// <0 means no folder is created and no data is stored
																						// 	0 means no data is recorded but a folder is still created for later use
																						//  >0 means data is recorded approximately every so many seconds
	Cycler cycler(c1, ID, verbose, timeCycleData);									// Make the Cycler
	double ahi, whi, timei;															// feedback variables we don't need


	// *********************************************************** 2 (dis)charge ***********************************************************************

	// We are going to start trying to simulate the (dis)charge with a normal time step (for time integration and ramping)
	// 	But smaller diffusion constants lead to decreased numerical stability and other discretiation errors
	// 	Therefore, these normal time steps might leads to errors in the code
	// So instead, there are loops which iteratively decrease the time steps to check if that solves the error.

	bool finished = false;															// boolean to indicate whether we have completed the discharge
	c1.getStates(s, &Iini);															// store the initial state

	// loop to decrease the integration time step (dt)
	for(int i=0;i<2 && !finished;i++){

		// Set the ramping time steps to their normal value
		Istep = 0.1;																// current step for ramping
		tstep = 0.001;																// time step for ramping

		// loop to decrease the time steps for ramping
		for(int j=0;j< 2 && !finished;j++){

			// Try to simulate the (dis)charge
			try{

				// restore the original battery state in case an error occurred earlier in the loop
				c1.setStates(s, Iini);
				c1.setRamping(Istep, tstep);
				cycler.setCyclingDataTimeResolution(0);								// don't collect cycling data during the charging

				// Bring the cell to the correct soc before simulating the (dis)charge
				if(Crate > 0){														// simulate a discharge
					// first fully charge the cell to the maximum voltage at 1C
					cycler.CC_V_CV_I(1, ocvfit.Vmax, Ccut2, dt, blockDegradation, &ahi, &whi, &timei);
					Vset = ocvfit.Vmin;
				}
				else{																// simulate a charge
					// first fully discharge the cell to the minimum voltage at 1C
					cycler.CC_V_CV_I(1, ocvfit.Vmin, Ccut2, dt, blockDegradation, &ahi, &whi, &timei);
					Vset = ocvfit.Vmax;
				}

				// simulate the CC CV (dis)charge
				cycler.setCyclingDataTimeResolution(dt);							// collect cycling data of every time step
				cycler.CC_V_CV_I(abs(Crate), Vset, Ccut, dt, blockDegradation, &ahi, &whi, &timei);

				// If we get here, no errors were thrown and we can leave the loop
				finished = true;													// indicate we have finished the simulation
				break;																// leave the loops
			}
			catch(int e){															// An error occurred while simulating the (dis)charge
				Istep = Istep / 10.0;												// reduce the ramping parameters by a factor of 10
				// now we are still in the loop which decreases the ramping time steps ('final' is still false), so you will try again
				// and the original battery state will be restored at the start of the loop, so the illegal battery state will be 'forgotten'
			}
		} // end loop to decrease the ramping time step

		// if we haven't finished the cycle, try again with a smaller time step
		if(!finished){
			dt = dt / 10.0;
			// now we are still in the loop which decreases the time step, so you will try again
			// and the original battery state will be restored at the start of the loop, so the illegal battery state will be 'forgotten'
		}
	}// end loop to decrease the integration time step


	// *********************************************************** 3 output ***********************************************************************

	// Throw an error if the simulation wasn't successful. This happens if 'finished' is still false or if ahi == 0 i.e. no charge could be discharged
	if(!finished || ahi == 0)
		throw 10002;

	// Get the cell voltage from the simulated (dis)charge from the Cycler
	double Ahi[nin], Vi[nin], Ti[nin];
	int nout2;
	try{
		cycler.returnCyclingData(nin, Ahi, Vi, Ti, &nout2);
	}
	catch(int e){
		cerr<<"Error in determineCharacterisation::CCCV when getting the cycling data "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}
	for(int i=0;i<nout2;i++){
		Vsim[i][0] = Ahi[i]; // always positive
		Vsim[i][1] = Vi[i];
		Tsim[i][0] = Ahi[i];
		Tsim[i][1] = Ti[i];
	}
	*nsim = nout2;
}

void fitDiffusionAndRate(int hierarchy, int ir, double R, int nDstep, double Dpstep, double Dnstep, bool logDstep, double Dpmin, double Dnmin,
		int nkstep, double kpstep, double knstep, bool logkstep, double kpmin, double knmin,
		int nCCCV, double weights[], string names[], int lengths[],
		double Crates[], double Ccuts[], double Tref, struct OCVparam ocvfit,
		double* err, double par[], double parindex[]){
	/*
	 * Function which goes trough the specified search space for Dp, Dn, kp and kn, for a constant value of the DC resistance R.
	 * The steps can be taken logarithmically (e.g. e-10, e-11, etc.) or linearly (e.g. 1e-10, 2e-10, etc.)
	 * For each combination of the 4 parameters, the voltage of the CCCV cycles is simulated and the error with the measured voltages is calculated.
	 * The best fit (i.e. 4 parameters with the lowest error on the OCV curve) is returned
	 *
	 * IN
	 * hierarchy level of the hierarchy in which we are
	 * ir 		index of where in the search space for r we are
	 * R 		DC resistance of the cell to be used in the search
	 * nDstep 	number of steps in the search space for diffusion constants
	 * Dpstep 	magnitude of one step in the search space for the cathode diffusion constant
	 * Dnstep 	magnitude of one step in the search space for the anode diffusion constant
	 * logDstep if true the steps in the search space for the diffusion constants are logarithmic, else they are linear
	 * Dpmin 	lowest value in the search space of Dp
	 * Dnmin 	lowest value n the search space of Dn
	 * nkstep 	number of steps in the search space for rate constants
	 * kpstep 	magnitude of one step in the search space for the cathode rate constant
	 * knstep 	magnitude of one step in the search space for the anode rate constant
	 * logkstep if true the steps in the search space for the rate constants are logarithmic, else they are linear
	 * kpmin 	lowest value in the search space for kp
	 * knmin 	lowest value in the search space for kn
	 * nCCCV 	number of CCCV cycles
	 * weights 	array with the weight that should be given to each CCCV curve when calculating the overall error
	 * 			the sum of all weights should be one
	 * names 	array with the names of the csv files with the measurements for the cell
	 * lengths 	length of the csv files
	 * Crates 	array with the C rates of the CC phases of each experiment, >0 for discharge, <0 for charge
	 * Ccuts 	array with the C rates of the current threshold for the CV phase of each experiment, >0 . (set to a very large value if you don't want a CV phase)
	 * Tref 	temperature at which the characterisation is done [K]
	 * ocvfit 	structure with the values of the OCV parameters determined by determineOCV::estimateOCVparam
	 *
	 * OUT
	 * err 		error of the best fit
	 * par 		values of r, Dp, Dn, kp and kn which achieved the best fit
	 * parindex indices in the search space for r, Dp, Dn, kp and kn which achieved the best fit
	 *
	 * THROWS
	 * 10000 	an input file is in the wrong format
	 * 10001 	the arrays for the data are too short, you have to increase the value of 'n'
	 */

	// *********************************************************** 1 variables ***********************************************************************

	// Read the measured voltage profiles
	int n = 10000;													// max length of the data
	for(int i=0;i<nCCCV;i++)										// ensure the arrays are long enough to store the data from all profiles
		n = max(n, lengths[i]);
	bool val;														// boolean indicating if the data is in the correct format
	bool checkRange = false;										// the first column has to start at 0, but not end at 1
	double Vdata[n][2];												// array with the data of one file
	double Vall[n][nCCCV], Ahall[n][nCCCV];							// matrix with the voltage of each experiment, and the charge throughput of each experiment
	for(int i=0;i<nCCCV;i++){										// loop to read the voltage profile of each cycle
		if(lengths[i] > n){
			cerr<<"ERROR in determineCharacterisation::fitDiffusionAndRate. The arrays to store the data are too short. They have a length of "<<n<<" but need a length of at least "<<lengths[i]<<". Throw an error."<<endl<<flush;
			throw 10001;
			// if you get this error, increase the value of 'n' defined in this function ('max length of the data')
		}
		loadCSV_2colMatrix(names[i], lengths[i], Vdata);			// read the csv file with the voltage profile

		// check the data is in the correct format
		val = validOCV(checkRange, lengths[i], Vdata);
		if(!val){
			cerr<<"ERROR in determineCharacterisation::fitDiffusionAndRate. Input file "<<names[i]<<" has the wrong format. throwing an error"<<endl<<flush;
			throw 10000;
		}

		// reorganise the data in separate matrices for charge and voltage
		for(int j=0;j<lengths[i];j++){
			Ahall[j][i] = Vdata[j][0];
			Vall[j][i] = Vdata[j][1];
		}
	}

	// Variables
	Model M;														// structure with the matrices for the spatial discretisation for the solid diffusion PDE
	Model_initialise(M);
	double dt = 2;													// time step to be used for the simulation
	int nin = 1/0.5 * 3600/dt * 10; 								// length of the arrays to store the simulation results. (for a 0.5 C rate CC phase, *10 for the CV phase and safety margin)
	double Vsim[nin][2], Tsim[nin][2]; 								// arrays to store the simulation results
	int nsim;														// number of data points in the simulated voltage curve
	double erri;													// error of this CCCV cycle with the given parameters
	double errcomb;													// error of all CCCV cycles with the given parameters
	double errmin = 10000000000;									// lowest error encountered so far
	double Dp, Dn, kp, kn;											// variables for the search space

	// *********************************************************** 2 loop through the search space ***********************************************************************

	// scan the search space for Dp
	for(int idp=0;idp<nDstep;idp++){
		if (logDstep)												// take steps logarithmically
			Dp = Dpmin * pow(Dpstep,(idp));
		else														// take steps linearly
			Dp = Dpmin + idp*Dpstep;

		// scan the search space for Dn
		for(int idn = 0;idn<nDstep;idn++){
			if (logDstep)											// take steps logarithmically
				Dn = Dnmin * pow(Dnstep,(idn));
			else													// take steps linearly
				Dn = Dnmin + idn*Dnstep;

			// scan the search space for kp
			for(int ikp=0;ikp<nkstep;ikp++){
				if (logkstep)										// take steps logarithmically
					kp = kpmin * pow(kpstep,(ikp));
				else												// take steps linearly
					kp = kpmin + ikp*kpstep;

				// scan the search space for kn
				for(int ikn=0;ikn<nkstep;ikn++){
					if (logkstep)									// take steps logarithmically
						kn = knmin * pow(knstep,(ikn));
					else											// take steps linearly
						kn = knmin + ikn*knstep;

					// Calculate the error for this set of parameters
					errcomb = 0;									// initialise the combined error of all CCCV experiments for this combination of Dp, Dn, kp and kn to 0
					for(int i=0;i<nCCCV;i++){						// loop through all CCCV cycles
						try{
							// Simulate this CCCV experiment
							CCCV(Crates[i], Ccuts[i], Tref, Dp, Dn, kp, kn, R, ocvfit, M, nin, Vsim, Tsim, &nsim);

							// Get the measured data from this CCCV experiment
							for(int j=0;j<lengths[i];j++){
								Vdata[j][0] = Ahall[j][i];
								Vdata[j][1] = Vall[j][i];
							}

							// calculate the error (function defined in determineOCV.cpp)
							erri = calculateError(false, lengths[i], nsim, Vdata, Vsim);
						}
						catch(int e){
							// An error occured while simulating. This means the parameters were infeasible.
							// Therefore, indicate this as a very large error.
							erri = 10000000000;
						}

						// calculate the total (weighted) error
						errcomb += abs(erri)*weights[i];
					} // loop for CCCV experiments

					// Store the minimum error
					if(errcomb < errmin){							// check if the error of this combination is better than the best fit so far
						par[0] = R;
						par[1] = Dp;
						par[2] = Dn;
						par[3] = kp;
						par[4] = kn;
						parindex[0] = ir;
						parindex[1] = idp;
						parindex[2] = idn;
						parindex[3] = ikp;
						parindex[4] = ikn;
						errmin = errcomb;
					}
				} // loop for kn
			} // loop for kp
		} // loop for Dn
	} // loop for Dp
	*err = errmin;													// return the lowest error
}

void oneLevelCharacterisationFit(int hierarchy, int nrstep, double rstep, double rmin,
		int nDstep, double Dpstep, double Dnstep, bool logDstep, double Dpmin, double Dnmin,
		int nkstep, double kpstep, double knstep, bool logkstep, double kpmin, double knmin,
		int nCCCV, double weights[], string names[], int lengths[],
		double Crates[], double Ccuts[], double Tref, struct OCVparam ocvfit,
		double* err, double par[], int parindex[]){
	/*
	 * Function to find the optimal points in one level of the hierarchical search for r, Dp, Dn, kp and kn
	 * I.e. it goes through the search space of 5 variables
	 * 	DC resistance
	 * 	four variables from fitDiffusionAndRate (diffusion and rate constants)
	 * This function uses three threads to accelerate the computation
	 * The best fit (i.e. four parameters with the lowest error on the voltage curve) is returned and a csv file with the optimal value of the parameters is written.
	 *
	 * IN
	 * hierarchy level of the hierarchy in which we are
	 * nrstep	number of steps in the search space for Rdc
	 * rstep 	magnitude of one step in the search space for Rdc (always linear)
	 * rmin		lowest value for Rdc
	 * nDstep 	number of steps in the search space for diffusion constants
	 * Dpstep 	magnitude of one step in the search space for the cathode diffusion constant
	 * Dnstep 	magnitude of one step in the search space for the anode diffusion constant
	 * logDstep if true the steps in the search space for the diffusion constants are logarithmic, else they are linear
	 * Dpmin 	lowest value in the search space of Dp
	 * Dnmin 	lowest value n the search space of Dn
	 * nkstep 	number of steps in the search space for rate constants
	 * kpstep 	magnitude of one step in the search space for the cathode rate constant
	 * knstep 	magnitude of one step in the search space for the anode rate constant
	 * logkstep if true the steps in the search space for the rate constants are logarithmic, else they are linear
	 * kpmin 	lowest value in the search space for kp
	 * knmin 	lowest value in the search space for kn
	 * nCCCV 	number of CCCV cycles
	 * weights 	array with the weight that should be given to each CCCV curve when calculating the overall error
	 * 			the sum of all weights should be one
	 * names 	array with the names of the csv files with the measurements for the cell
	 * lengths 	length of the csv files
	 * Crates 	array with the C rates of the CC phases of each experiment, >0 for discharge, <0 for charge
	 * Ccuts 	array with the Crate of the current threshold for the CV phase of each experiment, >0 . (set to a very large value if you don't want a CV phase)
	 * Tref 	temperature at which the characterisation is done [K]
	 * ocvfit 	structure with the values of the OCV parameters determined by determineOCV::estimateOCVparam
	 *
	 * OUT
	 * err 		error of the best fit
	 * par 		values of r, Dp, Dn, kp and kn which achieved the best fit
	 * parindex indices in the search space for R, Dp, Dn, kp, kn
	 */

	// variables
	double r1, r2, r3;												// resistance for each thread [Ohm]
	double err1, err2, err3;										// lowest error for that value of resistance
	double par1[5], par2[5], par3[5];								// parameters [R Dp Dn kp kn] giving the lowest error for that resistance
	double parindex1[5], parindex2[5], parindex3[5];				// indices for parameters [R Dp Dn kp kn] giving the lowest error for that resistance
	double errmin = 10000000000;									// lowest error encountered so far

	// loop through the search space for R, take 3 steps per iteration, each using one thread
	for(int ir=0;ir<nrstep/3+1;ir++){

		// calculate the resistance for the 3 steps
		r1 = rmin + (ir*3)*rstep;
		r2 = rmin + (ir*3+1)*rstep;
		r3 = rmin + (ir*3+2)*rstep;

		// For each resistance, search the best fit of the remaining 4 parameters (multi-threaded)
		thread t1(fitDiffusionAndRate, hierarchy,(ir*3), r1, nDstep, Dpstep, Dnstep, logDstep, Dpmin, Dnmin, nkstep, kpstep, knstep, logkstep, kpmin, knmin,
				nCCCV, weights, names, lengths, Crates, Ccuts, Tref, ocvfit, &err1, par1, parindex1);
		thread t2(fitDiffusionAndRate, hierarchy,(ir*3+1), r2, nDstep, Dpstep, Dnstep, logDstep, Dpmin, Dnmin, nkstep, kpstep, knstep, logkstep, kpmin, knmin,
				nCCCV, weights, names, lengths, Crates, Ccuts, Tref, ocvfit, &err2, par2, parindex2);
		thread t3(fitDiffusionAndRate, hierarchy,(ir*3+2), r3, nDstep, Dpstep, Dnstep, logDstep, Dpmin, Dnmin, nkstep, kpstep, knstep, logkstep, kpmin, knmin,
				nCCCV, weights, names, lengths, Crates, Ccuts, Tref, ocvfit, &err3, par3, parindex3);

		// join the threads again
		t1.join();
		t2.join();
		t3.join();

		// store the lowest error & the parameters leading to this error
		if(err1 < errmin){ 											// Check if the error from the first step is lower
			errmin = err1;
			for(int i=0;i<5;i++){
				par[i] = par1[i];
				parindex[i] = parindex1[i];
			}
		}
		if(err2 < errmin){ 											// Check if the error from the second step is lower
			errmin = err2;
			for(int i=0;i<5;i++){
				par[i] = par2[i];
				parindex[i] = parindex2[i];
			}
		}
		if(err3 < errmin){ 											// Check if the error from the third step is lower
			errmin = err3;
			for(int i=0;i<5;i++){
				par[i] = par3[i];
				parindex[i] = parindex3[i];
			}
		}
	}

	// Return the lowest error
	*err = errmin;

	// Print the best fit, and write in a CSV file
	double Rdc = par[0];
	double Dp = par[1];
	double Dn = par[2];
	double kp = par[3];
	double kn = par[4];
	cout<<"The best fit in hierarchy "<<hierarchy<<" is: Rdc = "<<Rdc<<", Dp = "<<Dp<<", Dn = "<<Dn<<", kp = "<<kp<<", kn = "<<kn<<endl<<flush;
	ofstream output;
	// write the parameters and the magnitude of the error
	string na = "characterisationFit_" + to_string(hierarchy) + "_param.csv";
	output.open(na);
	output<<"Rdc"<<","<<Rdc<<"\n";
	output<<"Dp"<<","<<Dp<<"\n";
	output<<"Dn"<<","<<Dn<<"\n";
	output<<"kp"<<","<<kp<<"\n";
	output<<"kn"<<","<<kn<<"\n";
	output<<"total RMSE"<<","<<errmin<<"\n";
	output.close();
}

void hierarchicalCharacterisationFit(int hmax, int nrstep, double Rstep, double Rmin,
		int nDstep, double Dpstep, double Dnstep, bool logDstep, double Dpmin, double Dnmin,
		int nkstep, double kpstep, double knstep, bool logkstep, double kpmin, double knmin,
		int nCCCV, double weights[], string names[], int lengths[],
		double Crates[], double Ccuts[], double Tref, struct OCVparam ocvfit,
		double* err, double par[]){
	/*
	 * Hierarchical search algorithm to converge on the best fit. For a convex problem, the optimal point is found.
	 * For a nonconvex problem, the algorithm might find a local minimum only.
	 *
	 * It iteratively 'zooms in' on the optimal point in the search space.
	 * I.e. you first call it with a large range and a large step size for each parameter
	 * 		it then finds the optimal combination (with this large step)
	 * 		the search space is then refined to a region around the optimal point, with a smaller range and a smaller step size
	 * 		the second iteration then finds the optimal point in this new (smaller) region (with the smaller step)
	 * 		the search space is again refined to the region around this 2nd optimal point, with an even smaller range, and even smaller steps
	 * 		etc.
	 * The region is always refined to the points before and after the optimal point for each parameter.
	 * E.g. if the optimal value is P and we were scanning with a step size of dp, then the new range in the next level is P-dp to P+dp
	 * And the number of steps remains the same, so the new step size is ( (P+dp) - (P-dp) ) / (number of steps - 1) = 2dp/(step-1)
	 *
	 * So if the function is called with a step number nstep and an initial step size of dp,
	 * then in hierarchical level n, the step size is dp* (2/(nstep-1)^(n-1)
	 * So after 'hmax' level, the accuracy is dp* (2/(nstep-1)^(hmax-1)
	 *
	 * IN
	 * hmax  	number of hierarchical steps we should take.
	 * 			the higher the value, the higher the accuracy of the fit will be but the longer the calculation will take
	 * nrstep	number of steps in the search space for Rdc
	 * Rstep 	magnitude of one step in the search space for Rdc (always linear)
	 * Rmin		lowest value for Rdc
	 * ir 		index of where in the search space for r we are
	 * R 		DC resistance of the cell to be used in the search
	 * nDstep 	number of steps in the search space for diffusion constants
	 * Dpstep 	magnitude of one step in the search space for the cathode diffusion constant
	 * Dnstep 	magnitude of one step in the search space for the anode diffusion constant
	 * logDstep if true the steps in the search space for the diffusion constants in the first level are logarithmic (e.g. e-10, e-11, etc.)
	 * 			else they are linear  (e.g. 1e-10, 2e-10, etc.)
	 * 			From the second level onwards, the steps are always linear.
	 * Dpmin 	lowest value in the search space of Dp
	 * Dnmin 	lowest value n the search space of Dn
	 * nkstep 	number of steps in the search space for rate constants
	 * kpstep 	magnitude of one step in the search space for the cathode rate constant
	 * knstep 	magnitude of one step in the search space for the anode rate constant
	 * logkstep if true the steps in the search space for the rate constants in the first level are logarithmic (e.g. e-10, e-11, etc.)
	 * 			else they are linear  (e.g. 1e-10, 2e-10, etc.)
	 * 			From the second level onwards, the steps are always linear.
	 * kpmin 	lowest value in the search space for kp
	 * knmin 	lowest value in the search space for kn
	 * nCCCV 	number of CCCV cycles
	 * weights 	array with the weight that should be given to each CCCV curve when calculating the overall error
	 * 			the sum of all weights should be one
	 * names 	array with the names of the csv files with the measurements for the cell
	 * lengths 	length of the csv files
	 * Crates 	array with the C rates of the CC phases of each experiment, >0 for discharge, <0 for charge
	 * Ccuts 	array with the Crate of the current threshold for the CV phase of each experiment, >0 . (set to a very large value if you don't want a CV phase)
	 * Tref 	temperature at which the characterisation is done [K]
	 * ocvfit 	structure with the values of the OCV parameters determined by determineOCV::estimateOCVparam
	 *
	 * OUT
	 * err 		error of the best fit
	 * par 		values of r, Dp, Dn, kp and kn which achieved the best fit
	 */

	// variables
	double erri; 								// error in this level of the hierarchy
	double pari[5];								// best fit parameters in this level of the hierarchy
	int parindex[5];							// indices of the best fit parameters in this level of the hierarchy

	// Calculate the maximum value for the parameters in the initial search space
	double Rmax = Rmin + nrstep*Rstep;			// maximum values of the DC resistance in the search space
	double Dpmax, Dnmax, kpmax, knmax;			// maximum value for the other parameters
	if (logDstep){
		Dpmax = Dpmin * pow(Dpstep,nDstep-1);
		Dnmax = Dnmin * pow(Dnstep,nDstep-1);
	}
	else{
		Dpmax = Dpmin + Dpstep * (nDstep-1);
		Dnmax = Dnmin + Dnstep * (nDstep-1);
	}
	if (logkstep){
		kpmax = kpmin * pow(kpstep,nkstep-1);
		knmax = knmin * pow(knstep,nkstep-1);
	}
	else{
		kpmax = kpmin + kpstep * (nkstep-1);
		knmax = knmin + knstep * (nkstep-1);
	}


	// Loop for each level in the search
	for (int h=0;h<hmax;h++){

		// print the search space of this level
		cout<<"Start hierarchy level "<<h<<" with the following search spaces: "<<endl;
		cout<<"R: from "<<Rmin<<" to "<<Rmax<<" in "<<nrstep<<" steps with magnitude "<<Rstep<<endl;
		cout<<"Dp: from "<<Dpmin<<" to "<<Dpmax<<" in "<<nDstep<<" steps with magnitude "<<Dpstep<<endl;
		cout<<"Dn: from "<<Dnmin<<" to "<<Dnmax<<" in "<<nDstep<<" steps with magnitude "<<Dnstep<<endl;
		cout<<"kp: from "<<kpmin<<" to "<<kpmax<<" in "<<nkstep<<" steps with magnitude "<<kpstep<<endl;
		cout<<"kn: from "<<knmin<<" to "<<knmax<<" in "<<nkstep<<" steps with magnitude "<<knstep<<endl;

		// Calculate the best fit in this level
		oneLevelCharacterisationFit(h, nrstep, Rstep, Rmin,
				nDstep, Dpstep, Dnstep, logDstep, Dpmin, Dnmin,
				nkstep, kpstep, knstep, logkstep, kpmin, knmin,
				nCCCV, weights, names, lengths, Crates, Ccuts, Tref, ocvfit, &erri, pari, parindex);

		// Update the search space
		// 	suppose the optimal value for a parameter p was in the search space at index i
		// 	then the new search space has as minimum value the value of p at i-1 and as maximum the value of p at i+1
		Rmin = Rmin + (parindex[0]-1)*Rstep;
		Rmax = Rmin + (parindex[0]+1)*Rstep;
		if (logDstep){
			Dpmax = Dpmin * pow(Dpstep,(parindex[1]+1));
			Dpmin = Dpmin * pow(Dpstep,(parindex[1]-1));
			Dnmax = Dnmin * pow(Dnstep,(parindex[2]+1));
			Dnmin = Dnmin * pow(Dnstep,(parindex[2]-1));
		}
		else{
			Dpmax = Dpmin + (parindex[1]+1)*Dpstep;
			Dpmin = Dpmin + (parindex[1]-1)*Dpstep;
			Dnmax = Dnmin + (parindex[2]+1)*Dnstep;
			Dnmin = Dnmin + (parindex[2]-1)*Dnstep;
		}
		if (logkstep){
			kpmax = kpmin * pow(kpstep,(parindex[3]+1));
			kpmin = kpmin * pow(kpstep,(parindex[3]-1));
			knmax = knmin * pow(knstep,(parindex[4]+1));
			knmin = knmin * pow(knstep,(parindex[4]-1));
		}
		else{
			kpmax = kpmin + (parindex[3]+1)*kpstep;
			kpmin = kpmin + (parindex[3]-1)*kpstep;
			knmax = knmin + (parindex[4]+1)*knstep;
			knmin = knmin + (parindex[4]-1)*knstep;
		}

		// only the first search level is logarithmic, afterwards the search has linear steps
		logDstep = false;
		logkstep = false;

		// 	the step size is chosen to produce the same number of steps as before
		Rstep = (Rmax - Rmin)/(nrstep-1);					// (nRstep-1) because for loop goes from 0 to nRstep-1
		Dpstep = (Dpmax - Dpmin)/(nrstep-1);
		Dnstep = (Dnmax - Dnmin)/(nrstep-1);
		kpstep = (kpmax - kpmin)/(nrstep-1);
		knstep = (knmax - knmin)/(nrstep-1);
	}

	// Make the output parameters
	*err = erri;
	for(int i=0;i<5;i++)
		par[i] = pari[i];
}


void estimateCharacterisation(){
	/*
	 * Function which will find the parameters for the diffusion constants, rate constants and DC resistance
	 * which best fit the measured voltage curve of the user.
	 *
	 * As input, the user has to supply csv files with voltage measurements of some full CCCV (dis)charges.
	 * E.g. starting from the maximum voltage, do a CC discharge at 1C followed by a CV discharge until the current is below some threshold.
	 * The first column has to have the charge throughput in Ah, starting at 0 and increasing as the cell is charged or discharged.
	 * The second column has to have the cell voltage at the corresponding point.
	 * The cycles have to be 'full', i.e. if it is a charge, than the first voltage must be Vmax and the last voltage must be Vmin (i.e. from fully charged to fully discharged) and vice versa.
	 *
	 * The user can specify in the code below which C rates were used for the CC phases, and what the cutoff current threshold is.
	 * Additionally, the user has to specify the OCV-related parameters for the cell, as calculated by determineOCV.
	 *
	 * The function finds optimal values for 5 parameters:
	 * 		diffusion constants (Dp and Dn)
	 * 		rate constants of the main li-insertion reaction (kp and kn)
	 * 		DC resistance of the cell (r)
	 *
	 * The values of the parameters which give the best fit are written in a CSV file.
	 * Also, the simulated voltage curves with these parameters are written to csv files so it can be compared with the measured ones.
	 */

	// *********************************************************** 1 USER INPUT ***********************************************************************

	// Specify the characterisation tests for which data is available
	double Tref = 273 + 25;												// Temperature at which the characterisation should be done [K]
	int nCCCV = 5;														// number of data sets for CCCV (dis)charges
	string names[nCCCV] = {"Characterisation_0.2C_CC_discharge.csv", "Characterisation_0.5C_CC_discharge.csv",
			"Characterisation_1C_CC_discharge.csv", "Characterisation_2C_CC_discharge.csv","Characterisation_3C_CC_discharge.csv"};// Name of the files with the voltage curve from the cell
	int lengths[nCCCV] = {320,129,66,37,27};							// the number of data points for each voltage curve
	double Crates[nCCCV] = {0.2, 0.5, 1, 2, 3};							// C rates of the CC phases for each voltage curve, <0 for charge, >0 for discharge
	double Ccuts[nCCCV] = {100, 100, 100, 100, 100};					// C rates of the current threshold for the CV phases [A], > 0
																		// 	a very high value (above the Crate of the CC) will avoid there is a CV phase at all, so there is only a CC phase
	double weights[nCCCV] = {0.2, 0.2, 0.2, 0.2, 0.2};					// array with the weight that should be attributed to each individual CCCV curve when calculating the overall error
																		// 	no restrictions are imposed, but it is recommended that all weights sum up to 1
																		// 	if a weight is set to 0, that curve is ignored. If a weight is negative, the error for that curve will be maximised so this is not recommended.

	// names of output files
	string nameparam = "characterisationFit_parameters.csv";			// name of the output csv file in which the optimal parameters will be written
	string nameCCCVfit = "characterisationFit_";						// prefix appended before the name of the output data files with the simulations for the best fit
																		// each file will have 3 columns: charge throughput, voltage and temperature

	// Specify the OCV parameters (calculated by determineOCV::estimateOCVparameters)
	OCVparam ocvfit;
		ocvfit.elec_surf = 0.0982;										// electrode surface
		ocvfit.ep = 0.5;												// volume fraction of active material in the cathode
		ocvfit.en = 0.5;												// volume fraction of active material in the anode
		ocvfit.thickp = 70*pow(10,-6);									// thickness of the cathode
		ocvfit.thickn = 73.5*pow(10,-6);								// thickness of the anode
		ocvfit.lifracpini = 0.6862;										// lithium fraction in the cathode at 50% soC
		ocvfit.lifracnini = 0.4843;										// lithium fraction in the anode at 50% SoC
		ocvfit.cmaxp = 51385;											// maximum lithium concentration in the cathode [mol m-3]
		ocvfit.cmaxn = 30555;											// maximum lithium concentration in the anode [mol m-3]
		ocvfit.cap = 2.7;												// the capacity of the cell [Ah]
		ocvfit.Vmax = 4.2;												// maximum voltage of the cell [V]
		ocvfit.Vmin = 2.7; 												// minimum voltage of the cell [V]

		ocvfit.namepos = "OCVfit_cathode.csv";							// name of the CSV file with the cathode OCV curve
		ocvfit.nameneg = "OCVfit_anode.csv";							// name of the CSV file with the anode OCV curve
		ocvfit.np = 49;													// number of points in the cathode OCV curve
		ocvfit.nn = 63;													// number of points in the anode OCV curve

	// ****************************************** 2 define the search space for fitting parameters ***********************************************************************

	// define the search space for the characterisation parameters at reference temperature
	//	 Dp		diffusion constant of the positive electrode at reference temperature
	//	 Dn		diffusion constant of the negative electrode at reference temperature
	//	 kp		rate constant of the lithium insertion at the positive electrode at reference temperature
	//	 kn		rate constant of the lithium insertion at the negative electrode at reference temperature
	//	 Rdc	DC resistance of the cell
	int nDstep = 10;													// number of steps in the search space for Dp and Dn
	double Dpmin = pow(10,-18);											// lowest value in the search space for Dp [m s-1]
	double Dnmin = pow(10,-18);											// lowest value in the search space for Dn [m s-1]
	double Dpstep = 13;													// step size in the search space for Dp
	double Dnstep = 13;													// step size in the search space for Dn
	bool logDstep = true;												// if true, the steps in the first level are logarithmically, i.e. D = Dmin * Dstep^i
																		// if false, the steps are linearly, i.e. D = Dmin + i*Dstep
																		// steps in the later search levels are always linearly
	int nkstep = 10;													// number of steps in the search space for kp and kn
	double kpmin = pow(10,-18);											// lowest value in the search space for kp
	double knmin = pow(10,-18);											// lowest value in the search space for kn
	double kpstep = 13;													// step size in the search space for kp
	double knstep = 13;													// step size in the search space for kn
	bool logkstep = true;												// if true, the steps in the first level are logarithmically, i.e. k = kmin * kstep^i
																		// if false, the steps are linearly, i.e. k = kmin + i*kstep
																		// steps in the later search levels are always linearly
	int nrstep = 6;														// number of steps in the search space for Rdc
	double rmin = 0.000001;												// minimum value of Rdc [Ohm]
	double rstep = 0.005;												// step size in the search space for Rdc [Ohm]

	// ***************************************************** 3 Fit the parameters ***********************************************************************

	// Call the hierarchical search algorithm, which does the fitting
	int hmax = 3;														// number of hierarchical levels to use. Increasing this number will improve the accuracy, but take longer to calculate
	double err;															// error in the best fit
	double par[5];														// parameters giving the lowest error [R Dp Dn kp kn]
	hierarchicalCharacterisationFit(hmax, nrstep, rstep, rmin,
			nDstep, Dpstep, Dnstep, logDstep, Dpmin, Dnmin,
			nkstep, kpstep, knstep, logkstep, kpmin, knmin,
			nCCCV, weights, names, lengths, Crates, Ccuts, Tref, ocvfit, &err, par);


	// ***************************************************** 4 write outputs ***********************************************************************

	// Print the best fit, and write in a CSV file
	double Rdc = par[0];
	double Dp = par[1];
	double Dn = par[2];
	double kp = par[3];
	double kn = par[4];
	cout<<"The best fit is: Rdc = "<<Rdc<<", Dp = "<<Dp<<", Dn = "<<Dn<<", kp = "<<kp<<", kn = "<<kn<<endl<<flush;

	// write the parameters in a csv file
	ofstream output;
	output.open(nameparam);
	output<<"Rdc"<<","<<Rdc<<"\n";
	output<<"Dp"<<","<<Dp<<"\n";
	output<<"Dn"<<","<<Dn<<"\n";
	output<<"kp"<<","<<kp<<"\n";
	output<<"kn"<<","<<kn<<"\n";
	output<<"temperature"<<","<<Tref<<"\n";
	output.close();

	// Simulate the voltage curves at the best fit
	Model M;															// structure with the matrices for the spatial discretisation for the solid diffusion PDE
	Model_initialise(M);
	double dt = 2;														// time step to be used for the simulation
	int nin = 1/0.5 * 3600/dt * 10; 									// length of the arrays to store the simulation results. (for a 0.5 C rate CC phase, *10 for the CV phase and safety margin)
	double Vsim[nin][2], Tsim[nin][2]; 									// arrays to store the simulation results
	int nsim;															// number of data points in the simulated voltage curve
	for(int i=0;i<nCCCV;i++){											// loop through all CCCV experiments
		CCCV(Crates[i], Ccuts[i], Tref, Dp, Dn, kp, kn, Rdc, ocvfit, M, nin, Vsim, Tsim, &nsim);

		// write the simulated voltages in a csv file
		output.open(nameCCCVfit + names[i]);
		for(int i=0;i<nsim;i++)
			output<<Vsim[i][0]<<","<<Vsim[i][1]<<","<<Tsim[i][1]<<"\n";
		output.close();
	}

	// Read the measured voltage profiles
	int n = 10000;														// max length of the data
	for(int i=0;i<nCCCV;i++)											// ensure the arrays are long enough to store the data from all profiles
		n = max(n, lengths[i]);
	bool val;															// boolean indicating if the data is in the correct format
	bool checkRange = false;											// the first column has to start at 0, but not end at 1
	double Vdata[n][2];													// array with the data of one file
	double Vall[n][nCCCV], Ahall[n][nCCCV];								// matrix with the voltage of each experiment, and the charge throughput of each experiment
	for(int i=0;i<nCCCV;i++){											// loop to read the voltage profile of each cycle
		if(lengths[i] > n){
			cerr<<"ERROR in determineCharacterisation::fitDiffusionAndRate. The arrays to store the data are too short. They have a length of "<<n<<" but need a length of at least "<<lengths[i]<<". Throw an error."<<endl<<flush;
			throw 10001;
			// if you get this error, increase the value of 'n' defined in this function ('max length of the data')
		}
		loadCSV_2colMatrix(names[i], lengths[i], Vdata);				// read the csv file with the voltage profile

		// check the data is in the correct format
		val = validOCV(checkRange, lengths[i], Vdata);
		if(!val){
			cerr<<"ERROR in determineCharacterisation::fitDiffusionAndRate. Input file "<<names[i]<<" has the wrong format. throwing an error"<<endl<<flush;
			throw 10000;
		}

		// reorganise the data in separate matrices for charge and voltage
		for(int j=0;j<lengths[i];j++){
			Ahall[j][i] = Vdata[j][0];
			Vall[j][i] = Vdata[j][1];
		}
	}

	// Calculate the error of the best fit
	double errcomb = 0;													// the combined error of all CCCV experiments for this combination of Dp, Dn, kp and kn
	double erri[nCCCV];													// the error for this CCCV cycle
	output.open(nameparam, std::ios_base::app);							// append the errors in the file with the parameter values
	for(int i=0;i<nCCCV;i++){											// loop through all CCCV experiments
		// Simulate this CCCV experiment
		CCCV(Crates[i], Ccuts[i], Tref, Dp, Dn, kp, kn, Rdc, ocvfit, M, nin, Vsim, Tsim, &nsim);

		// Get the measured data from this CCCV experiment
		for(int j=0;j<lengths[i];j++){
			Vdata[j][0] = Ahall[j][i];
			Vdata[j][1] = Vall[j][i];
		}

		// calculate the error (function defined in determineOCV.cpp)
		erri[i] = calculateError(false, lengths[i], nsim, Vdata, Vsim);
		errcomb += abs(erri[i])*weights[i];
	}
	output<<"combined RMSE"<<","<<errcomb<<"\n";

	// Write on which cycles this fit is based
	output<<"\n this was for the following cycles"<<"\n";
	output<<"name of the data file"<<","<<"length of the data file"<<","<<"C rate of the CC phase"<<","<<"C rate of the limit current for the CV phase"<<","<<"RMSE"<<","<<"weight"<<"\n";
	for(int i=0;i<nCCCV;i++){
		output<<names[i]<<","<<lengths[i]<<","<<Crates[i]<<","<<Ccuts[i]<<","<<erri[i]<<","<<weights[i] <<"\n";
	}
	output<<"\n";

	// then note down the settings used to get this fit
	output<<"\n";
	output<<"Below are the settings which produced this result"<<"\n";
	output<<"temperature"<<","<< Tref<<"\n";
	output<<"electrode surface"<<","<< ocvfit.elec_surf<<"\n";
	output<<"cathode active volume fraction"<<","<< ocvfit.ep<<"\n";
	output<<"anode active volume fraction"<<","<< ocvfit.en<<"\n";
	output<<"cathode thickness"<<","<< ocvfit.thickp<<"\n";
	output<<"anode thickness"<<","<< ocvfit.thickn<<"\n";
	output<<"cathode initial li-fraction at 50% SoC"<<","<< ocvfit.lifracpini<<"\n";
	output<<"anode initial li-fraction at 50% SoC"<<","<< ocvfit.lifracnini<<"\n";
	output<<"maximum voltage"<<","<< ocvfit.Vmax<<"\n";
	output<<"minimum voltage"<<","<< ocvfit.Vmin<<"\n";
	output<<"maximum li-concentration in the cathode"<<","<<ocvfit.cmaxp <<"\n";
	output<<"maximum li-concentration in the anode"<<","<<ocvfit.cmaxn <<"\n";
	output<<"name of the file with the cathode OCV curve"<<","<<ocvfit.namepos <<"\n";
	output<<"name of the file with the anode OCV curve"<<","<<ocvfit.nameneg <<"\n";

	// the search space
	output<<"\n";
	output<<"Below are the settings of the initial search space"<<"\n";
	output<<"number of steps in the search space for the diffusion constant"<<","<< nDstep<<"\n";
	output<<"minimum cathode diffusion constant"<<","<< Dpmin<<"\n";
	output<<"minimum anode diffusion constant"<<","<< Dnmin<<"\n";
	output<<"initial step size for the cathode diffusion constant"<<","<< Dpstep<<"\n";
	output<<"initial step size for the anode diffusion constant"<<","<< Dnstep<<"\n";
	output<<"logarithmic steps for the diffusion constant?"<<","<< logDstep<<"\n";
	output<<"number of steps in the search space for the rate constant"<<","<< nkstep<<"\n";
	output<<"minimum cathode rate constant"<<","<< kpmin<<"\n";
	output<<"minimum anode rate constant"<<","<< knmin<<"\n";
	output<<"initial step size for the cathode rate constant"<<","<< kpstep<<"\n";
	output<<"initial step size for the anode rate constant"<<","<< knstep<<"\n";
	output<<"logarithmic steps for the rate constant?"<<","<< logkstep<<"\n";
	output<<"number of steps in the search space for the DC resistance"<<","<< nrstep<<"\n";
	output<<"minimum DC resistance"<<","<< rmin<<"\n";
	output<<"initial step size in the search for the DC constant"<<","<< rstep<<"\n";
	output<<"number of levels in the search hierarchy"<<","<<hmax <<"\n";

	output.close();
}
