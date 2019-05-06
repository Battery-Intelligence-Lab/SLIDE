/*
 * detmineOCV.cpp
 *
 * The functions below can be used to find parameters relating to the OCV of a cell if a measured OCV curve is supplied.
 *
 * As input, the user has to supply the electrode OCV curves and the cell OCV curve.
 * See readOCVinput()
 *
 * Then the functions below calculate
 * 		the 'windows' of the half-cell OCV curve you want to use, i.e. the li-fractions of each electrode at the fully charged and fully discharged state
 * 		the amount of active material of each electrode
 *
 * These values can be used in the Cell-constructors:
 * 		the initial concentration (for a cell starting at 50% SoC)
 * 		the electrode surface, thickness, volume fraction, effective surface
 * 			changing the radius of the particles is not recommended because it means the spatial discretisation changes, and hence the model matrices have to be recalculated in Matlab
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <thread>
#include <cassert>
#include "determineOCV.h"
#include "ReadCSVfiles.h"
#include "Interpolation.h"

using namespace std;

bool validOCV(bool checkRange, int nin, double x[][2]){
	/*
	 * Function to check if the OCV curves supplied by the users have the correct format.
	 * the first column must have strictly increasing values starting at 0.
	 * 		the electrode OCV curves must end at 1 (they have the li-fractions, so the range must be from 0 to 1)
	 * 		the electrode OCV curve can end at any point
	 * the second column must have realistic voltages (0 <= V <= 10)
	 *
	 * If one of the conditions is not satisfied (i.e. the OCV curve is illegal)
	 * an error message is printed and the function returns false.
	 *
	 * IN
	 * checkRange 	boolean to indicate if we need to check that the range of the first column goes to 1
	 * 				set to true for the electrode OCV curves (where the first column gives the lithium fraction, so from 0 to 1)
	 * 				set to false for the cell OCV curve (where the first column gives the discharged charge, so from 0 to some positive value)
	 * nin 			number of rows in the curves
	 * x 			matrix with the OCV curve.
	 * 				the first column gives the 'x-values' (lithium fraction or discharged charge)
	 * 				the second column gives the 'y-values' (open circuit voltage)
	 *
	 * OUT
	 * bool 		true if the OCV curve has the correct format
	 */

	// Check if the range of the first column is valid
	bool range;										// boolean to indicate of the first column has the correct range
	if (checkRange){								// check the first column has the correct range, from 0 to 1
		range = abs(x[0][0]) < pow(10,-5) && abs(x[nin-1][0]- 1) < pow(10,-5); // allow a small error of e-5
		if(!range)
			cerr<<"illegal value in deterimineOCV::validOCV: the range of the first column is wrong. It goes from "<<x[0][0]<<" to "<<x[nin-1][0]<<" instead of from 0 to 1"<<endl<<flush;
	}
	else{											// check the first column has the correct starting point, 0
		range = abs(x[0][0]) < pow(10,-5);
		if(!range)
			cerr<<"illegal value in deterimineOCV::validOCV: the range of the first column is wrong. It starts at "<<x[0][0]<<" instead of starting at 0"<<endl<<flush;
	}

	// Loop through the curve to check that the values of the first column are increasing and the voltages in the second column are realistic.
	bool decreasei = false;							// boolean to indicate if this row has a decreasing value in the 1st column
	bool decrease = false;							// boolean to indicate if there are decreasing values in the 1st column
	bool unrealistici = false;						// boolean to indicate if this row an unrealistic value in the 2nd column
	bool unrealistic = false;						// boolean to indicate if there is an unrealistic value in the 2nd column
	for(int i=0;i<nin;i++){
		// check for increasing values
		if(i > 0)
			decreasei = (x[i][0] <= x[i-1][0]);
		if(decreasei){								// print an error message
			cerr<<"illegal value in deterimineOCV::validOCV: the first column has decreasing or double values at rows "<<i-1<<" and "<<i<<endl<<flush;
			cerr<<"values are "<<x[i-1][0]<<" and "<<x[i][0]<<endl<<flush;
			decrease = true;						// we have found a decreasing value
		}

		// check for realistic voltages
		unrealistici = (x[i][1] < 0) ||  (x[i][1] > 10);
		if(unrealistici){							// print an error message
			cerr<<"illegal value in deterimineOCV::validOCV: the second column has an illegal value of "<<x[i][1]<<" at row "<<i<<". The values have to be between 0 and 10"<<endl<<flush;
			unrealistic = true;						// we have found an unrealistic value
		}
	}

	// return true if the format is correct
	return range && !decrease && !unrealistic; 		// the range must be correct and we cannot have any decreasing or unrealistic values
}

void readOCVinput(string namepos, string nameneg, string namecell, int np, int nn, int ncell, double OCVp[][2], double OCVn[][2], double OCVcell[][2]){
	/*
	 * Function to read the csv files witht the OCV curves which the user has to supply as input.
	 *
	 * As input, the user has to supply the electrode OCV curves.
	 * These must come in the form of two csv files, each with two columns.
	 * The first column of each file must give the lithium fraction in increasing order (i.e. going from 0 to 1)
	 * The second column of each file must give the voltage vs li/li+ of the electrode at the corresponding li fraction
	 *
	 * Additionally, the user has to supply the cell OCV curve.
	 * This must be one csv file with 2 columns.
	 * The first column gives the discharged charge [Ah], starting from 0 (i.e. cell is fully charged) to c (i.e. the cell is fully discharged) [c is the capacity which can be discharged]
	 * The second column gives the OCV of the cell [V].
	 *
	 * IN
	 * namepos 	name of the CSV file with the OCV curve of the cathode
	 * nameneg 	name of the CSV file with the OCV curve of the anode
	 * namecell name of the CSV file with the OCV curve of the cell
	 * np 		number of data points in the file with the cathode OCV
	 * nn 		number of data points in the file with the anode OCV
	 * ncell 	number of data points in the file with the cell OCV
	 *
	 * OUT
	 * OCVp 	matrix with the cathode OCV
	 * OCVn 	matrix with the anode OCV
	 * OCVcell 	matrix with the cell OCV curve
	 *
	 * THROWS
	 * 10000 	one of the files has an illegal format
	 */

	// read the csv files
	try{
		loadCSV_2colMatrix(namepos, np, OCVp);
		loadCSV_2colMatrix(nameneg, nn, OCVn);
		loadCSV_2colMatrix(namecell, ncell, OCVcell);
	}
	catch(int e){
		cerr<<"ERROR in determineOCV::readOCVinput, an error "<<e<<" happened while reading the files. Throwing the error on."<<endl<<flush;
		throw e;
	}

	// Check the OCV curves are valid
	bool valp = validOCV(true, np, OCVp);						// check the cathode OCV curve. The first column must go from 0-1
	if(!valp){
		cerr<<"ERROR in determineOCV::readOCVinput, the file with the cathode OCV curve " <<namepos<< " has an illegal format"<<endl<<flush;
		throw 10000;
	}
	bool valn = validOCV(true, nn, OCVn);						// check the anode OCV curve. The first column must go from 0-1
	if(!valn){
		cerr<<"ERROR in determineOCV::readOCVinput, the file with the anode OCV curve " <<nameneg<< " has an illegal format"<<endl<<flush;
		throw 10000;
	}
	bool valcell = validOCV(false, ncell, OCVcell);				// check the cell OCV curve. The first column must start at 0 but can end anywhere
	if(!valcell){
		cerr<<"ERROR in determineOCV::readOCVinput, the file with the cell OCV curve " <<namecell<< " has an illegal format"<<endl<<flush;
		throw 10000;
	}
}

void discharge(int np, int nn, double OCVp[][2], double OCVn[][2], double cap, double AMp, double AMn, double cmaxp, double cmaxn, double sp, double sn,
		double Vend, int nin, int* nout, double OCV[][2], double OCVanode[][2], double OCVcathode[][2], double fp[], double fn[]){
	/*
	 * Function to simulate a CC discharge at a very low current, so we can ignore resistances and diffusion.
	 *
	 * It assumes the concentrations are uniform at all times, so we don't have to account for diffusion.
	 * Therefore, we don't need to use the full cell model.
	 * Rather, we can simulate two 'buckets of lithium' and we only need to exchange lithium between both.
	 *
	 * IN
	 * np		number of points on the cathode OCV curve
	 * nn 		number of points on the anode OCV curve
	 * OCVp 	cathode OCV curve
	 * OCVn 	anode OCV curve
	 * cap		capacity of the cell [Ah]
	 * AMp 		amount of active material on the cathode [m3], elec_surf * thicp * ep
	 * AMn 		amount of active material on the anode [m3], elec_surf * thickn * en
	 * cmaxp 	maximum li-concentration in the cathode [mol m-3]
	 * cmaxn	maximum li-concentration in the anode [mol m-3]
	 * sp 		start li-fraction at the cathode (between 0 and 1)
	 * sn 		start li-fraciton at the anode (between 0 and 1)
	 * Vend 	voltage when the OCV discharge should end [V]
	 * nin 		number of rows in the matrix OCV provided for output
	 *
	 * OUT
	 * nout 	number of point in the simulated OCV curve
	 * OCV 		OCV curve of the cell simulated with the given parameters
	 * 			the first column [i][0] gives the discharged charge in [Ah]
	 * 			the second column [i][1] gives the cell OCV in [V]
	 * OCVanode OCV of the anode
	 * 			the first column [i][0] gives the lithium fraction in the anode
	 * 			the second column [i][1] gives the voltage of the anode
	 * OCVcathode OCV of the cathode
	 * 			the first column [i][0] gives the lithium fraction in the cathode
	 * 			the second column [i][1] gives the voltage of the cathode
	 * fp		array with the lithium fraction in the cathode at 100% SoC, 50% SoC and 0% SoC
	 * fn	 	array with the lithium fraction in the anode at 100% SoC, 50% SoC and 0% SoC
	 *
	 * throws
	 * 10001	the output provided for output is too short.
	 */

	// *********************************************************** 1 variables ***********************************************************************

	// Variables
	double n = 1.0;						// number of electrons involved in the reaction
	double F = 96487.0;					// Faraday's constant
	double V;							// cell voltage in this time step
	double Ah = 0;						// discharged charge up to now
	double I = 0.1*cap;					// magnitude of the discharge current [A]
	double dt = 5;						// time step to use in the simulation [s]
	bool verbose = true;				// print error messages in the linear interpolation
	bool bound = true;					// check the surface concentration is within the allowed limits when doing linear interpolation. Throw an error if not
	bool end = false; 					// boolean to check if we have reached the end voltage
	bool rp;							// boolean to check if the li-fraction in the cathode is in the valid range (between 0 and 1)
	bool rn;							// boolean to check if the li-fraction in the anode is in the valid range (between 0 and 1)

	// If we store a data point for every time step, the output arrays will be very long.
	// Therefore, store only one value every 'nstore' time steps
	int nt = 0;							// number of time steps simulated so far
	double nstore = 100;				// store one data point every 100 time steps
	int i = 0; 							// number of stored values (index)

	// Store the lithium fractions at the start (=100% SoC)
	fp[0] = sp;
	fn[0] = sn;

	// *********************************************************** 2 discharge & measure voltage ***********************************************************************

	// a loop to keep discharging until we reach the end voltage
	double ocvpi, ocvni;
	while (!end){

		// get the open circuit voltage
		try{
			ocvpi = linInt_matrix(verbose, bound, OCVp, np, sp);
			ocvni = linInt_matrix(verbose, bound, OCVn, nn, sn);
			V = ocvpi - ocvni; 			// OCV = OCV_cathode - OCV_anode
		}
		catch(int e){
			cout<<"Error in deterimineOCV::discharge from linInt when getting the voltage. positive li-fraction is "<<sp<<" negaive li-fraction is "<<sn<<". Both should be between 0 and 1"<<endl<<flush;
			throw e;
		}

		// check if we want to store this point
		if(fmod(nt,nstore) == 0){																			// store if the number of time steps ('nt') is a multiple of 'nstore'
			// check the output arrays are long enough, throw an error if not
			if (nin <= nt){
				cerr<<"ERROR in deterimineOCV::discharge: the array provided for output has length "<<nin<<" but we need at least "<<nt+1<<" rows"<<endl<<flush;
				cerr<<"Voltage = "<<V<<", sp = "<<sp<<", sn = "<<sn<<endl<<flush;
				throw 10001;
			}
			// store the data point
			OCV[i][0] = Ah;																					// discharge charge
			OCV[i][1] = V;																					// open circuit voltage
			OCVanode[i][0] = sn;																			// anode lithium fraction
			OCVanode[i][1] = ocvni;																			// anode potential
			OCVcathode[i][0] = sp;																			// cathode lithium fraction
			OCVcathode[i][1] = ocvpi;																		// cathode potential
			i++;
		}

		// Discharge the cell for one more time step
			// the discharged charge is ahi =  I * dt / 3600 										[Ah]
			// the corresponding change in mol is mi = I * dt / (nF)								[mol]
			// the corresponding change in li-concentration in electrode e is ci_e = mi / AMe		[mol m-3]
			// the corresponding change in li-fraction in electrode is i si_e = ci_e / cmax_e		[-]
		Ah += I * dt / 3600;
		sp += I * dt / (n*F) / AMp / cmaxp;
		sn -= I * dt / (n*F) / AMn / cmaxn;
		nt ++;

		// stop if we have reached the end voltage, or if one of the li-fractions gets out of bounds
		rp = sp <= 0 || sp >= 1;																			// is the cathode lithium fraction is valid?
		rn = sn <= 0 || sn >= 1;																			// is the anode lithium fraction is valid?
		if ((V <= Vend)  || rp  || rn){
			// ensure we store this final value
			OCV[i][0] = Ah;
			OCV[i][1] = V;
			OCVanode[i][0] = sn;																			// anode lithium fraction
			OCVanode[i][1] = ocvni;																			// anode potential
			OCVcathode[i][0] = sp;																			// cathode lithium fraction
			OCVcathode[i][1] = ocvpi;																		// cathode potential
			i++;
			end = true;																						// indicate we have reached the end voltage
		}
	}

	// *********************************************************** 3 output parameters ***********************************************************************

	// Store the lithium fractions at the end (= 0% SoC)
	fp[2] = sp;
	fn[2] = sn;

	// Then the lithium fractions at 50% SoC are the average between the fractions at 100% and 0%
	// This is true because SoC is defined based on charge throughput, which is directly (linearly) linked with the lithium concentration
	fp[1] = (fp[0] + fp[2])/2.0;
	fn[1] = (fn[0] + fn[2])/2.0;

	*nout = i; 																								// the number of data points
}

double calculateError(bool bound, int ncell, int nsim, double OCVcell[][2], double OCVsim[][2]){
	/*
	 * Function to calculate the root mean square error between the OCV curve of the cell supplied by the user and the simulated OCV curve
	 *
	 * IN
	 * bound	boolean deciding what to do if the value of x is out of range of xdat for linear interpolation
	 * 			i.e. how to 'extend' OCVsim to the same capacity as OCVcell if OCVsim has a lower capacity
	 * 				if true, the value will be set to 0
	 * 				if false, the value will be set to the last point of OCVsim (e.g. 2.7)
	 * ncell	number of points in the OCV curve of the cell
	 * nsim 	number of points in the simulated OCV curve
	 * OCVcell 	OCV curve of the cell, 2 columns
	 * OCVsim 	simulated OCV curve, 2 columns
	 *
	 * OUT
	 * err 		RMSE between both curves
	 */

	double rmse = 0;			// root mean square error
	double ahi;					// discharged charge at this point
	double Vsimi;				// simulated OCV at this point

	bool verbose = false;		// don't print error messages since there will be errors if the simulated OCV curve is wrong

	// loop through all data points
	for(int i=0;i<ncell;i++){

		// Get the simulated OCV at the discharged charge of this point on the measured OCV curve of the cell
		ahi = OCVcell[i][0];
		try{					// interpolate the simulated OCV curve to find the simulated OCV at this discharged charge
			Vsimi = linInt_matrix(verbose, bound, OCVsim, nsim, ahi);
		}
		catch(int e){			// if bound is true, linInt throws an error if the x point is out of range (e.g. because the simulated OCV curve stops at 3Ah discharged while the OCV of the cell goes until 3.5Ah discharged
			Vsimi = 0;			// then set the simulated voltage to 0;
		}

		// Calculate the error
		rmse += pow(OCVcell[i][1] - Vsimi,2); 	// sum ( (Vcell[i] - Vsim[i])^2, i=0..ncell )
	}

	// Calculate the RMSE
	rmse = sqrt( rmse / ncell);
	return rmse;
}

void fitAMnAndStartingPoints(int hierarchy,int ap,double AMp,
		int nAMstep, double AMnstep, double AMnmin,
		int nsstep, double spstep, double snstep, double spmin, double snmin,
		int np, int nn, int ncell, string namepos, string nameneg, string namecell, double cmaxp, double cmaxn,
		double* err, double par[], int parindex[]){
	/*
	 * Function to scan a given search space for 3 parameters
	 * (amount of anode active material, starting point for the cathode and starting point for the anode).
	 * For each combination of the 3 parameters (and the given amount of cathode active material),
	 * the OCV curve is simulated and the error with the measured OCV curve is calculated.
	 * The best fit (i.e. three parameters with the lowest error on the OCV curve) is returned
	 *
	 * IN
	 * hierarchy 		level of the hierarchy in which we are
	 * ap 				index in the search space for the cathode active material
	 * AMp 				amount of cathode active material [m3]
	 * nAMstep 			number of steps in the search space for AMp and AMn
	 * AMnstep 			step size in the search space for AMn
	 * AMnmin			minimum value in the search space for AMn
	 * nsstep			number of steps in the search space for sp and sn
	 * spstep 			step size in the search space for sp
	 * snstep			step size in the search space for sn
	 * spmin 			minimum value in the search space for sp
	 * snmin 			minimum value in the search space for sn
	 * np 				number of data points in the cathode OCV curve
	 * nn 				number of data points in the anode OCV curve
	 * ncell			number of data points in the cell's OCV curve
	 * namepos 			name of the CSV file with the cathode OCV curve
	 * nameneg 			name of the CSV file with the anode OCV curve
	 * namecell			name of the CSV file with the cell's OCV curve
	 * cmaxp 			maximum lithium concentration in the cathode [mol m-3]
	 * cmaxn			maximum lithium concentration in the anode [mol m-3]
	 *
	 * OUT
	 * err 				lowest error for this amount of cathode active material
	 * par 				parameters giving the lowest error (AMp, AMn, sp, sn)
	 * parindex 		index indicating where in the search space the optimal parameters are (AMp, AMn, sp, sn)
	 */

	// *********************************************************** 1 variables ***********************************************************************

	// Read the OCV curves given as input (half-cell OCV curve and measured full cell OCV curve)
	double OCVp[np][2], OCVn[nn][2], OCVcell[ncell][2];
	try{
		readOCVinput(namepos, nameneg, namecell,np, nn, ncell, OCVp, OCVn, OCVcell);
	}
	catch(int e){
		cerr<<"Error in determineOCV::fitAMnAndStartingPoints, the input files have the wrong format."<<endl<<flush;
		return; 								// stop calculating because we can't do anything
	}

	// variables
	double AMn, sp, sn; 						// fitting parameters in this iteration
	double erri;								// error of this combination of parameters
	double errmin = 10000000000;				// lowest error encountered so far
	double Vend = OCVcell[ncell-1][1];			// minimum voltage of the measured OCV curve, i.e. the voltage until which the OCV curve should be simulated
	double cap = OCVcell[ncell-1][0];			// capacity of the cell [Ah] (the discharge Ah at the last point on the OCV curve)
	int nin = 7200*20;							// length of the array to store the simulated OCV curve
	double OCVsim[nin][2];						// array for the simulated OCV curve
	double ocvn[nin][2],ocvp[nin][2];			// array for the simulated voltage of each electrode
	int nsim;									// number of points in the simulated OCV curve
	double lifp[3], lifn[3];					// array to store the lithium fractions at 100% SoC, 50% SoC and 0% SoC (not needed here)


	// *********************************************************** 2 loop through the search space ***********************************************************************

	for(int an=0;an<nAMstep;an++){				// loop for the search space of AMn
		AMn = (AMnmin + AMnstep*an);
		for(int fp=0;fp<nsstep;fp++){			// loop through the search space of sp
			sp = spmin + spstep*fp;
			sp = max(0.0, sp);					// avoid that the starting points go out of range (the lithium fraction has to be between 0 and 1)
			sp = min(sp, 1.0);
			for(int fn=0;fn<nsstep;fn++){		// loop through the search space of sn
				sn = snmin + snstep*fn;
				sn = max(0.0, sn);				// avoid that the starting points go out of range (the lithium fraction has to be between 0 and 1)
				sn = min(sn, 1.0);

				// Simulate the OCV curve with these parameters
				try{
					discharge(np, nn, OCVp, OCVn, cap, AMp, AMn, cmaxp, cmaxn, sp, sn, Vend, nin, &nsim, OCVsim, ocvn, ocvp, lifp, lifn);
				}
				catch(int e){
					cout<<"Error in determineOCV::fitAMnAndStartingPoints: "<<e<<" for AMp = "<<AMp<<", AMn = "<<AMn<<", sp = "<<sp<<", sn = "<<sn<<endl<<flush;

					// just indicate this as a very large error, because it gives an infeasible combination of parameters
					erri = 10000000000;
				}

				// calculate the error between the simulated and measured OCV curve
				erri = calculateError(true, ncell, nsim, OCVcell, OCVsim);

				// Store the minimum error & parameters leading to this error
				if(erri < errmin){				// check if the error of this combination is lower than the best fit so far
					par[0] = AMp;
					par[1] = AMn;
					par[2] = sp;
					par[3] = sn;
					parindex[0] = ap;
					parindex[1] = an;
					parindex[2] = fp;
					parindex[3] = fn;
					errmin = erri;
				}
			}
		}
	}

	// Return the minimum error
	*err = errmin;
}

void oneLevelOCVfit(int hierarchy,
		int nAMstep, double AMpstep, double AMnstep, double AMpmin, double AMnmin,
		int nsstep, double spstep, double snstep, double spmin, double snmin,
		int np, int nn, int ncell, string namepos, string nameneg, string namecell, double cmaxp, double cmaxn,
		double* err, double par[], int parindex[]){
	/*
	 * Function to find the optimal points in one level of the hierarchical search for AMp, AMn, sp and sn.
	 * I.e. it goes through the search space of 4 variables
	 * 	amount of cathode active material
	 * 	three variables from fitAMnAndStartingPoint
	 * This function uses three threads to accelerate the computation
	 * The best fit (i.e. four parameters with the lowest error on the OCV curve) is returned and a csv file with the optimal value of the parameters is written.
	 *
	 * IN
	 * hierarchy 	level of the hierarchy in which we are
	 * nAMstep 		number of steps in the search space for AMp and AMn
	 * AMpstep		step size in the search space for AMp
	 * AMnstep 		step size in the search space for AMn
	 * AMpmin 		minimum value in the search space for AMp
	 * AMnmin		minimum value in the search space for AMn
	 * nsstep		number of steps in the search space for sp and sn
	 * spstep 		step size in the search space for sp
	 * snstep		step size in the search space for sn
	 * spmin 		minimum value in the search space for sp
	 * snmin 		minimum value in the search space for sn
	 * np 			number of data points in the cathode OCV curve
	 * nn 			number of data points in the anode OCV curve
	 * ncell		number of data points in the cell's OCV curve
	 * namepos 		ame of the CSV file with the cathode OCV curve
	 * nameneg 		name of the CSV file with the anode OCV curve
	 * namecell		name of the CSV file with the cell's OCV curve
	 * cmaxp 		maximum lithium concentration in the cathode [mol m-3]
	 * cmaxn		maximum lithium concentration in the anode [mol m-3]
	 *
	 * OUT
	 * err 			lowest error for this amount of cathode active material
	 * par 			parameters giving the lowest error (AMp, AMn, sp, sn)
	 * parindex 	indices indicating where in the search space the optimal parameters are (AMp, AMn, sp, sn)
	 */

	// Variables
	double errmin = 10000000000;				// lowest error encountered so far
	int ap1, ap2, ap3;							// indices for the search space of AMp
	double AMp1, AMp2, AMp3;					// amount of cathode active material for each
	double err1, err2, err3;					// lowest error for that amount of cathode active material
	double par1[4], par2[4], par3[4];			// parameters [AMp AMn sp sn] giving the lowest error for that amount of cathode active material
	int parindex1[4], parindex2[4], parindex3[4]; // indices indicating where in the search space the optimal parameters are

	// Loop through the search space of AMp, take 3 steps per iteration, each using one thread
	for(int ap=0;ap<nAMstep/3+1;ap++){			// loop for the search space of AMp in steps of 3

		// make the amount of cathode active material for the 3 steps
		ap1 = 3*ap;
		ap2 = 3*ap+1;
		ap3 = 3*ap+2;
		AMp1 = (AMpmin + AMpstep*ap1);
		AMp2 = (AMpmin + AMpstep*ap2);
		AMp3 = (AMpmin + AMpstep*ap3);

		// For each of the amounts of cathode material, search the best fit of the remaining 3 parameters (multi-threaded)
		thread t1(fitAMnAndStartingPoints, hierarchy, ap1, AMp1,nAMstep, AMnstep, AMnmin,nsstep, spstep, snstep, spmin, snmin,
				np, nn, ncell, namepos, nameneg, namecell, cmaxp, cmaxn,&err1, par1, parindex1);
		thread t2(fitAMnAndStartingPoints, hierarchy, ap2, AMp2,nAMstep, AMnstep, AMnmin,nsstep, spstep, snstep, spmin, snmin,
				np, nn, ncell, namepos, nameneg, namecell, cmaxp, cmaxn,&err2, par2, parindex2);
		thread t3(fitAMnAndStartingPoints, hierarchy, ap3, AMp3,nAMstep, AMnstep, AMnmin,nsstep, spstep, snstep, spmin, snmin,
				np, nn, ncell, namepos, nameneg, namecell, cmaxp, cmaxn,&err3, par3, parindex3);

		// join the threads again
		t1.join();
		t2.join();
		t3.join();

		// store the lowest error & the parameters leading to this error
		if(err1 < errmin){ 						// Check if the error from the first step is lower
			errmin = err1;
			for(int i=0;i<4;i++){
				par[i] = par1[i];
				parindex[i] = parindex1[i];
			}
		}
		if(err2 < errmin){ 						// Check if the error from the second step is lower
			errmin = err2;
			for(int i=0;i<4;i++){
				par[i] = par2[i];
				parindex[i] = parindex2[i];
			}
		}
		if(err3 < errmin){ 						// Check if the error from the third step is lower
			errmin = err3;
			for(int i=0;i<4;i++){
				par[i] = par3[i];
				parindex[i] = parindex3[i];
			}
		}
	}

	// Return the lowest error
	*err = errmin;

	// Print the best fit, and write in a CSV file
	cout<<"The best fit in hierarchy "<<hierarchy<<" is: AMp = "<<par[0]<<", AMn = "<<par[1]<<", sp = "<<par[2]<<", sn = "<<par[3]<<endl<<flush;
	ofstream output;
	// write the parameters
	string na = "OCVFit_" + to_string(hierarchy) + "_param.csv";
	output.open(na);
	output<<"AMp"<<","<<par[0]<<"\n";
	output<<"AMn"<<","<<par[1]<<"\n";
	output<<"start pos"<<","<<par[2]<<"\n";
	output<<"start neg"<<","<<par[3]<<"\n";
	output.close();
}

void hierarchicalOCVfit(int hmax, int nAMstep, double AMpstep, double AMnstep, double AMpmin, double AMnmin,
		int nsstep, double spstep, double snstep, double spmin, double snmin,
		int np, int nn, int ncell, string namepos, string nameneg, string namecell, double cmaxp, double cmaxn,
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
	 * then in hierarchical level n, the step size is dp* (2/(nstep-1))^(n-1)
	 * So after 'hmax' level, the accuracy is dp* (2/(nstep-1))^(hmax-1)
	 *
	 * IN
	 * hmax  		number of hierarchical steps we should take.
	 * 				the higher the value, the higher the accuracy of the fit will be but the longer the calculation will take
	 * nAMstep 		number of steps in the search space for AMp and AMn
	 * AMpstep		step size in the search space for AMp
	 * AMnstep 		step size in the search space for AMn
	 * AMpmin 		minimum value in the search space for AMp
	 * AMnmin		minimum value in the search space for AMn
	 * nsstep		number of steps in the search space for sp and sn
	 * spstep 		step size in the search space for sp
	 * snstep		step size in the search space for sn
	 * spmin 		minimum value in the search space for sp
	 * snmin 		minimum value in the search space for sn
	 * np 			number of data points in the cathode OCV curve
	 * nn 			number of data points in the anode OCV curve
	 * ncell		number of data points in the cell's OCV curve
	 * namepos 		ame of the CSV file with the cathode OCV curve
	 * nameneg 		name of the CSV file with the anode OCV curve
	 * namecell		name of the CSV file with the cell's OCV curve
	 * cmaxp 		maximum lithium concentration in the cathode [mol m-3]
	 * cmaxn		maximum lithium concentration in the anode [mol m-3]
	 *
	 * OUT
	 * err 			lowest error for this amount of cathode active material
	 * par 			parameters giving the lowest error (AMp, AMn, sp, sn)
	 */

	double erri; 									// error in this level of the hierarchy
	double pari[4];									// best fit parameters in this level of the hierarchy
	int parindexi[4];								// indices of the best fit parameters in this level of the hierarchy
	double AMpmax = AMpmin + AMpstep*(nAMstep-1);	// maximum values of the parameters in the search space
	double AMnmax = AMnmin + AMnstep*(nAMstep-1);
	double spmax = spmin + spstep*(nsstep-1);
	double snmax = snmin + snstep*(nsstep-1);

	// Loop for each level in the search
	for (int h=0;h<hmax;h++){

		// print the search space of this level
		cout<<"Start hierarchy level "<<h<<" with the following search spaces: "<<endl;
		cout<<"AMp: from "<<AMpmin<<" to "<<AMpmax<<" in "<<nAMstep<<" steps with magnitude "<<AMpstep<<endl;
		cout<<"AMn: from "<<AMnmin<<" to "<<AMnmax<<" in "<<nAMstep<<" steps with magnitude "<<AMnstep<<endl;
		cout<<"sp: from "<<spmin<<" to "<<spmax<<" in "<<nsstep<<" steps with magnitude "<<spstep<<endl;
		cout<<"sn: from "<<snmin<<" to "<<snmax<<" in "<<nsstep<<" steps with magnitude "<<snstep<<endl;

		// Calculate the best fit in this level
		oneLevelOCVfit(h, nAMstep, AMpstep, AMnstep, AMpmin, AMnmin,
			nsstep, spstep, snstep, spmin, snmin,
			np, nn, ncell, namepos, nameneg, namecell, cmaxp, cmaxn,
			&erri, pari, parindexi);

		// Update the search space
		// 	suppose the optimal value for a parameter p was in the search space at index i
		// 	then the new search space has as minimum value the value of p at i-1 and as maximum the value of p at i+1
		AMpmax = AMpmin + (parindexi[0]+1)*AMpstep;
		AMpmin = AMpmin + (parindexi[0]-1)*AMpstep;
		AMnmax = AMnmin + (parindexi[1]+1)*AMnstep;
		AMnmin = AMnmin + (parindexi[1]-1)*AMnstep;
		spmax = spmin + (parindexi[2]+1)*spstep;
		spmin = spmin + (parindexi[2]-1)*spstep;
		snmax = snmin + (parindexi[3]+1)*snstep;
		snmin = snmin + (parindexi[3]-1)*snstep;

		// 	the step size is chosen to produce the same number of steps as before
		AMpstep = (AMpmax - AMpmin)/(nAMstep-1);			// nAMstep-1 because for loop goes from 0 to nAMstep-1
		AMnstep = (AMnmax - AMnmin)/(nAMstep-1);
		spstep = (spmax - spmin)/(nsstep-1);
		snstep = (snmax - snmin)/(nsstep-1);
	}

	// Make the output parameters
	*err = erri;
	for(int i=0;i<4;i++)
		par[i] = pari[i];
}

void estimateOCVparameters(){
	/*
	 * Function which will find the parameters which best fit the OCV curve of the user.
	 *
	 * As input, the user has to supply the electrode OCV curves.
	 * These must come in the form of two csv files, each with two columns.
	 * The first column of each file must give the lithium fraction in increasing order (i.e. going from 0 to 1)
	 * The second column of each file must give the voltage vs li/li+ of the electrode at the corresponding li fraction
	 *
	 * Additionally, the user has to supply the measured cell OCV curve.
	 * This must be one csv file with 2 columns.
	 * The first column gives the discharged charge [Ah], starting from 0 (i.e. cell is fully charged) to c (i.e. the cell is fully discharged) [c is the capacity which can be discharged]
	 * The second column gives the OCV of the cell [V].
	 *
	 * The function finds optimal values for 4 parameters:
	 * 		cinip 		the initial lithium fraction for the cathode
	 * 		cinin 		the initial lithium fraction for the anode
	 * 		thickp		the thickness of the cathode
	 * 		thickn		the thickness of the anode
	 *
	 * The values for the other geometric parameters are assumed known
	 * 		elec_surf 	the geometric surface area of the electrode (assumed to be 0.982 m2, which is representative for an 18650 cell)
	 * 		ep			the volume fraction of active material on the cathode (assumed to be 50%)
	 * 		en 			the volume fraction of active material on the anode (assumed to be 50%)
	 * If you have different values, you have to change them in the code below (lines 791-793)
	 *
	 * The values of the parameters which give the best fit are written in a CSV file.
	 * Also, the simulated OCV curve with these parameters is written to a csv file so it can be compared with the measured one.
	 * 		The first column gives the simulated discharged charge
	 * 		The second column gives the simulated cell OCV
	 * 		The third column gives the simulated anode lithium fraction
	 * 		The fourth column gives the simulated anode OCV (vs Li/Li+)
	 * 		The fifth column gives the simulated cathode lithium fraction
	 * 		The sixth column gives the simulated cathode OCV (vs Li/Li+)
	 */

	// *********************************************************** 1 USER INPUT ***********************************************************************

	// input parameters
	string namepos = "OCVfit_cathode.csv";		// name of the file with the OCV curve of the cathode
	int np = 49;								// number of data points on the cathode OCV curve
	string nameneg = "OCVfit_anode.csv";		// name of the file with the OCV curve of the anode
	int nn = 63;								// number of data points on the anode OCV curve
	string namecell = "OCVfit_cell.csv";		// name of the file with the OCV curve of the cell
	int ncell = 74;								// number of data points on the cell OCV curve
	double cmaxp = 51385;						// maximum li-concentration in the cathode [mol m-3]
	double cmaxn = 30555;						// maximum li-concentration in the anode [mol m-3]

	// names of the output files
	string nameparam = "OCVfit_parameters.csv"; // name of the output file in which the parameters of the best fit will be written
	string nameOCV = "OCVfit_sim.csv";			// name of the output file in which the OCV curve with the best parameters will be written

	// Read the OCV curves
	double OCVp[np][2], OCVn[nn][2], OCVcell[ncell][2];
	try{
		readOCVinput(namepos, nameneg, namecell,np, nn, ncell, OCVp, OCVn, OCVcell);
	}
	catch(int e){
		cerr<<"Error in determineOCV::estimateOCVparameters, the input files have the wrong format."<<endl<<flush;
		return; 								// stop calculating because we can't do anything
	}

	// ****************************************** 2 define the search space for fitting parameters ***********************************************************************

	// We can play with 4 parameters:
		// AMp 	the amount of active material on the positive electrode
		// AMn 	the amount of active material on the negative electrode
		// sp	the starting point on the positive electrode OCV curve
		// sn 	the starting point on the negative electrode OCV curve

	// estimate the amount of active material needed to reach the cell capacity:
	double cap = OCVcell[ncell-1][0];			// capacity of the cell [Ah] (the discharge Ah at the last point on the OCV curve)
	double n = 1.0;								// number of electrons involved in the reaction
	double F = 96487.0;							// Faraday's constant
	double AMp_guess = cap*3600 / (n*F*cmaxp);	// the total charge in an electrode in Ah is given by; n * F * cmax * AM / 3600
	double AMn_guess = cap*3600 / (n*F*cmaxn);

	// Define the search space for the amount of active material on each electrode
	double step = 0.1;							// take steps of 10% of the guessed active material
	double AMmax = 5;							// the maximum amount of active material is 5 times the guessed amount
	double AMmin = 0;							// the minimum amount of active material is 0
	int nAMstep = (AMmax - AMmin)/step+1;		// the number of steps in the search space, +1 to include the end point
	double AMpmin = AMmin*AMp_guess;			// the lowest amount of AMp in the search space
	double AMnmin = AMmin*AMn_guess;			// the lowest amount of AMn in the search space
	double AMpstep = step*AMp_guess;			// the step size in the search space for AMp
	double AMnstep = step*AMn_guess;			// the step size in the search space for AMn

	// Define the search space for the initial lithium fractions at each electrode
	double df = 0.1;							// take steps of 10% of the lithium fraction
	double fmin = 0.0;							// the minimum lithium fraction is 0
	double fmax = 1.0;							// the maximum lithium fraction is 1
	int nsstep = (fmax - fmin)/df+1;			// number of steps in the search space
	double spmin = fmin;						// the lowest sp value in the search space
	double snmin = fmin;						// the lowest sn value in the search space
	double spstep = df;							// the step size in the search space for sp
	double snstep = df;							// the step size in the search space for sn

	// ***************************************************** 3 Fit the parameters ***********************************************************************

	// Call the hierarchical search algorithm, which does the fitting
	int hmax = 2;								// number of levels in the hierarchy to consider.
	double err;									// lowest error of the best fit
	double par[4];								// parameters of the best fit
	hierarchicalOCVfit(hmax, nAMstep, AMpstep, AMnstep, AMpmin, AMnmin,
			nsstep, spstep, snstep, spmin, snmin,
			np, nn, ncell, namepos, nameneg, namecell, cmaxp, cmaxn,
			&err, par);


	// ***************************************************** 4 write outputs ***********************************************************************

	// Print the best fit
	cout<<"The best fit is: AMp = "<<par[0]<<", AMn = "<<par[1]<<", sp = "<<par[2]<<", sn = "<<par[3]<<endl<<flush;
	ofstream output;

	// write the parameters in a csv file
	output.open(nameparam);
	output<<"AMp"<<","<<par[0]<<"\n";
	output<<"AMn"<<","<<par[1]<<"\n";
	output<<"start pos"<<","<<par[2]<<"\n";
	output<<"start neg"<<","<<par[3]<<"\n";
	output.close();

	// Simulate the best-fit OCV curve
	double Vend = OCVcell[ncell-1][1];			// minimum voltage of the OCV curve
	int nin = 7200*20;							// length of the array to store the simulated OCV curve
	double OCVsim[nin][2];						// array for the simulated OCV curve
	double OCVnsim[nin][2],OCVpsim[nin][2];		// array for the simulated voltage of each electrode
	int nsim;									// number of points in the simulated OCV curve
	double fp[3], fn[3];						// arrays to store the lithium fractions at 100%, 50% and 0% SoC
	discharge(np, nn, OCVp, OCVn, cap, par[0], par[1], cmaxp, cmaxn, par[2], par[3], Vend, nin, &nsim, OCVsim, OCVnsim, OCVpsim, fp, fn);

	// Write this best-fit OCV curve in a csv file so the user can check it using the matlab script readEstimateOCV.m
	output.open(nameOCV);
	for(int i=0;i<nsim;i++)
		output<<OCVsim[i][0]<<","<<OCVsim[i][1]<<","<<OCVnsim[i][0]<<","<<OCVnsim[i][1]<<","<<OCVpsim[i][0]<<","<<OCVpsim[i][1]<<"\n";
	output.close();

	// From the 4 fitted values, we need to determine the following parameters, needed by the single particle model implemented in Cell
		// elec_surf 	the geometric surface area of the electrode
		// thickp		the thickness of the cathode
		// thickn		the thickness of the anode
		// ep			the volume fraction of active material on the cathode
		// en 			the volume fraction of active material on the anode
		// cinip 		the initial lithium fraction for the cathode
		// cinin 		the initial lithium fraction for the anode

	// Assume the volume fractions are 50% and the electrode surface is 0.0982 m2.
	// Then we can calculate the thickness of the electrodes to give the desired amount of active material, which is given by:
	// 		AM = elec_surf * thick * e
	double ep = 0.5;
	double en = 0.5;
	double elec_surf = 0.0982;
	double thickp = par[0] / (elec_surf*ep);
	double thickn = par[1] / (elec_surf*en);

	// Append these parameters in the csv file where we had written the fitted parameters
	output.open(nameparam,std::ios_base::app);
	output<<"cathode volume fraction ep"<<","<<ep<<"\n";
	output<<"anode volume fraction en"<<","<<en<<"\n";
	output<<"electrode surface ele_surf"<<","<<elec_surf<<"\n";
	output<<"cathode thickness thickp"<<","<<thickp<<"\n";
	output<<"anode thickness thickn"<<","<<thickn<<"\n";
	output<<"cathode lithium fraction at 100% SoC"<<","<<fp[0]<<"\n";
	output<<"anode lithium fraction at 100% SoC"<<","<<fn[0]<<"\n";
	output<<"cathode lithium fraction at 50% SoC"<<","<<fp[1]<<"\n";
	output<<"anode lithium fraction at 50% SoC"<<","<<fn[1]<<"\n";
	output<<"cathode lithium fraction at 0% SoC"<<","<<fp[2]<<"\n";
	output<<"anode lithium fraction at 0% SoC"<<","<<fn[2]<<"\n";
	output<<"capacity of the cell in Ah"<<cap<<"\n";
	output<<"maximum voltage of the cell"<<OCVcell[0][1]<<"\n";
	output<<"minimum voltage of the cell"<<OCVcell[ncell-1][1]<<"\n";
	output<<"The error on the OCV curve with this fit is"<<","<<err<<"\n";

	// then note down the settings used to get this fit
	output<<"\n";
	output<<"Below are the settings which produced this result"<<"\n";
	output<<"maximum li-concentration in the cathode"<<","<<cmaxp <<"\n";
	output<<"maximum li-concentration in the anode"<<","<<cmaxn <<"\n";
	output<<"name of the file with the cathode OCV curve"<<","<<namepos <<"\n";
	output<<"name of the file with the anode OCV curve"<<","<<nameneg <<"\n";
	output<<"name of the file with the cell's OCV curve"<<","<<namecell <<"\n";
	output<<"capacity of the cell"<<","<<cap <<"\n";

	// the search space
	output<<"\n";
	output<<"Below are the settings of the initial search space"<<"\n";
	output<<"relative step size in the search for active material"<<","<<step<<"\n";
	output<<"relative minimum amount of active material"<<","<<AMmin <<"\n";
	output<<"relative maximum amount of active material"<<","<<AMmax <<"\n";
	output<<"step size in the search for the starting li-fraction"<<","<<df <<"\n";
	output<<"minimum li-fraction"<<","<<fmin <<"\n";
	output<<"maximum li-fraction"<<","<<fmax <<"\n";
	output<<"number of levels in the search hierarchy"<<","<<hmax <<"\n";

	output.close();
}
