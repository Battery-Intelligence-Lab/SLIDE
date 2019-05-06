/*
 * BasicCycler.hpp
 *
 * Header for the class implementing a basic cycler.
 * A basic cycler simulates a battery tester (and can be programmed similarly)
 * It offers functions to load a cell with a CC and/or CV.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_BASICCYCLER_HPP_
#define SRC_BASICCYCLER_HPP_

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <assert.h>
#include <memory>
#include <thread>

#include "Cell.hpp"

#define printCrit 1				// threshold of verbose of when to print error messages for critical errors
#define printNonCrit 2			// threshold of verbose of when to print error messages for noncritical errors
#define printCyclerFunctions 3	// threshold of verbose of when to print the start and end of functions of the BasicCycler
#define printCyclerHighLevel 4 	// threshold of verbose of when to print the high-level flow of the program in the BasicCycler
#define printCyclerDetail 5 	// threshold of verbose of when to print the low-level detailed flow of the program in the BasicCycler
#define printfindCVcurrentDetail 6 // threshold of verbose of when to print the details of how the current for the CV phase is found

using namespace std;

class BasicCycler {

protected:
	Cell c;						// cell connected to the basicCycler
	string ID;					// identification string for this basicCycler, will be included in the name of the subfolder in which data files are stored.
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

	// data of the cell which is being cycled
	int CyclingDataTimeInterval;// time resolution at which the cell data is stored [s], if 0 no data is stored
								// 		the higher the time resolution, the more accurate the stored data is
	int index;					// number of data points stored so far
	int fileIndex;				// number of data files written so far
	int maxLength; 				// length of the arrays in which the results are stored before they are written to a csv file
	double Iout[100000];		// current of the cell at every step [A]
	double Vout[100000];		// voltage of the cell at every step [V]
	double OCVpout[100000];		// cathode potential of the cell at every step [V]
	double OCVnout[100000];		// anode potential of the cell at every step [V]
	double Tout[100000];		// temperature of the cell at every step [K]
	double timeCha;				// cumulative time spent on charging since the start of this data collection [s]
	double timeDis;				// cumulative time spent on discharging since the start of this data collection [s]
	double timeRes;				// cumulative time spent on rest since the start of this data collection [s]
	double AhCha;				// cumulative charged throughput since the start of this data collection [A]
	double AhDis;				// cumulative discharged charge throughput since the start of this data collection [A]
	double WhCha;				// cumulative charged energy throughput since the start of this data collection [Wh]
	double WhDis;				// cumulative discharged energy throughput since the start of this data collection [Wh]
	double timeChaout[100000];	// cumulative time spent on charging since the start at every step [s]
	double AhChaout[100000];	// cumulative charged throughput since the start at every step [A]
	double WhChaout[100000];	// cumulative charged energy throughput since the start at every step [Wh]
	double timeDisout[100000];	// cumulative time spent on discharging since the start at every step [s]
	double AhDisout[100000];	// cumulative discharged charge throughput since the start at every step [A]
	double WhDisout[100000];	// cumulative discharged energy throughput since the start at every step [Wh]
	double timeResout[100000];	// cumulative time spent on rest since the start at every step [s]

	void storeResults(double I, double v, double ocvp, double ocvn, double tem); 														// store the cycling data of a cell
	int setCurrent(double I, double Vupp, double Vlow);																					// auxiliary function of CC_t_V to set the current
	void findCVcurrent_recursive(double Imin, double Imax, int sign, double Vset, double dt, bool blockDegradation, double* Il, double* Vl); // auxiliary function to solve the nonlinear equation to keep the voltage constant

public:
	BasicCycler(Cell& ci, string IDi, int verbose, int CyclingDataTimeIntervali);
	virtual ~BasicCycler();
	Cell& getCell();																													// returns (a reference to) the cell of the basicCycler
	void setCyclingDataTimeResolution(int timeResolution);																				// change the time resolution of the data collection

	// Functions to write the cycling data of the cell
	void writeCyclingData();																											// writes the data of the cell to a csv file with the standard name
	void writeCyclingData(string name, bool clear);																						// writes the data of the cell to the specified csv file
	void returnCyclingData(int nin, double Ah[], double V[], double T[], int* nout);													// return the voltage of the cell instead of writing it to the file

	// cycle battery at constant current
	virtual int CC_t_V(double I, double dt, bool blockDegradation, double time, double Vupp, double Vlow, double* ahi, double* whi, double* timei); // CC cycle with a time and two voltage end-condition
	virtual int CC_t(double I, double dt, bool blockDegradation, double time, double* ahi, double* whi, double* timei); 				// CC cycle for a fixed amount of time
	virtual int CC_V(double I, double dt, bool blockDegradation, double Vset, double* ahi, double* whi, double* timei); 				// CC cycle until a given voltage is reached
	virtual void CC_halfCell_full(double I, double dt, bool pos, int nin, double OCVi[], double* ahi, int* n);							// CC cycle only one electrode

	// cycle battery at constant voltage
	virtual void findCVcurrent(double Vset, double dt, bool blockDegradation, double* Il, double* Vl);									// find the current needed to keep the voltage constant at the specified value
	virtual int CV_t_I(double V, double dt, bool blockDegradation, double time, double Icut, double* ahi, double* whi, double* timei); 	// CV cycle with both a time and current limit
	virtual void CV_t(double V, double dt, bool blockDegradation, double time, double* ahi, double* whi, double* timei); 				// CV cycle for a fixed amount of time
	virtual void CV_I(double V, double dt, bool blockDegradation, double Icut, double* ahi, double* whi, double* timei); 				// CV cycle until the current is below a threshold value

	// hybrid CC CV
	virtual int CC_t_CV_t(double I, double dt, bool blockDegradation, double time, double Vupp, double Vlow, double* ahi, double* whi, double* timei);// load the cell for a given time, with a CC and switch automatically to a CV if a voltage limit is reached
	virtual void CC_V_CV_I(double Crate, double Vset, double Icut, double dt, bool blockDegradation, double* ahi, double* whi, double* timei); 	// bring the cell to a specified voltage, starting with a CC phase, followed by a CV phase

	// current profile
	virtual int followI(int nI, string nameI, bool blockDegradation, int limit, double Vupp, double Vlow, double* ahi, double* whi, double* timeio); // follow a predefined current pattern

};
#endif /* SRC_BASICCYCLER_HPP_ */
