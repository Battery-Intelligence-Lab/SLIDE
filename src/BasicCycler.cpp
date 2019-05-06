/*
 * BasicCycler.cpp
 *
 * Class to implement a 'basic cycler', which represents the battery tester (and can be used similarly)
 *
 * A basic cycler implements functions to load a cell with a (combination of) constant current (CC) and constant voltage (CV).
 * A basic cycler automatically handles the so-called 'cycling data', which are periodic measurements of the cell's voltage,  temperature and current.
 * This data is stored locally in arrays, and when the arrays are full (or when the user 'pushes' the data), the cycling data is written to a csv file.
 * The csv files are grouped in a subfolder of this project folder. The name of the subfolder is the identifier of this basic cycler.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "BasicCycler.hpp"
#include "ReadCSVfiles.h"
#include <direct.h>

using namespace std;

BasicCycler::BasicCycler(Cell& ci, string IDi, int verbosei, int CyclingDataTimeIntervali) {
	/*
	 * Constructor of a BasicCycler.
	 *
	 * IN
	 * ci			cell connected to the cycler (call by reference, i.e. the cell given to the cycler will be changed during execution)
	 * IDi 			identification string of for this BasicCycler.
	 * 				A new folder is made to store the data of this cell
	 * 				therefore, IDi must obey all restrictions for naming folders on your operating systems
	 * 				e.g. for Windows, the following are not allowed < > : " \ | / ? *
	 * 				https://docs.microsoft.com/en-us/windows/desktop/fileio/naming-a-file
	 * 				spaces are also not recommended
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
	 * CyclingDataTimeIntervali 	time interval at which the cycling data must be stored, >= 0 [s]
	 * 				if 0, no data is stored
	 * 				if < 0, no folder is even created to store results. The value is later set to 0 such that no data is stored
	 *
	 * Throws
	 * 1000 		the value of CyclingDataTimeIntervali is negative (it must be >= 0)
	 * 1019 		illegal identification string
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::BasicCycler starting"<<endl;

	c = ci;						// cell connected to the basicCycler
	ID = IDi;					// identification string for this basicCycler, also the name of the folder in which the data will be stored.
	verbose = verbosei;			// decide how verbose the output is

	// data of the cell which is being cycled
	CyclingDataTimeInterval = CyclingDataTimeIntervali;		// time resolution at which the cell data is stored [s], if 0 no data is stored
	index = 0;					// number of data points stored so far in the internal arrays
	fileIndex = 0;				// number of results files written so far
	maxLength = 100000; 		// length of the arrays in which the results are stored before they are written to a csv file
								// if you change this value, you also have to change the length of the arrays defined in BasicCycler.hpp

	// initialise all the data arrays at 0
	timeCha = 0;
	timeDis = 0;
	timeRes = 0;
	AhCha = 0;
	AhDis = 0;
	WhCha = 0;
	WhDis = 0;
	for(int i=0;i<maxLength;i++){
		timeChaout[i] = 0;
		AhChaout[i] = 0;
		WhChaout[i] = 0;
		timeDisout[i] = 0;
		AhDisout[i] = 0;
		WhDisout[i] = 0;
		timeResout[i] = 0;
		Iout[i] = 0;
		Vout[i] = 0;
		OCVpout[i] = 0;
		OCVnout[i] = 0;
		Tout[i] = 0;
	}

	// Make a subfolder for the data of this BasicCycler
	if(CyclingDataTimeIntervali >= 0){
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::BasicCycler is making a subfolder to store the results"<<endl;
		int s = mkdir((ID).c_str());
		if (s != 0){
			cerr<<"ERROR in BasicCycler::BasicCycler. The subfolder "<<ID<<" for the data could not be created. "<<endl;
			cout<<"This might have two reasons: either the directory already exists, or you are using some forbidden characters in the name."<<endl;
			cout<<"In the first case, the data from the folder would be overwritten so this is not allowed. Please delete the folder and try again"<<endl;
			cout<<"In the latter case, don't use spaces or special characters in the prefix or identification string. Underscores and dots are allowed"<<endl;
			throw 1019;
		}
	}
	else{
		CyclingDataTimeInterval = 0; // set to 0 to indicate we are not storing data
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::BasicCycler is NOT making a subfolder to store the results"<<endl;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::BasicCycler terminating"<<endl;
}

BasicCycler::~BasicCycler() {
}

Cell& BasicCycler::getCell(){
	return c;		// return a reference to the cell
}

void BasicCycler::setCyclingDataTimeResolution(int timeResolution){
	/*
	 * Change the data collection settings.
	 *
	 * IN
	 * timeResolution 	time interval in which the cycling data should be stored, >= 0 [s]
	 * 					if 0, no data is stored
	 * Throws
	 * 1000 			the value of timeResolution is negative (it must be >= 0)
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::setCyclingDataTimeResolution starting"<<endl;

	if (timeResolution<0){
		cerr<<"ERROR in BasicCycler::BasicCycler. The value for the data collection time interval is "<<timeResolution<<". It must be >= 0"<<endl<<flush;
		throw 1000;
	}

	CyclingDataTimeInterval = timeResolution;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::setCyclingDataTimeResolution terminating"<<endl;
}

void BasicCycler::writeCyclingData(){
	/*
	 * Function which will write the cycling data stored so far to a csv file
	 * The name of the CSV file is 'CyclingData_', followed by the index of the data file
	 * I.e. the name of the first data file is CyclingData_0.csv, the 2nd file is called CyclingData_1.csv, etc.
	 * The file is written in the subfolder of this BasicCycler, the name of the subfolder is the identification string of this BasicCycler.
	 *
	 * data written per column:
	 * 		total_time 		the total time in seconds since the start of data recording
	 * 		Ah_throughput 	the total charge throughput in [Ah] since the start of this data batch
	 * 		Wh_throughput 	the total energy throughput in [Wh] since the start of this data batch
	 * 		I 				the cell current in [A], positive for discharging, negative for charging
	 * 		V 				the cell voltage in [V]
	 * 		OCVp 			the cathode potential in [V]
	 * 		OCVn 			the anode potential in [V]
	 * 		Temperature 	the cell temperature in [K]
	 * 		charge_time 	the total time spend on charging in seconds since the start of this data batch
	 * 		charge_Ah 		the total charged charge in [Ah] since the start of this data batch
	 * 		charge_Wh 		the total charged energy in [Wh] since the start of this data batch
	 * 		discharge_time 	the total time spend on discharging in seconds since the start of this data batch
	 * 		discharge_Ah 	the total discharged charge in [Ah] since the start of this data batch
	 * 		discharge_Wh 	the total discharged energy in [Wh] since the start of this data batch
	 * 		rest_time		the total time spend on resting in seconds since the start of this data batch
	 *
	 * THROWS
	 * 1001	could not open file
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::writeCyclingData() starting"<<endl;

	// store data if we are collecting data and if there is new data available
	if(CyclingDataTimeInterval != 0 && index > 0){

		// Make the name of the results file. It is 'CyclingData_, followed by the index in the correct subfolder.
		string name = "CyclingData_" + to_string(fileIndex) + ".csv";			// name of the csv file
		string fullName = ".\\"+ID+"\\"+name;									// include the subfolder in the full name
		ofstream output(fullName);												// open the file

		if (output.is_open()){

			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::writeCyclingData() has opened "<<fullName<<" and is writing data"<<endl;

			// Write the data
			for(int i=0;i<index;i++){											// loop for each row (one row is one data point)
				output<<(timeChaout[i] + timeDisout[i] + timeResout[i])<<","; 	// total time [sec]
				output<<(AhChaout[i] + AhDisout[i])<<",";						// total charge throughput [Ah]
				output<<(WhChaout[i] + WhDisout[i])<<",";						// total energy throughput [Wh]
				output<<Iout[i]<<","<<Vout[i]<<","<<OCVpout[i]<<","<<OCVnout[i]<<","<<Tout[i]<<","; // I V OCV_pos OCV_neg T
				output<<timeChaout[i]<<","<<AhChaout[i]<<","<<WhChaout[i]<<",";	// time on charge, charged charge, charged energy
				output<<timeDisout[i]<<","<<AhDisout[i]<<","<<WhDisout[i]<<",";	// time on discharge, discharged charge, discharged energy
				output<<timeResout[i]<<"\n";									// time on rest
			}
			output.close();														// close the file
		}
		else{																	// the file could not be opened for some reason. E.g. because the subfolder doesn't exist
			cerr<<"ERROR on BasicCycler::writeCyclingStates(). File "<<fullName<<" could not be opened. Throwing an error"<<endl;
			throw 1001;
		}

		// clear the data stored in the array, reset the counter
		for(int i=0;i<maxLength;i++){
			timeChaout[i] = 0;
			AhChaout[i] = 0;
			WhChaout[i] = 0;
			timeDisout[i] = 0;
			AhDisout[i] = 0;
			WhDisout[i] = 0;
			timeResout[i] = 0;
			Iout[i] = 0;
			Vout[i] = 0;
			OCVpout[i] = 0;
			OCVnout[i] = 0;
			Tout[i] = 0;
		}
		index = 0;

		// reset the cumulative variables
		timeCha = 0;
		timeDis = 0;
		timeRes = 0;
		AhCha = 0;
		AhDis = 0;
		WhCha = 0;
		WhDis = 0;

		// Increase the counter for the number of data files we have written
		fileIndex ++;
	}
	else{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::writeCyclingData() didn't write data because we are not collecting data (time resolution is "<<CyclingDataTimeInterval<<") or because there is no new data (number of data points so far:  "<<index<<endl;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::writeCyclingData() terminating"<<endl;
}

void BasicCycler::writeCyclingData(string name, bool clear){
	/*
	 * write the cycling data to a csv file with the specified name (instead of to the standard name)
	 * data written per column:
	 * 		total_time 		the total time in seconds since the start of data recording
	 * 		Ah_throughput 	the total charge throughput in [Ah] since the start of this data batch
	 * 		Wh_throughput 	the total energy throughput in [Wh] since the start of this data batch
	 * 		I 				the cell current in [A], positive for discharging, negative for charging
	 * 		V 				the cell voltage in [V]
	 * 		OCVp 			the cathode potential in [V]
	 * 		OCVn 			the anode potential in [V]
	 * 		Temperature 	the cell temperature in [K]
	 * 		charge_time 	the total time spend on charging in seconds since the start of this data batch
	 * 		charge_Ah 		the total charged charge in [Ah] since the start of this data batch
	 * 		charge_Wh 		the total charged energy in [Wh] since the start of this data batch
	 * 		discharge_time 	the total time spend on discharging in seconds since the start of this data batch
	 * 		discharge_Ah 	the total discharged charge in [Ah] since the start of this data batch
	 * 		discharge_Wh 	the total discharged energy in [Wh] since the start of this data batch
	 * 		rest_time		the total time spend on resting in seconds since the start of this data batch
	 *
	 * IN
	 * name 	name of the CSV file in which the cycling data will be written
	 * clear 	if true, the data buffer is cleared after writing the file
	 * 			if false, the data is kept in the buffer so it will be written again if one of the writeCyclingData-functions is called.
	 *
	 * THROWS
	 * 1001 	could not open the file
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::writeCyclingData(string, bool) starting"<<endl;

	if(CyclingDataTimeInterval != 0 && index > 0){

		// Open the file
		string fullName = ".\\"+ID+"\\"+name; 									// append the subfolder to the name to write the file in the subfolder
		ofstream output(fullName);
		if (output.is_open()){

			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::writeCyclingData() has opened "<<fullName<<" and is writing data"<<endl;

			// Write the data
			for(int i=0;i<index;i++){											// loop for each row (one row is one data point)
				output<<(timeChaout[i] + timeDisout[i] + timeResout[i])<<","; 	// total time [sec]
				output<<(AhChaout[i] + AhDisout[i])<<",";						// total charge throughput [Ah]
				output<<(WhChaout[i] + WhDisout[i])<<",";						// total energy throughput [Wh]
				output<<Iout[i]<<","<<Vout[i]<<","<<OCVpout[i]<<","<<OCVnout[i]<<","<<Tout[i]<<","; // I V OCV_pos OCV_neg T
				output<<timeChaout[i]<<","<<AhChaout[i]<<","<<WhChaout[i]<<",";	// time on charge, charged charge, charged energy
				output<<timeDisout[i]<<","<<AhDisout[i]<<","<<WhDisout[i]<<",";	// time on discharge, discharged charge, discharged energy
				output<<timeResout[i]<<"\n";									// time on rest
			}
			output.close();														// close the file
		}
		else{																	// the file could not be opened for some reason. E.g. because the subfolder doesn't exist
			cerr<<"ERROR on BasicCycler::writeCyclingStates(string, bool). File "<<fullName<<" could not be opened. Throwing an error"<<endl;
			throw 1001;
		}

		// clear the data stored in the array, reset the counter
		if(clear){
			for(int i=0;i<maxLength;i++){
				timeChaout[i] = 0;
				AhChaout[i] = 0;
				WhChaout[i] = 0;
				timeDisout[i] = 0;
				AhDisout[i] = 0;
				WhDisout[i] = 0;
				timeResout[i] = 0;
				Iout[i] = 0;
				Vout[i] = 0;
				OCVpout[i] = 0;
				OCVnout[i] = 0;
				Tout[i] = 0;
			}
			index = 0;

			// reset the cumulative variables
			timeCha = 0;
			timeDis = 0;
			timeRes = 0;
			AhCha = 0;
			AhDis = 0;
			WhCha = 0;
			WhDis = 0;
		}

		// don't increase the file index because we haven't written another 'standard' cycling data file
	}
	else{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::writeCyclingData(string, bool) didn't write data because we are not collecting data (time resolution is "<<CyclingDataTimeInterval<<") or because there is no new data (number of data points so far:  "<<index<<endl;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::writeCyclingData(string, bool) terminating"<<endl;
}

void BasicCycler::returnCyclingData(int nin, double Ah[], double V[], double T[], int* nout){
	/*
	 * Function which returns the cycling data arrays.
	 * No csv files are written (and the data isn't cleared from the internal memory)
	 *
	 * IN
	 * nin 		length of the arrays provided
	 *
	 * OUT
	 * Ah 		array with the cumulative charge throughput
	 * V 		array with the voltage at every step
	 * T 		array with the temperature at every step
	 * nout 	number of data points written in the arrays
	 *
	 * THROWS
	 * 1002 	the arrays provided are too small
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::returnCyclingData starting"<<endl;

	if(nin < index){
		cerr<<"ERROR in BasicCycler::returnCyclingData. The arrays provided have a length "<<nin<<" but there are "<<index<<" data points. Throwing an error"<<endl<<flush;
		throw 1002;
	}

	for(int i=0;i<index;i++){
		Ah[i] = abs(AhChaout[i]) + abs(AhDisout[i]);
		V[i] = Vout[i];
		T[i] = Tout[i];
	}

	*nout = index;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::returnCyclingData terminating"<<endl;
}

void BasicCycler::storeResults(double I, double v, double ocvp, double ocvn, double tem){
	/*
	 * Function to store the cycling data to the internal memory.
	 *
	 * IN
	 * I 			current in this time step [A]
	 * v 			voltage in this time step [V]
	 * ocvp			cathode OCV in this time step [V]
	 * ocvn			anode OCV in this time step [V]
	 * tem			cell temperature in this time step [K]
	 */

	if(verbose >= printCyclerDetail)						// note: this function is called at a high time-granularity so only print its start and end if we want to print details
		cout<<"BasicCycler::storeResults starting"<<endl;

	// check if there is new data
	bool newdata;
	if (index == 0)
		newdata = true;										// this is the first data point
	else{
		bool newcha = timeChaout[index-1] != timeCha;		// we have spent time charging
		bool newdis = timeDisout[index-1] != timeDis;		// we have spent time discharging
		bool newres = timeResout[index-1] != timeRes;		// we have spent time resting
		newdata = newcha || newdis || newres;
	}

	// If CyclingDataTimeInterval is 0, or if no time has passed since the previous data point, no data is stored
	if(CyclingDataTimeInterval != 0 && newdata){

		// If the internal arrays are full, write the results to csv files & clear the internal arrays
		if(index > (maxLength-1)) 							// length -1 because you start at row [0]
			writeCyclingData();

		// add the data to the internal arrays
		Iout[index] = I;
		Vout[index] = v;
		OCVpout[index] = ocvp;
		OCVnout[index] = ocvn;
		Tout[index] = tem;
		timeChaout[index] = timeCha;
		AhChaout[index] = AhCha;
		WhChaout[index] = WhCha;
		timeDisout[index] = timeDis;
		AhDisout[index] = AhDis;
		WhDisout[index] = WhDis;
		timeResout[index] = timeRes;

		// increase the counter for the number of data points stored
		index++;
	}

	if(verbose >= printCyclerDetail)								// note: this function is called at a high time-granularity so only print its start and end if we want to print details
		cout<<"BasicCycler::storeResults terminating"<<endl;
}

int BasicCycler::setCurrent(double I, double Vupp, double Vlow){
	/*
	 * Auxiliary function of CC_t_V which tries to set the specified current to the cell without violating the voltage limits.
	 * It restores the initial state if the current can't be set
	 *
	 * IN
	 * I 			load current [A], positive = discharge, negative = charge
	 * Vupp			upper voltage limit, must be in the voltage limit of the cell (0 Cell.Vmin <= Vupp <= Cell.Vmax) [V]
	 * Vlow			lower voltage limit, must be in the voltage limit of the cell (0 Cell.Vmin <= Vlow <= Cell.Vmax) [V]
	 *
	 * OUT
	 * int 			which end condition was reached
	 * 					-3 	the minimum cell voltage was exceeded
	 * 					-2 	the maximum cell voltage was exceeded
	 * 					0 	an error occurred
	 * 					1 	the current could be set
	 * 					2 	the upper voltage limit was exceeded while charging the cell
	 * 					3 	the lower voltage limit was exceeded while discharging the cell
	 *
	 * THROWS
	 * 1004 		illegal voltage or current input
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::setCurrent with current = "<<I<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<" is starting"<<endl;

	// variables
	double v, ocvp, ocvn, etap, etan, rdrop, tem;
	State sini;												// initial state of the cell
	double Iini;											// initial current of the cell
	int endcriterion = 99;									// integer indicating why the function terminated
	bool check = false;										// don't check if the battery state is valid after setting a current in the underlying functions because we do the check here.

	// get the initial state
	c.getStates(sini, &Iini);

	// ************************************************ 1 check that the current has the correct sign ***********************************************************************
	// Check if the sign of the current and voltage limits are compatible.
	// i.e. if you want to charge the cell, the upper voltage limit must be above the OCV of the cell
	// and if you are discharging the cell, the lower cell voltage must be below the OCV of the cell
	// else you can never reach the specified cell voltage with the specified current
	// we use the OCV instead of the cell voltage because we want to allow the case where only the magnitude of the current is wrong
	// 	i.e. you need to charge to reach the specified voltage but at the specified current you overshoot the voltage limit (e.g. due to the resistive voltage drop)
	// 	and similar for discharge.
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::setCurrent is checking if the current has the correct sign"<<endl;

	// Get the OCV of the cell
	try{
		c.setI(verbose >= printCrit, check, 0);				// set the current to 0
		c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in BasicCycler::setCurrent when getting the OCV of the cell: "<<e<<". Throwing it on"<<endl;
		throw e;
	}

	// if we are discharging, check that Vlow is below the OCV of the cell
	if (I > 0){
		if(v < Vlow){										// throw an error if it isn't
			if(verbose >= printCrit)
				cerr<<"Error in BasicCycler::setCurrent. The current has the wrong sign. You are trying to discharge the cell at a current of "<<I
				<<" but the OCV of the cell "<<v<<" is already below the lower voltage limit of "<<Vlow
					<<" so you can never reach this lower voltage limit while discharging."<<endl<<flush;
			throw 1004;
		}
	}

	// if we are charging, check that Vupp is above the OCV of the cell
	if(I < 0){												// throw an error if it isn't
		if(v > Vupp){
			if(verbose >= printCrit)
				cerr<<"Error in BasicCycler::setCurrent. The current has the wrong sign. You are trying to charge the cell at a current of "<<I
				<<" but the OCV of the cell "<<v<<" is already above the upper voltage limit of "<<Vupp
					<<" so you can never reach this lower voltage limit while charging."<<endl<<flush;
			throw 1004;
		}
	}

	// ************************************************ 2 set the current ***********************************************************************

	// Try setting the cell current, and don't check if the resulting state is valid
	try{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::setCurrent is setting the cell current"<<endl;
		c.setI(verbose >= printNonCrit, check, I);			// check is false, so we don't throw an error if the resulting battery state is valid
	}
	catch(int e){
		// no error should happen since we didn't check the battery state for validity
		if(verbose >= printCrit)
			cout<<"error in BasicCycler::setCurrent when setting the current of "<<I<<"A, error: "<<e<<". Throwing it on"<<endl;
		throw e;
	}


	// ************************************************ 3 check the voltage limits ****************************************************************************

	endcriterion = 1;										// initially, no voltage limit is exceeded

	// Try getting the cell voltage
	try{
		c.getVoltage(verbose >= printNonCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);

		// maximum and upper voltage limits
		if(v > c.getVmax())										// above the maximum cell voltage
			endcriterion = -2;
		else if (v > Vupp)										// below the maximum cell voltage but above the upper voltage limit
			endcriterion = 2;

		// minimum and lower voltage limits
		if (v < c.getVmin())									// below the minimum cell voltage
			endcriterion = -3;
		else if (v < Vlow)										// above the minimum cell voltage but below the lower voltage limit
			endcriterion = 3;

	}
	catch(int e){
		if(verbose >= printNonCrit)
			cout<<"Error in BasicCycler::setCurrent when getting the voltage of the cell after setting the current: "<<endl;

		// If we get an error (e.g. concentration out of bounds), we can't set the current
		// based on the sign of the current we can guess which voltage limit would have been exceeded
		if (I < 0)
			endcriterion = -2;								// you are charging, so you probably exceeded Vmax
		if (I > 0)
			endcriterion = -3;								// you are discharging so you probably exceeded Vmin
		else
			endcriterion = 0;								// you got an error while resting the cell, which means the cell is in an illegal condition
	}

	// **************************************************** 4 output parameters *****************************************************************************

	// If we exceeded one of the voltage limits, restore the original battery state
	if(endcriterion != 1)
		c.setStates(sini, Iini);
	// else we have correctly set the current so we leave everything like it is

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::setCurrent with current = "<<I<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<" is terminating with value "<<endcriterion<<endl;

	return endcriterion;
}

int BasicCycler::CC_t_V(double I, double dt, bool blockDegradation, double time, double Vupp, double Vlow, double* ahi, double* whi, double* timei){
	/*
	 * function to load the battery at a constant current with tree end-conditions
	 * 		load for a given amount of time
	 * 		charge until the voltage is above an upper voltage limit
	 * 		discharge until the voltage is below a lower voltage limit
	 * The first condition satisfied terminates the CC
	 * 		i.e. you reach the full time, or the voltage is above the upper limit or the voltage is below the lower limit
	 * Note that the voltage limits are only triggered if the current has the correct sign.
	 * 		So if the cell is above the upper voltage limit but you are discharging, the cell keeps discharging until the cell reaches the lower voltage limit
	 * 		and if the cell is below the lower voltage limit but you are charging, the cell keeps charging until the cell reaches the upper voltage limit
	 *
	 * IN
	 * I 			load current [A], positive = discharge, negative = charge
	 * dt 			time step to be taken in the time integration [sec]
	 * 				should be small enough to ensure numerical stability. The stability depends on the cell current and the temperature
	 * 				the order of magnitude is 1 to 5 seconds.
	 * 				If the data-collection time interval is smaller than dt, the data-collection time interval is used instead of the specified value
	 * 				i.e. dt = min(dt,CyclingDataTimeInterval)
	 * blockDegradation if true, degradation is not accounted for during this CC
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 * 				if false, the cell is degraded during the CC phase.
	 * time 		total time for which this current should be applied [sec]
	 *					must be a multiple of dt
	 * Vupp			upper voltage limit, must be in the voltage limit of the cell (0 Cell.Vmin <= Vupp <= Cell.Vmax) [V]
	 * Vlow			lower voltage limit, must be in the voltage limit of the cell (0 Cell.Vmin <= Vlow <= Cell.Vmax) [V]
	 *
	 * OUT
	 * ahi 			the total discharged capacity [Ah]
	 * whi 			the total discharged energy [Wh]
	 * timei		the total time the cell has been loaded with the CC [sec]
	 * int 			which end condition was reached
	 * 					-3 	the minimum cell voltage was exceeded
	 * 					-2 	the maximum cell voltage was exceeded
	 * 					0 	an error occurred
	 * 					1 	the full time was completed
	 * 					2 	the upper voltage limit was reached while charging the cell
	 * 					3 	the lower voltage limit was reached while discharging the cell
	 *
	 * THROWS
	 * 1003 		the total time is not a multiple of dt
	 * 1004 		illegal voltage or current input
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_t_V with time = "<<time<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<", and current "<<I<<" is starting"<<endl;

	// Check that the total time is a multiple of the time step
	if (remainder(time, dt) > 0.01){
		cerr<<"Error in BasicCycler::CC_t_V. The total time "<<time<<" is not a multiple of the time step dt "<<dt<<endl<<flush;
		throw 1003;
	}

	// Check that the voltage limits are within the cell voltage limits
	bool uvmax = Vupp > c.getVmax(); 				// check if the upper voltage is below the cell maximum voltage
	if (uvmax)
		cerr<<"Error in BasicCycler::CC_t_V. The upper voltage "<<Vupp<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool uvmin = Vupp < c.getVmin();				// check if the upper voltage is above the cell minimum voltage
	if (uvmin)
		cerr<<"Error in BasicCycler::CC_t_V. The upper voltage "<<Vupp<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	bool lvmax = Vlow > c.getVmax(); 				// check if the lower voltage is below the cell maximum voltage
	if (lvmax)
		cerr<<"Error in BasicCycler::CC_t_V. The lower voltage "<<Vlow<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool lvmin = Vlow < c.getVmin();				// check if the lower voltage is above the cell minimum voltage
	if (lvmin)
		cerr<<"Error in BasicCycler::CC_t_V. The lower voltage "<<Vlow<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (uvmax || uvmin || lvmax || lvmin)
		throw 1004;


	// *********************************************************** 1 variables & settings ***********************************************************************

	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_t_V is making the variables"<<endl;

	// ensure the time steps is smaller than the data collection time resolution (if we are collecting data)
	double feedb = CyclingDataTimeInterval;
	if(feedb > 0)
		dt = min(dt, feedb);

	// check that the total time is still multiple of the time step, if not set the time step to 1sec (which always works)
	if (remainder(time, dt) > 0.01)
		dt = 1;

	// number of time steps
	int ttot = time / dt; 							// number of time steps needed in total
	int nstore = CyclingDataTimeInterval/dt;		// number of time steps between two data collection points

	// store initial states to restore them if needed
	State s, s2;
	double Iini, Iprev;
	c.getStates(s, &Iini);

	// variables
	double v, ocvp, ocvn, etap, etan, rdrop, tem;
	double ah = 0;									// discharged charge [Ah]
	double wh = 0;									// discharged energy [Wh]
	double tt = 0;									// time on load
	int ti = 0;										// number of time steps taken
	int endcriterion = 99;							// integer indicating why the function terminated
	bool er = false;								// indicated that an error occurred

	// *********************************************************** 2 set the current ****************************************************************************

	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_t_V is trying to set the current"<<endl;

	// store the initial battery state if we are storing data
	// If the cell is in an invalid condition, this will produce an error
	if(CyclingDataTimeInterval > 0){
		try{
			c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"error in BasicCycler::CC_t_V when getting the initial cell voltage. This means the cell is in an illegal state when this function is called: "<<e<<". Throwing on the error"<<endl;
			throw e;
		}
		timeRes += 0.000001;						// add a small amount to the rest time to ensure the new data point is different from the point before
		storeResults(c.getI(), v, ocvp, ocvn, tem); // store the initial data point
	}

	// try setting the current, this throws an error if the current has the wrong sign (1004) or if the battery is in an illegal state
	int end1;										// integer indicating if the current could sucessfuly be applied
	try{
		end1 = setCurrent(I, Vupp, Vlow);			// this underlying function tries setting the current and reports back if a voltage limit was reached
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"error in BasicCycler::CC_t_V when setting the cell current: "<<e<<". Throwing on the error"<<endl;
		throw e;
	}

	// If the current could not be set without violating the voltage limits, stop now
	if (end1 != 1){

		if(verbose >= printCyclerFunctions)
			cout<<"BasicCycler::CC_t_V with time = "<<time<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<", and current "<<I<<" is terminating because the current couldn't be set without hitting the voltage limitations: "<<end1<<endl;

		// restore the original battery state
		c.setStates(s, Iini);

		// Make the output parameters: we haven't been able to have any throughput
		*ahi = 0;
		*whi = 0;
		*timei = 0;

		// report back which voltage limit was hit
		return end1;
	}

	// store the battery state (after setting the current) if we are storing data
	if(CyclingDataTimeInterval > 0){
		c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
		timeRes += 0.0001;							// add a small amount to the rest time to ensure the new data point is different from the point before
		storeResults(c.getI(), v, ocvp, ocvn, tem); // store the data point
	}

	// *********************************************************** 3 loop for the total time ******************************************************************

	for (int t=0;t<ttot;t++){

		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::CC_t_V is applying a current of "<<I<<"A in iteration "<<t<<" with so far "<<tt<<" seconds done, cell voltage "<<v<<endl;

		// get the battery state to restore it if needed
		c.getStates(s2, &Iprev);

		// try to follow the current
		try{
			// follow the current
			c.ETI(verbose >= printCrit,dt,blockDegradation);

			// Get the voltage after the time step
			c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
		}
		// an error occurred while following the current or getting the voltage, therefore the current can't be sustained
		catch(int err){
			if(verbose >= printCrit)
				cout<<"Error in BasicCycler::CC_t_V while cycling, error "<<err<<" in time step "<<t<<" and the last voltage was "<<v<<endl;
			er = true;								// indicate the error happened
			v = 0;
			endcriterion = 0;
		}

		// Check the end criteria. Undo the last iteration if a limit is exceeded such that we stay within the limits at all times
		if(er){										// an error occurred
			c.setStates(s2, Iprev);
			endcriterion = 0;
			break;
		}
		else if(I < 0 && v > Vupp) {				// we are charging so stop if the voltage is above the upper limit
			c.setStates(s2, Iprev);
			endcriterion = 2;
			break;
		}
		else if(I > 0 && v < Vlow){					// we are discharging so stop if the voltage is below the limit
			c.setStates(s2, Iprev);
			endcriterion = 3;
			break;
		}
		else if (v > c.getVmax()){					// the cell has exceeded the maximum voltage limit
			c.setStates(s2, Iprev);
			endcriterion = -2;
			break;
		}
		else if (v < c.getVmin()){					// the cell has exceeded the minimum voltage limit
			c.setStates(s2, Iprev);
			endcriterion = -3;
			break;
		}
		else{										// this is a valid iteration and we can to store the results
			// update the time and discharged capacity and energy
			tt += dt;
			ah += I*dt/3600.0;
			wh += I * v *dt/3600.0;
			ti ++;

			// Update the data collection variables
			if (I > 0){								// discharging
				timeDis += dt;
				AhDis += abs(I)*dt/3600.0;
				WhDis += abs(I) * v *dt/3600.0;
			}
			else if (I < 0){						// charging
				timeCha += dt;
				AhCha += abs(I)*dt/3600.0;
				WhCha += abs(I) * v *dt/3600.0;
			}
			else{									// resting
				timeRes += dt;
			}

			// store the results at the specified time resolution
			if(remainder(t, nstore) == 0){
				if(verbose >= printCyclerDetail)
					cout<<"BasicCycler::CC_t_V is storing cycling data in time step "<<t<<endl;
				storeResults(I, v, ocvp, ocvn, tem);
			}
		} // end data storage for this time step

	} // end loop with time steps

	// Check if the loop ended because we have reached the time limit
	if(ti == ttot)
		endcriterion = 1;

	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_t_V has finished applying the current with end criterion "<<endcriterion<<endl;

	// *********************************************************** 4 output parameters ***********************************************************************

	// store the last data point if we are storing data
	if(CyclingDataTimeInterval > 0){
		try{
			c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in BasicCycler::CC_t_V when getting the voltage at the end, error "<<e<<". Return 0 "<<endl;
			endcriterion = 0;
		}
		storeResults(I, v, ocvp, ocvn, tem);
	}

	// Make the output parameters
	*ahi = ah;
	*whi = wh;
	*timei = tt;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_t_V with time = "<<time<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<", and current "<<I<<" is terminating with "<<endcriterion<<endl;

	return endcriterion;
}

int BasicCycler::CC_t(double I, double dt, bool blockDegradation, double time, double* ahi, double* whi, double* timei){
	/*
	 * function to load the battery at a constant current for a given amount of time.
	 *
	 * IN
	 * I 			load current [A], positive = discharge, negative = charge
	 * dt 			time step to be taken in the time integration [sec]
	 * 				should be small enough to ensure numerical stability. The stability depends on the current rating and the temperature
	 * 				the order of magnitude is 1 to 5 seconds.
	 * 				If the data-collection time interval is smaller than dt, the data-collection time interval is used instead of the specified value
	 * 				i.e. dt = min(dt,CyclingDataTimeInterval)
	 * blockDegradation if true, degradation is not accounted for during this CC
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 * time 		total time for which this current should be applied [sec]
	 *					must be a multiple of dt
	 *
	 * OUT
	 * ahi 			the total discharged capacity [Ah]
	 * whi 			the total discharged energy [Wh]
	 * timei		the total time it took [sec]
	 * int 			integer indicating why the CC phase stopped
	 * 					-3 	the minimum cell voltage was exceeded
	 * 					-2 	the maximum cell voltage was exceeded
	 * 					0 	an error occurred
	 * 					1 	the full time was completed
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_t with time = "<<time<<", and current "<<I<<" is starting"<<endl;

	// variables
	double ah = 0;									// discharged charge [Ah]
	double wh = 0;									// discharged energy [Wh]
	double tt = 0;									// time on load
	int endcr;										// integer indicating why the underlying function ended

	// Set the voltage limits wide enough so this is never a problem
	double Vupp = c.getVmax();
	double Vlow = c.getVmin();

	// Call CC_t_V with very wide voltage limits such that we always end because of the time-restrictions
	try{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::CC_t is calling CC_t_V with a current of "<<I<<", upper voltage limit "<<Vupp<<", lower voltage limit "<<Vlow<<" and time "<<time<<endl;
		endcr = CC_t_V(I, dt, blockDegradation, time, Vupp, Vlow, &ah, &wh, &tt);
	}
	catch(int err){
		if(verbose >= printCrit)
			cout<<"Error in BasicCycler::CC_t while cycling in the underlying function, error "<<err<<", throwing it on"<<endl;
		throw err;
	}

	// Make the output parameters
	*ahi = ah;
	*whi = wh;
	*timei = tt;
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_t is checking the reason why CC_t_V finished, which was "<<endcr<<endl;

	// double check that we have indeed completed the full time if the end criterion says so
	if(endcr == 1)
		assert(abs(tt - time) < 1.0);

	// there were no upper and lower voltage limits, so remove these end conditions and replace them with the cell voltage limits
	if (endcr == 2)
		endcr = -2;									// CC stopped because the maximum cell voltage was exceeded
	if (endcr == 3)
		endcr = -3;									// CC stopped because the minimum cell voltage was exceeded

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_t with time = "<<time<<", and current "<<I<<" is terminating with "<<endcr<<endl;

	return endcr;									// return why the CC phase stopped
}

int BasicCycler::CC_V(double I, double dt, bool blockDegradation, double Vset, double* ahi, double* whi, double* timei){
	/*
	 * function to load the battery at a constant current until a given voltage is reached.
	 * The function always slightly undershoots the voltage (i.e. it stops just before it hits the voltage rather than exceeding the voltage).
	 * This is to avoid illegal battery states (e.g. if you charge to Vlim = c.getVmax() and you overshoot, the battery state is illegal)
	 *
	 * IN
	 * I 			load current [A], positive = discharge, negative = charge
	 * 				must be compatible with the voltage limits.
	 * 					i.e. if the voltage limit is above the cell voltage, the current must be negative (charging)
	 * 					and vice versa (discharging if the cell voltage is above the voltage limit)
	 * dt 			time step to be taken in the time integration [sec]
	 * 				should be small enough to ensure numerical stability. The stability depends on the current rating and the temperature
	 * 				the order of magnitude is 1 to 5 seconds
	 * 				If the data-collection time interval is smaller than dt, the data-collection time interval is used instead of the specified value
	 * 				i.e. dt = min(dt,CyclingDataTimeInterval)
	 * blockDegradation if true, degradation is not accounted for during this CC
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 * Vset 		the voltage to which you want to (dis)charge the cell, Cell.Vmin <= Vset <= Cell.Vmax [V]
	 *
	 * OUT
	 * ahi 			the total discharged capacity [Ah]
	 * whi 			the total discharged energy [Wh]
	 * timei		the total time it took [sec]
	 * int 			integer indicating why the CC phase ended
	 * 					-3 	the minimum cell voltage was exceeded
	 * 					-2 	the maximum cell voltage was exceeded
	 * 					0 	an error occurred
	 * 					2 	the upper voltage limit was reached while charging the cell
	 * 					3 	the lower voltage limit was reached while discharging the cell
	 *
	 * THROWS
	 * 1004 		illegal voltage or current input
	 * 1016 		CC_t_V which did the CC phase ended because the time limit was reached, instead of the voltage limit.
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_V with voltage = "<<Vset<<", and current "<<I<<" is starting"<<endl;

	// *********************************************************** 1 variables & checks on input parameters ***********************************************************************

	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_V is making the variables and checking the voltage limits"<<endl;

	// variables
	double ah = 0;												// discharged charge [Ah]
	double wh = 0;												// discharged energy [Wh]
	double tt = 0;												// time on load
	double Vupp, Vlow;											// upper and lower voltage limits of the underlying function
	int endcr;													// integer indicating why the CC phase ended
	bool check = true;											// check if the battery state is valid after setting a current, throw an error if not
	double v;													// cell voltage
	double ocvp, ocvn, etap, etan, rdrop, tem;					// unneeded feedback variables

	// check if the voltage is allowed
	bool vmax = Vset > c.getVmax(); 							// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in BasicCycler::CC_V. The voltage "<<Vset<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vset < c.getVmin();								// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in BasicCycler::CC_V. The voltage "<<Vset<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (vmax || vmin)
		throw 1004;

	// Check if the current has the correct sign
	// get the OCV and if Vset > OCV then we must charge so I must be negative and vice versa
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_V is checking that the current has the correct sign"<<endl;
	try{
		c.setI(verbose >= printNonCrit, check, 0);				// set the current to 0
		c.getVoltage(verbose >= printNonCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
	}
	catch(int e){
		if(verbose >= printNonCrit)
			cout<<"Error in BasicCycler::CC_V when getting the OCV of the cell: "<<e<<". Throwing it on"<<endl;
		throw e;
	}
	bool sign = !( (v < Vset && I < 0) || (v > Vset && I > 0)); // check I has the correct sign: charge -> I<0 or discharge -> I > 0
	if(sign){
		if(verbose >= printCrit)
			cerr<<"Error in BasicCycler::CC_V. The current has the wrong sign. The cell voltage is "<<v<<" and you want to get to "<<Vset
			<<" but the current is "<<I<<". A current < 0 means the cell is charging, > 0 means it is discharging."<<endl<<flush;
		throw 1004;
	}

	// *********************************************************** 2 call CC_t_V to get the cell to the specified limit ***********************************************************************

	// Set the upper and lower voltage limit
	if (I > 0){
		// discharge the cell so set the lower voltage limit to the set voltage
		Vlow = Vset;
		Vupp = c.getVmax();
	}
	else{
		// charge the cell so set the upper voltage limit to the set voltage
		Vlow = c.getVmin();
		Vupp = Vset;
	}
	double time = 99999999;										// set the time limit very high so this is never the problem

	// Call CC_t_V with very long time limit such that we always end because of the voltage-restrictions
	try{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::CC_V is calling CC_t_V with a current of "<<I<<", upper voltage limit "<<Vupp<<", lower voltage limit "<<Vlow<<" and time "<<time<<endl;
		endcr = CC_t_V(I, dt, blockDegradation, time, Vupp, Vlow, &ah, &wh, &tt);
	}
	catch(int err){
		if(verbose >= printCrit)
			cout<<"Error in BasicCycler::CC_V while cycling in the underlying function, error "<<err<<", throwing it on"<<endl;
		throw err;
	}

	// Make the output parameters
	*ahi = ah;
	*whi = wh;
	*timei = tt;

	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_V is checking the reason why CC_t_V finished, which was "<<endcr<<endl;

	// ensure the CC phase didn't end because of the time limit
	if(endcr == 1){
		cerr<<"Error in BasicCycler::CC_V, the underlying function CC_t_V terminated because the time limit of the CC phase was reached instead of the voltage limit. Throwing an error"<<endl<<flush;
		throw 1016;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_V with voltage = "<<Vset<<", and current "<<I<<" is terminating with "<<endcr<<endl;

	return endcr;
}

void BasicCycler::CC_halfCell_full(double I, double dt, bool pos, int nin, double OCVi[], double* ahi, int* n){
	/*
	 * This function cycles one electrode with a constant current.
	 * It does a CC charge/discharge until the extreme concentration (0 or Cmax) of that electrode is reached.
	 * The other electrode is left as it was.
	 * The function ignores degradation while doing this.
	 * No cycling data is stored during this function (because the cell is in an illegal condition)
	 * Instead, the OCV of the electrode you are cycling is returned as an array
	 *
	 * USE THIS FUNCTION WITH CAUTION.
	 * It leaves the battery in an illegal state because only one electrode is cycled.
	 * i.e. it is as if the battery were cut in two and one electrode is cycled against a metallic li-electrode.
	 * The other electrode is left at its initial state
	 *
	 * IN
	 * I 			load current [A], positive = discharge, negative = charge
	 * dt 			time step to be taken in the time integration
	 * pos 			if true, the positive electrode is cycled
	 * 				if false, the negative electrode is cycled
	 * nin 			length of the arrays which the user supplied to store the output, i.e. length(OCVi)
	 *
	 * OUT
	 * OCVi			array in which the open circuit potential of the electrode to be considered will be written
	 * ahi 			the total discharged capacity [Ah], positive for discharge, negative for charge
	 * n 			number of data points written in the arrays, i.e. OCVi[0] until OCVi[n-1] contains the voltage
	 *
	 * throws
	 * 1002			the arrays provided to store the input are too short
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_halfCell_full is starting"<<endl;

	// variables
	State s;									// store initial states to restore them if needed
	double Iprev;								// cell current in the previous time step
	bool blockDegradation = true;				// bool to indicate we want to ignore degradation
	bool check = false;							// don't check if the battery is in a valid state after changing the current because we know it won't be
	double  V, OCVpi, OCVni, etapi, etani, Rdropi, tempi; // variables for unneeded feedback
	bool limit = true;							// boolean to indicate if we should continue cycling
	int nt = 0;									// number of steps until the extreme concentration is reached
	double ah = 0;								// total discharged capacity [Ah]

	// Get the initial cell voltage
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_halfCell_full is getting the initial cell voltage"<<endl;
	try{
		c.setI(verbose >= printCrit, check, I);			// set the current
		c.getVoltage(verbose >= printCrit, &V, &OCVpi, &OCVni, &etapi, &etani, &Rdropi, &tempi);
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in BasicCycler::CC_halfCell_full when setting the cell current "<<e<<". Throw it on"<<endl;
		throw e;
	}
	if (pos)									// we are cycling the positive electrode, store OCVpos
		OCVi[0] = OCVpi;
	else										// we are cycling the negative electrode, store OCVneg
		OCVi[0] = OCVni;
	nt ++;										// We have stored one data point

	// Cycle while you have not reached the extreme surface concentration
	while (limit) {
		c.getStates(s, &Iprev);					// battery state before this iteration

		try {
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::CC_halfCell_full is cycling electrode "<<pos<<" with cathode potential "<<OCVpi<<" and anode potential "<<OCVni<<endl;

			c.ETI_electr(verbose >= printNonCrit, I, dt, blockDegradation, pos); // cycle only one electrode
			c.getVoltage(verbose >= printNonCrit, &V, &OCVpi, &OCVni, &etapi, &etani, &Rdropi, &tempi); // get the voltage

			// check that the arrays provided for output are long enough to store the additional point
			if (nin <= nt){
				if(verbose >= printCrit)
					cerr<<"ERROR in BasicCycler::CC_halfCell_full the array provided for output has a length of "<<nin<<" but we need a length of at least "<<nt+1<<endl<<flush;
				throw 1002;
			}

			// Store the data point
			if (pos)							// store the OCV of the electrode you are interested in
				OCVi[nt] = OCVpi;
			else
				OCVi[nt] = OCVni;
			ah += I * dt/3600.0;
			nt ++;								// a new data point was stored
		}
		catch (int er) {						// getVoltage throws an error if the min or max concentration has been reached
			if (verbose >= printNonCrit)		// so this is a non critical error (we can expect it to happen because we want to (dis)charge an electrode to its min or max concentration
				cout<<"Error in BasicCycler::CC_halfCell_full when OCVp = "<<OCVpi<<" and OCVn = "<<OCVni<<endl;
			limit = false;						// indicate we want to stop cycling (the loop-conditions will be false so we stop (dis)charging)
			c.setStates(s, Iprev);				// undo the last iteration to end up in a valid battery state

			// if the error is 1002, the error was that the arrays for output are too short, rather than that we have reached the min or max concentration
			if(er == 1002){
				if(verbose >= printCrit)				// that is not supposed to happen so print this critical error message and throw the error on
					cout<<"Error in BasicCycler::CC_halfCell_full, the arrays were too short. throwing it on."<<endl;
				throw er;
			}
		}
	}

	// Make the output variables
	*n = nt;
	*ahi = ah;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_halfCell_full is terminating"<<endl;
}

void BasicCycler::findCVcurrent_recursive(double Imin, double Imax, int sign, double Vset, double dt, bool blockDegradation, double* Il, double* Vl){
	/*
	 * This is a recursive function which solves the nonlinear equation to find the current needed to reach a given voltage after a given time step.
	 *
	 * The current will produce a voltage which has a relative error below 0.1%, i.e. (Vl - Vset)/Vset < 0.001
	 * Additionally, it will always be an 'undershoot' (we always stop 'too soon'):
	 * 		if you are charging, Vl will be below Vset (by max 0.1%)
	 * 		if you are discharging, Vl will be above Vset (by max 0.1%)
	 * 		this is to avoid that you violate the voltage limits of the cell (e.g. if you do a CV at Vmax, an overshoot would produce V > Vmax), which would result in an error
	 * If no such current can be found (e.g. the step size becomes too small), the current returned is a best-guess.
	 * 		in most cases, this current will undershoot the set voltage but with a relative error above 0.1%
	 *
	 * The search algorithm recursively refines the lower and upper bound on the current which reaches the set voltage.
	 * Initially, this function must be called with some bounds.
	 * the algorithm then takes 11 discrete steps between these bounds and searches when the error changes sign
	 * 		i.e. step i gives a positive error while step i+1 gives a negative error or vice versa. then we know the sought current is between those two
	 * 		these two values will be the new bounds for the lower-level recursion (which will have a smaller range and therefore a smaller step size)
	 * 		this process continues until a current is found which satisfies the end-conditions
	 * In case none of the currents in the search space is small enough (no undershoot), the lower bound is halved in the next recursion level
	 * in case none of the currents in the search space is large enough (no overshoot), the upper bound is doubled in the next recursion level
	 * If the step size becomes too small, the algorithm is terminated with the best-guess so far in order to avoid a slow convergence or an an infinite loop.
	 * 		this best-guess will still produce an undershoot but with a larger relative error
	 *
	 * Note: function does not alter the battery state, i.e. at the end of the function, the original battery state is restored.
	 * It merely searches for the required current.
	 *
	 * IN
	 * Imin 	the lower bound on the absolute value of the current needed to reach the voltage [A], > 0
	 * Imax 	the upper bound on the absolute value of the current needed to reach the voltage [A], > 0
	 * sign 	the sign of the current (i.e. do we need to charge or discharge to reach the voltage), 1 for discharge, -1 for charge
	 * Vset 	the voltage we are trying to reach, [V]
	 * dt 		the time step after which we need to reach the voltage
	 * blockDegradation if true, degradation is not accounted for during this CV
	 * 			set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 *
	 * OUT
	 * Il 		the current which is necessary to reach this voltage after the given time step
	 *  			this current results in a voltage which has a relative error below 0.1% and is in undershoot (above Vset if discharging; below Vset if charging)
	 * 				in case on good current could be found (e.g. an error), Il is a best-guess is returned which will always undershoot the target voltage (but with a larger error)
	 * Vl 		voltage you can expect to get after applying Il for the specified time step
	 *
	 * THROWS
	 * 1005		the voltage Vset has an illegal value
	 * 1006 	the bounds have an illegal value: negative or the upper bound is below the lower bound (this should never happen)
	 * 1007 	the step size in the search algorithm has become too small before the end conditions have been met (Il becomes the best-guess)
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::findCVcurrent_recursive is starting with set voltage "<<Vset<<" and range "<<Imin<<" to "<<Imax<<endl;

	// check if the voltage is valid
	bool vmax = Vset > c.getVmax(); 							// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in BasicCycler::findCVcurrent_recursive. The voltage "<<Vset<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vset < c.getVmin();								// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in BasicCycler::findCVcurrent_recursive. The voltage "<<Vset<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (vmax || vmin)
		throw 1005;

	// check if the bounds on the current have legal values
	bool imin = Imin < 0;
	if (imin)
		cerr<<"ERROR in BasicCycler::findCVcurrent_recursive. The lower bound "<<Imin<<" is negative but it has to be positive. Throwing an error"<<endl;
	bool imax = Imax <= 0;
	if (imax)
		cerr<<"ERROR in BasicCycler::findCVcurrent_recursive. The upper bound "<<Imax<<" is negative but it has to be strictly positive. Throwing an error"<<endl;
	bool lim = Imax <= Imin;
	if(lim)
		cerr<<"ERROR in BasicCycler::findCVcurrent_recursive. The upper bound "<<Imax<<" is below the lower bound "<<Imin<<". Throwing an error"<<endl;
	if (imin || imax || lim)
		throw 1006;

	// *********************************************************** 1 variables & settings ***********************************************************************
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::findCVcurrent_recursive with set voltage "<<Vset<<" and range "<<Imin<<" to "<<Imax<<" is initialising"<<endl;

	// store the initial battery states to restore them at the end
	State s;
	double Iini;
	c.getStates(s, &Iini);

	// variables
	double ocvp, ocvn, etap, etan, rdrop, tem;					// unneeded feedback variables
	bool check = false;											// we don't want to check if the state is valid after setting a current, because we know that the search algorithm will occasionally exceed the voltage limit
	int nstep = 10;												// number of steps in the search space for the current
	double dI = sign*(Imax - Imin)/nstep;						// step size for the current [A]
	double Itest = sign*Imin;									// current in this iteration for the search for Il [A]
	double Vtest;												// voltage in this iteration for the search for Il
	double err;													// relative error between the actual cell voltage and the desired cell voltage [-]
	double iminnew = -1;										// new estimate for the lower bound on the absolute value of Il [A]
	double imaxnew = -1;										// new estimate for the upper bound on the absolute value of Il [A]
	bool found = false;											// boolean indicating if the current has been found

	// check that the step size is large enough in order to avoid an infinite loop or very slow convergence
	if (abs(dI) < pow(10,-6)){

		// our best guess is the lower bound, which should be an undershoot
		c.setI(verbose >= printCrit, check, Imin);
		c.ETI(verbose >= printCrit,dt,blockDegradation);
		c.getVoltage(verbose >= printCrit, &Vtest, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
		*Il = Imin;
		*Vl = Vtest;

		// print a warning and throw the error
		if(verbose >= printNonCrit)
			cerr<<"ERROR in BasicCycler::findCVcurrent_recursive. The current step "<<dI<<" is too small Throwing an error and returning best guess Il = "<<Imin<<" which produces voltage "<<Vtest<<" instead of "<<Vset<<endl;
		throw 1007;
	}

	// *********************************************************** 3 search for the required current ****************************************************************************
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::findCVcurrent_recursive is starting to scan the search space"<<endl;

	// loop through the search space, starting at the lower bound and gradually increasing the (absolute value of) the test current
	for(int i=0;i<nstep+1;i++){

		// restore the original battery state
		c.setStates(s, Iini);

		// calculate effect of the test current in this iteration
		try{
			c.setI(verbose >= printCrit, check, Itest);			// set the current
			c.ETI(verbose >= printCrit,dt,blockDegradation);	// take the time step
			c.getVoltage(verbose >= printCrit, &Vtest, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);// check the voltage

		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in BasicCycler::findCVcurrent_recursive. An error occurred when getting the effect of Itest "<<Itest<<"A, error: "<<e<<". The  voltage was "<<Vtest<<". Throwing the error on"<<endl;
			throw e;
		}

		// Calculate the error between the actual and the desired voltage
		err = sign*(Vtest-Vset)/Vset; 							//relative error, positive for undershoot, negative for overshoot

		if(verbose >= printfindCVcurrentDetail)
			cout<<"BasicCycler::findCVcurrent_recursive with set voltage "<<Vset<<" and range "<<Imin<<" to "<<Imax<<" has test current "<<Itest<<" which produces test voltage "<<Vtest<<", which is a relative error of "<<err*100<<"%"<<endl;

		// check if the current is good enough
		if (0<err && err < 0.001){								// the error in undershoot and it is below 1%
			*Il = Itest;
			*Vl = Vtest;
			found = true;										// indicate we have found the required current
			break;												// leave the loop, we don't need to check the other steps in the search space
		}

		// update the estimated lower and upper bounds
		if(err >= 0)											// we are undershooting the voltage, so this is a lower bound of the current
			iminnew = abs(Itest);								// we are increasing the (absolute value of) the current, so this will be updated with the largest bound found so far
		else{													// we are overshooting the voltage, so this is an upper bound of the current
			imaxnew = abs(Itest);								// we started at the lower bound so we know that the first time there is an overshoot, this will be the smallest upper bound
			break;												// we can stop (later iterations will just give larger upper bounds)
		}

		// update the test current
		Itest += dI;
	}

	// restore the original battery state
	c.setStates(s, Iini);

	// ****************************************************** 4 update the bounds & recurse to convergence ****************************************************************************

	// if one of the bounds hasn't been reached, set it to a 50% wider range then before
	if (iminnew == -1)											// the lower bound still has its original value of -1, so no current in the search space was an undershoot
		iminnew = Imin * 0.5;									// half the lower bound for the next level
	if (imaxnew == -1)											// the upper bound still has its original value of -1, so no current in the search space was an overshoot
		imaxnew = Imax * 2;										// double the upper bound for the next level

	// if we haven't found the current yet, recursively call the function again with the new bounds
	// 		if the desired current was between the lower and upper level, we call the function again with more narrow limits (higher lower limit, lower upper limit)
	// 		if the desired current was not between the lower and upper level, we call the function again with wider limits (lower lower limig or higher upper limit)
	// 		this converges around the desired current
	if (!found){
		try{
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::findCVcurrent_recursive with set voltage "<<Vset<<" and range "<<Imin<<" to "<<Imax<<" is recursively call itself with updated range "<<iminnew<<" to "<<imaxnew<<endl;

			findCVcurrent_recursive(iminnew, imaxnew, sign, Vset, dt, blockDegradation, Il, Vl);
		}
		catch(int e){
			if(verbose >= printNonCrit)
				cout<<"Error in findCVcurrent_recursive in a lower recursive level. This level had range "<<Imin<<" to "<<Imax<<" and set voltage "<<Vset<<". Terminating with Il = "<<*Il<<" and Vl "<<*Vl<<" and throwing the error on"<<endl;
			throw e;
		}
	}

	// the lowest level (the one in which the desired current was found), will return the output parameters to the higher level, which will in turn return them to the top level
	// so we can simply end the function, the output parameters have been set in the lowest level recursion.

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::findCVcurrent_recursive with set voltage "<<Vset<<" and range "<<Imin<<" to "<<Imax<<" is terminating with current "<<*Il<<" which should produce a voltage of "<<*Vl<<endl;
}

void BasicCycler::findCVcurrent(double Vset, double dt, bool blockDegradation, double* Il, double* Vl){
	/*
	 * Function to find the current necessary to get the cell voltage to Vset in the given time step.
	 * It is used to do a CV (dis)charge for a cell.
	 * The state space model takes the current as input and has the voltage as output.
	 * To keep the voltage constant at Vset, we need to solve the nonlinear equation: Vset = batter_model(Il)
	 *
	 * The accuracy of the function is such that the relative error on the voltage will be below 0.1%.
	 * The error is always on the side of the OCV of the cell ('undershoot').
	 * This means that on charge, Il will give a voltage V < Vset by max 0.1%.
	 * 					  discharge, Il will give a voltage V > Vset by max 0.1%
	 * If an error occurs and the search algorithm can't find a current which satisfies these conditions, a best-guess is returned.
	 * this best-guess current will still undershoot the set voltage (to avoid that you violate the voltage limits of the cell if you follow this current) but the error might be larger than 0.1%
	 *
	 * The search algorithm itself is implemented in findCVcurrent_recursive.
	 * See that function for a description on how it works.
	 *
	 * Note: function does not alter the battery state, i.e. at the end of the function, the original battery state is restored.
	 * It merely searches for the required current.
	 *
	 * IN
	 * Vset 	the voltage which the cell should achieve (+- 0.1%)
	 * 				Cell.Vmin <= V =< Cell.Vmax
	 * dt 		the time step of the time integration [s]
	 * 				dt should be a few seconds, e.g. 1 to 5 to ensure numerical stability
	 * blockDegradation if true, degradation is not accounted for during this CV
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 *
	 * OUT
	 * Il 		the current which is necessary to reach this voltage after the given time step
	 * Vl 		voltage you can expect to get after applying Il for the specified time step
	 *
	 * THROWS
	 * 1005		the voltage Vset has an illegal value
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::findCVcurrent is starting with set voltage "<<Vset<<endl;

	// check the input voltage
	bool vmax = Vset > c.getVmax(); 					// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in BasicCycler::findCVcurrent. The voltage "<<Vset<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vset < c.getVmin();						// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in BasicCycler::findCVcurrent. The voltage "<<Vset<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (vmax || vmin)
		throw 1005;

	// *********************************************************** 1 variables & settings ***********************************************************************
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::findCVcurrent is making the variables"<<endl;

	// store the initial battery states to restore them at the end
	State s;
	double Iini;
	c.getStates(s, &Iini);

	// variables
	double Itest;										// test current for the search for Il
	double Vtest;										// test votlage for the search for Il
	double ocvp, ocvn, etap, etan, rdrop, tem;			// unneeded feedback variables

	// Set the initial boundaries on the current (i.e. we think Imin < |Il| < Imax)
	double Imin;										// lower bound on the absolute value of the current needed to reach the set voltage
	double Imax;										// upper bound on the absolute value of the current needed to reach the set voltage
	if (c.getI() != 0){									// if there is a current in the cell, assume this current is a relatively good approximation
		Imin = abs(c.getI())*0.95;						// start the search space with 5% around the cell current, to ensure a fast convergence if the cell current is a good approximation
		Imax = abs(c.getI())*1.05;						// 		if Il is not in these bounds, the recursive algorithm will enlarge them  until Il is between them
	}
	else{												// if there is no current in the cell, assume 0 < |Il| < 1C
		Imin = 0;
		Imax = c.getNominalCap();
	}

	// *********************************************************** 2 initialisation ****************************************************************************
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::findCVcurrent is initialising"<<endl;

	// check what would happen if the current is set to 0 for one time step
	try{
		c.setI(verbose >= printCrit, true, 0);			// set the current to 0 and check if the resulting voltage is valid (which should always be the case)
		c.ETI(verbose >= printCrit,dt,blockDegradation);
		c.getVoltage(verbose >= printCrit, &Vtest, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);

		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::findCVcurrent the cell has an OCV of "<<Vtest<<", which is an error of "<<(Vtest-Vset)/Vset <<endl;
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in BasicCycler::findCVcurrent. An error occurred when getting the OCV at the start: "<<e<<". Throwing the error on"<<endl;
		throw e;
	}

	// determine whether we need to charge or discharge to reach the required voltage
	int sign;
	if (Vtest > Vset)									// discharge to get to Vset -> positive current
		sign = 1;
	else
		sign = -1; 										// charge to get to Vset -> negative current


	// ******************************************************* 3 search for the current ****************************************************************************

	// check if a current of 0 gets us to the required voltage
	double err = sign*(Vtest-Vset)/Vset; 				//relative error on the voltage, positive for undershoot, negative for overshoot
	if (0<err && err < 0.001){							// the error is an undershoot and it is below 0.1%
		Itest = 0;

		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::findCVcurrent has found that a current of 0 reaches the specified voltage"<<endl;
	}
	else{												// if not, call the recursive search algorithm
		try{
			// restore the original battery state
			c.setStates(s, Iini);

			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::findCVcurrent is searching for the current to maintain the voltage"<<endl;
			findCVcurrent_recursive(Imin, Imax, sign, Vset, dt, blockDegradation, &Itest, &Vtest);
		}
		catch(int e){
			if(verbose >= printNonCrit)
				cout<<"Error in BasicCycler::findCVcurrent. An error occurred in the findCVcurrent_recursive when looking for the current "<<e<<". Throwing the error on"<<endl;

			// set the output parameter to our best-guess (which was returned by findCVcurrent_recursive even if an error happened)
			*Il = Itest;
			*Vl = Vtest;

			throw e;
		}
	}

	// *********************************************************** 4 output parameters ***********************************************************************

	// set the output parameter
	*Il = Itest;
	*Vl = Vtest;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::findCVcurrent is terminating with current "<<*Il<<" which should produce a voltage of "<<*Vl<<endl;
}

int BasicCycler::CV_t_I(double Vset, double dt, bool blockDegradation, double time, double Icut, double* ahi, double* whi, double* timei){
	/*
	 * function to load the battery at a constant voltage until the current is below a given threshold or a time limit is reached.
	 * The first end condition reached terminates the CV phase
	 *
	 * IN
	 * Vset			the voltage at which the battery should be kept [V]
	 * dt 			time step to be taken in the time integration [sec]
	 * 				should be small enough to ensure numerical stability. The stability depends on the current rating and the temperature
	 * 				the order of magnitude is 1 to 5 seconds
	 * 				If the data-collection time interval is smaller than dt, the data-collection time interval is used instead of the specified value
	 * 				i.e. dt = min(dt,CyclingDataTimeInterval)
	 * blockDegradation if true, degradation is not accounted for during this CV
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 * time 		total time for which this voltage should be applied [sec]
	 *					must be a multiple of dt
	 * Icut 		the absolute value of the lowest current allowed [A]
	 * 					if a negative value is given, the current threshold is ignored
	 *
	 * OUT
	 * ahi 			the total discharged capacity [Ah]
	 * whi 			the total discharged energy [Wh]
	 * timei		the total time it took [sec]
	 * int 			integer indicating why the CV phase ended
	 * 					1 	the time limit was reached
	 * 					2 	the current limit was reached
	 *
	 * THROWS
	 * 1003 		the total time is not a multiple of dt
	 * 1005			illegal value for the voltage
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CV_t_I with time limit = "<<time<< ", and current limit "<<Icut<<"A, and set voltage "<<Vset<<" is starting"<<endl;

	// check if the voltage limit is allowed
	bool vmax = Vset > c.getVmax(); 							// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in BasicCycler::CV_t_I. The voltage "<<Vset<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vset < c.getVmin();								// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in BasicCycler::CV_t_I. The voltage "<<Vset<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (vmax || vmin)
		throw 1005;

	// *********************************************************** 1 variables & settings ***********************************************************************
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CV_t_I is making the variables"<<endl;

	// Check that the total time is a multiple of the time step
	if (remainder(time, dt) > 0.01){
		cerr<<"Error in BasicCycler::CV_t_I. The total time "<<time<<" is not a multiple of the time step dt "<<dt<<endl<<flush;
		throw 1003;
	}

	// ensure the time steps is smaller than the data collection time resolution (if we are collecting data)
	double feedb = CyclingDataTimeInterval;
	if(feedb > 0)
		dt = min(dt, feedb);

	// check that the total time is still multiple of the time step, if not set the time step to 1sec (which always works)
	if (remainder(time, dt) > 0.01)
		dt = 1;

	// variables
	double Il;													// current in this step needed to keep the voltage constant [A]
	double Vl;													// expected voltage when applying Il [V]
	bool found = false;											// boolean to indicate if the current is smaller than the cutoff current
	bool check = true;											// check if the battery state is valid after setting a current, throw an error if not
	double ah = 0;												// cumulative discharged charge up to now [Ah]
	double wh = 0;												// cumulative discharged energy up to now [Wh]
	double tt = 0;												// cumulative time up to now [sec]
	int t = 0;													// number of time steps
	int nstore = CyclingDataTimeInterval/dt;					// number of time steps it takes to cover 'CyclingDataTimeInterval' seconds
	double v;													// cell voltage in this time step [V]
	double ocvp, ocvn, tem, etap, etan, rdrop;					// feedback variables not needed
	int endcr = 99;												// integer indicating why the CV phase terminated

	// store the initial battery state if we are storing data
	// If the cell is in an invalid condition, this will produce an error
	if(CyclingDataTimeInterval > 0){
		try{
			c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"error in BasicCycler::CV_t_I when getting the initial cell voltage. This means the cell is in an illegal state when this function is called: "<<e<<". Throwing on the error"<<endl;
			throw e;
		}
		timeRes += 0.000001;									// add a small amount to the rest time to ensure the new data point is different from the point before
		storeResults(c.getI(), v, ocvp, ocvn, tem); 			// store the initial data point
	}

	// **************************** 2 loop until the current is below the cutoff current or the maximum time is reached ********************************************************************

	while (!found){

		 // find the current we need to keep the voltage constant in the next time step
		try{
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::CV_t_I is getting the current to maintain the voltage at "<<Vset<<"V in iteration "<<t<<" with so far "<<tt<<"seconds done"<<endl;
			findCVcurrent(Vset,dt, blockDegradation, &Il, &Vl);
		}
		catch(int e){
			if(verbose >= printNonCrit)
				cout<<"Error in BasicCycler::CV_t_I when finding the current "<<Il<<", would give voltage "<<Vl<<" instead of "<<Vset<<", error "<<e<<"."<<endl;

			// try to recover from the error:
			bool valid = Vl >= c.getVmin() && Vl <= c.getVmax();// is the expected voltage in the allowed range?
			bool err = abs(Vl - Vset) < 0.01;					// is the error on the voltage below 1%?
																// note that we always undershoot the voltage so there no danger of going outside the set voltage range
			if(valid && err){
				// We can accept the current, so print a message to the user that we recovered from the error
				if(verbose >= printNonCrit)
					cout<<"BasicCycler::CV_t_I accepted the CV current because the estimated voltage is in the valid region and the error is still below 1%"<<endl;
			}
			else{
				// the current is not acceptable, write an error message and throw the error
				if(verbose >= printCrit)
					cerr<<"Error in BasicCycler::CV_t_I when finding the CV current, best guess = "<<Il<<", which would give voltage "<<Vl<<" instead of "<<Vset<<", error "<<e<<". Throwing it on"<<endl;
				throw e;
			}
		}

		// check if we've reached the cut off current
		if (abs(Il) <= Icut){									// if Icut is negative, this never happens
			found = true;										// the current is smaller than the cutoff current -> stop searching
			endcr = 2;
		}

		// else take one time step at this current (i.e. keep the voltage constant for one more time step)
		else {

			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::CV_t_I is applying a current of "<<Il<<"A in iteration "<<t<<" with so far "<<tt<<" seconds done, cell current "<<Il<<endl;

			try{
				c.setI(verbose >= printCrit, check, Il);		// set the current
				c.ETI(verbose >= printCrit,dt,blockDegradation);// apply the current for one time step
				c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem); // get the cell voltage
			}
			catch(int e){
				if(verbose >= printCrit)
					cout<<"Error in BasicCycler::CV_t_I when applying the current "<<Il<<", error "<<e<<". Throwing it on"<<endl;
				throw e;
			}

			// update the throughput parameters
			ah += Il*dt/3600.0;									// discharged charge in Ah
			wh += Il*v*dt/3600.0;								// discharged energy in Wh
			tt += dt;											// time in sec

			// Update the data collection variables
			if (Il > 0){										// discharging
				timeDis += dt;
				AhDis += abs(Il)*dt/3600.0;
				WhDis += abs(Il) * v *dt/3600.0;
			}
			else if (Il < 0){									// charging
				timeCha += dt;
				AhCha += abs(Il)*dt/3600.0;
				WhCha += abs(Il) * v *dt/3600.0;
			}
			else{												// resting
				timeRes += dt;
			}

			// store the results at the specified time resolution
			if(remainder(t, nstore) == 0)
				storeResults(Il, v, ocvp, ocvn, tem);
			t++;
		}

		// check if we have reached the time limit
		if(tt >= time){
			found = true;
			endcr = 1;
		}

	} // end loop until the current is below the threshold or the time limit is reached

	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CV_t_I with time limit = "<<time<< ", and current limit "<<Icut<<"A has stopped time stepping with current"<<Il<<"A with "<<tt<<" seconds done and is now storing data"<<endl;

	// *********************************************************** 3 output parameters ***********************************************************************

	// store what happened since the last time we stored data if we are storing data
	if(CyclingDataTimeInterval > 0){
		c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);
		storeResults(Il, v, ocvp, ocvn, tem);
	}

	// make the output variables
	*ahi = ah;
	*whi = wh;
	*timei = tt;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CV_t_I with time limit = "<<time<< ", and current limit "<<Icut<<"A, and set voltage "<<Vset<<" is terminating with "<<endcr<<endl;

	return endcr;
}

void BasicCycler::CV_t(double Vset, double dt, bool blockDegradation, double time, double* ahi, double* whi, double* timei){
	/*
	 * function to load the battery at a constant voltage for a given time
	 *
	 * IN
	 * Vset			the voltage at which the battery should be kept [V]
	 * dt 			time step to be taken in the time integration [sec]
	 * 				should be small enough to ensure numerical stability. The stability depends on the current rating and the temperature
	 * 				the order of magnitude is 1 to 5 seconds
	 * 				If the data-collection time interval is smaller than dt, the data-collection time interval is used instead of the specified value
	 * 				i.e. dt = min(dt,CyclingDataTimeInterval)
	 * blockDegradation if true, degradation is not accounted for during this CV
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 * time 		total time for which this voltage should be applied [sec]
	 *					must be a multiple of dt
	 *
	 * OUT
	 * ahi 			the total discharged capacity [Ah]
	 * whi 			the total discharged energy [Wh]
	 * timei		the total time it took [sec]
	 *
	 * THROWS
	 * 1003 		the total time is not a multiple of dt
	 * 1005 		the voltage has an illegal value
	 * 1017 		the underlying function which did the CV loading terminated because of the current limit instead of the time limit
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CV_t with time limit = "<<time<< ", and set voltage "<<Vset<<" is starting"<<endl;

	// check if the voltage is allowed
	bool vmax = Vset > c.getVmax(); 				// check if the maximum voltage is below the cell's maximum voltage
	if (vmax)
		cerr<<"Error in BasicCycler::CV_t. The voltage "<<Vset<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vset < c.getVmin();					// check if the minimum voltage is above the cell's minimum voltage
	if (vmin)
		cerr<<"Error in BasicCycler::CV_t. The voltage "<<Vset<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (vmax || vmin)
		throw 1005;

	// Set a negative current threshold so that it is always the time limit which is reached
	double Ilim = -1;

	// variables
	double ah = 0;									// discharged charge [Ah]
	double wh = 0;									// discharged energy [Wh]
	double tt = 0;									// time on load
	int endcr;										// integer indicating why the underlying function ended

	// Call CV_t_I to load the cell with the CV
	try{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::CV_t is calling CV_t_I with a set voltage of "<<Vset<<", current limit "<<Ilim<<", time limit "<<time<<endl;
		endcr = CV_t_I(Vset, dt, blockDegradation, time, Ilim, &ah, &wh, &tt);
	}
	catch(int err){
		if(verbose >= printCrit)
			cout<<"Error in BasicCycler::CV_t while cycling in the underlying function, error "<<err<<", throwing it on"<<endl;
		throw err;
	}

	// Make the output parameters
	*ahi = ah;
	*whi = wh;
	*timei = tt;

	// double check that we have indeed completed the full time
	if(endcr == 2){
		cerr<<"ERROR in BasicCycler::CV_t because the underlying function CV_t_I finished because the current threshold was reached which means we didn't complete the full time. "<<endl<<flush;
		throw 1017;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CV_t with time limit = "<<time<< ", and set voltage "<<Vset<<" is terminating"<<endl;
}

void BasicCycler::CV_I(double Vset, double dt, bool blockDegradation, double Icut, double* ahi, double* whi, double* timei){
	/*
	 * function to load the battery at a constant voltage until the current is below a given threshold
	 *
	 * IN
	 * Vset			the voltage at which the battery should be kept [V]
	 * dt 			time step to be taken in the time integration [sec]
	 * 				should be small enough to ensure numerical stability. The stability depends on the current rating and the temperature
	 * 				the order of magnitude is 1 to 5 seconds
	 * 				If the data-collection time interval is smaller than dt, the data-collection time interval is used instead of the specified value
	 * 				i.e. dt = min(dt,CyclingDataTimeInterval)
	 * blockDegradation if true, degradation is not accounted for during this CV
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 * Icut 		the absolute value of the lowest current allowed [A]
	 * 					> 0
	 *
	 * OUT
	 * ahi 			the total discharged capacity [Ah]
	 * whi 			the total discharged energy [Wh]
	 * timei		the total time it took [sec]
	 *
	 * THROWS
	 * 1005			illegal value for the voltage
	 * 1008 		illegal value for the cut off current
	 * 1018 		the underlying function which did the CV loading terminated because of the time limit instead of the current limit
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CV_I with current limit "<<Icut<<"A, and set voltage "<<Vset<<" is starting"<<endl;

	// check if the voltage and cutoff current are allowed
	bool vmax = Vset > c.getVmax(); 				// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in BasicCycler::CV_I. The voltage "<<Vset<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vset < c.getVmin();					// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in BasicCycler::CV_I. The voltage "<<Vset<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (vmax || vmin)
		throw 1005;
	bool curr = Icut < 0;
	if(curr){
		cerr<<"Error in BasicCycler::CV_I. The cutoff current "<<Icut<<" is negative. It must be positive."<<endl<<flush;
		throw 1008;
	}

	// Set a very large time limit so that it is always the current limit which is reached
	double timelim = 99999999;

	// variables
	double ah = 0;									// discharged charge [Ah]
	double wh = 0;									// discharged energy [Wh]
	double tt = 0;									// time on load
	int endcr;										// integer indicating why the underlying function ended

	// Call CV_t_I to load the cell with the CV
	try{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::CV_I is calling CV_t_I with a set voltage of "<<Vset<<", current limit "<<Icut<<", time limit "<<timelim<<endl;
		endcr = CV_t_I(Vset, dt, blockDegradation, timelim, Icut, &ah, &wh, &tt);
	}
	catch(int err){
		if(verbose >= printCrit)
			cout<<"Error in BasicCycler::CV_I while cycling in the underlying function, error "<<err<<", throwing it on"<<endl;
		throw err;
	}

	// Make the output parameters
	*ahi = ah;
	*whi = wh;
	*timei = tt;

	// double check that we have indeed reached the current threshold
	if(endcr == 1){
		cerr<<"ERROR in BasicCycler::CV_I because the underlying function CV_t_I finished because the time threshold was reached which means we didn't reach the current threshold. "<<endl<<flush;
		throw 1018;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CV_I with current limit "<<Icut<<"A, and set voltage "<<Vset<<" is terminating"<<endl;
}

int BasicCycler::CC_t_CV_t(double I, double dt, bool blockDegradation, double time, double Vupp, double Vlow, double* ahi, double* whi, double* timei){
	/*
	 * Function to load a cell for a given time with a constant current (CC).
	 * If a voltage limit is hit in the meantime, the cell is kept at constant voltage (CV) for the rest of the time.
	 *
	 * IN
	 * I 			load current [A], positive = discharge, negative = charge
	 * dt 			time step to be taken in the time integration [sec]
	 * blockDegradation if true, degradation is not accounted for during this (dis)charge
	 * time 		time for which this current should be applied [sec]
	 *				must be a multiple of dt
	 * 				If the data-collection time interval is smaller than dt, the data-collection time interval is used instead of the specified value
	 * 				i.e. dt = min(dt,CyclingDataTimeInterval)
	 * Vupp			upper voltage limit, must be in the voltage range of the cell, Cell.Vmin <= Vupp <= Cell.Vmax, [V]
	 * Vlow			lower voltage limit, nust be in the voltage range of the cell, Cell.Vmin <= Vlow <= Cell.Vmax, [V]
	 *
	 * OUT
	 * ahi 			the total discharged capacity [Ah]
	 * whi 			the total discharged energy [Wh]
	 * timei		the total time the cell has been loaded with the CC [sec]
	 * int 			integer indicating why the CC phase finished
	 * 					-3 	the minimum cell voltage was exceeded
	 * 					-2 	the maximum cell voltage was exceeded
	 * 					0 	an error occurred
	 * 					1 	the full time was completed (i.e. no CV was needed because no voltage limit was reached)
	 * 					2 	the upper voltage limit was reached
	 * 					3 	the lower voltage limit was reached
	 *
	 * THROWS
	 * 1003 		the total time is not a multiple of dt
	 * 1009 		the CC phase finished with an error and we don't know at which voltage limit this happened.
	 * 				so we don't know which voltage we should keep during the CV phase and we can't complete the full time
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_t_CV_t with time = "<<time<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<", and current "<<I<<" is starting"<<endl;

	// *********************************************************** 1 variables & settings ***********************************************************************

	// Check that the total time is a multiple of the time step
	if (remainder(time, dt) > 0.01){
		cerr<<"Error in BasicCycler::CC_t_CV_t. The total time "<<time<<" is not a multiple of the time step dt "<<dt<<endl<<flush;
		throw 1003;
	}

	// ensure the time steps is smaller than the data collection time resolution (if we are collecting data)
	double feedb = CyclingDataTimeInterval;
	if(feedb > 0)
		dt = min(dt, feedb);

	// check that the total time is still multiple of the time step, if not set the time step to 1sec (which always works)
	if (remainder(time, dt) > 0.01)
		dt = 1;

	// variables for output
	double ah1, wh1, tt1; 										// feedback variables for the CC_t function
	int end1;													// integer indicating why the CC phase stopped
	double ah2, wh2, tt2; 										// feedback variables for the CV_t function
	double ah = 0;												// total discharged charge [Ah]
	double wh = 0;												// total discharged energy [Wh]
	double tt = 0;												// total time [sec]

	// *********************************************************** 2 apply a constant current ***********************************************************************

	// First, try to do a CC_t. It will stop cycling as soon as it hit a voltage limit.
	try{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::CC_t_CV_t is starting the CC phase for a time of "<<time<<endl;
		end1 = CC_t_V(I, dt, blockDegradation, time, Vupp, Vlow, &ah1, &wh1, &tt1);
	}
	catch(int e){
		if(verbose >= printCrit)								// this was considered a critical error in CC_t_V so print that in this case, it wasn't a critical error and we might recover from it
			cout<<"Error in a subfunction in BasicCycler::CC_t_CV_t while loading the cell with a constant current: "<<e<<". Skipping the CC phase and hoping this solves the problem."<<endl<<flush;
		ah1 = 0;
		wh1 = 0;
		tt1 = 0;

		// try to check the cell voltage to see which voltage limit was violated
		double v, ocvp, ocvn, etap, etan, rdrop, T;
		try{
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::CC_t_CV_t with is checking the voltage of the cell after an error in the CC phase"<<endl;
			c.getVoltage(verbose >= printCrit, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &T);
		}
		catch(int e2){
			if(verbose >= printCrit)
				cout<<"Error in a subfunction in BasicCycler::CC_t_CV_t when getting the cell voltage while trying to recover from the error in CC phase: "<<e<<". giving up and throwing on the error"<<endl<<flush;
			throw e;
		}
		if(v <= Vlow)
			end1 = 3;
		else if (v >= Vupp)
			end1 = 2;
		else
			end1 = 0;											// we don't know -> this will lead to errors later in the code
	}

	// copy the output parameters
	ah = ah1;
	wh = wh1;
	tt = tt1;

	// *********************************************************** 3 apply a constant voltage if needed ***********************************************************************

	// check if we need to do a CV part or not
	if (end1 == 1){												// no voltage limit was hit, so we should have completed the full cycle

		//todo
		if (!(abs(tt1- time)<pow(10,-3))){
			cout<<"cc_t_cv_t error with time. CC time was "<<tt1<<" \t while total time was "<<time<<" giving an error of "<<abs(tt1- time)<<endl;
		}


		assert(abs(tt1- time)<pow(10,-3));						// double check we have indeed completed the full period, allow a small margin of error (this should never fail, if it does, there is a mistake in the code of CC_t)
	}
	else{														// we haven't completed the full time, so we need to do a CV for the remaining time
		double time2 = time - tt1;								// the remaining time to complete the full time [sec]

		// find which voltage limit was hit by the CC function
		double Vset;
		if(end1 == 2 || end1 == -2)
			Vset = Vupp;
		else if(end1 == 3 || end1 == -3)
			Vset = Vlow;
		else {
			if(verbose >= printCrit)							// we don't know which voltage limit was hit
				cerr<<"Error in BasicCycler::CC_t_CV_t the CC phase finished with an error so we don't know at which voltage to do a CV. Throwing an error."<<endl<<flush;
			throw 1009;
		}

		// apply the constant voltage
		try{
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::CC_t_CV_t with is starting the CV phase for a further "<<time2<<" seconds at a voltage of "<<Vset<<endl;
			CV_t(Vset, dt, blockDegradation, time2, &ah2, &wh2, &tt2);

		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in a subfunction in BasicCycler::CC_t_CV_t while loading the cell with a constant voltage: "<<e<<". Throwing it on."<<endl<<flush;
			throw e;
		}

		// copy the output parameters
		ah += ah2;
		wh += wh2;
		tt += tt2;
	}

	// *********************************************************** 4 output parameters ***********************************************************************

	*ahi = ah;
	*whi = wh;
	*timei = tt;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_t_CV_t with time = "<<time<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<", and current "<<I<<" is terminating with "<<end1<<endl;

	// Return an integer indicating which voltage limit was reached during the CC phase
	return end1;
}

void BasicCycler::CC_V_CV_I(double Crate, double Vset, double Ccut, double dt, bool blockDegradation, double* ahi, double* whi, double* timei){
	/*
	 * Function to get a cell to a specified voltage doing first a CC (dis)charge to the specified voltage,
	 * followed by a CV (dis)charge until the current is below the cutoff current
	 *
	 * IN
	 * Crate 	Crate of the current in the CC phase [-] >0
	 * Vset		the voltage to which the cell should be (dis)charged [V], c.getVmin() <= Vset <= c.getVmax()
	 * Ccut 	the Crate of the cutoff current for the CV (dis)charge [-], >0
	 * 				if Ccut > Crate, no CV phase is done
	 * dt 		the time step to be used
	 * blockDegradation if true, degradation is not accounted for during this CV
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 *
	 * OUT
	 * ahi 		the total discharged capacity [Ah]
	 * whi 		the total discharged energy [Wh]
	 * timei	the total time it took [sec]
	 *
	 * THROWS
	 * 1005		illegal value for the voltage which should be reached
	 * 1008 	illegal value for the current threshold for the CV phase
	 * 1010		illegal value for the C rate for the CC phase
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_V_CV_I with voltage = "<<Vset<< ", CC current Crate "<<Crate<<"C and CV current cutoff Crate "<< Ccut<<"C, is starting"<<endl;

	// check the voltage which should be reached
	bool vmax = Vset > c.getVmax(); 							// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in BasicCycler::CC_V_CV_I. The voltage "<<Vset<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vset < c.getVmin();								// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in BasicCycler::CC_V_CV_I. The voltage "<<Vset<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	if (vmax || vmin)
		throw 1005;

	// check the current threshold for the CV phase
	bool currCV = Ccut < 0;
	if(currCV){
		cerr<<"Error in BasicCycler::CC_V_CV_I. The cutoff C rate "<<Ccut<<" is negative. It must be positive."<<endl<<flush;
		throw 1008;
	}

	// check the C rate for the current during the CC phase
	bool currCC = Crate < 0;
	if(currCC){
		cerr<<"Error in BasicCycler::CC_V_CV_I. The Crate "<<Crate<<" is negative. It must be positive."<<endl<<flush;
		throw 1010;
	}

	// *********************************************************** 1 variables & settings ***********************************************************************

	// variables
	double v;													// voltage of the cell
	double ocvp, ocvn, etap, etan, rdrop, tem;					// unneeded feedback variables
	double ah1 = 0;												// cumulative discharged charge in the CC phase [Ah]
	double wh1= 0;												// cumulative discharged energy in the CC phase [Wh]
	double tt1 = 0;												// cumulative time in the CC phase [sec]
	double ah2 = 0;												// cumulative discharged charge in the CV phase [Ah]
	double wh2= 0;												// cumulative discharged energy in the CV phase [Wh]
	double tt2 = 0;												// cumulative time in the CV phase [sec]
	bool check = true;											// check if the battery state is valid after setting a current, throw an error if not

	// check whether we need to charge or discharge in the CC phase
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::CC_V_CV_I is determining whether to charge or discharge in the CC phase"<<endl;
	int sign;													// integer deciding whether we need to charge or discharge
	try{
		c.setI(verbose >= printCrit, check, 0);					// set the cell current to 0
		c.getVoltage(verbose >= printCrit,&v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem);// get the OCV
		if (v < Vset)
			sign = -1;											// we need to charge (the OCV is lower than the voltage we want to achieve)
		else
			sign = 1;											// we need to discharge discharge
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in CC_V_CV_I when checking if we need to charge or discharge: "<<e<<". Throwing it on"<<endl;
		throw e;
	}

	// *********************************************************** 2 CC phase ***********************************************************************

	// Do the CC phase.
	// the CC phase is terminated if:
	// 		the full time has been reached
	// 		a voltage limit is reached (which might happen instantaneously due to the resistive voltage drop when the current is applied)
	// 		an error occurs
	// In the latter two cases, a CV will be done for the remaining time
	try{
		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::CC_V_CV_I is starting a CC phase with I = "<<sign*Crate*c.getNominalCap()<<"A until the set voltage of "<<Vset<<endl;
		CC_V(sign*Crate*c.getNominalCap(), dt, blockDegradation, Vset, &ah1, &wh1, &tt1); // CC (dis)charge at the given Crate

	}
	catch(int e){
		if(verbose >= printNonCrit)
			cout<<"Error in CC_V_CV_I: error during the CC phase. The OCV was "<<v<<" and we were trying to apply a current of "<<sign*Crate*c.getNominalCap()<<". Skip and go to the CV phase"<<endl;
	} // error during the CC phase so go immediately to the CV

	// *********************************************************** 3 CV phase ***********************************************************************

	// do a CV phase until the given cutoff current if needed
	if(Ccut < Crate){											// cutoff current is larger than the CC current, so we don't have to do a CV phase
		double Icut = Ccut * c.getNominalCap();
		try{
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::CC_V_CV_I with is starting the CV phase"<<endl;
			CV_I(Vset, dt, blockDegradation, Icut,&ah2, &wh2, &tt2);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in CC_V_CV_I: error during the CV phase with at cutoff current of "<<Icut<<"A. error "<<e<<". Throwing it on."<<endl<<flush;
			throw e;
		}
	}

	// *********************************************************** 4 output variables ***********************************************************************

	// Make the output variables
	*ahi = ah1 + ah2;
	*whi = wh1 + wh2;
	*timei = tt1+tt2;

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::CC_V_CV_I with voltage = "<<Vset<< ", CC current Crate "<<Crate<<"C and CV current cutoff Crate "<< Ccut<<"C, is terminating"<<endl;
}

int BasicCycler::followI(int nI, string nameI, bool blockDegradation, int limit, double Vupp, double Vlow, double* ahi, double* whi, double* timei){
	/*
	 * function to follow a certain current pattern.
	 * The cell tries to follow the current every time step.
	 *
	 * IN
	 * nI 		length of the current profile (i.e. number of rows in the csv file)
	 * 				if n is too small, only the first n steps of the profile are simulated
	 * 				if n is too large (i.e. you think there are more rows than there actually are in the csv file)
	 * 					the behaviour is not guaranteed. Probably, the current in the additional steps will be 0 and the function will work
	 * 					but the function might also crash (an error might happen in when the program tries to read the additional rows)
	 * nameI 	name of the CSV-file with the current profile
	 * 				the first column contains the current in [A], positive for discharge, negative for charge
	 * 				the second column contains the time in [sec] the current should be maintained
	 * blockDegradation if true, degradation is not accounted for during this current profile
	 * 				set this to 'true' if you want to ignore degradation for now (e.g. if you're characterising a cell)
	 * limit 	integer describing what to do if the current can't be maintained because a voltage limit is reached
	 * 				0 	immediately go to the next current step of the profile (i.e. reduce the time of this step)
	 * 				1 	keep the voltage constant for the rest of this step of the profile (i.e. reduce the current for the rest of this step)
	 * Vupp		upper voltage limit, Cell.Vmin <= Vupp <= Cell.Vmax, [V]
	 * Vlow		lower voltage limit, Cell.Vmin <= Vlow <= Cell.Vmax, [V]
	 *
	 * OUT
	 * ahi		charge throughput while following the profile [Ah]
	 * whi		energy throughput while following the profile [Wh]
	 * timei	time spent while following the profile [sec]
	 * int 		integer indicating if a voltage limit was hit while following the current profile
	 * 				-1 if the lower voltage limit was hit while following the original profile
	 * 				0 if the profile could be followed, i.e. no voltage limit was hit while following this profile
	 * 				1 if the upper voltage limit was hit while following the original profile
	 * 				10 if both the lower and upper voltage limits were hit while following the original profile
	 * 				100 if we don't know which voltage limit was hit while following the original profile
	 *
	 * throws
	 * 1011 	limit has an illegal value (it is not 0, 1 or 2)
	 * 1012 	if the cell is in an illegal state (i.e. even a current of 0A still violates the voltage constraints of the cell)
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::followI with profile = "<<nameI<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<", is starting"<<endl;

	if (limit <0 || limit > 1){
		cerr<<"ERROR in BasicCycler::followI, illegal value for the limit setting: "<<limit<<", only the values 0 and 1 are allowed. Throwing an error"<<endl<<flush;
		throw 1011;
	}

	// *********************************************************** 1 read the current profile & variables ***********************************************************************
	if(verbose >= printCyclerDetail)
		cout<<"BasicCycler::followI is reading the current profile"<<endl;

	// Read the current profile
	double I[nI], T[nI];
	try{
		loadCSV_2col(nameI, nI, I, T);								// read the file
	}
	catch(int e){
		cout<<"error in BasicCycler::followI when reading the file with the current profile called "<<nameI<<", error "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// variables
	double dt; 														// time step to be used for this step in the profile [sec]
	bool Tlow = c.getTenv() < (273+45);								// boolean indicating if the environmental temperature is below 45 degrees
	bool Ilow;														// boolean indicating if the magnitude of the current in this step of the profile is below 1.5C
	bool Imed;														// boolean indicating if the magnitude of the current in this step of the profile is below 3C
	State s;														// state of the battery
	double Iprev;													// current in the previous time step
	double ah;														// capacity discharged during this step in the profile [Ah]
	double wh;														// energy discharged during this step in the profile [Wh]
	double tt;														// time spent during this step in the profile [sec]
	int vlim;														// integer indicating why the CC phase finished
	double ahtot = 0;												// charge throughput up to this step in the profile [Ah]
	double whtot = 0;												// energy throughput up to this step in the profile [Wh]
	double tttot = 0; 												// cumulative time up to this step in the profile [sec]
	bool vminlim = false;											// boolean to indicate if the minimum voltage limit was hit
	bool vmaxlim = false;											// boolean to indicate if the maximum voltage limit was hit
	bool verr = false;												// boolean to indicate if an unknown error occurred

	// ****************************************************** 2 loop through the profile ***********************************************************************

	for(int i=0;i<nI;i++){

		if(verbose >= printCyclerDetail)
			cout<<"BasicCycler::followI is in step "<<i<<" with a current of "<<I[i]<<" and time of "<<T[i]<<" seconds."<<endl;

		// Determine the time step to be used for the time integration in this step of the profile
		Ilow = abs(I[i]) < 1.5*c.getNominalCap();					// is the current below 1.5C?
		Imed = abs(I[i]) < 3*c.getNominalCap();						// is the current below 3C?
		if (fmod(T[i], 3) == 0 && Tlow && Ilow)						// if the temperature is low, the current is low, and the time of the step is a multiple of 3 sec
			dt = 3;													// then use 3 seconds as time step
		else if (fmod(T[i], 2) == 0 && Imed)						// the current is medium, and the time of the step is a multiple of 2 sec
			dt = 2;													// then use 2 seconds as time step
		else														// else use the lower value of 1 second or the time step
			dt = min(1.0, T[i]);
		// note: if the data collection time interval is smaller than dt, this will be corrected in the underlying functions (in CC_t_V)

		// print a warning if the current is large since this might cause numerical problems
		if(!Imed)
			cout<<"Warning in BasicCycler::followI. The current in step "<<i<<" is "<<I[i]<<", which is above 3C. Very large currents might lead to errors in the model"<<endl<<flush;
			// large currents are a problem because then charge throughput in one time step is too high, and the concentration difference between 2 time steps is too large
			// such that it is possible that in time step t, all is fine (e.g. v = 3.9), and in time step t+1, the li-concentration is larger than the maximum concentration (and the cell voltage would be 8V)

		// store initial states such they can be restored if needed
		c.getStates(s, &Iprev);

		// Follow this step of the profile
		try{
			if(limit == 0)											// follow the profile as long as we can, skip the rest of the step if we hit a voltage limit
				vlim = CC_t_V(I[i], dt, blockDegradation, T[i], Vupp, Vlow, &ah, &wh, &tt);
			else if (limit == 1) 									// keep the voltage constant for the rest of the step if we hit a voltage limit
				vlim = CC_t_CV_t(I[i], dt, blockDegradation, T[i], Vupp, Vlow, &ah, &wh, &tt);
		}
		catch(int e){
			cout<<"Error in a subfunction of BasicCycler::followI when following step "<<i<<" of the profile, which has current "<<I[i]
				<<" A for a duration of "<<T[i]<<" seconds. Error"<< e<<". Throwing it on "<<endl<<flush;
			throw e;
		}

		// check if a voltage limit was hit
		if (vlim == 2 || vlim == -2){								// the upper voltage limit (or the maximum cell voltage) was reached
			vmaxlim = true;
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::followI has hit the upper voltage limit during step "<<i<<endl;
		}
		else if (vlim == 3 || vlim == -3){							// the lower voltage limit (or the minimum cell voltage) was reached
			vminlim = true;
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::followI has hit the lower voltage limit during step "<<i<<endl;
		}
		else if (vlim == 1){}										// no voltage limit was reached
		else{
			verr = true;											// an error occurred and we don't know which voltage limit was reached
			if(verbose >= printCyclerDetail)
				cout<<"BasicCycler::followI has encountered an error so we don't know if a voltage limit was reached during step "<<i<<endl;
		}

		// update the cumulative throughput
		ahtot += abs(ah);
		whtot += abs(wh);
		tttot += tt;
	}

	// *********************************************************** 3 output parameters ***********************************************************************

	// return the throughput
	*ahi = ahtot;
	*whi = whtot;
	*timei = tttot;

	// set the value of the return-integer to indicate which voltage limit was hit
	int endvalue;
	if (verr)														// an unknown voltage limit was hit while following the profile
		endvalue = 100;
	else if (vminlim && vmaxlim)									// both lower and upper voltage limits were hit
		endvalue = 10;
	else if (vminlim)												// lower voltage limit was hit
		endvalue = -1;
	else if (vmaxlim)												// upper voltage limit was hit
		endvalue = 1;
	else															// no voltage limit was hit
		endvalue = 0;
	if(verbose >= printCyclerFunctions)
		cout<<"BasicCycler::followI with profile = "<<nameI<< ", and voltage limits "<<Vupp<<" to "<< Vlow<<", is terminating with "<<endvalue<<endl;

	return endvalue;
}

