/*
 * Cycling.cpp
 *
 * Script which groups some top-level functions to simulate how a cell is cycled.
 * The voltage and temperatures while cycling are written to csv files, which can be read by Matlab scripts.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

// Include header files
#include "Cycling.h"
#include "Degradation.h"
#include "Cycler.hpp"
#include "Cell_user.hpp"

void CCCV(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose){
	/*
	 * Function to load a cell with a few CCCV cycles.
	 * The voltage and temperature are measured every few seconds, and the results are written to a csv file in a subfolder.
	 * the csv file contains one row per data entry, and one column per measured property.
	 * The columns are as follows:
	 * 		1			2				3				4	5	6		7		8				9			10			11			12				13				14				15
	 * 		total_time 	Ah_throughput 	Wh_throughput 	I 	V 	OCVp 	OCVn 	Temperature 	charge_time charge_Ah 	charge_Wh 	discharge_time 	discharge_Ah 	discharge_Wh 	rest_time
	 * This csv file can be read by the Matlab function 'readCCCV.m'
	 *
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * pref 		string with which the name of the subfolder in which the results should be written, will begin
	 * 				e.g. if pref = "1", then the subfolder will be called 1_xxxxxx (with xxxx the degradation identifier)
	 * 				use this as an identifier for different simulations, e.g. the next time you simulate, change pref to '2'
	 * 				this avoids that you overwrite previous results
	 * 				avoid spaces and special characters in the string
	 * ageingID 	struct with degradation settings (which degradation models to be used)
	 * cellType 	integer deciding which cell to use for the simulation
	 *  				0 	Kokam cell (high power Kokam NMC)
	 *  				1 	Panasonic cell (high energy LGChem NMC)
	 *  				2 	user cell
	 * verbose 		integer indicating how verbose the simulation has to be.
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
	 */

	// *********************************************************** 1 Make the cell ***********************************************************************

	// Make a cell, the type of the cell depending on the value of 'cellType'
	Cell c1;
	if (cellType ==0)
		c1 = Cell_KokamNMC (M, degid, verbose);					// a high power NMC cell made by Kokam
	else if (cellType ==1)
		c1 = Cell_LGChemNMC (M, degid, verbose);				// a high energy NMC cell made by LG Chem
	else
		c1 = Cell_user(M, degid, verbose);						// a user-defined cell

	// Settings of the cycles
	double T = 273+25;											// temperature at which the cycling should be done [K]
	double Crate = 0.1;											// C rate of the current for the CC phase of the initial charge.
	double Vmax = c1.getVmax();									// voltage to which the cell should be charged, use the maximum cell voltage [V]
	double Vmin = c1.getVmin();									// voltage to which to discharge the cell, use the minimum cell voltage [V]
	double Ccut_cha = 0.05;										// use a current threshold of 0.05C for the CV phases on charging
	double Ccut_dis = 10;										// use a current threshold of 10C for the CV phases on charging
																// 	if the threshold is larger than the C rate on discharge, no CV will be done (because the terminal condition is reached immediately)
	bool blockDegradation = true;								// don't account for degradation while doing the cycles
	double ahi, whi, timei;										// discharged charge, energy and time spend during the CCCV cycle

	// Set the temperature of the cell
	c1.setTenv(T);												// set the environmental temperature
	c1.setT(T);													// set the cell temperature

	// *********************************************************** 2 Make the cycler ***********************************************************************
	// A BasicCyler represents a battery tester (i.e. it has functions to load a cell with a CC or CV)
	// A Cycler represents the programs one would write on the battery tester

	// settings of the cycler
	string ID = "CCCV";											// identification string for the data of this cell (will also be used to name the folder so must have the same restrictions as the prefix, i.e. no spaces, no special characters, etc.)
	int timeCycleData = 4;										// time interval at which cycling data has to be recorded [s]
																// 	0 means no data is recorded
																//  if not 0, data is recorded approximately every so many seconds
	double dt = 2;												// use a time step of 2 seconds for time integration

	// append the ageing identifiers to the prefix
	pref = pref + "_";
	pref += print_DEG_ID(degid);
	pref += "_";

	// Make a Cycler
	Cycler cycler(c1, pref + ID, verbose, timeCycleData);


	// *********************************************************** 3 simulate the cycles ***********************************************************************
	// The functionality is implemented by the functions in the BasicCycler (which is part of the Cycler)
	// the function CC_V_CV_I first does a CC phase at the specified C rate to the specified voltage, and then does a CV at that voltage until the current falls below the specified limit
	// It automatically collects the data while the cell is being loaded.

	// Charge the cell
	cycler.CC_V_CV_I(Crate, Vmax, Ccut_cha, dt, blockDegradation, &ahi, &whi, &timei);

	// Do a full 0.5C cycle
	Crate = 0.5;
	cycler.CC_V_CV_I(Crate, Vmin, Ccut_cha, dt, blockDegradation, &ahi, &whi, &timei);
	cycler.CC_V_CV_I(Crate, Vmax, Ccut_dis, dt, blockDegradation, &ahi, &whi, &timei);

	// Do a full 1C cycle
	Crate = 1;
	cycler.CC_V_CV_I(Crate, Vmin, Ccut_cha, dt, blockDegradation, &ahi, &whi, &timei);
	cycler.CC_V_CV_I(Crate, Vmax, Ccut_dis, dt, blockDegradation, &ahi, &whi, &timei);

	// Do a full 2C cycle
	Crate = 2;
	cycler.CC_V_CV_I(Crate, Vmin, Ccut_cha, dt, blockDegradation, &ahi, &whi, &timei);
	cycler.CC_V_CV_I(Crate, Vmax, Ccut_dis, dt, blockDegradation, &ahi, &whi, &timei);

	// *********************************************************** 4 write the data to a file ***********************************************************************
	// The cycling data of a cell is stored automatically by the Cycler.
	// A cycler has a certain 'buffer' of data to store (100,000 data points).
	// When the buffer is full, the data is written to a file and the buffer is emptied, so that new data points can be stored.
	// The user can force a write of the stored data by calling the function 'writeCyclingData()' to write results before the buffer is full (which is the case here)

	// The csv file with the data is put in a subfolder of the cycler (each cycler has its own folder)
	// The name of the subfolder is determined by the values of 'pref', degid and 'ID'.
	// 		it is called ppp_degid_iii (where 'ppp', 'degid' and 'iii' are the strings stored in the variables 'pref', 'degid' and 'ID')
	// 		the function print_DEG_ID converts the degradation identifiers to a string representation. It is located in Degradation.cpp
	// The name of the csv file is CyclingData_i.csv where i is a number indicating which data file it is.
	// 		the file with the first data batch is called CyclingData_0.csv
	// 		the file with the second data batch is called CyclingData_1.csv
	// 		etc.
	// Alternatively, the user can specify the name of the csv file in which the data should be written
	// 		then you have to call writeCyclingData(name), with 'name' a string with the name of the csv file.
	// 		This is implemented in the two lines below. Comment line 141 which writes the data to a file with the standard name, and enable the line 144-145

	// write to a csv file with the standard naming convention (CyclingData_0.csv)
	cycler.writeCyclingData();

	// write to a specific file (comment in line 123 and uncomment the lines below)
	// string name = "CCCV.csv";
	// cycler.writeCyclingData(name);

	/*
	 * the written csv file has one row per data entry.
	 * Each column represents one measured property.
	 * The columns are as follows:
	 * 		1			2				3				4	5	6		7		8				9			10			11			12				13				14				15
	 * 		total_time 	Ah_throughput 	Wh_throughput 	I 	V 	OCVp 	OCVn 	Temperature 	charge_time charge_Ah 	charge_Wh 	discharge_time 	discharge_Ah 	discharge_Wh 	rest_time
	 */
}

void FollowCurrent(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose){
	/*
	 * Function to let a cell has follow a specified current profile.
	 * The voltage and temperature are measured every few seconds, and the results are written to a csv file in a subfolder.
	 * the csv file contains one row per data entry, and one column per measured property.
	 * The columns are as follows:
	 * 		1			2				3				4	5	6		7		8				9			10			11			12				13				14				15
	 * 		total_time 	Ah_throughput 	Wh_throughput 	I 	V 	OCVp 	OCVn 	Temperature 	charge_time charge_Ah 	charge_Wh 	discharge_time 	discharge_Ah 	discharge_Wh 	rest_time
	 * This csv file can be read by the Matlab function 'readFollowCurrent.m'
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * pref 		string with which the name of the subfolder in which the results should be written, will begin
	 * 				e.g. if pref = '1', then the subfolder will be called 1_xxxxxx (with xxxx the degradation identifier)
	 * 				use this as an identifier for different simulations, e.g. the next time you simulate, change pref to '2'
	 * 				this avoids that you overwrite previous results
	 * 				avoid spaces and special characters
	 * ageingID 	struct with degradation settings (which degradation models to be used)
	 * cellType 	integer deciding which cell to use for the simulation
	 *  				0 	Kokam cell (high power Kokam NMC)
	 *  				1 	Panasonic cell (high energy LGChem NMC)
	 *  				2 	user cell
	 * verbose 		integer indicating how verbose the simulation has to be.
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
	 */

	// *********************************************************** 1 the current profile ***********************************************************************
	// Give the current profile which should be followed
		// the first column must give the current [A] which should be followed
		// 		>0 is discharge
		//		<0 is charge
		// the second column should give the time [sec] for which this current should be maintained
	// 5 example current profiles are provided:
	// 	name										length		description
	// 	Current Profile random.csv 					100			a random current profile with a maximum current of 5A, each current step takes maximum 1000 seconds
	//	Current Profile drive cycle HWFET.csv		766			HWFET drive cycle with a maximum current of 8.1A (around 3C), each current step takes 1 second
	//	Current Profile drive cycle NYCC.csv		599			NYCC drive cycle with a maximum current of 8.1A (around 3C), each current step takes 1 second
	//	Current Profile drive cycle UDDS.csv		1370		UDDS drive cycle with a maximum current of 8.1A (around 3C), each current step takes 1 second
	// 	Current Profile drive cycle US06.csv		601			US06 drive cycle with a maximum current of 8.1A (around 3C), each current step takes 1 second

	string profile = "Current Profile drive cycle UDDS.csv";	// name of the csv file with the current profile
	int length = 1370;											// length of the profile (number of rows in the csv file)

	int limit = 1;												// what to do if the voltage limits are reached while following the profile:
												  					// 0 immediately go to the next current step of the profile (i.e. reduce the time of the step)
												  					// 1 keep the voltage constant for the rest of this step of the profile (i.e. reduce the current for the rest of the step)

	// *********************************************************** 2 Make the cell ***********************************************************************

	// Make a cell, the type of the cell depending on the value of 'cellType'
	Cell c1;
	if (cellType ==0)
		c1 = Cell_KokamNMC (M, degid, verbose);					// a high power NMC cell made by Kokam
	else if (cellType ==1)
		c1 = Cell_LGChemNMC (M, degid, verbose);				// a high energy NMC cell made by LG Chem
	else
		c1 = Cell_user(M, degid, verbose);						// a user-defined cell

	// Settings of the cycles
	double T = 273+25;											// temperature at which the cycling should be done [K]
	double Vmax = c1.getVmax();									// voltage to which the cell should be charged, now use the maximum cell voltage [V]
	double Vmin = c1.getVmin();									// voltage to which to discharge the cell, now use the minimum cell voltage [V]
	bool blockDegradation = true;								// don't account for degradation while doing the cycles
	double ahi, whi, timei;										// discharge charge, energy and time spend during the current profile

	// Set the temperature and voltage limits of the cell
	c1.setTenv(T);												// set the environmental temperature
	c1.setT(T);													// set the cell temperature

	// *********************************************************** 3 Make the cycler ***********************************************************************
	// A BasicCyler represents a battery tester (i.e. it has functions to load a cell with a CC or CV)
	// A Cycler represents the programs one would write on the battery tester

	// settings of the cycler
	string ID = "followCurrent";								// identification string for the data of this cell (will also be used to name the folder so must have the same restrictions as the prefix, i.e. no spaces, no special characters, etc.)
	int timeCycleData = 2;										// time interval at which cycling data has to be recorded [s]
																// 	0 means no data is recorded
																//  if not 0, data is recorded approximately every so many seconds
																// 	but it is ensured that there is at least one data point per current step
																// 		so if a given current has to be maintained for 1 second, and 'timeCycleData' == 10
																// 		then there will still be a data point for that 1 second

	// append the ageing identifiers to the prefix
	pref = pref + "_";
	pref += print_DEG_ID(degid);
	pref += "_";

	// Make a Cycler
	Cycler cycler(c1, pref + ID, verbose, timeCycleData);

	// *********************************************************** 4 follow the current profile ***********************************************************************
	// The functionality is implemented by the functions in the BasicCycler (which is part of the Cycler)
	// followI follows a current profile by doing a series of CC (dis)charges.
	// The cycling data is automatically stored

	cycler.followI(length, profile, blockDegradation, limit, Vmax, Vmin, &ahi, &whi, &timei);


	// *********************************************************** 5 write the data to a file ***********************************************************************
	// The cycling data of a cell is stored automatically by the Cycler.
	// A cycler has a certain 'buffer' of data to store (100,000 data points).
	// When the buffer is full, the data is written to a file and the buffer is emptied, so that new data points can be stored.
	// The user can force a write of the stored data by calling the function 'writeCyclingData()' to write results before the buffer is full (which is the case here)

	// The csv file with the data is put in a subfolder of the cycler (each cycler has its own folder)
	// The name of the subfolder is determined by the values of 'pref', degid and 'ID'.
	// 		it is called ppp_degid_iii (where 'ppp', 'degid' and 'iii' are the strings stored in the variables 'pref', 'degid' and 'ID')
	// 		the function print_DEG_ID converts the degradation identifiers to a string representation. It is located in Degradation.cpp
	// The name of the csv file is CyclingData_i.csv where i is a number indicating which data file it is.
	// 		the file with the first data batch is called CyclingData_0.csv
	// 		the file with the second data batch is called CyclingData_1.csv
	// 		etc.
	// Alternatively, the user can specify the name of the csv file in which the data should be written
	// 		then you have to call writeCyclingData(name), with 'name' a string with the name of the csv file.
	// 		This is implemented in the two lines below. Comment line 284 which writes the data to a file with the standard name, and enable the line 287-288

	// write to a csv file with the standard naming convention (CyclingData_0.csv)
	cycler.writeCyclingData();

	// write to a specific file (comment in line 123 and uncomment the lines below)
	// string name = "followCurrent.csv";
	// cycler.writeCyclingData(name);

	/*
	 * the written csv file has one row per data entry.
	 * Each column represents one measured property.
	 * The columns are as follows:
	 * 		1			2				3				4	5	6		7		8				9			10			11			12				13				14				15
	 * 		total_time 	Ah_throughput 	Wh_throughput 	I 	V 	OCVp 	OCVn 	Temperature 	charge_time charge_Ah 	charge_Wh 	discharge_time 	discharge_Ah 	discharge_Wh 	rest_time
	 */
}


