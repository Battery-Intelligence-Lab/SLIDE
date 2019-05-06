/*
 * Degradation.cpp
 *
 * Implements simulations for degradation experiments where multiple cells undergo similar degradation experiments but with different parameters (e.g. different temperatures, voltage windows, C rates, etc).
 * The exact degradation procedures are defined in the Cycler.
 * In the functions defined here, the functions from the Cycler are called various times with slightly different parameters.
 * As such, we can simulate the effect the different parameters have on battery degradation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

// Include header files
#include "Degradation.h"
#include "Cycler.hpp"
#include "Cell_user.hpp"


using namespace std;

string print_DEG_ID(DEG_ID degid){
	/*
	 * Function to get a string representation of the struct with the degradation settings
	 * This string is part of the names of the subfolders in which results are written.
	 *
	 * IN
	 * degid 	struct with the degradation identifiers
	 *
	 * OUT
	 * string 	string of the degradation identifiers
	 * 			identifiers of the same mechanism are separated by -
	 * 			identifiers of different mechanisms are separated by _
	 * 			e.g. if we use SEI model 1, no SEI porosity effect, no surface cracks, LAM model 2 and LAM model 3 and lithium plating model 1:
	 * 				2-0_0-0_2-3_1
	 * 				2 		SEI model 1
	 * 				0		no SEI porosity
	 * 				0		no surface cracks
	 * 				0 		don't decrease the diffusion due to surface cracks
	 * 				2		LAM model 2
	 * 				3		LAM model 3
	 * 				1		lithium plating model 1
	 *
	 */

	// start with an empty string
	string id = "";

	// print SEI models and SEI_porosity (decreasing the active volume fraction due to SEI growth), separated by -
	for(int i=0;i<degid.SEI_n;i++){
		id += to_string(degid.SEI_id[i]);
		id += "-";
	}
	id += to_string(degid.SEI_porosity);

	// mechanism separator
	id += "_";

	// print CS models (surface cracking) and CS_diffusion (reducing the diffusion constant because of cracks), separated by -
	for(int i=0;i<degid.CS_n;i++){
		id += to_string(degid.CS_id[i]);
		id += "-";
	}
	id += to_string(degid.CS_diffusion);

	// mechanism separator
	id += "_";

	// print LAM models separated by -
	for(int i=0;i<degid.LAM_n;i++){
		id += to_string(degid.LAM_id[i]);
		if(i < degid.LAM_n-1)					// print separator - only between LAM models
			id += "-";							// not after the last LAM model
	}

	// mechanism separator
	id += "_";

	// print plating model
	id += to_string(degid.pl_id);

	// output
	return id;
}

void Calendar_one(const struct Model& M, const struct DEG_ID& degid, int cellType, int verbose,
		double V, double Ti, int Time, int mode, int timeCycleData, int timeCheck, struct checkUpProcedure proc, string name){
	/*
	 * Function which simulates one calendar ageing regime.
	 * It does little more than calling the corresponding function from Cycler.cpp
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * degid 		struct with degradation settings (which degradation models to be used)
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
	 * V 			the voltage at which the battery has to rest [V]
	 * Ti 			the temperature at which the battery has to rest [K]
	 * Time 		the time for which the battery has to rest [days]
	 * mode 		integer deciding how often to recharge the cell to the specified voltage.
	 * 				When a cell rests, the voltage will slip a bit due to degradation, so the question is how often to recharge.
	 * 				0 	recharge only after a check-up
	 * 				1 	recharge every day
	 * 				2 	float the cell at the specified voltage, instead of letting it rest
	 * 					this option is not recommended, it takes ages to calculate
	 * 				Any other value will produce an error
	 * timeCycleData the time interval at which cycle data (such as the voltage of the cell) should be stored [s]
	 * 				if 0, no cycle data is stored
	 * timeCheck 	the time after which a check-up has to be done [days]
	 * 				the total rest time is floored to the nearest multiple of timeCheck
	 * proc 		structure with the parameters of the check-up procedure with the following fields:
	 * 		blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
	 * 								and at the end of the check-up, the exact battery state from the start of the check-up is restored
	 * 								if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
	 * 								and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
	 * 		capCheck 			boolean indicating if the capacity should be checked
	 * 		OCVCheck			boolean indicating if the half-cell OCV curves should be checked
	 * 		CCCVCheck			boolean indicating if some CCCV cycles should be done as part of the check-up procedure
	 * 		pulseCheck			boolean indicating if a pulse discharge test should be done as part of the check-up procedure
	 * 		includeCycleData	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
	 * 		nCycles				number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
	 * 		Crates				array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
	 * 		Ccut				C rate of the cutoff current for the CV phase for the CCCV check-up, must be positive
	 * 		profileName 		name of the csv file which contains the current profile for the pulse test
	 * 								the first column contains the current in [A] (positive for discharge, negative for charge)
	 * 								the second column contains the time in [sec] the current should be maintained
	 * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
	 * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
	 * name 		the name of the subfolder in which all the data for this simulation is written, must obey the naming convention for folders
	 * 				avoid special characters or spaces
	 */

	// settings of the cycler
	double dt = 2;												// use a time step of 2 seconds
	if (Ti < 40)
		dt = 5;													// a lower temperature allows a larger time step without numerical problems

	// Make a cell, the type of the cell depending on the value of 'cellType'
	Cell c1;
	if (cellType ==0)
		c1 = Cell_KokamNMC (M, degid, verbose);					// a high power NMC cell made by Kokam
	else if (cellType ==1)
		c1 = Cell_LGChemNMC (M, degid, verbose);				// a high energy NMC cell made by LG Chem
	else
		c1 = Cell_user(M, degid, verbose);						// a user-defined cell

	// Make the cycler
	Cycler cycler(c1, name, verbose, timeCycleData);

	// Call the Calendar-function of the cycler. Wrap it in a try-catch to avoid fatal errors
	try{
		cycler.calendarAgeing(dt, V, Ti, Time, timeCheck, mode, proc);
	}
	catch(int err){
		cout<<"Calendar_one experienced error "<<err<<" during execution of "<<name<<", abort this test"<<endl << std::flush;
		if(err == 15){
			cout<<"Error 15 means that the cell had degraded too much to continue simulating. \n"
				"This can be due to too much SEI growth, too much loss of lithium, too much loss of active material (thin electrodes, low volume fraction, or low effective surface)  \n"
				"too much surface cracking, too low diffusion constants, too high resistance, or too much lithium plating. \n"
				"Continue simulating might lead to errors in the code so the simulation of this degradation test is stopped. \n"
				"The results which have been written are all valid and you can ignore the error messages above which explained where this error comes from."<<endl<<flush;
		}
	}
}

void Cycle_one(const struct Model& M, const struct DEG_ID& degid, int cellType, int verbose, double Vma, double Vmi,
		double Ccha, bool CVcha, double Ccutcha, double Cdis, bool CVdis, double Ccutdis, double Ti, int timeCycleData, int nrCycles, int nrCap, struct checkUpProcedure proc, string name){
	/*
	 * Function which simulates one cycle ageing regime.
	 * It does little more than calling the corresponding function from Cycler.cpp
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * degid 		struct with degradation settings (which degradation models to be used)
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
	 * Vma 			maximum voltage of this cycle, below the maximum voltage of the cell [V]
	 * Vmi 			minimum voltage of this cycle, above the minimum voltage of the cell [V]
	 * Ccha 		C rate at which the battery should be charged during the CC phase[-]
	 * CVcha 		boolean indicating if a CV charge should be done after the CC charge (true) or not (false)
	 * 				if true, charging is CC CV
	 * 				if false, charging is CC only
	 * Ccutcha 		cutoff C rate for the CV charge [A], > 0
	 * Cdis 		C rate at which the battery should be discharged [-]
	 * CVdis 		boolean indicating if a CV discharge should be done after the CC discharge (true) or not (false)
	 * 				if true, discharging is CC CV
	 * 				if false, discharging is CC only
	 * Ccutdis 		cutoff C rate for the CV discharge [A], > 0
	 * Ti 			temperature at which the cell has to be cycled [K]
	 * timeCycleData the time interval at which cycling data (e.g. the voltage of the cell) should be stored [s]
	 * 				if 0, no cycle data is stored
	 * nrCycles 	number of cycles to be simulated in total [-]
	 * nrCap 		the number of cycles between consecutive check-ups [-]
	 * proc 		structure with the parameters of the check-up procedure with the following fields:
	 * 		blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
	 * 								and at the end of the check-up, the exact battery state from the start of the check-up is restored
	 * 								if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
	 * 								and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
	 * 		capCheck 			boolean indicating if the capacity should be checked
	 * 		OCVCheck			boolean indicating if the half-cell OCV curves should be checked
	 * 		CCCVCheck			boolean indicating if some CCCV cycles should be done as part of the check-up procedure
	 * 		pulseCheck			boolean indicating if a pulse discharge test should be done as part of the check-up procedure
	 * 		includeCycleData	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
	 * 		nCycles				number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
	 * 		Crates				array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
	 * 		Ccut				C rate of the cutoff current for the CV phase for the CCCV check-up, must be positive
	 * 		profileName 		name of the csv file which contains the current profile for the pulse test
	 * 								the first column contains the current in [A] (positive for discharge, negative for charge)
	 * 								the second column contains the time in [sec] the current should be maintained
	 * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
	 * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
	 * name 		the name of the subfolder in which all the data for this simulation is written, must obey the naming convention for folders
	 * 				avoid special characters or spaces
	 */

	// settings of the cycler
	double dt = 2;												// use a time step of 2 seconds to ensure numerical stability
	if (Ti < 40)
		dt = 3;													// a lower temperature allows a larger time step without numerical problems

	// Make a cell, the type of the cell depending on the value of 'cellType'
	Cell c1;
	if (cellType ==0)
		c1 = Cell_KokamNMC (M, degid, verbose);					// a high power NMC cell made by Kokam
	else if (cellType ==1)
		c1 = Cell_LGChemNMC (M, degid, verbose);				// a high energy NMC cell made by LG Chem
	else
		c1 = Cell_user(M, degid, verbose);						// a user-defined cell

	// Make the cycler
	Cycler cycler(c1, name, verbose, timeCycleData);

	// Call the cycle ageing function from the cycler
	try{
		cycler.cycleAgeing(dt, Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, nrCycles, nrCap, proc);
	}
	catch(int err){
		cout<<"Cycle_one experienced error "<<err<<" during execution of "<<name<<", abort this test"<<endl << std::flush;
		if(err == 15){
			cout<<"Error 15 means that the cell had degraded too much to continue simulating. \n"
				"This can be due to too much SEI growth, too much loss of lithium, too much loss of active material (thin electrodes, low volume fraction, or low effective surface)  \n"
				"too much surface cracking, too low diffusion constants, too high resistance, or too much lithium plating. \n"
					"Continue simulating might lead to errors in the code so the simulation of this degradation test is stopped. \n"
				"The results which have been written are all valid and you can ignore the error messages above which explained where this error comes from."<<endl<<flush;
		}
	}
}

void Profile_one(const struct Model& M, const struct DEG_ID& degid, int cellType, int verbose, string profName, int n, int limit,
		double Vma, double Vmi, double Ti, int timeCycleData, int nrProfiles, int nrCap, struct checkUpProcedure proc, string name){
	/*
	 * Calls the ProfileAgeing() function of a Cycler
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * degid	 	struct with degradation settings (which degradation models to be used)
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
	 * profName 	name of the csv file with the current profile
	 * n			length of the profile (number of rows in the csv file)
	 * limit		integer describing what to do if the current can't be maintained because a voltage limit is reached
	 * 					0 	immediately go to the next current step of the profile (i.e. reduce the time of the step)
	 * 					1 	keep the voltage constant for the rest of this step of the profile (i.e. reduce the current of the step)
	 * Vma 			maximum voltage under which the cell should stay, below the maximum voltage of the cell [V]
	 * Vmi 			minimum voltage above which the cell should stay, above the minimum voltage of the cell [V]
	 * Ti 			temperature at which the cell has to be cycled [K]
	 * timeCycleData the time interval at which cycling data (e.g. the voltage of the cell) should be stored [s]
	 * 				if 0, no cycle data is stored
	 * nrProfiles 	number of profiles to be simulated in total [-]
	 * nrCap 		number of profile repetitions between consecutive check-ups [-]
	 * proc 		structure with the parameters of the check-up procedure with the following fields:
	 * 		blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
	 * 								and at the end of the check-up, the exact battery state from the start of the check-up is restored
	 * 								if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
	 * 								and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
	 * 		capCheck 			boolean indicating if the capacity should be checked
	 * 		OCVCheck			boolean indicating if the half-cell OCV curves should be checked
	 * 		CCCVCheck			boolean indicating if some CCCV cycles should be done as part of the check-up procedure
	 * 		pulseCheck			boolean indicating if a pulse discharge test should be done as part of the check-up procedure
	 * 		includeCycleData	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
	 * 		nCycles				number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
	 * 		Crates				array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
	 * 		Ccut				C rate of the cutoff current for the CV phase for the CCCV check-up, must be positive
	 * 		profileName 		name of the csv file which contains the current profile for the pulse test
	 * 								the first column contains the current in [A] (positive for discharge, negative for charge)
	 * 								the second column contains the time in [sec] the current should be maintained
	 * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
	 * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
	 * name 		the name of the subfolder in which all the data for this simulation is written, must obey the naming convention for folders
	 * 				avoid special characters or spaces
	 */

	// Make a cell, the type of the cell depending on the value of 'cellType'
	Cell c1;
	if (cellType ==0)
		c1 = Cell_KokamNMC (M, degid, verbose);					// a high power NMC cell made by Kokam
	else if (cellType ==1)
		c1 = Cell_LGChemNMC (M, degid, verbose);				// a high energy NMC cell made by LG Chem
	else
		c1 = Cell_user(M, degid, verbose);						// a user-defined cell

	// Make the cycler
	Cycler cycler(c1, name, verbose, timeCycleData);

	// Print a warning if you want to store cycling data
	// In profileAgeing, you are guaranteed to get one point per step in the profile
	// so if the steps are very short (e.g. 1sec), you are storing a huge amount of data (e.g. every second)
	// 	This seriously slowing down the code (taking hours instead of minutes) and produces huge amounts of data (several GB)
	if (timeCycleData != 0)
		cout<<"Warning for profile ageing: the cycling data of the cell is going to be stored, which will lead to much slower calculation (several hours) and a huge amount of data (several GB)"<<endl;

	// Call the profile ageing function from the Cycler
	try{
		cycler.profileAgeing(n, profName, limit, Vma, Vmi, Ti, nrProfiles, nrCap, proc);
	}
	catch(int err){
		cout<<"Profile_one experienced error "<<err<<" during execution of "<<name<<", abort this test"<<endl << std::flush;
		if(err == 15){
			cout<<"Error 15 means that the cell had degraded too much to continue simulating. \n"
				"This can be due to too much SEI growth, too much loss of lithium, too much loss of active material (thin electrodes, low volume fraction, or low effective surface)  \n"
				"too much surface cracking, too low diffusion constants, too high resistance, or too much lithium plating. \n"
					"Continue simulating might lead to errors in the code so the simulation of this degradation test is stopped. \n"
				"The results which have been written are all valid and you can ignore the error messages above which explained where this error comes from."<<endl<<flush;
		}
	}
}

void CycleAgeing(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose){
	/*
	 * Function to simulate a selection of cycle ageing experiments.
	 * The cells are cycled continuously with the same cycle (with occasional check-up procedures)
	 * Various simulations are done, with different settings for the cycles (temperatures, voltage windows, currents, etc.)
	 * The results from the check-up procedures and the cycling data from cells are written in one subfolder per simulation.
	 *
	 * The function uses multi-threaded code to accelerate the calculation.
	 * The optimal number of threads to use, is the total number of logical cores on your computer minus one (to avoid flooding the CPU)
	 * Currently, this function is programmed to use 3 threads, which is optimal for a dual-core PC.
	 * The function prints a warning if you PC has more or fewer cores, and a different number of threads is better
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * pref 		string with which the name of the subfolder in which the results should be written, will begin
	 * 				e.g. if pref = '1', then the subfolder will be called 1_xxxxxx (with xxxx the degradation identifier)
	 * 				use this as an identifier for different simulations, e.g. the next time you simulate, change pref to '2'
	 * 				this avoids that you overwrite previous results
	 * 				must obey the naming convention for folders, avoid special characters or spaces
	 * degid	 	struct with degradation settings (which degradation models to be used)
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
	 *
	 * OUT
	 * The function writes many csv files with the results of the check-up tests and cycling data of the cell.
	 * The following files are written in a subfolder called pref_degid
	 * (with pref the string in the variable pref and degid the string representing the degradation identifier, see print_DEG_ID):
	 *
	 * The cycling data (voltage and temperature on cycling) is written in files called CyclingData_x.csv
	 * 			(where x = 0 for the first data batch, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in BasicCycler::writeCyclingData or Cycling::CCCV or in the word documents
	 * The capacity measurement are written in a csv file called DegradationData_batteryState.csv
	 * 			exact description of the file can be found in Cycler::checkUp_batteryStates or in the word document
	 * The half-cell OCV curves are written in a csv file called DegradationData_OCV.csv
	 * 			exact description of the file can be found in Cycler::checkUp_OCV
	 * The voltage measurements during the CCCV cycles from the check-ups are written in csv files called DegradationData_CheckupCycle_x.csv
	 * 			(x = 0 for the first check-up, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in Cycler::checkUp_CCCV
	 * The voltage measurements during the pulse discharge test from the check-ups are written in csv files called DegradationData_CheckupPulse_x.csv
	 * 			(x = 0 for the first check-up, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in Cycler::checkUp_pulse
	 */

	// *********************************************************** 1 variables ***********************************************************************

	// append the ageing identifiers to the prefix
	pref = pref + "_";
	pref += print_DEG_ID(degid);
	pref += "_";

	// find the number of parallel threads that are optimal to use
	unsigned int Ncor = std::thread::hardware_concurrency();	// Ncor is the number of (logical) cores, 0 if c++ can't identify it
	if (Ncor != 0 && Ncor != 4)
		cout<<"The CycleAgeing function uses a non-optimal number of threads. Ideally, you should use "<<Ncor-1<<" cores."<<endl;
		// if you see this warning, you should reorganise the code below.
		// every 'std::thread xxxx'-command starts a new thread.
		// every 'xxxx.join()'-command terminates a thread.
		// the number of threads started before they are all joined is the number of threads used simultaneously.
		// to change this, simply cut-and-paste the code.

	// Make variables to describe the cycling regimes
		string name;												// identification strings for each experiment (will also be used to name the folder so must have the same restrictions as the prefix, i.e. no spaces, no special characters, etc.)
		double Ti;													// environmental temperature [K]
		double Vma, Vmi;											// maximum and minimum voltages of the cycle [V]
		double Ccha, Cdis;											// Crate of the CC charge and CC discharge
		bool CVcha = true;											// we want to have a CC CV charge (if false, then charge has only a CC phase)
		double Ccutcha = 0.05;										// Crate of the cutoff current for the CV phase of the charge [-]
		bool CVdis = false;											// we want to have a CC discharge (if true, then discharge has both a CC and CV phase)
		double Ccutdis = 1.0;										// Crate of the cutoff current for the CV phase of the discharge [-]
		int nrCycles = 3000;										// the number of cycles which has to be simulated
		int nrCap = 500;											// the number of cycles between check-ups
		int timeCycleData = 60;										// time interval at which cycling data (voltage and temperature) has to be recorded [s]
																	// 	0 means no data is recorded
																	//  if not 0, data is recorded approximately every so many seconds

	// *********************************************************** 2 check-up procedure ******************************************************************

	// Make a struct to describe the check-up procedure
		struct checkUpProcedure proc;
		proc.blockDegradation = true; 								// boolean indicating if degradation is accounted for during the check-up, [RECOMMENDED: TRUE]
		proc.capCheck = true;										// boolean indicating if the capacity should be checked
		proc.OCVCheck = true;										// boolean indicating if the half-cell OCV curves should be checked
		proc.CCCVCheck = true;										// boolean indicating if some CCCV cycles should be done as part of the check-up procedure
		proc.pulseCheck = true;										// boolean indicating if a pulse discharge test should be done as part of the check-up procedure
		proc.includeCycleData = true;								// boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
		proc.nCycles = 3;											// number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
																	// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[0] = 0.5;										// do a 0.5C cycle as part of the CCCV check-up, must be positive
																	// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[1] = 1.0;										// do a 1C cycle as part of the CCCV check-up, must be positive
																	// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[2] = 2.0;										// do a 2C cycle as part of the CCCV check-up, must be positive
																	// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Ccut_cha = 0.05;										// C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
		proc.Ccut_dis = 100;										// C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
																	// 		if the cutoff is larger than the value in the CC phase, no CV phase is done
		proc.profileName = "CheckupPulseProfile.csv";				// name of the csv file which contains the current profile for the pulse test
																	//	the first column contains the current in [A] (positive for discharge, negative for charge)
																	//	the second column contains the time in [sec] the current should be maintained
																	//	the profile must be a net discharge, i.e. sum (I*dt) > 0
		proc.profileLength = 13;									// length of the current profiles for the pulse test (number of rows in the csv file)

		// *********************************************************** 3 simulations ******************************************************************

		// The cycling-window needs the voltage windows between which it should cycle the battery.
		// Users typically find it easier to think in terms of state of charge windows.
		// For the Kokam cell used here, the conversion is as follows (for other cells the user has to derive the conversion from the OCV curve)
			// 0%		2.7V
			// 10%		3.42V
			// 20%		3.49V
			// 50% 		3.65V
			// 80%		3.98V
			// 90%		4.08V
			// 100%		4.2V
		// The voltages will be slightly different for the LG Chem cell because it has different OCV curves

	// Print a message that we are starting the first batch of simulations
	cout<<pref<<"\t Batch 1"<<endl;

		// cycle at environmental temperature T=45, 1C charge, 1C discharge, SoC windows 0-100%
		cout<<"Cycle 1 T45, 1C1D, 0-100%"<<endl;					// print a message to the user saying which cycle we are simulating
		Vma = 4.2;
		Vmi = 2.7;
		Ti = 273+45;
		Ccha = 1;
		Cdis = 1;
		name = pref + "T45_1C1D_SoC0-100";
		//Cycle_one(M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);					// simulate this (without using multi-threaded computation
		std::thread cyc1(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name); 	// make a new thread and simulate this

		// T=45, 2C charge, 1C discharge , 0-100%
		cout<<"Cycle 2, T=45, 2C1D, 0-100%"<<endl;
		Vma = 4.2;
		Vmi = 2.7;
		Ti = 273+45;
		Ccha = 2;
		Cdis = 1;
		name = pref + "T45_2C1D_SoC0-100";
		std::thread cyc2(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=45, 3C, 0-100%
		cout<<"Cycle 3 T=45, 3C1D, 0-100%"<<endl;
		Vma = 4.2;
		Vmi = 2.7;
		Ti = 273+45;
		Ccha = 3;
		Cdis = 1;
		name = pref + "T45_3C1D_SoC0-100";
		std::thread cyc3(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// We have now started three threads (called cyc1, cyc2 and cyc3).
		// We want to have 3 running at the same time, so we need to join them here.
		// If you want to have more threads running simultaneously (e.g. 7):
			// cut the three join() commands below, and paste them below cyc7 (when you will have started 7 threads)
			// also cut the three join() commands for threads cyc4, cyc5 and cyc6 to that location.
			// and finally cut the cyc7.join() there too.
			// then you will have started 7 threads before you join them, so there will be 7 running simultaneously.
			// do the same for the other threads.
	cyc1.join();
	cyc2.join();
	cyc3.join();

	cout<<pref<<"\t Batch 2"<<endl;

		// T=45, 1C, 0-80%
		cout<<"Cycle 4 T=45, 1C1D, 0-80%"<<endl;
		Vma = 3.98;
		Vmi = 2.7;
		Ti = 273+45;
		Ccha = 1;
		Cdis = 1;
		name = pref + "T45_1C1D_SoC0-80";
		std::thread cyc4(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=45, 2C, 0-80%
		cout<<"Cycle 5 T=45, 2C1D, 0-80%"<<endl;
		Vma = 3.98;
		Vmi = 2.7;
		Ti = 273+45;
		Ccha = 2;
		Cdis = 1;
		name = pref + "T45_2C1D_SoC0-80";
		std::thread cyc5(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=45, 3C, 0-80%
		cout<<"Cycle 6 T=45, 3C1D, 0-80%"<<endl;
		Vma = 3.98;
		Vmi = 2.7;
		Ti = 273+45;
		Ccha = 3;
		Cdis = 1;
		name = pref + "T45_3C1D_SoC0-80";
		std::thread cyc6(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

	// if you want to have 7 threads, cut these three commands to below std::thread cyc7
	cyc4.join();
	cyc5.join();
	cyc6.join();

	cout<<pref<<"\t Batch 3"<<endl;

		// T=45, 1C, 20-100%
		cout<<"Cycle 7 T=45, 3C1D, 20-100%"<<endl;
		Vma = 4.2;
		Vmi = 3.49;
		Ti = 273+45;
		Ccha = 1;
		Cdis = 1;
		name = pref + "T45_1C1D_SoC20-100";
		std::thread cyc7(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// If you want to have 7 threads running simultaneously, cut-paste the 6 previous join() commands here
		// also cut-paste the cyc7.join() command here

		// T=45, 2C, 20-100%
		cout<<"Cycle 8 T=45, 2C1D, 20-100%"<<endl;
		Vma = 4.2;
		Vmi = 3.49;
		Ti = 273+45;
		Ccha = 2;
		Cdis = 1;
		name = pref + "T45_2C1D_SoC20-100";
		std::thread cyc8(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=45, 3C, 20-100%
		cout<<"Cycle 9 T=45, 3C1D, 20-100%"<<endl;
		Vma = 4.2;
		Vmi = 3.49;
		Ti = 273+45;
		Ccha = 3;
		Cdis = 1;
		name = pref + "T45_3C1D_SoC20-100";
		std::thread cyc9(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

	cyc7.join(); 			// this line has to be cut-pasted to just after you start thread cyc7 to have 7 threads running simultaneously
	cyc8.join();
	cyc9.join();

	cout<<pref<<"\t Batch 4"<<endl;

		// Cycle ageing 1
		// T=45, 1C, 10-90%
		cout<<"Cycle 10"<<"\t"<<"T=45, 1C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+45;
		Ccha = 1;
		Cdis = 1;
		name = pref + "T45_1C1D_SoC10-90";
		std::thread cyc10(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=45, 2C, 10-90%
		cout<<"Cycle 11"<<"\t"<<"T=45, 2C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+45;
		Ccha = 2;
		Cdis = 1;
		name = pref + "T45_2C1D_SoC10-90";
		std::thread cyc11(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=45, 3C, 10-90%
		cout<<"Cycle 12"<<"\t"<<"T=45, 3C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+45;
		Ccha = 3;
		Cdis = 1;
		name = pref + "T45_3C1D_SoC10-90";
		std::thread cyc12(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

	cyc10.join();
	cyc11.join();
	cyc12.join();

	cout<<pref<<"\t Batch 5"<<endl;

		// T=5, 1C, 10-90%
		cout<<"Cycle 13"<<"\t"<<"T=5, 1C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+5;
		Ccha = 1;
		Cdis = 1;
		name = pref + "T5_1C1D_SoC10-90";
		std::thread cyc13(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=5, 2C, 10-90%
		cout<<"Cycle 14"<<"\t"<<"T=5, 2C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+5;
		Ccha = 2;
		Cdis = 1;
		name = pref + "T5_2C1D_SoC10-90";
		std::thread cyc14(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=5, 3C, 10-90%
		cout<<"Cycle 15"<<"\t"<<"T=5, 3C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+5;
		Ccha = 3;
		Cdis = 1;
		name = pref + "T5_3C1D_SoC10-90";
		std::thread cyc15(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

	cyc13.join();
	cyc14.join();
	cyc15.join();

	cout<<pref<<"\t Batch 6"<<endl;

		// T=25, 1C, 10-90%
		cout<<"Cycle 16"<<"\t"<<"T=25, 1C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+25;
		Ccha = 1;
		Cdis = 1;
		name = pref + "T25_1C1D_SoC10-90";
		std::thread cyc16(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=25, 2C, 10-90%
		cout<<"Cycle 17"<<"\t"<<"T=25, 2C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+25;
		Ccha = 2;
		Cdis = 1;
		name = pref + "T25_2C1D_SoC10-90";
		std::thread cyc17(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=25, 3C, 10-90%
		cout<<"Cycle 18"<<"\t"<<"T=25, 3C1D, 10-90%"<<endl;
		Vma = 4.08;
		Vmi = 3.42;
		Ti = 273+25;
		Ccha = 3;
		Cdis = 1;
		name = pref + "T25_3C1D_SoC10-90";
		std::thread cyc18(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

	cyc16.join();
	cyc17.join();
	cyc18.join();

	cout<<pref<<"\t Batch 7"<<endl;

		// T=25, 1C, 10-90%
		cout<<"Cycle 19"<<"\t"<<"T=25, 1C1D, 0-100%"<<endl;
		Vma = 4.2;
		Vmi = 2.7;
		Ti = 273+25;
		Ccha = 1;
		Cdis = 1;
		name = pref + "T25_1C1D_SoC0-100";
		std::thread cyc19(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=25, 2C, 0-100%
		cout<<"Cycle 20"<<"\t"<<"T=25, 2C1D, 0-100%"<<endl;
		Vma = 4.2;
		Vmi = 2.7;
		Ti = 273+25;
		Ccha = 2;
		Cdis = 1;
		name = pref + "T25_2C1D_SoC0-100";
		std::thread cyc20(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

		// T=25, 3C, 0-100%
		cout<<"Cycle 21"<<"\t"<<"T=25, 3C1D, 0-100%"<<endl;
		Vma = 4.2;
		Vmi = 2.7;
		Ti = 273+25;
		Ccha = 3;
		Cdis = 1;
		name = pref + "T25_3C1D_SoC0-100";
		std::thread cyc21(Cycle_one,M, degid, cellType, verbose,Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, timeCycleData, nrCycles, nrCap, proc, name);

	cyc19.join();
	cyc20.join();
	cyc21.join();
}

void CalendarAgeig(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose){
	/*
	 * Function to simulate a selection of calendar ageing experiments.
	 * The cells are rested for long times (with occasional check-up procedures)
	 * Various simulations are done, with different settings for the resting (temperatures, state of charge)
	 * The results from the check-up procedures and the cycling data from cells are written in one subfolder per simulation.
	 *
	 * The function uses multi-threaded code to accelerate the calculation.
	 * The optimal number of threads to use, is the total number of logical cores on your computer minus one (to avoid flooding the CPU)
	 * Currently, this function is programmed to use 3 threads, which is optimal to a dual-core PC.
	 * The function prints a warning if you PC has more or fewer cores, and a different number of threads is better
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * pref 		string with which the name of the subfolder in which the results should be written, will begin
	 * 				e.g. if pref = '1', then the subfolder will be called 1_xxxxxx (with xxxx the degradation identifier)
	 * 				use this as an identifier for different simulations, e.g. the next time you simulate, change pref to '2'
	 * 				this avoids that you overwrite previous results
	 * 				must obey the naming convention for folders, avoid special characters or spaces
	 * degid	 	struct with degradation settings (which degradation models to be used)
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
	 *
	 * OUT
	 * The function writes many csv files with the results of the check-up tests and cycling data of the cell.
	 * The following files are written in a subfolder called pref_degid
	 * (with pref the string in the variable pref and degid the string representing the degradation identifier, see print_DEG_ID):
	 *
	 * The cycling data (voltage and temperature on cycling) is written in files called CyclingData_x.csv
	 * 			(where x = 0 for the first data batch, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in BasicCycler::writeCyclingData or Cycling::CCCV
	 * The capacity measurement are written in a csv file called DegradationData_batteryState.csv
	 * 			exact description of the file can be found in Cycler::checkUp_batteryStates
	 * The half-cell OCV curves are written in a csv file called DegradationData_OCV.csv
	 * 			exact description of the file can be found in Cycler::checkUp_OCV
	 * The voltage measurements during the CCCV cycles from the check-ups are written in csv files called DegradationData_CheckupCycle_x.csv
	 * 			(x = 0 for the first check-up, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in Cycler::checkUp_CCCV
	 * The voltage measurements during the pulse discharge test from the check-ups are written in csv files called DegradationData_CheckupPulse_x.csv
	 * 			(x = 0 for the first check-up, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in Cycler::checkUp_pulse
	 */

	// *********************************************************** 1 variables ***********************************************************************

	// append the ageing identifiers to the prefix
	pref = pref + "_";
	pref += print_DEG_ID(degid);
	pref += "_";

	// find the number of parallel threads we can have:
	unsigned int Ncor = std::thread::hardware_concurrency();	// Ncor is the number of (logical) cores, 0 if c++ can't identify it
	if (Ncor != 0 && Ncor != 4)
		cout<<"The CalendarAgeing function uses a non-optimal number of threads. Ideally, you should use "<<Ncor-1<<" cores."<<endl;
		// if you see this warning, you should reorganise the code below.
		// every 'std::thread xxxx'-command starts a new thread.
		// every 'xxxx.join()'-command terminates a thread.
		// the number of threads started before they are all joined is the number of threads used now.
		// to change this, simply cut-and-paste the code.

	// Make variables to describe the cycling regimes
		string name;											// identification strings for each experiment (will also be used to name the folder so must have the same restrictions as the prefix, i.e. no spaces, no special characters, etc.)
		int mode = 0;											// integer deciding how often to recharge the cells (due to degradation, the voltage will decrease over time. But no self-discharge is simulated)
																// 0 means we recharge to the specified voltage only after a check-up
																// 1 means we recharge to the specified voltage every day
																// 2 means we float the cell at a constant voltage, rather than resting it
																// 		floating at constant voltage takes very long to simulate
		double V;												// voltage at which the cell has to rest [V] above the minimum and below the maximum voltage of the cell
		double Ti;												// temperature at which the cell has to rest [K]
		int Time = 20*30;										// time the cell has to rest [days]
		int timeCheck = 30;										// time between consecutive check-ups [days]
		int timeCycleData = 3600;								// time interval at which cycling data (voltage and temperature) has to be recorded [s]
																// 	0 means no data is recorded
																//  if not 0, data is recorded approximately every so many seconds

	// *********************************************************** 2 check-up procedure ******************************************************************

	// Make a struct to describe the check-up procedure
		struct checkUpProcedure proc;
		proc.blockDegradation = true; 							// boolean indicating if degradation is accounted for during the check-up, [RECOMMENDED: TRUE]
		proc.capCheck = true;									// boolean indicating if the capacity should be checked
		proc.OCVCheck = true;									// boolean indicating if the half-cell OCV curves should be checked
		proc.CCCVCheck = true;									// boolean indicating if some CCCV cycles should be done as part of the check-up procedure
		proc.pulseCheck = true;									// boolean indicating if a pulse discharge test should be done as part of the check-up procedure
		proc.includeCycleData = true;							// boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
		proc.nCycles = 3;										// number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates), maximum 100
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[0] = 0.5;									// do a 0.5C cycle as part of the CCCV check-up, must be positive
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[1] = 1.0;									// do a 1C cycle as part of the CCCV check-up, must be positive
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[2] = 2.0;									// do a 2C cycle as part of the CCCV check-up, must be positive
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Ccut_cha = 0.05;									// C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
		proc.Ccut_dis = 100;									// C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
																// 		if the cutoff is larger than the value in the CC phase, no CV phase is done
		proc.profileName = "CheckupPulseProfile.csv";			// name of the csv file which contains the current profile for the pulse test
																//	the first column contains the current in [A] (positive for discharge, negative for charge)
																//	the second column contains the time in [sec] the current should be maintained
																//	the profile must be a net discharge, i.e. sum (I*dt) > 0
		proc.profileLength = 13;								// length of the current profiles for the pulse test (number of rows in the csv file)

	// *********************************************************** 3 simulations ******************************************************************

		// The cycling-window needs the voltage windows between which it should cycle the battery.
		// Users typically find it easier to think in terms of state of charge windows.
		// For the Kokam cell used here, the conversion is as follows (for other cells the user has to derive the conversion from the OCV curve)
			// 0%		2.7V
			// 10%		3.42V
			// 20%		3.49V
			// 50% 		3.65V
			// 80%		3.98V
			// 90%		4.08V
			// 100%		4.2V
		// The voltages will be slightly different for the LG Chem cell because it has different OCV curves


	// Print a message that we are starting the first batch of simulations
	cout<<pref<<"\t Batch 1"<<endl;

	// T=5, 100%
		cout<<"calendar 1, T=5, 100%"<<endl;
		V = 4.2;
		Ti = 273+5;
		name = pref + "Cal-T5_SoC100";
		//Calendar_one(M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);					// regular (not multi-threaded) call to simulate this
		std::thread cal1 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);	// make a new thread and simulate this

		// T=5, 90%
		cout<<"calendar 2, T=5, 90%"<<endl;
		V = 4.08;
		Ti = 273+5;
		name = pref + "Cal-T5_SoC90";
		std::thread cal2 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

		// T=5, 50%
		cout<<"calendar 3, T=5, 50%"<<endl;
		V = 3.65;
		Ti = 273+5;
		name = pref + "Cal-T5_SoC50";
		std::thread cal3 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

		// We have now started three threads (cal1, cal2, cal3).
		// We want to have 3 running at the same time, so we need to join them here.
		// If you want to have more threads running simultaneously (e.g. 7):
			// cut the three join() commands below, and paste them below cal7 (when you will have started 7 threads)
			// also cut the three join() commands for threads cal4, cal5, cal6 to that location.
			// and finally cut the cal7.join() there too.
			// then you will have started 7 threads before you join them, so there will be 7 running simultaneously.
			// do the same for the other threads.
	cal1.join();
	cal2.join();
	cal3.join();

	cout<<pref<<"\t Batch 2"<<endl; 								// print to let the user know

		// T=25, 100%
		cout<<"calendar 4, T=25, 100%"<<endl;
		V = 4.2;
		Ti = 273+25;
		name = pref + "Cal-T25_SoC100";
		std::thread cal4 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

		// T=25, 90%
		cout<<"calendar 5, T=25, 90%"<<endl;
		V = 4.08;
		Ti = 273+25;
		name = pref + "Cal-T25_SoC90";
		std::thread cal5 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

		// T=25, 50%
		cout<<"calendar 6, T=25, 50%"<<endl;
		V = 3.65;
		Ti = 273+25;
		name = pref + "Cal-T25_SoC50";
		std::thread cal6 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

	// if you want to have 7 threads, cut these three commands to below std::thread cal7
	cal4.join();
	cal5.join();
	cal6.join();

	cout<<pref<<"\t Batch 3"<<endl; 								// print to let the user know

		// T=45, 100%
		cout<<"calendar 7, T=45, 100%"<<endl;
		V = 4.2;
		Ti = 273+45;
		name = pref + "Cal-T45_SoC100";
		std::thread cal7 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

		// If you want to have 7 threads running simultaneously, cut-paste the 6 previous join() commands here
		// also cut-paste the cal7.join() command here


		// T=45, 90%
		cout<<"calendar 8, T=45, 90%"<<endl;
		V = 4.08;
		Ti = 273+45;
		name = pref + "Cal-T45_SoC90";
		std::thread cal8 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

		// T=45, 50%
		cout<<"calendar 9, T=45, 50%"<<endl;
		V = 3.65;
		Ti = 273+45;
		name = pref + "Cal-T45_SoC50";
		std::thread cal9 (Calendar_one,M, degid, cellType, verbose, V, Ti, Time, mode, timeCycleData, timeCheck, proc, name);

	cal7.join();				// this line has to be cut-pasted to just after you start thread cal7 to have 7 threads running simultaneously
	cal8.join();
	cal9.join();
}

void ProfileAgeing(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose){
	/*
	 * Function to simulate a selection of drive cycle ageing experiments.
	 * The cells are cycled continuously with the same current profile (with occasional check-up procedures)
	 * Various simulations are done, with different settings for the cycles (temperatures, voltage windows, currents, etc.)
	 * The results from the check-up procedures and the cycling data from cells are written in one subfolder per simulation.
	 *
	 * The function uses multi-threaded code to accelerate the calculation.
	 * The optimal number of threads to use, is the total number of logical cores on your computer minus one (to avoid flooding the CPU)
	 * Currently, this function is programmed to use 3 threads, which is optimal to a dual-core PC.
	 * The function prints a warning if you PC has more or fewer cores, and a different number of threads is better
	 *
	 * IN
	 * M 			matrices of the spatial discretisation for the solid diffusion PDE
	 * pref 		string with which the name of the subfolder in which the results should be written, will begin
	 * 				e.g. if pref = '1', then the subfolder will be called 1_xxxxxx (with xxxx the degradation identifier)
	 * 				use this as an identifier for different simulations, e.g. the next time you simulate, change pref to '2'
	 * 				this avoids that you overwrite previous results
	 * 				must obey the naming convention for folders, avoid special characters or spaces
	 * degid	 	struct with degradation settings (which degradation models to be used)
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
	 *
	 * OUT
	 * The function writes many csv files with the results of the check-up tests and cycling data of the cell.
	 * The following files are written in a subfolder called pref_degid
	 * (with pref the string in the variable pref and degid the string representing the degradation identifier, see print_DEG_ID):
	 *
	 * The cycling data (voltage and temperature on cycling) is written in files called CyclingData_x.csv
	 * 			(where x = 0 for the first data batch, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in BasicCycler::writeCyclingData or Cycling::CCCV
	 * The capacity measurement are written in a csv file called DegradationData_batteryState.csv
	 * 			exact description of the file can be found in Cycler::checkUp_batteryStates
	 * The half-cell OCV curves are written in a csv file called DegradationData_OCV.csv
	 * 			exact description of the file can be found in Cycler::checkUp_OCV
	 * The voltage measurements during the CCCV cycles from the check-ups are written in csv files called DegradationData_CheckupCycle_x.csv
	 * 			(x = 0 for the first check-up, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in Cycler::checkUp_CCCV
	 * The voltage measurements during the pulse discharge test from the check-ups are written in csv files called DegradationData_CheckupPulse_x.csv
	 * 			(x = 0 for the first check-up, x = 1 for the second, etc.)
	 * 			exact description of the file can be found in Cycler::checkUp_pulse
	 */

	// *********************************************************** 1 variables ***********************************************************************


	// append the ageing identifiers to the prefix
	pref = pref + "_";
	pref += print_DEG_ID(degid);
	pref += "_";

	// find the number of parallel threads that are optimal to use
	unsigned int Ncor = std::thread::hardware_concurrency();	// Ncor is the number of (logical) cores, 0 if c++ can't identify it
	if (Ncor != 0 && Ncor != 4)
		cout<<"The ProfileAgeing function uses a non-optimal number of threads. Ideally, you should use "<<Ncor-1<<" cores."<<endl;
		// if you see this warning, you should reorganise the code below.
		// every 'std::thread xxxx'-command starts a new thread.
		// every 'xxxx.join()'-command terminates a thread.
		// the number of threads started before they are all joined is the number of threads used now.
		// to change this, simply cut-and-paste the code.

	// Give the current profile which should be followed
		// the first column must give the current [A] which should be followed
		// 		>0 is discharge
		//		<0 is charge
		// the second column should give the time [sec] for which this current should be maintained
		//		values for time are floored to the integer value (ie. the numbers behind the comma are ignored)
	// 5 example current profiles are provided:
	// 	name										length		description
	// 	Current Profile random.csv 					100			a random current profile with a maximum current of 5A, each current step takes maximum 1000 seconds
	//	Current Profile drive cycle HWFET.csv		766			HWFET drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
	//	Current Profile drive cycle NYCC.csv		599			NYCC drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
	//	Current Profile drive cycle UDDS.csv		1370		UDDS drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
	// 	Current Profile drive cycle US06.csv		601			US06 drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
	string profile;											// name of the csv file with the current profile
	int length;												// length of the profile (number of rows)
	int limit = 0;											// what to do if the voltage limits are reached while following the profile:
												  				// 0 immediately go to the next current step of the profile (i.e. reduce the time of this step)
												  				// 1 keep the voltage constant for the rest of this step of the profile (i.e. reduce the current of this step)
	double Vma;												// maximum voltage under which the cell should stay while following the current profile [V]
	double Vmi;												// minimum voltage above which the cell should stay while following the current profile [V]
	double Tenv;											// environmental temperature while the cell is following the profile [K]
	string name;											// identification strings for each experiment (will also be used to name the folder so must have the same restrictions as the prefix, i.e. no spaces, no special characters, etc.)
	int nrProfiles = 10000;									// number of times the current profile should be repeated [-]
	int nrCap = 1000;										// number of times the current profile should be repeated between consecutive check-ups [-]
	int timeCycleData = 0;									// time interval at which cycling data (voltage and temperature) has to be recorded [s]
															// 	0 means no data is recorded
															//  if not 0, data is recorded approximately every so many seconds
															// For this function (profileAgeing) it is highly recommended to keep timeCycleData at 0.
															// 	Otherwise, a huge amount of data is generated (many GB) and consequently the code slows down dramatically
															// 	This is because you are guaranteed to get one data point per step in the current profile if timeCycleData is bigger than 0,
															// 	independent of the time interval at which you are otherwise collecting data (i.e. even if you set it to 100, you will still get a data point every step)
															// 	As many drive cycles will have a step per second (i.e. the current changes every second), this means you store a data point every second (even if timeCycleData = 100)
															// 	which is obviously a huge amount of data if you simulate long term battery usage

	// *********************************************************** 2 check-up procedure ******************************************************************

	// Make a struct to describe the check-up procedure
		struct checkUpProcedure proc;
		proc.blockDegradation = true; 							// boolean indicating if degradation is accounted for during the check-up, [RECOMMENDED: TRUE]
		proc.capCheck = true;									// boolean indicating if the capacity should be checked
		proc.OCVCheck = true;									// boolean indicating if the half-cell OCV curves should be checked
		proc.CCCVCheck = true;									// boolean indicating if some CCCV cycles should be done as part of the check-up procedure
		proc.pulseCheck = true;									// boolean indicating if a pulse discharge test should be done as part of the check-up procedure
		proc.includeCycleData = true;							// boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
		proc.nCycles = 3;										// number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates).
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[0] = 0.5;									// do a 0.5C cycle as part of the CCCV check-up, must be positive
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[1] = 1.0;									// do a 1C cycle as part of the CCCV check-up, must be positive
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Crates[2] = 2.0;									// do a 2C cycle as part of the CCCV check-up, must be positive
																// If you change this variable, also change it in the mablab script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
		proc.Ccut_cha = 0.05;									// C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
		proc.Ccut_dis = 100;									// C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
																// 		if the cutoff is larger than the value in the CC phase, no CV phase is done
		proc.profileName = "CheckupPulseProfile.csv";			// name of the csv file which contains the current profile for the pulse test
																//	the first column contains the current in [A] (positive for discharge, negative for charge)
																//	the second column contains the time in [sec] the current should be maintained
																//	the profile must be a net discharge, i.e. sum (I*dt) > 0
		proc.profileLength = 13;								// length of the current profiles for the pulse test (number of rows in the csv file)

	// *********************************************************** 3 simulations ******************************************************************

		// The cycling-window needs the voltage windows between which it should cycle the battery.
		// Users typically find it easier to think in terms of state of charge windows.
		// For the Kokam cell used here, the conversion is as follows (for other cells the user has to derive the conversion from the OCV curve)
			// 0%		2.7V
			// 10%		3.42V
			// 20%		3.49V
			// 50% 		3.65V
			// 80%		3.98V
			// 90%		4.08V
			// 100%		4.2V
		// The voltages will be slightly different for the LG Chem cell because it has different OCV curves


	// Print a message that we are starting the first batch of simulations
	cout<<pref<<"\t Batch 1"<<endl;

	// HWFET cycle, 0-100% SoC, T = 25
	profile = "Current Profile drive cycle HWFET.csv";
	length = 766;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+25;
	name = pref + "prof-HWFET-T25_SoC0-100";
//	Profile_one(M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);					// regular (not multi-threaded) call to simulate this experiment
 	thread p1(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);		// start a new thread to simulate this experiment

	// HWFET cycle, 0-100% SoC, T = 45
	profile = "Current Profile drive cycle HWFET.csv";
	length = 766;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+45;
	name = pref + "prof-HWFET-T45_SoC0-100";
	thread p2(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// HWFET cycle, 10-90% SoC, T = 25
	profile = "Current Profile drive cycle HWFET.csv";
	length = 766;
	Vma = 4.08;
	Vmi = 3.42;
	Tenv = 273+25;
	name = pref + "prof-HWFET-T25_SoC10-90";
	thread p3(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// We have now started three threads (p1, p2, p3).
	// We want to have 3 running at the same time, so we need to join them here.
	// If you want to have more threads running simultaneously (e.g. 7):
		// cut the three join() commands below, and paste them below p7 (when you will have started 7 threads)
		// also cut the three join() commands for threads p4, p5, p6 to that location.
		// and finally cut the p7.join() there too.
		// then you will have started 7 threads before you join them, so there will be 7 running simultaneously.
		// do the same for the other threads.
	p1.join();
	p2.join();
	p3.join();

	// Start the first batch of simulations
	cout<<pref<<"\t Batch 2"<<endl; 								// print to let the user know

	// NYCC cycle, 0-100% SoC, T = 25
	profile = "Current Profile drive cycle NYCC.csv";
	length = 599;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+25;
	name = pref + "prof-NYCC-T25_SoC0-100";
	thread p4(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// NYCC cycle, 0-100% SoC, T = 45
	profile = "Current Profile drive cycle NYCC.csv";
	length = 599;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+45;
	name = pref + "prof-NYCC-T45_SoC0-100";
	thread p5(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// NYCC cycle, 10-90% SoC, T = 25
	profile = "Current Profile drive cycle NYCC.csv";
	length = 599;
	Vma = 4.08;
	Vmi = 3.42;
	Tenv = 273+25;
	name = pref + "prof-NYCC-T25_SoC10-90";
	thread p6(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// if you want to have 7 threads, cut these three commands to below std::thread p7
	p4.join();
	p5.join();
	p6.join();

	// Start the first batch of simulations
	cout<<pref<<"\t Batch 3"<<endl; 								// print to let the user know

	// UDDS cycle, 0-100% SoC, T = 25
	profile = "Current Profile drive cycle UDDS.csv";
	length = 1370;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+25;
	name = pref + "prof-UDDS-T25_SoC0-100";
	thread p7(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// If you want to have 7 threads running simultaneously, cut-paste the 6 previous join() commands here
	// also cut-paste the p7.join() command here

	// UDDS cycle, 0-100% SoC, T = 45
	profile = "Current Profile drive cycle UDDS.csv";
	length = 1370;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+45;
	name = pref + "prof-UDDS-T45_SoC0-100";
	thread p8(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// UDDS cycle, 10-90% SoC, T = 25
	profile = "Current Profile drive cycle UDDS.csv";
	length = 1370;
	Vma = 4.08;
	Vmi = 3.42;
	Tenv = 273+25;
	name = pref + "prof-UDDS-T25_SoC10-90";
	thread p9(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	p7.join();														// this line has to be cut-pasted to just after you start thread cal7 to have 7 threads running simultaneously
	p8.join();
	p9.join();

	// Start the first batch of simulations
	cout<<pref<<"\t Batch 4"<<endl; 								// print to let the user know

	// US06 cycle, 0-100% SoC, T = 25
	profile = "Current Profile drive cycle US06.csv";
	length = 601;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+25;
	name = pref + "prof-US06-T25_SoC0-100";
	thread p10(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// US06 cycle, 0-100% SoC, T = 45
	profile = "Current Profile drive cycle US06.csv";
	length = 601;
	Vma = 4.2;
	Vmi = 2.7;
	Tenv = 273+45;
	name = pref + "prof-US06-T45_SoC0-100";
	thread p11(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	// US06 cycle, 10-90% SoC, T = 25
	profile = "Current Profile drive cycle US06.csv";
	length = 601;
	Vma = 4.08;
	Vmi = 3.42;
	Tenv = 273+25;
	name = pref + "prof-US06-T25_SoC10-90";
	thread p12(Profile_one,M, degid, cellType, verbose, profile, length, limit, Vma, Vmi, Tenv, timeCycleData, nrProfiles, nrCap, proc, name);

	p10.join();
	p11.join();
	p12.join();
}
