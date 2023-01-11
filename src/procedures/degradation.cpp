/*
 * degradation.cpp
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

//!< Include header files
#include "degradation.hpp"
#include "cycler.hpp"
#include "cell_user.hpp"
#include "../utility/utility.hpp"

void Cycle_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose, //!< simulate one cycle ageing experiment
               const struct CycleAgeingConfig &cycAgConfig, bool CVcha, double Icutcha, bool CVdis, double Icutdis, int timeCycleData, int nrCycles, int nrCap, struct checkUpProcedure &proc, const std::string &pref)
{
  Cycle_one(M, degid, cellType, verbose, cycAgConfig.Vma, cycAgConfig.Vmi, //!< simulate one cycle ageing experiment
            cycAgConfig.Ccha,
            CVcha,
            Icutcha,
            cycAgConfig.Cdis,
            CVdis,
            Icutdis,
            cycAgConfig.Ti(),
            timeCycleData,
            nrCycles,
            nrCap,
            proc,
            cycAgConfig.get_name(pref));
}

void Calendar_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose,
                  double V, double Ti, int Time, int mode, int timeCycleData, int timeCheck, struct checkUpProcedure &proc, std::string name)
{
  /*
   * Function which simulates one calendar ageing regime.
   * It does little more than calling the corresponding function from cycler.cpp
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

  //!< settings of the cycler
  double dt = 2; //!< use a time step of 2 seconds
  if (Ti < 40)
    dt = 5; //!< a lower temperature allows a larger time step without numerical problems

  //!< Make a cell, the type of the cell depending on the value of 'cellType' #TODO  -> There is a cast.
  auto createCell = [&] {
    if (cellType == 0)
      return (Cell)Cell_KokamNMC(M, degid, verbose); //!< a high power NMC cell made by Kokam
    else if (cellType == 1)
      return (Cell)Cell_LGChemNMC(M, degid, verbose); //!< a high energy NMC cell made by LG Chem
    else
      return (Cell)slide::Cell_user(M, degid, verbose); //!< a user-defined cell
  };

  Cell c1 = createCell();

  //!< Make the cycler
  Cycler cycler(c1, name, verbose, timeCycleData);

  //!< Call the Calendar-function of the cycler. Wrap it in a try-catch to avoid fatal errors
  try {
    cycler.calendarAgeing(dt, V, Ti, Time, timeCheck, mode, proc);
  } catch (int err) {
    //!< std::cout << "Throw test: " << 75 << '\n';
    std::cout << "Calendar_one experienced error " << err << " during execution of " << name << ", abort this test.\n";
    if (err == 15) //!< #TODO this is not valid anymore.
    {
      std::cout << "Error 15 means that the cell had degraded too much to continue simulating.\n"
                   "This can be due to too much SEI growth, too much loss of lithium, too much loss of active material (thin electrodes, low volume fraction, or low effective surface)\n"
                   "too much surface cracking, too low diffusion constants, too high resistance, or too much lithium plating.\n"
                   "Continue simulating might lead to errors in the code so the simulation of this degradation test is stopped.\n"
                   "The results which have been written are all valid and you can ignore the error messages above which explained where this error comes from.\n";
    }
  }
}

void Cycle_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose, double Vma, double Vmi,
               double Ccha, bool CVcha, double Ccutcha, double Cdis, bool CVdis, double Ccutdis, double Ti, int timeCycleData,
               int nrCycles, int nrCap, struct checkUpProcedure &proc, std::string name)
{
  /*
   * Function which simulates one cycle ageing regime.
   * It does little more than calling the corresponding function from cycler.cpp
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

  //!< settings of the cycler
  double dt = 2; //!< use a time step of 2 seconds to ensure numerical stability
  if (Ti < 40)
    dt = 3; //!< a lower temperature allows a larger time step without numerical problems

  //!< Make a cell, the type of the cell depending on the value of 'cellType' #TODO  -> There is a cast.
  auto createCell = [&] {
    if (cellType == 0)
      return (Cell)Cell_KokamNMC(M, degid, verbose); //!< a high power NMC cell made by Kokam
    else if (cellType == 1)
      return (Cell)Cell_LGChemNMC(M, degid, verbose); //!< a high energy NMC cell made by LG Chem
    else
      return (Cell)slide::Cell_user(M, degid, verbose); //!< a user-defined cell
  };

  Cell c1 = createCell();

  //!< Make the cycler
  Cycler cycler(c1, name, verbose, timeCycleData);

  //!< Call the cycle ageing function from the cycler
  try {
    cycler.cycleAgeing(dt, Vma, Vmi, Ccha, CVcha, Ccutcha, Cdis, CVdis, Ccutdis, Ti, nrCycles, nrCap, proc);
  } catch (int err) {
    //!< std::cout << "Throw test: " << 76 << '\n';
    std::cout << "Cycle_one experienced error " << err << " during execution of " << name << ", abort this test.\n";
    if (err == 15) //!< #TODO this is not valid anymore.
    {
      std::cout << "Error 15 means that the cell had degraded too much to continue simulating.\n"
                   "This can be due to too much SEI growth, too much loss of lithium, too much loss of active material (thin electrodes, low volume fraction, or low effective surface) \n"
                   "too much surface cracking, too low diffusion constants, too high resistance, or too much lithium plating.\n"
                   "Continue simulating might lead to errors in the code so the simulation of this degradation test is stopped.\n"
                   "The results which have been written are all valid and you can ignore the error messages above which explained where this error comes from.\n";
    }
  }
}

void Profile_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose, std::string profName, int n, int limit,
                 double Vma, double Vmi, double Ti, int timeCycleData, int nrProfiles, int nrCap, struct checkUpProcedure &proc, std::string name)
{
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

  //!< Make a cell, the type of the cell depending on the value of 'cellType' #TODO  -> There is a cast.
  auto createCell = [&] {
    if (cellType == 0)
      return (Cell)Cell_KokamNMC(M, degid, verbose); //!< a high power NMC cell made by Kokam
    else if (cellType == 1)
      return (Cell)Cell_LGChemNMC(M, degid, verbose); //!< a high energy NMC cell made by LG Chem
    else
      return (Cell)slide::Cell_user(M, degid, verbose); //!< a user-defined cell
  };

  Cell c1 = createCell();

  //!< Make the cycler
  Cycler cycler(c1, name, verbose, timeCycleData);

  //!< Print a warning if you want to store cycling data
  //!< In profileAgeing, you are guaranteed to get one point per step in the profile
  //!< so if the steps are very short (e.g. 1sec), you are storing a huge amount of data (e.g. every second)
  //!< 	This seriously slowing down the code (taking hours instead of minutes) and produces huge amounts of data (several GB)
  if (timeCycleData != 0)
    std::cout << "Warning for profile ageing: the cycling data of the cell is going to be stored, which will lead to much slower "
                 "calculation (several hours) and a huge amount of data (several GB).\n";
  //!< Call the profile ageing function from the Cycler
  try {
    cycler.profileAgeing(profName, limit, Vma, Vmi, Ti, nrProfiles, nrCap, proc);
  } catch (int err) {
    //!< std::cout << "Throw test: " << 77 << '\n';
    std::cout << "Profile_one experienced error " << err << " during execution of " << name << ", abort this test.\n";
    if (err == 15) //!< #TODO this is not valid anymore.
    {
      std::cout << "Error 15 means that the cell had degraded too much to continue simulating.\n"
                   "This can be due to too much SEI growth, too much loss of lithium, too much loss of active material (thin electrodes, low volume fraction, or low effective surface)\n"
                   "too much surface cracking, too low diffusion constants, too high resistance, or too much lithium plating.\n"
                   "Continue simulating might lead to errors in the code so the simulation of this degradation test is stopped.\n"
                   "The results which have been written are all valid and you can ignore the error messages above which explained where this error comes from.\n";
    }
  }
}

void CycleAgeing(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, int verbose)
{
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
   * pref 		std::string with which the name of the subfolder in which the results should be written, will begin
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

  //!< *********************************************************** 1 variables ***********************************************************************

  //!< append the ageing identifiers to the prefix
  pref += +"_" + degid.print() + "_";

  //!< Make variables to describe the cycling regimes
  bool CVcha = true;      //!< we want to have a CC CV charge (if false, then charge has only a CC phase)
  double Ccutcha = 0.05;  //!< Crate of the cutoff current for the CV phase of the charge [-]
  bool CVdis = false;     //!< we want to have a CC discharge (if true, then discharge has both a CC and CV phase)
  double Ccutdis = 1.0;   //!< Crate of the cutoff current for the CV phase of the discharge [-]
  int nrCycles = 3000;    //!< the number of cycles which has to be simulated
  int nrCap = 500;        //!< the number of cycles between check-ups
  int timeCycleData = 60; //!< time interval at which cycling data (voltage and temperature) has to be recorded [s]
                          //!< 	0 means no data is recorded
                          //!<  if not 0, data is recorded approximately every so many seconds

  //!< *********************************************************** 2 check-up procedure ******************************************************************

  //!< Make a struct to describe the check-up procedure
  struct checkUpProcedure proc;
  proc.blockDegradation = true;                    //!< boolean indicating if degradation is accounted for during the check-up, [RECOMMENDED: TRUE]
  proc.capCheck = true;                            //!< boolean indicating if the capacity should be checked
  proc.OCVCheck = true;                            //!< boolean indicating if the half-cell OCV curves should be checked
  proc.CCCVCheck = true;                           //!< boolean indicating if some CCCV cycles should be done as part of the check-up procedure
  proc.pulseCheck = true;                          //!< boolean indicating if a pulse discharge test should be done as part of the check-up procedure
  proc.includeCycleData = true;                    //!< boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
  proc.nCycles = 3;                                //!< number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[0] = 0.5;                            //!< do a 0.5C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[1] = 1.0;                            //!< do a 1C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[2] = 2.0;                            //!< do a 2C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Ccut_cha = 0.05;                            //!< C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
  proc.Ccut_dis = 100;                             //!< C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
                                                   //!< if the cutoff is larger than the value in the CC phase, no CV phase is done
  proc.set_profileName("CheckupPulseProfile.csv"); //!< name of the csv file which contains the current profile for the pulse test
                                                   //	the first column contains the current in [A] (positive for discharge, negative for charge)
                                                   //	the second column contains the time in [sec] the current should be maintained
                                                   //	the profile must be a net discharge, i.e. sum (I*dt) > 0
  proc.profileLength = 13;                         //!< length of the current profiles for the pulse test (number of rows in the csv file)

  //!< *********************************************************** 3 simulations ******************************************************************

  //!< The cycling-window needs the voltage windows between which it should cycle the battery.
  //!< Users typically find it easier to think in terms of state of charge windows.
  //!< For the Kokam cell used here, the conversion is as follows (for other cells the user has to derive the conversion from the OCV curve)
  //!< 0%		2.7V
  //!< 10%		3.42V
  //!< 20%		3.49V
  //!< 50% 		3.65V
  //!< 80%		3.98V
  //!< 90%		4.08V
  //!< 100%		4.2V
  //!< The voltages will be slightly different for the LG Chem cell because it has different OCV curves

  //!< Use for loops to create a cycle ageing experiment configuration.
  std::vector<CycleAgeingConfig> cycleAgConfigVec;
  cycleAgConfigVec.reserve(25);

  //!< cycleAgConfigVec.emplace_back(4.08, 3.42, 45, 2, 1, 90, 10); //!< For debugging purposes.

  for (double Ccha : { 1, 2, 3 }) //!< Crate of the CC charge
  {
    double Cdis{ 1 }; //!< Crate of the CC discharge

    //!< Corresponds to SOC window 100% -- 0%
    double Vma{ 4.2 }, Vmi{ 2.7 }, SOCma{ 100 }, SOCmi{ 0 }; //!< maximum and minimum voltages of the cycle [V] with corresponding SOC values [%].
    for (double Tc : { 45, 25 })                             //!< 45 and 25 Celsius degrees of environmental temperature.
      cycleAgConfigVec.emplace_back(Vma, Vmi, Tc, Ccha, Cdis, SOCma, SOCmi);

    //!< Corresponds to SOC window 80% -- 0%
    Vma = 3.98, Vmi = 2.7, SOCma = 80, SOCmi = 0;
    for (double Tc : { 45 }) //!< 45 Celsius degrees of environmental temperature.
      cycleAgConfigVec.emplace_back(Vma, Vmi, Tc, Ccha, Cdis, SOCma, SOCmi);

    //!< Corresponds to SOC window 100% -- 20%
    Vma = 4.2, Vmi = 3.49, SOCma = 100, SOCmi = 20;
    for (double Tc : { 45 }) //!< 45 Celsius degrees of environmental temperature.
      cycleAgConfigVec.emplace_back(Vma, Vmi, Tc, Ccha, Cdis, SOCma, SOCmi);

    //!< Corresponds to SOC window 90% -- 10%
    Vma = 4.08, Vmi = 3.42, SOCma = 90, SOCmi = 10;
    for (double Tc : { 45, 25, 5 }) //!< 45, 25, and 5 Celsius degrees of environmental temperature.
      cycleAgConfigVec.emplace_back(Vma, Vmi, Tc, Ccha, Cdis, SOCma, SOCmi);
  }

  auto task_indv = [&](int i_begin) {
    //!< simulate one cycle ageing experiment
    Cycle_one(M, degid, cellType, verbose, cycleAgConfigVec[i_begin], CVcha, Ccutcha, CVdis, Ccutdis, timeCycleData, nrCycles, nrCap, proc, pref);
  };

  //!< Print a message that we are starting the simulations
  std::cout << "\t Cycle ageing experiments are started.\n";
  slide::run(task_indv, cycleAgConfigVec.size()); //!< Runs individual simulation in parallel or sequential depending on settings.
}

void CalendarAgeing(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, int verbose)
{
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
   * pref 		std::string with which the name of the subfolder in which the results should be written, will begin
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

  //!< *********************************************************** 1 variables ***********************************************************************

  //!< append the ageing identifiers to the prefix
  pref += "_" + degid.print() + "_";

  //!< Make variables to describe the cycling regimes
  constexpr int mode = 0;             //!< integer deciding how often to recharge the cells (due to degradation, the voltage will decrease over time. But no self-discharge is simulated)
                                      //!< 0 means we recharge to the specified voltage only after a check-up
                                      //!< 1 means we recharge to the specified voltage every day
                                      //!< 2 means we float the cell at a constant voltage, rather than resting it
                                      //!< 		floating at constant voltage takes very long to simulate
  constexpr int Time = 20 * 30;       //!< time the cell has to rest [days]
  constexpr int timeCheck = 30;       //!< time between consecutive check-ups [days]
  constexpr int timeCycleData = 3600; //!< time interval at which cycling data (voltage and temperature) has to be recorded [s]
                                      //!< 	0 means no data is recorded
                                      //!<  if not 0, data is recorded approximately every so many seconds

  //!< *********************************************************** 2 check-up procedure ******************************************************************

  //!< Make a struct to describe the check-up procedure
  struct checkUpProcedure proc;
  proc.blockDegradation = true;                    //!< boolean indicating if degradation is accounted for during the check-up, [RECOMMENDED: TRUE]
  proc.capCheck = true;                            //!< boolean indicating if the capacity should be checked
  proc.OCVCheck = true;                            //!< boolean indicating if the half-cell OCV curves should be checked
  proc.CCCVCheck = true;                           //!< boolean indicating if some CCCV cycles should be done as part of the check-up procedure
  proc.pulseCheck = true;                          //!< boolean indicating if a pulse discharge test should be done as part of the check-up procedure
  proc.includeCycleData = true;                    //!< boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
  proc.nCycles = 3;                                //!< number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates), maximum 100
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[0] = 0.5;                            //!< do a 0.5C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[1] = 1.0;                            //!< do a 1C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[2] = 2.0;                            //!< do a 2C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Ccut_cha = 0.05;                            //!< C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
  proc.Ccut_dis = 100;                             //!< C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
                                                   //!< if the cutoff is larger than the value in the CC phase, no CV phase is done
  proc.set_profileName("CheckupPulseProfile.csv"); //!< name of the csv file which contains the current profile for the pulse test
                                                   //	the first column contains the current in [A] (positive for discharge, negative for charge)
                                                   //	the second column contains the time in [sec] the current should be maintained
                                                   //	the profile must be a net discharge, i.e. sum (I*dt) > 0
  proc.profileLength = 13;                         //!< length of the current profiles for the pulse test (number of rows in the csv file)

  //!< *********************************************************** Simulations ******************************************************************

  //!< The cycling-window needs the voltage windows between which it should cycle the battery.
  //!< Users typically find it easier to think in terms of state of charge windows.
  //!< For the Kokam cell used here, the conversion is as follows (for other cells the user has to derive the conversion from the OCV curve)
  //!< 0%		2.7V
  //!< 10%		3.42V
  //!< 20%		3.49V
  //!< 50% 		3.65V
  //!< 80%		3.98V
  //!< 90%		4.08V
  //!< 100%		4.2V
  //!< The voltages will be slightly different for the LG Chem cell because it has different OCV curves

  std::vector<CalendarAgeingConfig> calAgConfig;
  std::array<double, 3> V_arr{ 4.2, 4.08, 3.65 }, SOC_arr{ 100, 90, 50 };

  for (double Tc : { 5, 25, 45 })             //!< temperature at which the cell has to rest [oC]
    for (size_t i = 0; i < V_arr.size(); i++) //!< voltage at which the cell has to rest [V] above the minimum and below the maximum voltage of the cell
      calAgConfig.emplace_back(V_arr[i], Tc, SOC_arr[i]);

  auto task_indv = [&](size_t i) {
    Calendar_one(M, degid, cellType, verbose, calAgConfig[i].V, calAgConfig[i].Ti(), Time, mode, timeCycleData, timeCheck, proc, calAgConfig[i].get_name(pref));
  };

  //!< Print a message that we are starting the simulations
  std::cout << "\t Calendar ageing experiments are started.\n";
  slide::run(task_indv, calAgConfig.size());
}

void ProfileAgeing(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, int verbose)
{
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
   * pref 		std::string with which the name of the subfolder in which the results should be written, will begin
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

  //!< *********************************************************** 1 variables ***********************************************************************

  //!< append the ageing identifiers to the prefix
  pref += "_" + degid.print() + "_";

  //!< find the number of parallel threads that are optimal to use
  //!< unsigned int Ncor = std::thread::hardware_concurrency(); //!< Ncor is the number of (logical) cores, 0 if c++ can't identify it

  //!< Give the current profile which should be followed
  //!< the first column must give the current [A] which should be followed
  //!< 		>0 is discharge
  //		<0 is charge
  //!< the second column should give the time [sec] for which this current should be maintained
  //		values for time are floored to the integer value (ie. the numbers behind the comma are ignored)
  //!< 5 example current profiles are provided:
  //!< 	name										length		description
  //!< 	Current Profile random.csv 					100			a random current profile with a maximum current of 5A, each current step takes maximum 1000 seconds
  //	Current Profile drive cycle HWFET.csv		766			HWFET drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
  //	Current Profile drive cycle NYCC.csv		599			NYCC drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
  //	Current Profile drive cycle UDDS.csv		1370		UDDS drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
  //!< 	Current Profile drive cycle US06.csv		601			US06 drive cycle with a maximum current of 8.1A (= 3C), each current step takes 1 second
  int length;    //!< length of the profile (number of rows)
  int limit = 0; //!< what to do if the voltage limits are reached while following the profile:
                 //!< 0 immediately go to the next current step of the profile (i.e. reduce the time of this step)
                 //!< 1 keep the voltage constant for the rest of this step of the profile (i.e. reduce the current of this step)

  int nrProfiles = 10000; //!< number of times the current profile should be repeated [-]
  int nrCap = 1000;       //!< number of times the current profile should be repeated between consecutive check-ups [-]
  int timeCycleData = 0;  //!< time interval at which cycling data (voltage and temperature) has to be recorded [s]
                          //!< 	0 means no data is recorded
                          //!<  if not 0, data is recorded approximately every so many seconds
                          //!< For this function (profileAgeing) it is highly recommended to keep timeCycleData at 0.
                          //!< 	Otherwise, a huge amount of data is generated (many GB) and consequently the code slows down dramatically
                          //!< 	This is because you are guaranteed to get one data point per step in the current profile if timeCycleData is bigger than 0,
                          //!< 	independent of the time interval at which you are otherwise collecting data (i.e. even if you set it to 100, you will still get a data point every step)
                          //!< 	As many drive cycles will have a step per second (i.e. the current changes every second), this means you store a data point every second (even if timeCycleData = 100)
                          //!< 	which is obviously a huge amount of data if you simulate long term battery usage

  //!< *********************************************************** 2 check-up procedure ******************************************************************

  //!< Make a struct to describe the check-up procedure
  struct checkUpProcedure proc;
  proc.blockDegradation = true;                    //!< boolean indicating if degradation is accounted for during the check-up, [RECOMMENDED: TRUE]
  proc.capCheck = true;                            //!< boolean indicating if the capacity should be checked
  proc.OCVCheck = true;                            //!< boolean indicating if the half-cell OCV curves should be checked
  proc.CCCVCheck = true;                           //!< boolean indicating if some CCCV cycles should be done as part of the check-up procedure
  proc.pulseCheck = true;                          //!< boolean indicating if a pulse discharge test should be done as part of the check-up procedure
  proc.includeCycleData = true;                    //!< boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
  proc.nCycles = 3;                                //!< number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates).
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[0] = 0.5;                            //!< do a 0.5C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[1] = 1.0;                            //!< do a 1C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Crates[2] = 2.0;                            //!< do a 2C cycle as part of the CCCV check-up, must be positive
                                                   //!< If you change this variable, also change it in the MATLAB script which reads the results from the check-up (the variable Crates in readAgeing_CCCV.m)
  proc.Ccut_cha = 0.05;                            //!< C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
  proc.Ccut_dis = 100;                             //!< C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
                                                   //!< if the cutoff is larger than the value in the CC phase, no CV phase is done
  proc.set_profileName("CheckupPulseProfile.csv"); //!< name of the csv file which contains the current profile for the pulse test
                                                   //	the first column contains the current in [A] (positive for discharge, negative for charge)
                                                   //	the second column contains the time in [sec] the current should be maintained
                                                   //	the profile must be a net discharge, i.e. sum (I*dt) > 0
  proc.profileLength = 13;                         //!< length of the current profiles for the pulse test (number of rows in the csv file)

  //!< *********************************************************** 3 simulations ******************************************************************

  //!< The cycling-window needs the voltage windows between which it should cycle the battery.
  //!< Users typically find it easier to think in terms of state of charge windows.
  //!< For the Kokam cell used here, the conversion is as follows (for other cells the user has to derive the conversion from the OCV curve)
  //!< 0%		2.7V
  //!< 10%		3.42V
  //!< 20%		3.49V
  //!< 50% 		3.65V
  //!< 80%		3.98V
  //!< 90%		4.08V
  //!< 100%		4.2V
  //!< The voltages will be slightly different for the LG Chem cell because it has different OCV curves

  //!< Use for loops to create a profile ageing experiment configuration.
  std::vector<ProfileAgeingConfig> profAgConfigVec;
  profAgConfigVec.reserve(25);

  //!< names of the csv file with the current profile:
  std::array<std::string, 4> profile_arr = { "Current Profile drive cycle HWFET.csv",
                                             "Current Profile drive cycle NYCC.csv",
                                             "Current Profile drive cycle UDDS.csv",
                                             "Current Profile drive cycle US06.csv" };

  //!< identification strings for each experiment (will also be used to name the folder so must have the same restrictions as the prefix, i.e. no spaces, no special characters, etc.)
  std::array<std::string, 4> prefName_arr = { "prof-HWFET", "prof-NYCC", "prof-UDDS", "prof-US06" }; //!< Should be same size as profile_arr.

  for (size_t i{}; i < profile_arr.size(); i++) //!< Profile file names.
  {
    //!< Corresponds to SOC window 100% -- 0%
    double Vma{ 4.2 }, Vmi{ 2.7 }, SOCma{ 100 }, SOCmi{ 0 }; //!< maximum and minimum voltages which the cell should stay while following the current profile [V]
                                                             //!< with corresponding SOC values [%].
    for (double Tc : { 45, 25 })                             //!< 45 and 25 Celsius degrees of environmental temperature while the cell is following the profile [oC]
      profAgConfigVec.emplace_back(Vma, Vmi, Tc, SOCma, SOCmi, profile_arr[i], prefName_arr[i]);

    //!< Corresponds to SOC window 90% -- 10%
    Vma = 4.08, Vmi = 3.42, SOCma = 90, SOCmi = 10;
    for (double Tc : { 25 }) //!< 25 Celsius degree of environmental temperature.
      profAgConfigVec.emplace_back(Vma, Vmi, Tc, SOCma, SOCmi, profile_arr[i], prefName_arr[i]);
  }

  auto task_indv = [&](size_t i) {
    //!< simulate one profile ageing experiment
    auto &conf = profAgConfigVec[i];
    Profile_one(M, degid, cellType, verbose, conf.csvName, length, limit, conf.Vma, conf.Vmi, conf.Ti(), timeCycleData, nrProfiles, nrCap, proc, conf.get_name(pref));
  };

  //!< Print a message that we are starting the simulations
  std::cout << "\t Profile ageing experiments are started.\n";
  slide::run(task_indv, profAgConfigVec.size());
}