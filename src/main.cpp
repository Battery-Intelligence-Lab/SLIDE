/*
 * main.cpp
 *
 * This file implements the main-function.
 * Here, the user has to choose what to simulate by commenting things in (by placing \\ at the start of the line) and out (by removing the \\)
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */


// Include header files
#include "determineOCV.h"
#include "determineCharacterisation.h"
#include "Cycling.h"
#include "Degradation.h"

using namespace std;

int main(int argv, char* argc[]){
	/*
	 * This function decides what will be simulated.
	 * In the code below, you have to uncomment the line which you want to execute (to uncomment, remove the // in front of the line)
	 * and comment all the other lines in the blocks with FUNCTION CALLS (to comment in, add // in front of the lines)
	 */

	// print that you start simulations
	cout<<"start simulations "<<endl;

	// Make a clock to measure how long the simulation takes
	std::clock_t tstart;
	double duration;
	tstart = std::clock();

	// Read the values for the matrices for the spatial discretisation of the solid diffusion PDE.
	// these values are calculated by the Matlab-script 'modelSetup.m', which writes them to csv files.
	// Here, we only need to read these csv files
	Model M;
	Model_initialise(M);

	// Choose how much messages should be printed to the terminal
	int verbose = 1;				// integer deciding how verbose the simulation should be
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

	// Choose the prefix, a string with which the names of all subfolders with the outputfiles will start (the folders will be called prefix_xxxx).
	// This is typically a number, and the reason why we use it is that results of previous simulations are not overwritten.
	// i.e. for one simulation, use pref "0", then the next simulation could have pref "1".
	// this ensures the results of the second simulation are written in different files and folders, rather than overwriting the files from the first simulation
	string pref = "0";//todo				// will be used to name a folder with the results, so must comply with the naming convention on your operating system
									// letters, numbers, underscores and dots are fine
									// don't use spaces
									// don't use forbidden characters, e.g. for Windows the following are not allowed < > : " \ | / ? *
									// https://docs.microsoft.com/en-us/windows/desktop/fileio/naming-a-file

	// Choose which cell should be used for the simulations
	int cellType = 0; 				// which cell to use for the simulation.
									// 0 	high power Kokam NMC cell (18650)
									// 1 	high energy LG Chem NMC cell (18650)
	   								// 2 	user cell (template class provided for where the user can define his own parameters)

	// Choose which degradation models to use for the simulations
	// The models are identified by a number, see DEG_ID in Cell.hpp and below
		/* SEI_id 		identification for which model(s) to use for SEI growth
		 * 					0	no SEI growth
		 * 					1 	kinetic model only (Tafel kinetics)
		 * 							ref: Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004
		 * 					2 	Pinson&Bazant model: linear diffusion + Tafel kinetics
		 * 							ref: Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
		 * 					3	Christensen and Newman model: diffusion and kinetics model
		 * 							ref: Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
		 * SEI_porosity	integer deciding whether we reduce the volume fraction of active electrode material due to SEI growth
		 *  				0	don't subtract anything
		 * 					1	use correlation from Ashwin
		 * 							ref: Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
		 * CS_id	identification for which model(s) to use for surface crack growth
		 *  				0 	no surface cracking
		 * 					1 	Laresgoiti's stress + crack growth model
		 * 							ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
		 * 					2 	Dai stress model + Laresgoiti crack growth
		 * 							ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
		 * 								 Dai, Cai, White, Journal of Power sources 247, 2014
		 * 					3 	model based on Deshpande and Bernardi
		 * 							ref: Deshpande & Bernardi,Journal of the Electrochemical Society 164 (2), 2017
		 * 					4 	model from Barai
		 * 							ref: Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
		 * 					5 	model from Ekstrom
		 * 							ref: Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
		 * CS_diffusion	integer deciding whether to reduce the negative diffusion constant due to surface cracks
		 *  				0 	don't decrease diffusion
		 * 					1	decrease according to Barai
		 * 							ref: Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
		 * LAM_id  identification for which model(s) to use for loss of active material
		 *  				0 	no LAM
		 * 					1	Dai's stress model and Laresgoiti's correlation to relate the stress to LAM
		 * 							ref: Dai, Cai, White, Journal of Power sources 247, 2014
		 * 								 Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
		 * 					2	delacourt's correlation between abs(current density) and porosity
		 * 							ref: Delacourt & Safari, Journal of the Electrochemical Society 159 (8), 2012
		 * 					3 	Kindermann's model for cathode dissolution: tafel kinetics for increased porosity
		 * 							ref: Kindermann, Keil, Frank, Jossen, Journal of the Electrochemical Society 164 (12), 2017
		 * 					4 	Narayanrao's correlation between effective surface area and current density
		 * 							ref: Narayanrao, Joglekar, Inguva, Journal of the Electrochemical Society 160 (1), 2012
		 * pl_id 	identification for which model to use for li-plating
		 *  				0 	no plating
		 * 					1	Yang et al's model for thermodynamic plating (Tafel kinetics)
		 * 							ref: Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
		 */
	// For now, assume we want to include 'Pinson&Bazant'-SEI growth, 'Delacourt'-LAM, 'Kindermann'-LAM and 'Yang'-lithium plating
	DEG_ID deg;
		deg.SEI_n = 1;										// there is 1 SEI model
		deg.SEI_id[0] = 2;									// Pinson & Bazant SEI growth
		deg.SEI_porosity = 0;								// don't decrease the porosity (set to 1 if you do want to decrease the porosity)

		deg.CS_n = 1;										// there is 1 model (that there are no cracks)
		deg.CS_id[0] = 0;									// no surface cracks
		deg.CS_diffusion = 0;								// don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

		deg.LAM_n = 2;										// there are 2 LAM models (Delacourt and Kindermann)
		deg.LAM_id[0] = 2;									// Delacourt LAM
		deg.LAM_id[1] = 3;									// Kindermann LAM

		deg.pl_id = 1;										// Yang litihium plating

	// Then the user has to choose what is simulated.
	// In the code below, uncomment the line which calls the function you want to execute (uncommenting means removing the // in front of the line)
	// and comment all the other lines (commenting means putting // in front of the line)

	// *********************************************** PARAMETRISATION FUNCTION CALLS *********************************************************************

//	estimateOCVparameters();								// OCV parametrisation
//	estimateCharacterisation(); 							// parametrisation of diffusion constant, rate constant and DC resistance

	// *********************************************** CYCLING FUNCTION CALLS ********************************************************
	CCCV(M, pref, deg, cellType, verbose);				// a cell does a few CCCV cycles
//	FollowCurrent(M, pref, deg, cellType, verbose); 		// a cell follows the current profile specified in a csv file

	// *********************************************** DEGRADATION FUNCTION CALLS ********************************************************

//	CalendarAgeig(M, pref, deg,cellType, verbose);		// simulates a bunch of calendar degradation experiments
//	CycleAgeing(M, pref, deg, cellType, verbose);			// simulates a bunch of cycle degradation experiments
//	ProfileAgeing(M, pref, deg, cellType, verbose);		// simulates a bunch of drive cycle degradation experiments

	// *********************************************** END ********************************************************
	// Now all the simulations have finished. Print this message, as well as how long it took to do the simulations
	duration = ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC;
	cout<<"finished all simulations in "<<floor(duration/60)<<":"<<duration - floor(duration/60)*60<<endl;
}
