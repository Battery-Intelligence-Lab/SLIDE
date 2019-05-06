/*
 * Cycler.cpp
 *
 * Class to implement a 'cycler', which extends a 'basic cycler'.
 * A Cycler consists of the programs one would write on battery cycles.
 * I.e. it defines the procedures for loading a cell in a degradation experiment.
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Cycler.hpp"
#include "Interpolation.h"
#include "ReadCSVfiles.h"

using namespace std;

Cycler::Cycler(Cell& ci, string IDi, int verbosei, int feedbacki)
: BasicCycler(ci, IDi, verbosei, feedbacki)  {
	/*
	 * Constructor of a Cycler.
	 *
	 * IN
	 * ci			cell connected to the cycler (call by reference, i.e. the cell given to the cycler will be changed during execution)
	 * IDi 			identification string of for this cycler.
	 * 				A new folder is made to store the data of this cell
	 * 				therefore, IDi must obey all restrictions for naming folders on your operating systems
	 * 				e.g. for Windows, the following are not allowed < > : " \ | / ? *
	 * 				https://docs.microsoft.com/en-us/windows/desktop/fileio/naming-a-file
	 * 				spaces are also not recommended
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
	 * CyclingDataTimeIntervali 	time interval at which the cycling data must be stored, >= 0 [s]
	 * 				if 0, no data is stored
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::Cycler is starting"<<endl;

	indexdegr = 0;						// index number of the check-up (how many check-ups have we done so far)

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::Cycler is terminating"<<endl;
}

Cycler::~Cycler() {
}

double Cycler::getCapacity(bool blockDegradation){
	/*
	 * Function to measure the capacity of a cell.
	 * It does a full charge / discharge cycle (C/25, CV until a current threshold of 0.005C) to the voltage limits of the cell.
	 * The capacity of the cell is the charge that can be extracted during the full discharge.
	 * After the capacity is measured, the original states are restored
	 *
	 * IN
	 * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this capacity measurement,
	 * 							and at the end of the capacity measurement, the exact battery state from the start of the capacity measurement is restored
	 * 						if false [NOT RECOMMENDED], degradation is accounted for during this capacity measurement,
	 * 							and at the end of the capacity measurement, the battery state is left as it is
	 * OUT
	 * cap 					capacity which could be discharged between the voltage limits of the cell [Ah]
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::getCapacity is starting"<<endl;


	// *********************************************************** 1 variables & settings ***********************************************************************

	// variables
	State s;																			// initial battery state
	double Iold;																		// initial battery current
	double dt = 2;																		// use a 2 second time step for accuracy (probably 5 would be fine as well)
	double crate = 1.0/25; 																// C rate we will be using
	double Ccut = 0.005;																// C rate of the cut-off current during CV charge
	double I = crate*c.getNominalCap(); 												// charging & discharging current during the CC part [A]
	double Icut = Ccut*c.getNominalCap();												// cut-off current during the CV discharge [A]
	double CCcap, CVcap, Cen, Ven;														// charge and energy capacity during CC and CV part
	double timei;																		// time spent during the cc and cv part

	// copy all the states to restore them later
	double Tenvi, Trefi;
	c.getStates(s, &Iold);

	// Set the cell temperature to the reference temperature
	c.getTemperatures(&Tenvi, &Trefi);
	c.setTenv(Trefi);
	c.setT(Trefi);

	// *********************************************************** 2 full charge / discharge cycle ***********************************************************************

	// fully charge and discharge the cell
	try{
		// fully charge battery to its specified maximum voltage (Cell::getVmax())
		if(verbose >= printCyclerHighLevel)												// print high level flow of the simulation if desired
			cout<<"Cycler::getCapacity going to fully charge the cell"<<endl;
		CC_V_CV_I(crate, c.getVmax(), Ccut, dt, blockDegradation, &CCcap, &Cen, &timei );// do a CC charge until the maximum voltage, followed by a CV at this voltage

		// fully discharge battery in separate CC and CV phases so we can get the CC and CV capacity separately (in the future we might be interested in them)
		if(verbose >= printCyclerHighLevel)												// print high level flow of the simulation if desired
			cout<<"Cycler::getCapacity going to CC discharge the cell"<<endl;
		CC_V(I, dt, blockDegradation, c.getVmin(), &CCcap, &Cen, &timei);				// discharge at a constant current until the minimum voltage is reached
		if(verbose >= printCyclerHighLevel)												// print high level flow of the simulation if desired
			cout<<"Cycler::getCapacity going to CV discharge the cell"<<endl;
		CV_I(c.getVmin(), dt, blockDegradation, Icut, &CVcap, &Ven, &timei);			// do a CV discharge until the current is below the cut-off current
	}
	catch(int e){
		if (verbose >= printCrit)
			cout<<"Error in a subfunction of Cycler::getCapacity "<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// *********************************************************** 3 output variables ***********************************************************************

	// restore the original battery state if we are not accounting for degradation
	if(blockDegradation)
		c.setStates(s, Iold);

	// always reset the original environmental temperature
	c.setTenv(Tenvi);

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::getCapacity is terminating with a capacity of "<<CCcap + CVcap<<"Ah"<<endl;

	// Return the capacity (for now we don't care about the CC or CV capacity so just return the total capacity
	return CCcap + CVcap;																// capacity discharged during the CC and CV discharge
}

void Cycler::getOCV(int nin, double Ah[], double OCVp[], double OCVn[], int* nout){
	/*
	 * function to calculate the half cell OCV curves.
	 * Each curve is fully simulated (i.e. from a lithium fraction of 0 to a lithium fraction of 1).
	 * Both curves are aligned as they are in the actual cell, and the common x-axis gives the discharged charge [Ah].
	 * The discharged charge is '0' when the cell is charged to its maximum voltage, and negative for all points with a higher cell voltage.
	 * When one electrode has reached the full lithium concentration (but the other electrode can still further (dis)charge, its potential is kept constant.
	 *
	 * Degradation is never accounted for while recording the OCV curves because the code does not support this
	 * After the OCV curves are recorded, the original states are restored
	 *
	 * IN
	 * nin 		length of the arrays provided for feedback, i.e. length(OCVp)
	 *
	 * OUT
	 * Ah 		discharged charge, 0 when the cell OCV is at its maximum voltage, negative when the cell OCV is above its maximum [Ah]
	 * OCVp		open circuit potential of the positive electrode, constant after we have reached the extreme points of the cathode OCV curve [V]
	 * OCVn		open circuit potential of the negative electrode, constant after we have reached the extreme points of the cathode OCV curve [V]
	 * nout 	total length of the output arrays
	 *
	 * THROWS
	 * 1002		the arrays provided for feedback are too short (ninocv < nout)
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::getOCV is starting"<<endl;

	// *********************************************************** 1 variables & settings ***********************************************************************

	// variables
	State sini, sCharged;																	// initial & fully charged battery state
	double Iini, ICharged;																	// initial & fully charged battery current
	double Tenvi, Trefi;																	// initial environmental & reference temperatures
	double Crate = 0.04;																	// Crate for the CC phase.
	double Ccutcha = 0.005;																	// C rate of the cut-off current during CV charge
	double Icha = -c.getNominalCap()*Crate; 												// charging current [A]
	double Idis = c.getNominalCap()*Crate; 													// discharging current [A]
	double dt = 2;																			// time step for time integration
	bool blockDegradation = true;															// don't account for degradation while we measure the OCV curves
	double ah, wh, tt;																		// unneeded feedback variables
	double ah_neg_p;																		// charge on the cathode OCV curve with a negative x-value (i.e. the points with a higher OCVp than the OCVp at the fully charged condition)
	double ah_neg_n;																		// charge on the anode OCV curve with a negative x-value (i.e. the points with a lower OCVn than the OCVn at the fully charged condition)
	bool bound = false;																		// in linear interpolation return the closest point if you are out of range
																							// 	this ensures the electrode OCV remains constant after it has reached the full lithium concentration

	// Data collection: we don't want to store a point every time step because that leads to very large arrays.
	int nStore = 750;																		// store a value every 750 time steps
	int nin2 = 1/Crate*3600/dt*2; 															// expected duration 1/Crate [hours] * 3600/dt [time steps], *2 such that the arrays are for sure long enough
	double OCVni[nin2], OCVpi[nin2], Ahpi[nin2], Ahni[nin2]; 								// make the arrays to store the half cell OCV curves
	int np, nn;																				// actual number of data points in the full discharge of each electrode

	// Store the initial battery states so we can restore them later
	c.getStates(sini, &Iini);

	// Set the cell and environmental temperature to the reference temperature
	c.getTemperatures(&Tenvi, &Trefi);
	c.setTenv(Trefi);
	c.setT(Trefi);

	// *********************************************************** 2 get the half-cell curves ***********************************************************************

	try{
		// fully charge battery to its specified maximum voltage (Cell::getVmax())
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::getOCV is fully charging the cell to the maximum cell voltage"<<endl;
		CC_V_CV_I(Crate, c.getVmax(), Ccutcha, dt, blockDegradation, &ah, &wh, &tt); 		// do a CC charge until the maximum voltage, followed by a CV at this voltage

		// Store this battery state
		c.getStates(sCharged, &ICharged);

		// We can cycle each electrode separately using the function CC_halfCell_full.
		// This function returns the electrode OCV (the resistive voltage drop is just ignored)
		// and the electrode is cycled until the lithium is exhausted (i.e. li-concentration is 0 or the maximum concentration)

		// fully cycle the cathode until it reaches its extreme lithium concentrations
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::getOCV is fully charging the cathode"<<endl;
		CC_halfCell_full(Icha, dt, true, nin2, OCVpi, &ah_neg_p, &np);						// ah_neg_p contains the capacity we charged additional (starting from the 'fully charged' state) on the positive electrode
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::getOCV is fully discharging the cathode"<<endl;
		CC_halfCell_full(Idis, dt, true, nin2, OCVpi, &ah, &np); 							// now OCVpi contains the cathode OCV for the full lithium concentration range (from 0 to the maximum concentration)

		// restore the fully charged battery state (undo the cycling of the cathode)
		c.setStates(sCharged, ICharged);

		// fully cycle the anode until it reaches its extreme lithium concentrations
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::getOCV is fully charging the anode"<<endl;
		CC_halfCell_full(Icha, dt, false, nin2, OCVni, & ah_neg_n, &nn);					// ah_neg_n contains the capacity we charged additional (starting from the 'fully charged' state) on the negative electrode
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::getOCV is fully discharging the anode"<<endl;
		CC_halfCell_full(Idis, dt, false, nin2, OCVni, & ah, &nn); 							// now OCVni contains the anode OCV for the full lithium concentration range (from maximum concentration to 0)
	}
	catch(int e){
		if (verbose >= printCrit)
			cout<<"Error in a subfunction of Cycler::getOCV when half-cycling the cells"<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// *********************************************************** 3 make a common x-axis ***********************************************************************
	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::getOCV is making the x-axis"<<endl;
	// the x-axis we want gives discharged charge [Ah]
	// We want the discharge charge to be 0 when the cell OCV is 4.2V (OCVp - OCVn = 4.2) [where 4.2 is the maximum cell voltage]
	// The points 'to the left', i.e. with a higher value for OCVp - OCVn must have a negative discharged charge
	// The points 'to the right', i.e. with a lower value for OCVp - OCVn must have a positive discharged charge

	// make the x-axis for the cathode
	// 		the variable 'ah_neg_p' has the charge we charged additionally starting from the 'fully charged' position (i.e. how far we went 'to the left')
	// 		so the 'leftmost' point of the cathode (with the highest OCVp) must have an x-value of ah_neg_p (which will be negative)
	// 		then the point on the OCVp curve where the cell voltage is 4.2V will have an x-value of 0
	// From that leftmost point, the x-axis gives the discharged charge in Ah.
	// 		the current during the discharge was constant (the variable Idis), so we just need to multiply the current by the time until now
	// 		every time step takes 'dt' seconds, so the time [h] in step i is given by 'i*dt/3600'
	for(int i=0;i<np;i++)
		Ahpi[i] = ah_neg_p + i*dt/3600.0*Idis;

	// make the x-axis for the anode, exactly the same as for the cathode
	for(int i=0;i<nn;i++)
		Ahni[i] = ah_neg_n + i*dt/3600.0*Idis;

	// Then we need to align both half-cell curves on a common x-axis
	// 		the most negative point on the x-axis must be the leftmost point of both curves: xmin = min (Ahpi[0], Ahni[0])
	// 		the most positive point on the x-axis must be the rightmost point of both curves: xmax = max (Ahpi[np-1], Ahni[nn-1])
	// 			where np is the number of points on the cathode, i.e. Ahpi[np-1] is the 'rightmost' point of the cathode OCV. Similar for the anode
	double Ahmin = min(Ahni[0], Ahpi[0]);													// the most negative point of the half-cell curves
	double Ahmax = max(Ahni[nn-1], Ahpi[np-1]);												// the most positive point of the half-cell curves

	// between these two extreme points, we make a data point every 'nstore' points.
	// 		i.e. the second point on the x-axis is Ahmin + nstore*dt/3600*Idis
	// 			 the thirds point on the x-axis is Ahmin + 2* nstore*dt/3600*Idis
	// 			 etc.
	// 		and we stop when the discharged charge is larger than Ahmax
	bool reachedEnd = false;																// boolean indicating if we have right the 'rightmost' point
	int nfull = 0;																			// number of data points on the common x-axis
	for(int i=0;i<nin;i++){
		Ah[i] = Ahmin + i*nStore*dt/3600.0*Idis;											// discharged charge
		if (Ah[i] > Ahmax){																	// stop if we've reached the most positive point of the OCV curve
			reachedEnd = true;																// we've reached the rightmost points
			nfull = i+1;																	// the number of points we've added (+1 because i starts at 0)
			break;
		}
	}

	// Check the array provided for output was long enough
	if (!reachedEnd){																		// we've reached the end of the loop but not the rightmost point
		int nguess = (Ahmax - Ahmin) / (nStore*dt/3600.0*Idis);								// the estimated number of points needed: the total Ah from left to right divided by the charge done every step
		cerr<<"ERROR in Cycler::getOCV, the arrays provided for input are too short. Their length was "<<nin<<" and the estimated required length is "<<nguess<<endl<<flush;
		throw 1002;
	}

	// *********************************************************** 4 align the half-cell curves on the common x-axis ***********************************************************************
	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::getOCV is aligning the half-cell curves"<<endl;
	// Now we have points of the electrodes each with their own x-axis
	// We can use linear interpolation to get the points on the common x-axis
	double ocvn, ocvp;																		// electrode OCVs for this data point on the common x-axis
	try{
		for(int i=0;i<nfull;i++){
			ocvp = linInt(verbose >= printCrit, bound, Ahpi, OCVpi, np, Ah[i]); 			// interpolate the negative half cell curve at the x-point we're at
			ocvn = linInt(verbose >= printCrit, bound, Ahni, OCVni, nn, Ah[i]);				// interpolate the positive half-cell curve at the x-point we're at
			OCVn[i] = ocvn;
			OCVp[i] = ocvp;
		}
	}
	catch(int e){
		if (verbose >= printCrit)
			cout<<"Error in a subfunction of Cycler::getOCV when getting the half-cell curves to the common x-axis"<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// *********************************************************** 5 output variables ***********************************************************************
	// restore the original battery state
	c.setStates(sini, Iini);																// reset the states and current
	c.setTenv(Tenvi);																		// reset the environmental temperature

	*nout = nfull;																		// return the number of data points in the arrays

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::getOCV terminating"<<endl;
}


double Cycler::checkUp_batteryStates(bool blockDegradation, bool checkCap, int cumCycle, double cumTime, double cumAh, double cumWh){
	/*
	 * Function to do a capacity check and write the battery states as part of a check-up.
	 * It will add one row of data in the csv file with the results (DegradationData_batteryState.csv in the subfolder of this Cycler)
	 * the row has the following entries:
	 * 		number of cycles until now
	 * 		time the cell has been cycled until now [h]
	 * 		cumulative Ah throughput up to now [Ah]
	 * 		cumulative Wh throughput up to now [Wh]
	 * 		the remaining capacity of the cell [Ah]
	 * 		the transformed li concentration at the positive inner nodes of the positive particle (nch values)
	 * 		the transformed li concentration at the positive inner nodes of the negative particle (nch values)
	 * 		the cell temperature [K]
	 * 		the thickness of the SEI layer [m]
	 * 		the lost lithium [As]
	 * 		the thickness of the cathode [m]
	 * 		the thickness of the anode [m]
	 * 		the volume fraction of active material in the cathode [-]
	 * 		the volume fraction of active material in the anode [-]
	 * 		the effective surface area of the cathode [m2 m-3]
	 * 		the effective surface area of the anode [m2 m-3]
	 * 		the surface area of the cracks at the surface of the negative particle [m2]
	 * 		the diffusion constant at reference temperature of the cathode [m s-1]
	 * 		the diffusion constant at reference temperature of the anode [m s-1]
	 * 		the specific resistance of the combined electrodes [Ohm m2]
	 * 		the thickness of the plated lithium layer [m]
	 * 		the DC resistance of the cell [Ohm]
	 * 		the active surface area of the anode [m2]
	 *
	 * IN
	 * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
	 * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
	 * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
	 * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
	 * checkCap 			if true, the capacity is measured.
	 * 						if false, the capacity is set to 0 and only the battery states are written
	 * cumCycle				number of cycles up to now [-]
	 * cumTime				time this cell has been cycled up to now [hour]
	 * cumAh				cumulative Ah throughput up to now [Ah]
	 * cumWh				cumulative Wh throughput up to now [Wh]
	 *
	 * OUT
	 * cap 					the remaining capacity of the cell [Ah]
	 *
	 * THROWS
	 * 1001 				the file in which to write the results couldn't be opened
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_batteryStates is starting"<<endl;

	// variables
	string nameState = "DegradationData_batteryState.csv";							// names of the csv file in which the results of this part of the check-up will be written
	State si;																		// initial battery state
	double Ii;																		// initial battery current [A]
	double states[ns];																// array with the battery's state

	// Get the capacity and the states
	double R = c.getR();															// DC resistance [Ohm]
	double An = c.getAnodeSurface();												// active anode surface area excluding cracks [m2]
	c.getStates(si, &Ii);															// State-object with the cell's state
	si.getStates(ns, states);														// array with the battery states
	double cap;																		// the remaining capacity [Ah]
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_batteryStates is getting the capacity"<<endl;
		if(checkCap)
			cap = getCapacity(blockDegradation);
		else
			cap = 0; // capacity is not measured
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkUp_batteryStates when measuring the cell capacity: "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// Write the results to the csv file
	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::checkUp_batteryStates is writing the capacity to the csv file"<<endl;
	string fullName = ".\\"+ID+"\\"+nameState;										// we want to write the file in a subfolder, so append the name of the subfolder before the name of the csv file
	ofstream output(fullName,std::ios_base::app);									// open the csv file and append the new data after the already existing data
	if (output.is_open()){
		output<<cumCycle<<","<<cumTime<<","<<cumAh<<","<<cumWh; 					// write the data points for where the cell is in it's life
		output<<","<<cap;															// write the cell capacity
		for(int s=0;s<ns;s++)
			output<<","<<states[s];													// write the state variables of the cell
		output<<","<<R;																// write the total cell resistance
		output<<","<<An;															// write the active anode surface area an*thickn*elec_surf
		output<<"\n";																// write an end-line (we have written everything we want)
		output.close();																// close the file
	}
	else{																			// we couldn't open the file
		if(verbose >= printCrit)
			cerr<<"ERROR in Cycler::checkUp_batteryStates. File "<<fullName<<" could not be opened. Throwing an error"<<endl;
		throw 1001;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_batteryStates terminating"<<endl;

	return cap;
}

void Cycler::checkUp_OCVcurves(bool blockDegradation, double ocvpini, double ocvnini){
	/*
	 * Function to measure and write the half-cell OCV curves as part of a check-up.
	 * This function will add 4 rows to the csv file with the results (DegradationData_OCV.csv in the subfolder of this Cycler)
	 * 		the first row contains the 'x-axis' of the half cell curves with the discharged charge [Ah]. It is 0 when the OCV of the cell is the maximum voltage, negative for points where the cell OCV is larger and positive for points where the cell OCV is smaller
	 *		the second row contains the cathode potential at the corresponding points on the x-axis
	 *		the third row contains the anode potential at the corresponding point on the x-axis
	 *		the fourth row is an empty row to indicate the previous 3 rows are one 'data set'
	 * At the end of the third row, also the 'operating point' of the cell is written, i.e. the potentials of the electrodes when the check-up procedure was started
	 *		this is especially interesting when simulating calendar ageing because it allows to check the anode potential at which the cell is actually resting
	 *
	 * IN
	 * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
	 * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
	 * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
	 * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
	 * ocvpini				cathode potential at the operating point when the check-up is called [V]
	 * ocvnini				anode potential at the operating point when the check-up is called [V]
	 *
	 * THROWS
	 * 1001 				the file in which to write the results couldn't be opened
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_OCVcurves is starting"<<endl;

	// variables
	string nameOCV = "DegradationData_OCV.csv";										// name of the csv file in which the OCV curves of the check-up will be written
	int Nocv;																		// number of actual data points of the OCV curves
	int Ninocv = 25.0*3600.0/2.0/750.0*5.0;											// length of the arrays with the OCV curves (0.04C -> 25 hours, 2s time steps, store one in every 750 points, *5 to ensure they are long enough)
	double ocvAh[Ninocv], ocvp[Ninocv], ocvn[Ninocv]; 								// arrays to contain the OCV curves

	// record the OCV cycles
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_OCVcurves is getting the OCV curves"<<endl;
		getOCV(Ninocv, ocvAh, ocvp, ocvn, &Nocv);									// measure the OCV curves
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkUp_OCVcurves when measuring the OCV curves: "<<e<<". Trying to recover by making longer arrays."<<endl;

		int Ninocv2 = 10.0*Nocv;													// make 10 times longer arrays
		double ocvAh2[Ninocv2], ocvp2[Ninocv2], ocvn2[Ninocv2];
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp_OCVcurves is getting the OCV curves with larger arrays"<<endl;
			getOCV(Ninocv2, ocvAh2, ocvp2, ocvn2, &Nocv);							// try again with these larger arrays
			if(verbose >= printCrit)
				cout<<"we have recovered from the error in Cycler::checkUp_OCVcurves when measuring the OCV curves by making longer arrays."<<endl;
		}
		catch(int e2){
			if(verbose >= printCrit)
				cout<<"Error in Cycler::checkUp_OCVcurves when measuring the OCV curves with longer arrays: "<<e<<". Throwing the error on."<<endl;
			throw e2;
		}
	}

	// Write the curves to the csv file
	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::checkUp_OCVcurves is writing the OCV curves to a csv file"<<endl;
	string fullName = ".\\"+ID+"\\"+nameOCV;										// we want to write the file in a subfolder, so append the name of the subfolder before the name of the csv file
	ofstream output(fullName,std::ios_base::app);									// open the csv file and append the new data after the already existing data
	if (output.is_open()){

		// write the x-axis, discharged Ah
		for (int i=0;i<Nocv;i++)
			output<<ocvAh[i]<<",";
		output<<"\n";

		// write the cathode OCV
		for (int i=0;i<Nocv;i++)
			output<<ocvp[i]<<",";
		output<<"\n";

		// writhe the anode OCV
		for (int i=0;i<Nocv;i++)
			output<<ocvn[i]<<",";

		// write the voltages of the operating point of the cell next to the anode OCV
		output<<","<<","<<ocvpini<<","<<ocvnini<<",";
		output<<"\n";

		// leave one row blank to indicate these three lines were one measurement
		output<<"\n";
		output.close();																	// close the file
	}
	else{
		if(verbose >= printCrit)
			cerr<<"ERROR on Cycler::checkUp_OCVcurves. File "<<fullName<<" could not be opened. Throwing an error"<<endl;
		throw 1001;
	}

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_OCVcurves terminating"<<endl;
}

void Cycler::checkUp_CCCV(bool blockDegradation, int nCycles, double Crates[], double Ccut_cha, double Ccut_dis, bool includeCycleData){
	/*
	 * Function to simulate and write some CCCV cycles as part of a check-up
	 * The voltage and temperature while cycling are written in a csv file in the subfolder of this Cycler.
	 * The data from each check-up is written in a new csv file, with as name DegradationData_CheckupCycle_x.csv
	 * 		where 'x' is 1 (for the first check-up), 2 (for the second check-up), etc.
	 * In the file, a data entry is recorded for every 2 seconds. Each data entry is written on one row of the csv file.
	 * Each row has 15 columns (with the values for that entry):
	 * 	 	total_time 		the total (cumulative) time since the start of this function, i.e. it is 0 for the first row [s]
	 * 	 	Ah_throughput 	the total (cumulative) charge throughput since the start of this function, i.e. it is 0 for the first row [Ah]
	 * 	 	Wh_throughput 	the total (cumulative) energy throughput since the start of this function, i.e. it is 0 for the first row [Wh]
	 * 	 	I 				the battery current at this point in time [A], > 0 for discharge, < 0 for charge
	 * 	 	V  				the battery voltage at this point in time [V]
	 * 	 	OCVp 			the cathode potential at this point in time [V]
	 * 	 	OCVn  			the anode potential at this point in time [V]
	 * 	 	Temperature 	the battery temperature at this point in time [K]
	 * 	 	charge_time  	the cumulative time spent on charging since the start of this function, i.e. it is 0 for the first row [s]
	 * 	 	charge_Ah  		the cumulative charged charge since the start of this function, i.e. it is 0 for the first row [Ah]
	 * 	 	charge_Wh  		the cumulative charged energy since the start of this function, i.e. it is 0 for the first row [Wh]
	 * 	 	discharge_time 	the cumulative time spent on discharging since the start of this function, i.e. it is 0 for the first row [s]
	 * 	 	discharge_Ah 	the cumulative discharged charge since the start of this function, i.e. it is 0 for the first row [Ah]
	 * 	 	discharge_Wh 	the cumulative discharged energy since the start of this function, i.e. it is 0 for the first row [Wh]
	 * 	 	rest_time 		the cumulative time spent on resting since the start of this function, i.e. it is 0 for the first row [s]
	 *
	 * IN
	 * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
	 * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
	 * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
	 * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
	 * nCycles				number of different cycles to be simulated (i.e. the length of the array crates)
	 * Crates				array with the Crates of the CC phase of the different cycles to be checked, must be all positive
	 * Ccut_cha				C rate of the cutoff current for the CV phase in the charges, must be positive
	 * Ccut_dis				C rate of the cutoff current for the CV phase in the discharges, must be positive
	 * includeCycleData	 	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
	 * 							if false, the cycling data from the CCCV curves is only stored in a separate data file (and not part of the 'regular cycling data files': 'cyclingData_x.csv')
	 * 							if true, the cycling data from the CCCV curves is stored both in a separate document and in the regular cycling data files ('cyclingData_x.csv')
	 *
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_CCCV starting"<<endl;

	// variables
	string nameCCCV = "DegradationData_CheckupCycle_" + to_string(indexdegr) + ".csv";	// name of the csv file in which the cycling data will be written
	double ahi, whi, timei;																// unneeded feedback variables (charge, energy and time spent in underlying functions)
	double dt = 2;																		// use time steps of 2 seconds
	int feedb_old = CyclingDataTimeInterval;											// the original data collection time interval
	int dataTimeInterval = 2;															// time resolution at which we want to store the cycling data from the CCCV curves
	CyclingDataTimeInterval = dataTimeInterval;											// update the time resolution at which cycling data is stored for the CCCV cycles

	// fully charge the cell
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_CCCV is charging the cell before doing the CCCV cycles"<<endl;
		CC_V_CV_I(1.0, c.getVmax(), Ccut_cha,dt, blockDegradation, &ahi, &whi, &timei); // 1C charge to the specified cutoff Crate
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkup_CCCV when initially charging the cell: "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// flush the cycling data of the cell which was still stored in the arrays so we can start with an empty data buffer to store the results of the check-up
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_CCCV is flushing the previously stored cycling data"<<endl;
		writeCyclingData();
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkup_CCCV when writing the original cycling data: "<<e<<endl<<flush;
		throw e;
	}

	// do the cycles
	for(int i=0;i<nCycles; i++){
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp_CCCV is simulating a CCCV cycle at Crate = "<<Crates[i]<<endl;
			CC_V_CV_I(Crates[i], c.getVmin(), Ccut_dis, dt, blockDegradation, &ahi, &whi, &timei);	// discharge to the minimum cell voltage (CC at the given C rate and CV to the given cutoff C rate)
			CC_V_CV_I(Crates[i], c.getVmax(), Ccut_cha, dt, blockDegradation, &ahi, &whi, &timei);	// charge to the maximum cell voltage (CC at the given C rate and CV to the given cutoff C rate)
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in Cycler::checkup_CCCV when following CCCV cycle at Crate: "<<Crates[i]<<". error "<<e<<". Throwing it on"<<endl<<flush;
			throw e;
		}
	}

	// write the cycling data from the CCCV cycles
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_CCCV is writing the cycling data from the CCCV cycles"<<endl;
		writeCyclingData(nameCCCV, !includeCycleData);												// if we don't want to include the cycling data, clear the buffer
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkup_CCCV when writing the cycling data from the CCCV cycles: "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// restore the old data-collection setting
	CyclingDataTimeInterval = feedb_old;

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_CCCV terminating"<<endl;
}

void Cycler::checkUp_pulse(bool blockDegradation, string profileName, int profileLength, bool includeCycleData){
	/*
	 * Function to simulate and write a pulse test as part of a check-up.
	 * The voltage and temperature while cycling are written in a csv file in the subfolder of this Cycler.
	 * The data from each check-up is written in a new csv file, with as name DegradationData_CheckupPulse_x.csv
	 * 		where 'x' is 1 (for the first check-up), 2 (for the second check-up), etc.
	 * In the file, a data entry is recorded for every 2 seconds. Each data entry is written on one row of the csv file.
	 * Each row has 15 columns (with the values for that entry):
	 * 	 	total_time 		the total (cumulative) time since the start of this function, i.e. it is 0 for the first row [s]
	 * 	 	Ah_throughput 	the total (cumulative) charge throughput since the start of this function, i.e. it is 0 for the first row [Ah]
	 * 	 	Wh_throughput 	the total (cumulative) energy throughput since the start of this function, i.e. it is 0 for the first row [Wh]
	 * 	 	I 				the battery current at this point in time [A], > 0 for discharge, < 0 for charge
	 * 	 	V  				the battery voltage at this point in time [V]
	 * 	 	OCVp 			the cathode potential at this point in time [V]
	 * 	 	OCVn  			the anode potential at this point in time [V]
	 * 	 	Temperature 	the battery temperature at this point in time [K]
	 * 	 	charge_time  	the cumulative time spent on charging since the start of this function, i.e. it is 0 for the first row [s]
	 * 	 	charge_Ah  		the cumulative charged charge since the start of this function, i.e. it is 0 for the first row [Ah]
	 * 	 	charge_Wh  		the cumulative charged energy since the start of this function, i.e. it is 0 for the first row [Wh]
	 * 	 	discharge_time 	the cumulative time spent on discharging since the start of this function, i.e. it is 0 for the first row [s]
	 * 	 	discharge_Ah 	the cumulative discharged charge since the start of this function, i.e. it is 0 for the first row [Ah]
	 * 	 	discharge_Wh 	the cumulative discharged energy since the start of this function, i.e. it is 0 for the first row [Wh]
	 * 	 	rest_time 		the cumulative time spent on resting since the start of this function, i.e. it is 0 for the first row [s]
	 *
	 * IN
	 * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
	 * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
	 * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
	 * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
	 * profileName 			name of the csv file which contains the current profile for the pulse test
	 * 							the first column contains the current in [A] (positive for discharge, negative for charge)
	 * 							the second column contains the time in [sec] the current should be maintained
	 * 							the profile must be a net discharge, i.e. sum (I*dt) > 0
	 * profileLength		length of the current profiles (number of rows in the csv file)
	 * includeCycleData	 	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
	 * 							if false, the cycling data from the pulse discharge is only stored in a separate data file (and not part of the 'regular cycling data files': 'cyclingData_x.csv')
	 * 							if true, the cycling data from the pulse discharge is stored both in a separate document and in the regular cycling data files ('cyclingData_x.csv')
	 *
	 * THROWS
	 * 1013 				the pulse profile is not a net discharge
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_pulse starting"<<endl;

	// variables
	string namePulse = "DegradationData_CheckupPulse_" + to_string(indexdegr) + ".csv";	// name of the csv file in which the cycling data will be written
	double dt = 2;																		// take time steps of 2 seconds
	int limit = 1;																		// if a voltage limit is reached during the pulse discharge profile, do a CV at that voltage for the rest of the time in that step
	int feedb_old = CyclingDataTimeInterval;											// the old data collection time interval
	int dataTimeInterval = 2;															// time resolution at which we want to store the cycling data from the pulse discharge
	CyclingDataTimeInterval = dataTimeInterval;											// update the time resolution of the cycling data collection to the new value for the pulse discharge
	double ahi, whi, timei;																// unneeded feedback variables
	double Ccut = 0.05;																	// C rate of the cutoff current for the CV phase

	// *********************************************************** 1 get the pulse profile ***********************************************************************
	// read the pulse profile
	double I[profileLength], Tprof[profileLength];										// arrays for the profile
	int T[profileLength];																// array for the duration of each step as integers (when reading the file, they become doubles)
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_pulse is reading the current profile"<<endl;
		loadCSV_2col(profileName, profileLength, I, Tprof);								// read the file
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkUp_pulse when reading the pulse profile: "<<e<<". Throwing it on"<<endl;
		throw e;
	}

	// convert the duration of each step to integers
	for(int i=0;i<profileLength;i++)
		T[i] = Tprof[i];

	// check that the pulse profile is a discharge
	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::checkUp_pulse is checking that the current profile is a discharge"<<endl;
	double aht = 0;																		// total charge throughput of the profile [As]
	for(int i=0;i<profileLength;i++)
		aht += I[i] * T[i];
	if (aht <= 0){																		// if the total charge throughput is negative, the profile is a net charge
		if(verbose >= printCrit)
			cerr<<"ERROR in Cycler::checkUp_pulse: the pulse profile is not a net discharge. "
				"The total charge throughput was "<<aht/3600.0<<"Ah, where a negative number indicates a net charge."<<endl<<flush;
		// The total charge throughput MUST BE positive, i.e. there must be more 'discharge' than 'charge' in the profile.
		// Else, the user has to change the function 'checkUp_pulse' to deal with a net charge
		throw 1013;
	}


	// *********************************************************** 2 charge the cell ***********************************************************************
	// fully charge the cell
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_pulse is charging the cell"<<endl;
		CC_V_CV_I(1.0, c.getVmax(), Ccut,dt, blockDegradation, &ahi, &whi, &timei); 	// fully charge the cell to the maximum voltage (1C CC, CV until 0.05C cut off current)
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkup_pulse when initially charging the cell: "<<e<<endl<<flush;
	}

	// flush the cycling data of the cell which was still stored in the arrays so we can start with an empty data buffer to store the results of the check-up
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_pulse is flushing out the previously stored cycling data"<<endl;
		writeCyclingData();
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkup_pulse when writing the original cycling data: "<<e<<endl<<flush;
		throw e;
	}


	// *********************************************************** 3 do the pulse discharge until the lower voltage is reached ***********************************************************************
	// The current pulses must be in a file with the name given by 'CheckPulseProfile'
	// 	and length given by 'length_checkPulseProfile'
	// 	the first column gives the current [A]
	// 	the second column gives the time each current should be maintained [s]
	// the cell starts the pulse tests at the maximum cell voltage (SoC = 100%). Keep doing the pulses until the lower voltage is reached
	bool Vmin = false;																	// has the lower voltage limit been reached?
	double Vupp = c.getVmax();															// the upper voltage limit is the maximum voltage of the cell
	double Vlow = c.getVmin();															// the lower voltage limit is the minimum voltage of the cell
	int Vlim;																			// return-integer of followI indicating which voltage limit was hit
	try{
		// a loop to repeatedly apply the pulse profile until the lower voltage is reached
		while(!Vmin){
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp_pulse is applying the pulse profile and has so far discharged "<<ahi<<"Ah"<<endl;
			Vlim = followI(profileLength, profileName, blockDegradation, limit, Vupp, Vlow, &ahi, &whi, &timei); // follow the current pulses, the return-integer indicates which voltage limit was hit while following the profile
			Vmin = (Vlim == -1) || (Vlim==10);											// a return-integer of -1 or 10 means the lower voltage was hit, stop when it is
		}
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkup_pulse when following the pulses: "<<e<<endl<<flush;
	}

	// write the cycling data from the pulse cycles
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp_pulse is writing the cycling data from the pulse profile"<<endl;
		writeCyclingData(namePulse, !includeCycleData);									// if we don't want to include the cycling data, clear the buffer
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkup_pulse when writing the cycling data from the pulse discharge: "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// restore the old data-collection setting
	CyclingDataTimeInterval = feedb_old;

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp_pulse terminating"<<endl;
}

double Cycler::checkUp(struct checkUpProcedure proc, int cumCycle, double cumTime, double cumAh, double cumWh){
	/*
	 * Function to do a check-up during a degradation experiment
	 * a check-up can consists of:
	 * 		a capacity measurement
	 * 		a half-cell OCV measurement
	 * 		some CC CV cycles (CC discharge, CV discharge, CC charge, CV charge) at 0.5, 1 and 2 C
	 * 		a pulse discharge test
	 *
	 * The results are written in csv files, in the subfolder of this Cycler.
	 * 		DegradationData_batteryState.csv				one new row of data appended at the end of the existing csv file
	 * 		DegradationData_OCV.csv							4 new rows of data appended at the end of the existing csv file
	 * 		DegradationData_CheckupCycle_x.csv 				a new file for the data of this check-up, x is the index of the check-up (1 for the first check-up, 2 for the second, etc.)
	 * 		DegradationData_CheckupPulse_x.csv				a new file for the data of this check-up, x is the index of the check-up (1 for the first check-up, 2 for the second, etc.)
	 * See the individual functions (checkUp_yyy) for an exact description of what is in each file.
	 *
	 * IN
	 * proc 				structure with the parameters of the check-up procedure with the following fields:
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
	 * 		Ccut_cha			C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
	 * 		Ccut_dis			C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
	 * 		profileName 		name of the csv file which contains the current profile for the pulse test
	 * 								the first column contains the current in [A] (positive for discharge, negative for charge)
	 * 								the second column contains the time in [sec] the current should be maintained
	 * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
	 * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
	 * cumCycle				number of cycles up to now [-]
	 * cumTime				time this cell has been cycled up to now [hour]
	 * cumAh				cumulative Ah throughput up to now [Ah]
	 * cumWh				cumulative Wh throughput up to now [Wh]
	 *
	 * OUT
	 * cap 					the remaining capacity of the cell [Ah]
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp starting"<<endl;

	// variables
	State sini;																// initial cell state
	double Iini;															// initial cell current
	double v;																// initial cell voltage
	double tcellini, Tenvini, Trefi;										// initial cell, environmental and reference temperature
	double dt = 2;															// time step when cycling and resting [s]
	double cap = 0;															// capacity of the cell [Ah]
	double ocvpini, ocvnini, etap, etan, rdrop;								// unneeded feedback variables
	double ahi, whi, timei;													// unneeded feedback variables

	// Store the initial battery state, OCV and temperature from before the check-up
	c.getStates(sini, &Iini);												// initial cell state and current
	c.getVoltage(verbose >= printCrit, &v, &ocvpini, &ocvnini, &etap, &etan, &rdrop, &tcellini); // initial cell voltage and temperature
	c.getTemperatures(&Tenvini, &Trefi);									// initial environmental and reference temperature
	int feedbackini = CyclingDataTimeInterval;								// initial data collection time interval

	// if we don't want to include the cycling data from the check-up in the cycling data from the cell
	// push the cycling data and set the data collection to 0 to stop recording
	if (!proc.includeCycleData){
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp is flushing the previously stored cycling data"<<endl;
			writeCyclingData();												// write the cycling data which was still stored
			CyclingDataTimeInterval = 0;									// don't record cycling data from the check-up
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in Cycler::checkUp when writing the previously stored cycling data: "<<e<<endl<<flush;
			throw e;
		}
	}


	// *********************************************************** 1 precharge the cell & set the temperature ***********************************************************************

	// Get the cell to 0.2V below the maximum voltage.
	// This is needed because later we are changing the temperature, which might affect the voltage a little bit.
	// If the cell is fully charged and the voltage changes a little bit, you can get an illegal voltage (e.g. Vmax + 0.001V), which results in an error
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::checkUp is bringing the cell to just below its maximum voltage"<<endl;
		CC_V_CV_I(1.0,  c.getVmax()-0.2, 0.05, dt, proc.blockDegradation, &ahi, &whi, &timei); // 1C current, 0.05C cutoff current
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::checkUp when setting the cell voltage to Vmax-0.2V: "<<e<<". Throwing it on"<<endl;
		throw e;
	}

	// do the check-up at reference temperature
	c.setTenv(Trefi);														// set the environmental temperature to the reference temperature
	c.setT(Trefi);															// set the cell temperature to the reference temperature

	// *********************************************************** 2 do the different check-up tests ***********************************************************************
	// Measure the capacity
	if(proc.capCheck){
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp is starting a capacity check"<<endl;
			cap = checkUp_batteryStates(proc.blockDegradation, true, cumCycle, cumTime, cumAh, cumWh);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in Cycler::checkUp when measuring the cell capacity: "<<e<<". Skip the capacity measurement"<<endl;
		}
	}

	// measure the OCV curves
	if (proc.OCVCheck){
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp is starting an OCV check"<<endl;
			checkUp_OCVcurves(proc.blockDegradation, ocvpini, ocvnini);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in Cycler::checkUp when measuring the OCV curves: "<<e<<". Skip the OCV measurements"<<endl;
			throw e;
		}
	}

	// measure the CCCV cycles
	if(proc.CCCVCheck){
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp is starting a CCCV cycle check"<<endl;
			checkUp_CCCV(proc.blockDegradation, proc.nCycles, proc.Crates, proc.Ccut_cha, proc.Ccut_dis, proc.includeCycleData);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in Cycler::checkUp when measuring the CCCV cycles: "<<e<<". Skip the CCCV measurement"<<endl;
		}
	}

	// measure the discharge pulses
	if(proc.pulseCheck){
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::checkUp is starting a pulse profile check"<<endl;
			checkUp_pulse(proc.blockDegradation, proc.profileName, proc.profileLength, proc.includeCycleData);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in Cycler::checkUp when measuring the pulse test: "<<e<<". Skip the CCCV measurement"<<endl;
		}
	}

	// increase the counter of the number of check-ups we have done
	indexdegr++;

	// *********************************************************** 3 reset the original battery state ***********************************************************************

	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::checkUp is restoring the initial battery state"<<endl;

	c.setTenv(Tenvini);														// restore the environmental temperature
	CyclingDataTimeInterval = feedbackini;									// restore the initial data-collection setting
	if(proc.blockDegradation)
		c.setStates(sini, Iini);											// restore the exact initial battery state
	else{
		double Trest = 3600;												// rest time is 1 hour
		double Ccut = 0.05;													// C rate of cutoff current from the CV phase
		CC_V_CV_I(1, ocvpini-ocvnini, Ccut, dt, proc.blockDegradation, &ahi, &whi, &timei); // bring the cell back to the original OCV (OCV_cell = OCV_p - OCV_n)
		CC_t(0.0, dt, proc.blockDegradation, Trest, &ahi, &whi, &timei);					// rest such that the cell temperatures goes back to the environmental temperature
	}

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::checkUp terminating with capacity "<<cap<<"Ah"<<endl;

	// return the remaining capacity
	return cap;
}

void Cycler::cycleAgeing(double dt, double Vma, double Vmi, double Ccha, bool CVcha, double Ccutcha,
		double Cdis, bool CVdis, double Ccutdis, double Ti, int nrCycles, int nrCap, struct checkUpProcedure proc){
	/*
	 * Function to simulate a cycle degradation experiment where a cell is continuously cycled at constant current and/or constant voltage.
	 * The parameters of the cycling regime are set by the inputs.
	 * After a set number of cycles, a check-up is done where the cell capacity, OCV curves, etc. are measured.
	 * These results are written to csv files.
	 *
	 * IN
	 * dt 		the time step to be used in the cycling, small enough to ensure stability and accuracy (1-5 seconds) [s]
	 * Vma 		maximum voltage of this cycle, must be below the maximum voltage of the cell [V]
	 * Vmi 		minimum voltage of this cycle, must be above the minimum voltage of the cell [V]
	 * Ccha 	C-rate at which the battery should be charged during the CC phase[-]
	 * CVcha 	boolean indicating if a CV charge should be done after the CC charge (true) or not (false)
	 * 			if true, charging is CC CV
	 * 			if false, charging is CC only
	 * Ccutcha 	C rate of the cutoff current for the CV charge [-], > 0
	 * Cdis 	C rate at which the battery should be discharged [-]
	 * CVdis 	boolean indicating if a CV discharge should be done after the CC discharge (true) or not (false)
	 * 			if true, discharging is CC CV
	 * 			if false, discharging is CC only
	 * Ccutdis 	C rate of the cutoff current for the CV discharge [-], > 0
	 * Ti 		environmental temperature at which the test should be performed [K]
	 * nrCycles	the number of cycles that has to be simulated [-]
	 * nrCap 	the number of cycles between consecutive check-ups [-]
	 * proc 	structure with the parameters of the check-up procedure.
	 * 			The struct is defined in Cycler.hpp and has the following fields:
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
	 * 								the second column contains the time in [sec] the current should be maintained for
	 * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
	 * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
	 *
	 *
	 * throws
	 * 1014		the input parameters describing the cycling regime are invalid
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::cycleAgeing starting"<<endl;

	// Check the input parameters
	bool vmax = Vma > c.getVmax(); 													// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in Cycler::cycleAgeing. The maximum voltage "<<Vma<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vmi < c.getVmin();													// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in Cycler::cycleAgeing. The minimum voltage "<<Vmi<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	bool Temin = Ti < TMIN;															// check the temperature is above 0 degrees, TMIN is defined in State.hpp
	if (Temin)
		cerr<<"Error in Cycler::cycleAgeing. The temperature "<<Ti<<"K is too low. The minimum value is 273"<<endl<<flush;
	bool Temax = Ti > TMAX;															// check the temperature is below 60 degrees, TMAX is defined in State.hpp
	if (Temax)
		cerr<<"Error in Cycler::cycleAgeing. The temperature "<<Ti<<" is too high. The maximum value is (273+60)"<<endl<<flush;
	bool cchar = Ccha <= 0;															// check the charging Crate is positive
	if (cchar)
		cerr<<"Error in Cycler::cycleAgeing. The charging Crate "<<Ccha<<" is negative. It must be positive"<<endl<<flush;
	bool cdis = Cdis <= 0;															// check the discharging Crate is positive
	if (cdis)
		cerr<<"Error in Cycler::cycleAgeing. The discharging Crate "<<Cdis<<" is negative. It must be positive"<<endl<<flush;
	bool ccvchar = Ccutcha <= 0;													// check the charging cutoff current is positive
	if (ccvchar)
		cerr<<"Error in Cycler::cycleAgeing. The Crate of the cutoff current for the CV charge "<<Ccutcha<<" is negative. It must be positive"<<endl<<flush;
	bool ccvdis = Ccutdis <= 0;														// check the discharging cutoff current is positive
	if (ccvdis)
		cerr<<"Error in Cycler::cycleAgeing. The Crate of the cutoff current for the CV discharge "<<Ccutdis<<" is negative. It must be positive"<<endl<<flush;
	bool cycles = nrCycles <= nrCap;												// check the number of cycles between consecutive check-ups is lower than the total number of cycles
	if (cycles)
		cerr<<"Error in Cycler::cycleAgeing. The number of cycles between two check ups "<<nrCap<<" is higher than the total number of cycles "<<nrCycles<<endl<<flush;
	if(vmax || vmin || Temin || Temax || cchar || cdis || ccvchar || ccvdis || cycles)
		throw 1014;


	// *********************************************************** 1 variables & settings ***********************************************************************

	// variables
	double ahi;																		// discharged charge in this charge or discharge [Ah]
	double whi;																		// discharged energy in this charge or discharge [Wh]
	double ti;																		// time spent in this charge or discharge [sec]
	double Ahtot = 0;																// charge throughput from the cycling regime until now [Ah]
	double Whtot = 0;																// energy throughput from the cycling regime until now [Wh]
	double timetot = 0;																// time spent while following the cycling regime until now [hour]
	double cap;																		// capacity of the cell at this point in time [Ah]
	bool blockDegradation = false; 													// account for degradation while we cycle
	bool final = true;																// boolean to indicate if a check-up at the end of the cycling regime is needed
	double capnom = c.getNominalCap();												// nominal cell capacity [Ah] to convert Crate to Amperes


	// *********************************************************** 2 cell initialisation ***********************************************************************

	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::cycleAgeing is initialising the cell"<<endl;

	// set the temperature
	c.setT(Ti);																		// set the cell temperature
	c.setTenv(Ti);																	// set the environmental temperature

	// Get the battery to Vma so the cycling can start with a discharge
	try{
		if (CVcha)
			CC_V_CV_I(Ccha, Vma, Ccutcha, dt, blockDegradation, &ahi, &whi, &ti);	// CC and CV charge
		else
			CC_V(-Ccha*capnom, dt, blockDegradation, Vma, &ahi, &whi, &ti);			// CC charge (must have current as input, not C rate)
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in a subfunction of Cycler::cycleAgeing when getting the cell to Vma, error "<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// do an initial check up
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::cycleAgeing is doing an initial check-up"<<endl;
		cap = checkUp(proc, 0, timetot, Ahtot, Whtot);
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in the initial check-up Cycler::cycleAgeing "<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// *********************************************************** 3 cycle age the cell ***********************************************************************

	for (int i=0;i<nrCycles;i++){													// loop for the specified number of cycles
		try{

			// discharge
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::cycleAgeing is discharging the cell in cycle number "<<i<<endl;
			if (CVdis)
				CC_V_CV_I(Cdis, Vmi, Ccutdis, dt, blockDegradation, &ahi, &whi, &ti);// CC and CV discharge
			else
				CC_V(Cdis*capnom, dt, blockDegradation, Vmi, &ahi, &whi, &ti);		// CC discharge (must have current as input, not C rate)
			Ahtot += abs(ahi);														// increase the charge throughput with the throughput of this discharge
			Whtot += abs(whi);														// increase the energy throughput with the throughput of this discharge
			timetot += (ti/3600);													// increase the total time the time of this discharge

			// charge
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::cycleAgeing is charging the cell in cycle number "<<i<<endl;
			if (CVcha)
				CC_V_CV_I(Ccha, Vma, Ccutcha, dt, blockDegradation, &ahi, &whi, &ti);// CC and CV charge
			else
				CC_V(-Ccha*capnom, dt, blockDegradation, Vma, &ahi, &whi, &ti);		// CC charge (must have current as input, not C rate)
			Ahtot += abs(ahi);														// increase the charge throughput with the throughput of this charge
			Whtot += abs(whi);														// increase the energy throughput with the throughput of this charge
			timetot += (ti/3600);													// increase the total time the time of this charge

			// do a check-up every nrCap cycles
			// 	i is the cycle number, so when it is a multiple of nrCap we need to do a check-up
			//  do i+1 to avoid doing a check-up in the first cycle
			if (fmod(i+1,nrCap) == 0){												// remainder is 0 -> it is a multiple of nrCap
				if(verbose >= printCyclerHighLevel)
					cout<<"Cycler::cycleAgeing is doing a check-up in cycle number "<<i<<endl;
				cap = checkUp(proc, i+1, timetot, Ahtot, Whtot);					// do the check-up procedure

				// End the experiment if the cell capacity has decreased too much
				if(cap<capnom/2.0){													// end simulation if the cell has only 50% of capacity left
					cout<<"Cycler::cycleAgeing has finished cycling regime "<<ID<<" early because the cell has already lost 50% of its capacity.";
					cout<<" We have done "<<i<<" cycles instead of "<<nrCycles<<" and the remaining capacity now is "<<cap<<" [Ah]"<<endl<<flush;
					final = false;													// skip the final check-up because we just did one
					break;															// stop cycling
				}
			}
		} // end try block

		// Catch an error which occurred while cycling the cell (or during the check-up procedure)
		catch(int e){

			// print an error message (this is a critical error)
			if(verbose >= printCrit){
				cout<<"Error in Cycler::cycleAgeing while cycling the cell according to  cycling regime "<<ID<<". Error encountered is "<<e<<". Stop cycling now.";
				cout<<" We have done "<<i<<" cycles instead of "<<nrCycles<<" and the capacity last measured is "<<cap<<" [Ah]"<<endl<<flush;
			}

			// we probably cannot do a full check-up procedure because the cell is in an illegal state.
			// Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
			checkUp_batteryStates(proc.blockDegradation, false, i+1, timetot, Ahtot, Whtot);

			final = false;															// skip the final check-up
			break;																	// stop cycling
		} // end try-catch block
	} // end loop to cycle the cell

	// *********************************************************** 4 final check-up ***********************************************************************


	// Don't do the check-up if the simulation ended due to an error (or shortage in capacity)
	if (final){
		try{
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::cycleAgeing is doing a final check-up"<<endl;
			checkUp(proc, nrCycles, timetot, Ahtot, Whtot);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in a in the final check-up of Cycler::cycleAgeing "<< e<<". Throwing it on "<<endl;
			throw e;
		}
	}

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::cycleAgeing terminating"<<endl;
}

void Cycler::calendarAgeing(double dt, double V, double Ti, int Time, int timeCheck, int mode, struct checkUpProcedure proc){
	/*
	 * Function to simulate a calendar degradation experiment where a cell is rested at a given voltage.
	 * The parameters of the calendar regime are set by the inputs.
	 * After a given amount of time, a check-up is done where the cell capacity, OCV curves, etc. are measured.
	 * These results are written to csv files.
	 *
	 * IN
	 * dt 		the time step to be used [sec]
	 * V 		the voltage at which the battery has to rest [V]
	 * Ti 		the temperature at which the battery has to rest [K]
	 * Time 	the time for which the battery has to rest [days]
	 * 			must be a multiple of timeCheck
	 * timeCheck the time after which a check-up has to be done [days]
	 * mode 	integer deciding how often to recharge to the specified voltage.
	 * 			When a cell rests, the voltage will slip a bit due to degradation, so the question is how often to recharge.
	 * 			0 	recharge only after a check-up
	 * 			1 	recharge every day
	 * 			2 	float the cell at the specified voltage, instead of letting it rest
	 * 				this option is not recommended, it takes ages to calculate
	 * 			Any other value will produce an error
	 * proc 	structure with the parameters of the check-up procedure with the following fields:
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
	 *
	 * throws
	 * 1014		the input parameters describing the calendar regime are invalid
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::calendarAgeing starting"<<endl;

	// Check the input parameters
	bool vmax = V > c.getVmax(); 				// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in Cycler::CalendarAgeing. The voltage "<<V<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = V < c.getVmin();				// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in Cycler::CalendarAgeing. The voltage "<<V<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	bool Temin = Ti < TMIN;						// check the temperature is above 0 degrees, TMIN is defined in State.hpp
	if (Temin)
		cerr<<"Error in Cycler::CalendarAgeing. The temperature "<<Ti<<"K is too low. The minimum value is 273"<<endl<<flush;
	bool Temax = Ti > TMAX;						// check the temperature is below 60 degrees, TMAX is defined in State.hpp
	if (Temax)
		cerr<<"Error in Cycler::CalendarAgeing. The temperature "<<Ti<<" is too high. The maximum value is (273+60)"<<endl<<flush;
	bool mod = (mode < 0 || mode > 2);			// check the setting in 'mode' is allowed
	if (mod)
		cerr<<"Error in Cycler::CalendarAgeing. The resting mode "<<mode<<" is illegal. Only values of 0, 1 or 2 are allowed"<<endl<<flush;
	bool cycles = (fmod(Time,timeCheck) > 1);	// check the total resting time is a multiple of the the time between two check-ups
	if (cycles)									// allow a margin of 1 day (for numerical errors)
		cerr<<"Error in Cycler::CalendarAgeing. The total resting time "<<Time<<" is not a multiple of the time between two check ups "<<timeCheck<<endl<<flush;
	if(vmax || vmin || Temin || Temax || mod || cycles)
		throw 1014;

	// print a warning if mode is 2 because it will take very long to simulate
	if(mode == 2)
		cout<<"Warning in Cycler::CalendarAgeing. You have chosen to float the cell at a constant voltage for a long period. "
				"It will take long to simulate this."<<endl<<flush;


	// *********************************************************** 1 variables & settings ***********************************************************************

	// Variables
	double ahi, ahi2;																// discharged charge in this charge or discharge [Ah]
	double whi, whi2;																// discharged energy in this charge or discharge [Wh]
	double ti, ti2;																	// time spent in this charge or discharge [sec]
	double Ahtot = 0;																// cumulative charge throughput until now [Ah]
	double Whtot = 0;																// cumulative energy throughput until now [Wh]
	double timetot = 0;																// cumulative time until now [hour]
	double cap;																		// cell capacity [Ah]
	double capnom = c.getNominalCap();												// nominal cell capacity
	double Ccut = 0.005;															// Crate of the cutoff current for CV phases [-]
	bool blockDegradation = false; 													// do account for degradation during calendar

	// time at which a check-up has to be done
	int nrdt = Time/timeCheck; 														// number of check-ups to be done
	double trest = timeCheck * (24.0*3600.0);										// resting time between consecutive check-ups in seconds


	// *********************************************************** 2 cell initialisation ***********************************************************************

	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::calendarAgeing is initialising the cell"<<endl;

	// set temperature
	c.setT(Ti);																		// set the Cell temperature
	c.setTenv(Ti);																	// set the environmental temperature

	// Get the battery to the resting voltage
	try{
		CC_V_CV_I(1.0, V, Ccut, dt, blockDegradation, &ahi, &whi, &ti); 			// do a 1C CC charge until the voltage, followed by a CV at this voltage
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in a subfunction of Cycler::CalendarAgeing when getting the cell to the resting voltage "<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// initial check-up
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::calendarAgeing is doing an initial check-up"<<endl;
		cap = checkUp(proc, 0, timetot, Ahtot, Whtot);
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in a subfunction of Cycler::CalendarAgeing in the initial check up "<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// *********************************************************** 3 calendar age the cell ***********************************************************************

	// loop for each resting period between two check-ups
	for (int i=0;i<nrdt;i++){
		try{
			// Calendar ageing type depending on 'mode'
			// rest the cell for the full period between consecutive check-ups and recharge at the end
			if (mode == 0){
				if(verbose >= printCyclerHighLevel)
					cout<<"Cycler::calendarAgeing is resting the cell in period "<<i<<endl;
				CC_t(0, dt, blockDegradation, trest, &ahi, &whi, &ti);				// rest (= CC charge at 0A)
				if(verbose >= printCyclerHighLevel)
					cout<<"Cycler::calendarAgeing is recharging the cell in period "<<i<<endl;
				CV_I(V, dt, blockDegradation, Ccut, &ahi2, &whi2, &ti2);			// recharge to the specified voltage
				timetot += (ti+ti2)/3600.0;											// number of hours we have rested additionally
			}
			// recharge every day
			else if (mode == 1){
				for(int j=0;j<timeCheck;i++){										// a loop for every day
					if(verbose >= printCyclerHighLevel)
						cout<<"Cycler::calendarAgeing is resting the cell in day "<<j<<" of period "<<i<<endl;
					CC_t(0, dt, blockDegradation, 24.0*3600.0, &ahi, &whi, &ti);	// rest for one day
					if(verbose >= printCyclerHighLevel)
						cout<<"Cycler::calendarAgeing is recharging the cell in day "<<j<<" of period "<<i<<endl;
					CV_I(V, dt, blockDegradation, Ccut, &ahi2, &whi2, &ti2);		// recharge to the specified voltage
					timetot += (ti+ti2)/3600.0;										// number of hours we have rested additionally
				}
			}
			// float at the set voltage (takes very long to calculate)
			else if (mode == 2){
				if(verbose >= printCyclerHighLevel)
					cout<<"Cycler::calendarAgeing is floating the cell in period "<<i<<endl;
				CV_t(V, dt, blockDegradation, trest, &ahi, &whi, &ti);				// do a CV (dis)charge for the entire trest period. Takes ages to compute
				timetot += (ti)/3600.0;												// number of hours we have rested additionally
			}
			else
				assert(false);														// not allowed, we have checked at the start that this can't happen

			// do a check-up
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::calendarAgeing is doing a check-up in period "<<i<<endl;
			cap = checkUp(proc, 0, timetot, Ahtot, Whtot);

			// End the experiment if the cell capacity has decreased too much
			if(cap<capnom/2.0){														// end simulation if the cell has only 50% of capacity left
				cout<<"Cycler::CalendarAgeing has finished calendar regime "<<ID<<" early do to too little remaining capacity.";
				cout<<" We have rested "<<i*timeCheck<<" days instead of "<<Time<<" and the remaining capacity now is "<<cap<<" [Ah]"<<endl<<flush;
				break;
			}
		} // end try block
		catch(int e){
			cout<<"Error in Cycler::CalendarAgeing while resting the cell according to calendar regime "<<ID<<". Error encountered is "<<e<<". Stop resting now.";
			cout<<" We have rested "<<i*timeCheck<<" days instead of "<<Time<<" and the capacity last measured is "<<cap<<endl<<flush;
			// we probably cannot do a full check-up because the cell is in an illegal state.
			// Therefore, only write the BatteryStates
			checkUp_batteryStates(proc.blockDegradation, false, 0, timetot, Ahtot, Whtot);
			break;
		} // end try-catch block

	} // end loop to rest and check-up

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::calendarAgeing terminating"<<endl;
}

void Cycler::profileAgeing(int length, string nameI, int limit, double Vma, double Vmi, double Ti, int nrProfiles, int nrCap, struct checkUpProcedure proc){
	/*
	 * Function to simulate degradation by continuously cycling a cell with a certain current profile.
	 * the profile is repeated until the cell is fully charged/discharged, then the cell is discharged/charged at 1C CC CV, and the profile is repeated again.
	 * The parameters of the cycling regime are set by the inputs.
	 * After a set number of profile repetitions, a check-up is done where the cell capacity, OCV curves, etc. are measured.
	 * These results are written to csv files.
	 *
	 * IN
	 * length	length of the current profile (i.e. number of rows in the csv file)
	 * nameI 	name of the CSV-file with the current profile
	 * 				the first column contains the current in [A], positive for discharge, negative for charge
	 * 				the second column contains the time in [sec] the current should be maintained
	 * limit 	integer describing what to do if the current can't be maintained because a voltage limit is reached
	 * 				0 	immediately go to the next current step of the profile
	 * 				1 	keep the voltage constant for the rest of this step of the profile
	 * Vma 		maximum voltage of this degradation experiment, must be below the maximum voltage of the cell [V]
	 * Vmi 		minimum voltage of this degradation experiment, must be above the minimum voltage of the cell [V]
	 * Ti 		temperature at which the test should be performed [K]
	 * nrProfiles the number of times the the current profile should be repeated
	 * nrCap 	the number of times the current profile is repeated approximately between consecutive check-ups
	 * 				a check-up is always done after a voltage limit was hit and we have re(dis)charged the cell
	 * 				e.g. you want to do a check-up every 20 cycles. Suppose you can do 15 cycles in the given voltage window (before you hit a voltage limit and you have to re(dis)charge the cell to the other voltage limit)
	 * 				then the check-up will be done after 30 cycles.
	 * proc 	structure with the parameters of the check-up procedure with the following fields:
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
	 *
	 *
	 * throws
	 * 1014 		the input parameters are invalid
	 * 1015 		the profile is invalid or the voltage limits are too small: the full voltage window was 'traversed' by just one repetition of the current profile
	 * 				i.e.: 	starting from a cell at the lower voltage, you still hit the upper voltage while following the profile
	 * 						or starting from a cell at the upper voltage, you still hit the lower voltage while following the profile
	 * 				reduce the absolute value of the currents in the current profile (or reduce the durations) to produce a valid current profile
	 */

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::profileAgeing starting"<<endl;

	// Check the input parameters
	bool vmax = Vma > c.getVmax(); 					// check if the maximum voltage is below the cell maximum voltage
	if (vmax)
		cerr<<"Error in Cycler::profileAgeing. The maximum voltage "<<Vma<<" is too high. The maximum value is "<<c.getVmax()<<endl<<flush;
	bool vmin = Vmi < c.getVmin();					// check if the minimum voltage is above the cell minimum voltage
	if (vmin)
		cerr<<"Error in Cycler::profileAgeing. The minimum voltage "<<Vmi<<" is too low. The minimum value is "<<c.getVmin()<<endl<<flush;
	bool Temin = Ti < TMIN;							// check the temperature is above 0 degrees, TMIN is defined in State.hpp
	if (Temin)
		cerr<<"Error in Cycler::profileAgeing. The temperature "<<Ti<<"K is too low. The minimum value is 273"<<endl<<flush;
	bool Temax = Ti > TMAX;							// check the temperature is below 60 degrees, TMIN is defined in State.hpp
	if (Temax)
		cerr<<"Error in Cycler::profileAgeing. The temperature "<<Ti<<" is too high. The maximum value is (273+60)"<<endl<<flush;
	bool cycles = nrProfiles <= nrCap;				// check the number of cycles between consecutive check-ups is lower than the total number of cycles
	if (cycles)
		cerr<<"Error in Cycler::profileAgeing. The number of cycles between two check ups "<<nrCap<<" is lower than the total number of cycles "<<nrProfiles<<endl<<flush;
	if(vmax || vmin || Temin || Temax || cycles)
		throw 1014;


	// *********************************************************** 1 variables & settings ***********************************************************************

	// Variables
	double capi;																		// capacity at this step [Ah]
	double timei;																		// time spent in this step [sec]
	double ahi;																			// charge discharged in this step [Ah]
	double whi;																			// energy discharged in this step [Wh]
	double Ahtot = 0;																	// cumulative charge throughput until now [Ah]
	double Whtot = 0;																	// cumulative energy throughput until now [Wh]
	double timetot = 0;																	// cumulative time until now [hour]
	int nrep = 0;																		// number of times the profile was repeated before a voltage limit was hit
	int nreptot = 0;																	// total number of times the profile was repeated
	int sign;																			// is the current profile a net charge (1) or discharge (-1)
	int vlim;																			// integer to indicate which voltage limit was hit while we were following the current profile
	bool Vlimhit = false;																// boolean to indicate if we need to re(dis)charge the cell to continue following the profile
	bool check = false;																	// boolean to indicate if we need to do a check-up after re(dis)charging the cell
	bool blockDegradation = false;														// account for degradation while following the current profile
	bool final = true;																	// boolean to check if a final check-up is needed
	double Ccc = 1;																		// C rate for the CC phase of the recharge or redischarge between profiles [-]
	double Ccut = 0.05;																	// C rate of the cutoff current for the CV phases  of the recharge or redischarge between profiles [-]

	// Read the current profile
	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::profileAgeing is reading the current profile"<<endl;
	double I[length], T[length];														// arrays to store the current profile as doubles
	try{
		loadCSV_2col(nameI, length, I, T);												// read the file
	}
	catch(int e){
		cout<<"error in Cycler::profileAgeing when reading the file with the current profile called "<<nameI<<", error "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// Determine if the profile is a net charge or a net discharge
	double aht = 0;																		// charge throughput of the profile
	for(int i=0;i<length;i++)
		aht += I[i] * T[i];
	if (aht > 0)																		// the profile is a net discharge
		sign = -1;
	else																				// the profile is net charge
		sign = 	1;

	// *********************************************************** 2 cell initialisation ***********************************************************************

	if(verbose >= printCyclerHighLevel)
		cout<<"Cycler::profileAgeing is initialising the cell"<<endl;

	// set the temperature
	c.setT(Ti);																			// set the cell temperature
	c.setTenv(Ti);																		// set the environmental temperature

	// (dis)charge the cell
	try{
	if (sign == -1)																		// we need to recharge to Vma because the profile is a discharge
		CC_V_CV_I(Ccc, Vma, Ccut, 2, blockDegradation, &ahi, &whi, &timei); 			// CC charge, followed by CV at a time step of 2 seconds
	else																				// we need to discharge to Vmi because the profile is a charge
		CC_V_CV_I(Ccc, Vmi, Ccut, 2, blockDegradation, &ahi, &whi, &timei); 			// CC discharge, followed by CV at a time step of 2 seconds
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in Cycler::profileAgeing when bringing the cell to the initial voltage "<<e<<". Throwing it on."<<endl<<flush;
		throw e;
	}

	// initial check up
	try{
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::profileAgeing is doing an initial check-up"<<endl;
		checkUp(proc, 0, timetot, Ahtot, Whtot);
	}
	catch(int e){
		if(verbose >= printCrit)
			cout<<"Error in a subfunction of Cycler::profileAgeing in the initial check-up "<< e<<". Throwing it on "<<endl;
		throw e;
	}

	// *********************************************************** 3 age the cell by continuously following the profile and re(dis)charging ***********************************************************************

	// Cycle the battery by following the profile by repeating the 3 steps
	// 		keep applying the profile until you hit a voltage limit
	// 		re(dis)charge the cell
	// 		do a check-up if needed
	while (nreptot < nrProfiles){														// loop to follow the profile the set number of times

		// *********************************************************** 3A  follow the profile until you hit a voltage limit ***********************************************************************
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::profileAgeing has applied the profile "<<nreptot<<" times and is starting the next series of profile repetitions"<<endl;
		try{
			Vlimhit = false;															// boolean to check when we've hit the voltage limit
			nrep = 0;																	// counter for how often we can follow the profile before hitting the voltage limit
			check = false;																// boolean to check if we need to do a check-up
			while(!Vlimhit){															// loop to keep applying the profile until you hit a voltage limit

				// follow the profile
				vlim = followI(length, nameI, blockDegradation, limit, Vma, Vmi, &ahi, &whi, &timei);

				// update the throughput
				Ahtot += abs(ahi);
				Whtot += abs(whi);
				timetot += timei/3600.0;
				nrep++;																	// increase the counter of repetitions before hitting the voltage limit
				nreptot++;																// increase the total counter
				if (fmod(nreptot,nrCap) == 0)											// the total number of repetitions is a multiple of the nr profiles between a check, so store that we need to do a check-up
					check = true;

				// check if we have hit the voltage limit
				if(sign == -1)															// the profile is a net discharge, so stop if you hit the lower voltage limit
					Vlimhit = (vlim == -1) || (vlim==10);								// a return-integer of -1 or 10 means the lower voltage was hit
				else
					Vlimhit = (vlim == 1) || (vlim==10);								// a return-integer of 1 or 10 means the upper voltage was hit
			}
		}
		catch(int e){
			if(verbose >= printCrit){
				cout<<"Error in Cycler::profileAgeing when cycling the cell with the current profile according to profile regime "<<ID<<". Error encountered is "<<e<<". Stop cycling now.";
				cout<<" We have done "<<nreptot<<" repetitions instead of "<<nrProfiles<<" and the capacity last measured is "<<capi<<endl<<flush;
			}

			// we probably cannot do a full check-up procedure because the cell is in an illegal state.
			// Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
			checkUp_batteryStates(proc.blockDegradation, false, nreptot, timetot, Ahtot, Whtot);

			final = false;																// skip the final check-up
			break;																		// stop cycling
		}

		// ensure we could repeat the profile multiple times before hitting the voltage limit
		if(nrep <= 1){
			// this means that you can never fully follow the current profile in the specified voltage window.
			// this function does not allow this sort of behaviour, and an error is thrown.
			// the user has to decrease the currents in the profile, or enlarge the voltage window (if possible)
			cerr<<"ERROR in Cycler::profileAgeing: the profile "<<ID<<" could only be repeated once before a voltage limit was hit."<<endl<<flush;
			// throw the error
			throw 1015;
		}

		// *********************************************************** 3B re(dis)charge ***********************************************************************
		// re(dis)charge the cell to the other voltage limit from the one you have hit
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::profileAgeing has applied the profile "<<nreptot<<" times and is going to re(dis)charge the cell"<<endl;
		try{
			if (sign == -1)																// we need to recharge to Vma
				CC_V_CV_I(Ccc, Vma, Ccut, 2, blockDegradation, &ahi, &whi, &timei); 	// CC charge, followed by CV at a time step of 2 seconds
			else																		// we need to discharge to Vmi
				CC_V_CV_I(Ccc, Vmi, Ccut, 2, blockDegradation, &ahi, &whi, &timei); 	// CC discharge, followed by CV at a time step of 2 seconds
			Ahtot += abs(ahi);
			Whtot += abs(whi);
			timetot += timei/3600.0;
		}
		catch(int e){
			if(verbose >= printCrit){
				cout<<"Error in Cycler::profileAgeing when re(dis)charging the cell according to profile regime "<<ID<<". Error encountered is "<<e<<". Stop cycling now.";
				cout<<" We have done "<<nreptot<<" repetitions instead of "<<nrProfiles<<" and the capacity last measured is "<<capi<<endl<<flush;
			}

			// we probably cannot do a full check-up procedure because the cell is in an illegal state.
			// Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
			checkUp_batteryStates(proc.blockDegradation, false, nreptot, timetot, Ahtot, Whtot);

			final = false;																// skip the final check-up
			break;																		// stop cycling
		}

		// *********************************************************** 3C check-up ***********************************************************************
		// do a check-up if needed
		if(check){
			if(verbose >= printCyclerHighLevel)
				cout<<"Cycler::profileAgeing has applied the profile "<<nreptot<<" times and is going to do a check-up"<<endl;
			try{
				capi = checkUp(proc, nreptot, timetot, Ahtot, Whtot);
				check = false;
			}
			catch(int e){
				if(verbose >= printCrit){
					cout<<"Error in Cycler::profileAgeing when doing a check-up according to profile regime "<<ID<<". Error encountered is "<<e<<". Stop cycling now.";
					cout<<" We have done "<<nreptot<<" repetitions instead of "<<nrProfiles<<" and the capacity last measured is "<<capi<<endl<<flush;
				}

				// we probably cannot do a full check-up procedure because the cell is in an illegal state.
				// Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
				checkUp_batteryStates(proc.blockDegradation, false, nreptot, timetot, Ahtot, Whtot);

				final = false;															// skip the final check-up
				break;																	// stop cycling
			}

			// End the experiment if the cell capacity has decreased too much
			if(capi<c.getNominalCap()/2.0){												// end simulation if the cell has only 50% of capacity left
				cout<<"Cycler::ProfileAgeing has finished cycling regime "<<ID<<" early do to too little remaining capacity.";
				cout<<" We have done "<<nreptot<<" repetitions of the profile instead of "<<nrProfiles<<" and the remaining capacity now is "<<capi<<endl<<flush;
				final = false;															// skip the final check-up because we just did one
				break;																	// stop cycling
			}
		}

		// keep repeating these 3 steps (repeat profile until voltage limit, re(dis)charge, check-up) until you have done enough repetitions

	} // end loop of profile ageing


	// *********************************************************** 4 final check-up ***********************************************************************
	// Don't do the check-up if the simulation ended due to an error (or shortage in capacity)
	if (final){
		if(verbose >= printCyclerHighLevel)
			cout<<"Cycler::profileAgeing has applied the profile "<<nreptot<<" times and is going to do a final check-up"<<endl;
		try{
			checkUp(proc, nreptot, timetot, Ahtot, Whtot);
		}
		catch(int e){
			if(verbose >= printCrit)
				cout<<"Error in a subfunction of Cycler::profileAgeing in the final check-up "<< e<<". Throwing it on "<<endl;
			throw e;
		}
	}

	if(verbose >= printCyclerFunctions)
		cout<<"Cycler::profileAgeing terminating"<<endl;
}
