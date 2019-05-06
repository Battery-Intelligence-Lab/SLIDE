/*
 * State.hpp
 *
 * Defines a class State which defines the state-variables of a cell for the state-space model formulation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#ifndef SRC_STATE_HPP_
#define SRC_STATE_HPP_

#define nch 5			// number of points in the spatial discretisation of the solid diffusion equation *** DON'T CHANGE THE VALUE ***
						// this is the number of positive inner Chebyshev nodes
						// 		the full chebyshev interval is from x = -1 to x = 1
						// 		the positive points go from x = 0 to x = 1
						// 		the inner positive points are the positive points excluding the point at x=0 and at x=1
						// 		so nch is the number of Chebyshev points with 0 < x < 1
						// do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied Matlab scripts. See the word document '2 overview of the code', section 'Matlab setup before running the C++ code'
#define ns (2*nch+14) 	// number of state variables
						// the parenthesis are important because this code segment is copy-pasted
						// so if you have somewhere in the code '1/ns' this becomes '1/(2*nch+14)'
						// while without the parenthesis it becomes '1/2*nch + 14'
#define TMIN (273+0) 	// the minimum temperature allowed in the simulation [K]
#define TMAX (273+60) 	// the maximum temperature allowed in the simulation [K]

namespace std {

class State {
public:
	State();											// constructor which initialises all state variables to 0
	virtual ~State();

	virtual void initialise(int nin, double zpi[], double zni[], double Ti, double deltai, double LLIi,
			double thickpi, double thickni, double epi, double eni, double api, double ani,
			double CSi, double Dpi, double Dni, double Ri, double deltalii); // initialises all state variables to the given values

	virtual void getStates(int nin, double states[]);	// get the battery states in an array
	virtual void getStates(State& si);					// get the battery states in another State-object
	void getZ(int nin, double zpi[], double zni[]);		// get the transformed concentrations
	double getT();										// get the temperature
	double getDelta();									// get the SEI thickness
	double getLLI();									// get the lost lithium
	double getThickp();									// get the thickness of the cathode
	double getThickn();									// get the thickness of the anode
	double getEp();										// get the volume fraction of active material of the cathode
	double getEn();										// get the volume fraction of active material of the anode
	double getAp();										// get the effective surface of the cathode
	double getAn();										// get the effective surface of the anode
	double getCS();										// get the crack surface area
	double getDp();										// get the diffusion constant of the cathode
	double getDn();										// get the diffusion constant of the anode
	double getr();										// get the specific resistance
	double getR(double elec_surf);						// get the total DC resistance
	double getDelta_pl();								// get the thickness of the plated lithium layer
	void getIniStates(int nin, double si[]);			// get the initial states

	virtual void setT(double Ti);						// set the temperature
	virtual void setZ(int nin, double zpi[], double zni[]);	// set the transformed concentration
	virtual void setStates(int nin, double states[]);	// set the states to the values in the array
	virtual void setStates(State si);					// set the states and initial states to the values of the given State-object
	virtual void setIniStates(int nin, double si[]);	// set the initial states to the values in the array
	virtual void overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani); // overwrite the states related to the geometry of a cell
	virtual void overwriteCharacterisationStates(double Dpi, double Dni, double ri); // overwrite the states related to the characterisation of a cell

	virtual void validState();							// check if the states are valid

private:
	// battery states
	// zp[nch] zn[nch] T delta LLI thickp thickn ep en ap an CS Dp Dn R delta_liPlating

	double zp[nch];		// twice transformed concentration at positive electrode
	double zn[nch];		// twice transformed concentration at negative electrode

	double T;			// battery temperature [K]

	double Dp; 			// diffusion constant of positive electrode at reference temperature [m/s]
	double Dn; 			// diffusion constant of negative electrode at reference temperature [m/s]

	double r;			// resistance times real surface area of the combined electrodes [Ohm m2]
						// often, you don't have data for the resistance of individual electrodes
						// so this model uses the 'average' resistance of both electrodes, which can be obtained from the measured cell DC resistance

	double delta;		// SEI thickness [m]
	double delta_pl;	// thickness of the lithium-plating layer [m]
	double LLI;			// lost Li [As]

	double thickp; 		// thickness of positive electrode [m]
	double thickn; 		// thickness of negative electrode [m]
	double ep; 			// volume fraction of active material in the positive electrode
	double en; 			// volume fraction of active material in the negative electrode
	double ap; 			// effective surface area of positive electrode
	double an; 			// effective surface area of negative electrode

	double CS;			// crack surface [m2]

	double sini[ns];	// array with the initial battery states
};

} /* namespace std */

#endif /* SRC_STATE_HPP_ */
