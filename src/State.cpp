/*
 * State.cpp
 *
 * Implements a class State which defines the state-variables of a cell for the state-space model formulation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "State.hpp"

#include <cassert>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <math.h>

namespace std {

State::State(){
	/*
	 * Default constructor which DOESN'T initialise the states.
	 * All states are set to 0
	 */

	for(int i=0;i<nch;i++){
		zp[i] = 0;
		zn[i] = 0;
	}
	T		= 0;
	delta 	= 0;
	LLI		= 0;
	thickp	= 0;
	thickn 	= 0;
	ep 		= 0;
	en 		= 0;
	ap		= 0;
	an 		= 0;
	CS		= 0;
	Dp 		= 0;
	Dn 		= 0;
	r		= 0;
	delta_pl = 0;
	for(int i=0;i<ns;i++)
		sini[i] = 0;
}

void State::initialise(int nin, double zpi[], double zni[], double Ti, double deltai, double LLIi,
			double thickpi, double thickni, double epi, double eni, double api, double ani,
			double CSi, double Dpi, double Dni, double ri, double deltalii){
	/*
	 * Initialise the state variables to the given values.
	 * Use this function only once, immediately after calling the constructor.
	 * Calling it again at a later point will generate an error (11)
	 *
	 * IN
	 * nin 		length of the arrays with the (transformed) concentration
	 * zpi 		transformed concentration at the positive inner Chebyshev nodes of the positive particle
	 * zni		transformed concentration at the positive inner Chebyshev nodes of the negative particle
	 * Ti 		cell temperature [K]
	 * deltai 	thickness of the SEI layer [m]
	 * LLIi 	lost lithium inventory [As]
	 * thickpi 	thickness of the cathode [m]
	 * thickni 	thickness of the anode [m]
	 * epi 		volume fraction of active material in the cathode [-]
	 * eni 		volume fraction of active material in the anode [-]
	 * api 		effective surface area of the porous cathode [m2 m-3]
	 * ani 		effective surface area of the porous anode [m2 m-3]
	 * CSi 		surface area of the cracks at the surface of the negative particle [m2]
	 * Dpi 		diffusion constant of the cathode at reference temperature [m s-1]
	 * Dni 		diffusion constant of the anode at reference temperature [m s-1]
	 * ri 		specific resistance of both electrodes combined [Ohm m2]
	 * deltalii thickness of the plated lithium layer [m]
	 *
	 * Note on ri:
	 * this is the resistance times the electrode surface, averaged between both electrodes.
	 * The total cell resistance is (see Cell::getR() ): r /( (thickp*ap*elec_surf + thickn*an*elec_surf)/2 )
	 * 	with r the resistance times the average electrode surface
	 * 		 thicki the thickness of electrode i
	 * 		 ai the specific surface area of electrode i
	 * 		 elec_surf the geometric surface area of the electrode (height of the electrode * width of the electrode)
	 * so if the measured DC resistance of the cell is R, the value of r can be calculated using:
	 * 		 ri = R * ( (thickp*ap*elec_surf + thickn*an*elec_surf)/2 )
	 *
	 * THROWS
	 * 10		the arrays have the wrong length
	 * 11 		this function is called when the initial states have already been initialised
	 * 12 		the state suggested is illegal
	 */

	// Check that the input arrays have the correct length
	if (nin!=nch){
		cerr<<"ERROR in State::initialise, the input arrays have a length of "<<nin<<" instead of "<<nch<<". Throwing an error"<<endl<<flush;
		throw 10;
	}

	// Check that the initial states haven't been stored yet
	bool first = false;
	for (int i=0;i<ns;i++)
		first = first || (sini[i] != 0);	// the constructor has set the initial states to 0, so if one of them is not 0, the initial states have been set
	if(first){
		cerr<<"ERROR in State::initialise, the initial states had already been set so you can't set them again"<<endl<<flush;
		throw 11;
	}

	// Set the state variables
	for(int i=0;i<nch;i++){
		zp[i] = zpi[i];
		zn[i] = zni[i];
	}
	T		= Ti;
	delta 	= deltai;
	LLI		= LLIi;
	thickp	= thickpi;
	thickn 	= thickni;
	ep 		= epi;
	en 		= eni;
	ap		= api;
	an 		= ani;
	CS		= CSi;
	Dp 		= Dpi;
	Dn 		= Dni;
	r		= ri;
	delta_pl = deltalii;

	// Store the initial states in a separate variable which won't ever be changed
	double inistates[ns];
	getStates(ns, inistates);			// get the values we just gave to the State
	setIniStates(ns, inistates);		// store them as initial states

	// Check if this was a valid state
	try{
		validState();
	}
	catch(int e){
		cout<<"Error in State::initialise, one of the states has an illegal value, throwing an error"<<endl<<flush;
		throw 12;
	}
}

void State::getStates(int nin, double states[]){
	/*
	 * Returns an array with all the states.
	 *
	 * IN
	 * nin 		length of the array provided, must be ns
	 *
	 * OUT
	 * states	array in which the battery states will be put
	 * 				zp 			the transformed li concentration at the positive inner nodes of the positive particle (nch values)
	 * 				zn			the transformed li concentration at the positive inner nodes of the negative particle (nch values)
	 * 				T 			the cell temperature [K]
	 * 				delta 		the thickness of the SEI layer [m]
	 * 				LLI 		the lost lithium [As]
	 * 				thickp 		the thickness of the cathode [m]
	 * 				thickn		the thickness of the anode [m]
	 * 				ep 			the volume fraction of active material in the cathode [-]
	 * 				en 			the volume fraction of active material in the anode [-]
	 * 				ap 			the effective surface area of the cathode [m2 m-3]
	 * 				an			the effective surface area of the anode [m2 m-3]
	 * 				CS			the surface area of the cracks at the surface of the negative particle [m2]
	 * 				Dp			the diffusion constant at reference temperature of the cathode [m s-1]
	 * 				Dn			the diffusion constant at reference temperature of the anode [m s-1]
	 * 				r 			the specific resistance of the combined electrodes [Ohm m2]
	 * 				delta_pl 	the thickness of the plated lithium layer [m]
	 *
	 * THROWS
	 * 10		the array has the wrong length
	 */

	// Check that the input arrays have the correct length
	if (nin!=ns){
		cerr<<"ERROR in State::getStates, the input arrays have a length of "<<nin<<" instead of "<<ns<<". Throwing an error"<<endl<<flush;
		throw 10;
	}

	// Put the states in the array
	for (int j=0;j<nch;j++){
		states[j] = zp[j];
		states[nch+j] = zn[j];
	}
	states[2*nch + 0] = T;
	states[2*nch + 1] = delta;
	states[2*nch + 2] = LLI;
	states[2*nch + 3] = thickp;
	states[2*nch + 4] = thickn;
	states[2*nch + 5] = ep;
	states[2*nch + 6] = en;
	states[2*nch + 7] = ap;
	states[2*nch + 8] = an;
	states[2*nch + 9] = CS;
	states[2*nch + 10] = Dp;
	states[2*nch + 11] = Dn;
	states[2*nch + 12] = r;
	states[2*nch + 13] = delta_pl;

	assert(ns == 2*nch + 14); 	// check that ns has the correct value. We have copied (2*nch + 14) values so that must be the value of ns
	// if this assertion fails, the user has changed something in the code at some point, without accounting for this change somewhere else.
	// e.g. if you add an extra state-variable, you have to increase the value of 'ns' (defined in State.hpp), and add it in all functions in State.
}

void State::getStates(State& si){
	/*
	 * Copy all the battery states to the given State object.
	 * Also the values of the initial states are copied (which are not included in the array representation)
	 *
	 * OUT
	 * si 	State-object (called by reference) to which the states will be copied
	 * 		the initial states of si will be overwritten by this function.
	 * 		this means that si must satisfy the conditions of setIniStates:
	 * 			si must be a newly made State (the initial states are 0)
	 * 			or the initial states of si must have the same value as the initial states of this state
	 * 				This latter conditions allows to reuse the same State-object to get the states of one Cell-object
	 */

	// copy the initial battery states
	try{
		si.setIniStates(ns, sini);			// copy the initial states from this State to si
	}
	catch(int e){
		cout<<"Error in getStates(State& si) when setting the initial states in the new state object: "<<e<<". Throwing it on"<<endl;
		throw e;
	}

	// copy the states using the array representation
	double states[ns];
	try{
		getStates(ns, states);
		si.setStates(ns, states);			// if you haven't set the initial states, this will produce an error because the state is not valid
	}
	catch(int e){
		cout<<"Error in State::getStates(State& si) when setting the states in the new object: "<<e<<". Throwing it on"<<endl;
		throw e;
	}
}

void State::getZ(int nin, double zpi[], double zni[]){
	/*
	 * Function to get the transformed lithium concentrations at the positive inner nodes of the particles
	 *
	 * IN
	 * nin 	length of the arrays, must be nch
	 *
	 * OUT
	 * zpi	transformed concentration at the positive inner nodes of the positive particles
	 * zni	transformed concentration at the positive inner nodes of the negative particles
	 *
	 * THROWS
	 * 10		the arrays have the wrong length
	 */

	// Check that the input arrays have the correct length
	if (nin!=nch){
		cerr<<"ERROR in State::getZ, the input arrays have a length of "<<nin<<" instead of "<<nch<<". Throwing an error"<<endl<<flush;
		throw 10;
	}

	// copy the states to the arrays
	for(int i=0;i<nch;i++){
		zpi[i] = zp[i];
		zni[i] = zn[i];
	}
}
double State::getT(){
	return T;
}
double State::getDelta(){
	return delta;
}
double State::getLLI(){
	return LLI;
}
double State::getThickp(){
	return thickp;
}
double State::getThickn(){
	return thickn;
}
double State::getEp(){
	return ep;
}
double State::getEn(){
	return en;
}
double State::getAp(){
	return ap;
}
double State::getAn(){
	return an;
}
double State::getCS(){
	return CS;
}
double State::getDp(){
	return Dp;
}
double State::getDn(){
	return Dn;
}
double State::getr(){
	return r;
}
double State::getR(double elec_surf){
	/*
	 * Returns the total DC resistance of the electrodes.
	 *
	 * IN
	 * elec_surf 	geometric surface area of the electrodes [m2]
	 *
	 * OUT
	 * resistance of the electrodes [Ohm]
	 */

	double surfp = ap*thickp*elec_surf; 	// real surface area of the positive electrode
	double surfn = an*thickn*elec_surf;		// real surface area of the negative electrode
	double surf = (surfp + surfn)/2;		// mean surface area
	return r / surf;
}
double State::getDelta_pl(){
	return delta_pl;
}
void State::getIniStates(int nin, double si[]){
	/*
	 * returns the initial states
	 *
	 * IN
	 * nin 	length of the array provided, must be ns
	 *
	 * OUT
	 * si 	array with the initial states of this State-object
	 *
	 * THROWS
	 * 10 	if the arrays have the wrong length
	 */

	// Check that the input arrays have the correct length
	if (nin!=ns){
		cerr<<"ERROR in State::getIniStates, the input arrays have a length of "<<nin<<" instead of "<<ns<<". Throwing an error"<<endl<<flush;
		throw 10;
	}

	// copy the initial states
	for(int i=0;i<ns;i++)
		si[i] = sini[i];
}

void State::setT(double Ti){
	/*
	 * Sets the temperature to the given value
	 *
	 * IN
	 * Ti		temperature [K], must be between 0 and 60 degrees, so 273 and 273+60
	 *
	 * THROWS
	 * 13 		illegal value of T
	 */

	// the temperature limits are defined in State.hpp
	bool Tmin = Ti < TMIN;							// check the temperature is above 0 degrees
	if (Tmin)
		cerr<<"Error in State::setT. The temperature "<<Ti<<"K is too low. The minimum value is "<<TMIN<<endl<<flush;
	bool Tmax = Ti > TMAX;							// check the temperature is below 60 degrees
	if (Tmax)
		cerr<<"Error in State::setT. The temperature "<<Ti<<"K is too high. The maximum value is "<<TMAX<<endl<<flush;
	if(Tmin || Tmax)
		throw 13;

	T = Ti;
}
void State::setZ(int nin, double zpi[], double zni[]){
	/*
	 * Sets the transformed concentration of the particles.
	 *
	 * IN
	 * nin	length of the arrays, must be nch
	 * zpi 	transformed concentrations at the positive inner nodes in the positive particle
	 * zni 	transformed concentrations at the positive inner nodes in the negative particle
	 *
	 * THROWS
	 * 10 	if the arrays have the wrong length
	 */

	// Check that the input arrays have the correct length
	if (nin!=nch){
		cerr<<"ERROR in State::setZ, the input arrays have a length of "<<nin<<" instead of "<<nch<<". Throwing an error"<<endl<<flush;
		throw 10;
	}

	// copy the concentrations
	for (int j=0;j<nch;j++){
		zp[j] = zpi[j];
		zn[j] = zni[j];
	}
}
void State::setStates(int nin, double states[]){
	/*
	 * sets the states to the values given in the array.
	 * At the end, it is checked if the state is allowed, and an error (15) is thrown up if not
	 *
	 * IN
	 * nin 		length of the array, must be ns
	 * states 	array with the states
	 *
	 * THROWS
	 * 10 		if the arrays have the wrong length
	 */

	// Check that the input arrays have the correct length
	if (nin!=ns){
		cerr<<"ERROR in State::setStates, the input arrays have a length of "<<nin<<" instead of "<<ns<<". Throwing an error"<<endl<<flush;
		throw 10;
	}

	// set the state variables
	for (int j=0;j<nch;j++){
		zp[j] = states[j];
		zn[j] = states[nch+j];
	}
	T 		= states[2*nch + 0];
	delta 	= states[2*nch + 1];
	LLI 	= states[2*nch + 2];
	thickp 	= states[2*nch + 3];
	thickn 	= states[2*nch + 4];
	ep 		= states[2*nch + 5];
	en 		= states[2*nch + 6];
	ap 		= states[2*nch + 7];
	an 		= states[2*nch + 8];
	CS 		= states[2*nch + 9];
	Dp 		= states[2*nch + 10];
	Dn 		= states[2*nch + 11];
	r 		= states[2*nch + 12];
	delta_pl= states[2*nch + 13];

	assert(ns == 2*nch + 14); 	// check that ns has the correct value. We have copied (2*nch + 13 + 1) values so that must be the value of ns
	// if this assertion fails, the user has changed something in the code at some point, without accounting for this change somewhere else.
	// e.g. if you add an extra state-variable, you have to increase the value of 'ns' (defined in State.hpp), and add it in all functions in State.

	try{
		validState(); 				// throw 15 if the state is illegal
	}
	catch(int e){
		cout<<"Error in State::setStates(double states[]), the suggested state is illegal: "<<e<<". throwing it on"<<endl<<flush;
		throw e;
	}
}

void State::setStates(State si){
	/*
	 * Copies all the battery states from the given one to the this battery state.
	 * The initial states are not changed, since they can only be set once using the function iniStates
	 *
	 * IN
	 * si 	State-object (called by reference) from which the states will be copied to this State object
	 */

	// get the states which require an array
	double zpi[nch], zni[nch];
	try{
		si.getZ(nch, zpi, zni); 		// transformed concentrations
	}
	catch(int e){
		cout<<"Error in State::setStates(State& si) when getting the concentrations and the initial states: "<<e<<". Throwing it on"<<endl;
		throw e;
	}

	// copy the main battery states
	for (int j=0;j<nch;j++){
		zp[j] = zpi[j];
		zn[j] = zni[j];
	}
	T 		= si.getT();
	delta 	= si.getDelta();
	LLI 	= si.getLLI();
	thickp 	= si.getThickp();
	thickn 	= si.getThickn();
	ep 		= si.getEp();
	en 		= si.getEn();
	ap 		= si.getAp();
	an 		= si.getAn();
	CS 		= si.getCS();
	Dp 		= si.getDp();
	Dn 		= si.getDn();
	r 		= si.getr();
	delta_pl= si.getDelta_pl();

	assert(ns == 2*nch + 14); 	// check that ns has the correct value. We have copied (2*nch + 13 + 1) values so that must be the value of ns
	// if this assertion fails, the user has changed something in the code at some point, without accounting for this change somewhere else.
	// e.g. if you add an extra state-variable, you have to increase the value of 'ns' (defined in State.hpp), and add it in all functions in State.
}

void State::setIniStates(int nin, double si[]){
	/*
	 * Function to set the initial states of this State object
	 * Setting the initial states is only allowed in two cases:
	 * 		this State-object must be a newly made State (the initial states are 0)
	 * 		or the initial states of this State-object must have the same value as the initial states of si
	 * 			i.e. the initial states of this State-object are not changed by this function
	 *
	 * IN
	 * nin 	length of the array
	 * si 	array with the initial states
	 *
	 * THROWS
	 * 10 	if the arrays have the wrong length
	 * 14 	this State-object does not satisfy either of the conditions:
	 * 				this state object is not a newly made state (i.e. the initial states have already been set)
	 * 				and the initial states of this state object are different from the ones in the array
	 * 16 	At least one of the states (except the concentration) is invalid, ie has a negative value
	 */

	// Check that the input arrays have the correct length
	if (nin!=ns){
		cerr<<"ERROR in State::setIniStates, the input arrays have a length of "<<nin<<" instead of "<<ns<<". Throwing an error"<<endl<<flush;
		throw 10;
	}

	// Check if setting the initial states is allowed, i.e. if they haven't been set yet or if the values don't change
	bool illegal = false;									// boolean indicating whether changing the initial states is illegal
	bool newi, samei;
	for (int i=0;i<ns;i++){
		newi = (sini[i] == 0);								// initial state was 0, i.e. it hadn't been set yet
		samei = ( abs((sini[i]-si[i])/sini[i]) < pow(10,-6)); // the initial state had the same value (with a e-6 relative tolerance) so won't be changed
		illegal = illegal || !(newi || samei);				// to be legal, the state must be new or identical
	}
	if(illegal){
		cerr<<"ERROR in State::setIniStates, the initial states had already been set so you can't set them again"<<endl<<flush;
		throw 14;
	}

	// Check the initial states are allowed
	// We can't say much about the (transformed) concentration. But all other states have to be positive at the very least
	bool negi;
	bool negtot = false;
	for (int i=0;i<14;i++){
		negi = si[2*nch + i] < 0;
		negtot = negtot | negi;
		if (negi)
			cerr<<"Error in State::setIniStates. State "<<2*nch+i<<" has as value "<<si[2*nch + i]  <<". Negative values are not allowed"<<endl<<flush;
	}
	if(negtot)
		throw 16;

	// copy the initial states
	for(int i=0;i<ns;i++)
		sini[i] = si[i];
}

void State::overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani){
	/*
	 * Function to overwrite the geometric parameters of the state.
	 * It also overwrites the initial states, so use it with extreme caution.
	 * It should only be called when you are parametrising a cell (determineCharacterisation.cpp), never while cycling a cell.
	 *
	 * IN
	 * thickpi 	thickness of the cathode [m]
	 * thickni 	thickness of the anode [m]
	 * epi 		volume fraction of active material in the cathode [-]
	 * eni 		volume fraction of active material in the anode [-]
	 * api 		effective surface area of the porous cathode [m2 m-3]
	 * ani 		effective surface area of the porous anode [m2 m-3]
	 */

	// set the states
	thickp = thickpi;
	thickn = thickni;
	ep = epi;
	en = eni;
	ap = api;
	an = ani;

	// overwrite the corresponding initial states
	sini[2*nch + 3] = thickp;
	sini[2*nch + 4] = thickn;
	sini[2*nch + 5] = ep;
	sini[2*nch + 6] = en;
	sini[2*nch + 7] = ap;
	sini[2*nch + 8] = an;

	assert(ns == 2*nch + 14); 	// check that ns has the correct value. Else we might have copied states to the wrong location
}

void State::overwriteCharacterisationStates(double Dpi, double Dni, double ri){
	/*
	 * Function to overwrite the parameters related to the characterisation of the cell.
	 * The states and initial states are overwritten so use this function with caution.
	 * It should only be called when you are parametrising a cell (determineCharacterisation.cpp), never while cycling a cell.
	 *
	 * IN
	 * Dpi	diffusion constant of the cathode at rate temperature [m s-1]
	 * Dni 	diffusion constant of the anode at rate temperature [m s-1]
	 * r 	specific resistance of the combined electrodes [Ohm m2]
	 */

	// Set the states
	Dp = Dpi;
	Dn = Dni;
	r = ri;

	// overwrite the corresponding initial states
	sini[2*nch + 10] = Dp;
	sini[2*nch + 11] = Dn;
	sini[2*nch + 12] = r;

	assert(ns == 2*nch + 14); 	// check that ns has the correct value. Else we might have copied states to the wrong location
}

void State::validState(){
	/*
	 * Check if this State object has valid parameters.
	 * Errors are thrown if they are not valid.
	 *
	 * Note to the user: most of these limits are purely intended to indicate something strange is going on.
	 * I.e. exceeding the limits below normally doesn't lead to errors in the code.
	 * They just indicate things which should not be happening (e.g. increasing the amount of active material).
	 * If a limit is 'critical' (i.e. exceeding this limit will give errors in the code), this is indicated in comments below
	 *
	 * The main purpose of the limits on degradation parameters (e.g. too low amount of active material)
	 * is to avoid a very long, or even infinitely long, calculation (e.g. when the cell simply can't reach a certain voltage because there is not enough active material)
	 * or to avoid problems with the time discretisation
	 * (e.g. the lower amount of active material, the larger the concentration difference for constant current and time step;
	 * at some point, the concentration difference over one time step will become too large, and lead to numerical instability or other discretisation errors)
	 *
	 * THROWS
	 * 15 	illegal state. The error message will detail which state has an illegal value
	 */

	// There are no limits on the transformed concentration, because this is the (twice) transformed concentration.
	// The limits are on the real concentrations, which must be calculated first.
	// This can't be done here because it requires parameters of the Cell itself.
	// see Cell::getC or Cell::getCsurf

	// temperature of the cell (in Kelvin), the temperature limits are defined in State.hpp
	bool Tmin = T < TMIN;							// check the temperature is above 0 degrees
	if (Tmin)
		cerr<<"Error in State::validState. The temperature "<<T<<"K is too low. The minimum value is "<<TMIN<<endl<<flush;
	bool Tmax = T > TMAX;							// check the temperature is below 60 degrees
	if (Tmax)
		cerr<<"Error in State::validState. The temperature "<<T<<" is too high. The maximum value is "<<TMAX<<endl<<flush;

	// thickness of the SEI layer
	bool del = delta <= 0;
	if (del)
		cerr<<"Error in State::validState. The SEI thickness delta is "<<delta<<", which is too low. Only strictly positive values are allowed"<<endl<<flush;
		// a value of 0 gives problems in some equations, which have a term 1/delta, which would give nan or inf if delta = 0
		// a negative value might lead to a further decrease in SEI thickness, so it will keep getting more and more negative

	// lost lithium
	bool li = LLI < 0;
	if (li)
		cerr<<"Error in State::validState. The lost lithium LLI is "<<LLI<<", which is too low. Only non-negative values are allowed"<<endl<<flush;
		// a value of 0 is allowed (it won't lead to errors in the code)

	// thickness of the electrodes
	bool tpmin = thickp <= sini[2*nch + 3] / 5;
	if (tpmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the thickness of the positive electrode is "<<thickp<<", which is too low. The minimum is "<<sini[2*nch + 3] / 5 <<", 1/5 of the original thickness"<<endl<<flush;
		// errors will happen if the thickness becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool tpmax = thickp > sini[2*nch + 3]*1.001; // leave a 0.1% margin for numerical errors
	if (tpmax)
		cerr<<"Error in State::validState. The thickness of the positive electrode is "<<thickp<<", which is too high. The maximum is "<<sini[2*nch + 3] <<", the original thickness"<<endl<<flush;
	bool tnmin = thickn <= sini[2*nch + 4] / 5;
	if (tnmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the thickness of the negative electrode is "<<thickn<<", which is too low. The minimum is "<<sini[2*nch + 4] / 5 <<", 1/5 of the original thickness"<<endl<<flush;
		// errors will happen if the thickness becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool tnmax = thickn > sini[2*nch + 4]*1.001;
	if (tnmax)
		cerr<<"Error in State::validState. The thickness of the negative electrode is "<<thickn<<", which is too high. The maximum is "<<sini[2*nch + 4] <<", the original thickness"<<endl<<flush;

	// volume fraction of the active material in the electrodes
	bool epmin = ep <= sini[2*nch + 5] / 5;
	if (epmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the volume fraction of the positive electrode is "<<ep<<", which is too low. The minimum is "<<sini[2*nch + 5] / 5 <<", 1/5 of the original volume fraction"<<endl<<flush;
		// errors will happen if the volume fraction becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool epmax = ep > sini[2*nch + 5]*1.001;
	if (epmax)
		cerr<<"Error in State::validState. The volume fraction of the positive electrode is "<<ep<<", which is too high. The maximum is "<<sini[2*nch + 5] <<", the original volume fraction"<<endl<<flush;
	bool enmin = en <= sini[2*nch + 6] / 5;
	if (enmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the volume fraction of the negative electrode is "<<en<<", which is too low. The minimum is "<<sini[2*nch + 6] / 5 <<", 1/5 of the original volume fraction"<<endl<<flush;
		// errors will happen if the volume fraction becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool enmax = en > sini[2*nch + 6]*1.001;
	if (enmax)
		cerr<<"Error in State::validState. The volume fraction of the negative electrode is "<<en<<", which is too high. The maximum is "<<sini[2*nch + 6] <<", the original volume fraction"<<endl<<flush;

	// effective surface area of the porous electrodes
	bool apmin = ap <= sini[2*nch + 7] / 5;
	if (apmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the effective surface of the positive electrode is "<<ap<<", which is too low. The minimum is "<<sini[2*nch + 7] / 5 <<", 1/5 of the original effective surface"<<endl<<flush;
		// errors will happen if the effective surface becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool apmax = ap > sini[2*nch + 7]*1.001;
	if (apmax)
		cerr<<"Error in State::validState. The effective surface of the positive electrode is "<<ap<<", which is too high. The maximum is "<<sini[2*nch + 7] <<", the original effective surface"<<endl<<flush;
	bool anmin = an <= sini[2*nch + 8] / 5;
	if (anmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the effective surface of the negative electrode is "<<an<<", which is too low. The minimum is "<<sini[2*nch + 8] / 5 <<", 1/5 of the original effective surface"<<endl<<flush;
		// errors will happen if the effective surface becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool anmax = an > sini[2*nch + 8]*1.001;
	if (anmax)
		cerr<<"Error in State::validState. The effective surface of the negative electrode is "<<an<<", which is too high. The maximum is "<<sini[2*nch + 8] <<", the original effective surface"<<endl<<flush;

	// surface area of the cracks growing at the particle surface
	bool csmin = CS <= 0;
	if(csmin)
		cerr<<"Error in State::validState. The crack surface area is "<<CS<<", which is too low. It must be strictly positive"<<endl<<flush;
		// don't allow 0 because it might give problems in some of the equations, which might grow CS proportional to the existing CS (so a value of 0 gives 0 increase)
		// a negative value might lead to a further decrease in CS, so it will keep getting more and more negative
	bool csmax = CS > pow(10,4)*sini[2*nch+9];
	if (csmax)
		cerr<<"Error in State::validState. The cell has degraded too much and the crack surface area is "<<CS<<", which is too high. The maximum is "<<pow(10,4)*sini[2*nch + 9] <<", 10,000 times the original crack surface area"<<endl<<flush;
		// normally, the initial value is 1% of the total real electrode surface area, so 10,000*initial value = 100 * total electrode surface area
		// but in theory no maximum value will give errors in the code

	// diffusion constant at reference temperature for the electrodes
	bool dpmin = Dp <= sini[2*nch + 10] / 5;
	if (dpmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the diffusion constant of the positive electrode is "<<Dp<<", which is too low. The minimum is "<<sini[2*nch + 10] / 5 <<", 1/5 of the original diffusion constant"<<endl<<flush;
		// errors will happen if the diffusion constant becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too bad diffusion -> too large surface concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool dpmax = Dp > sini[2*nch + 10]*1.001;
	if (dpmax)
		cerr<<"Error in State::validState. The diffusion constant of the positive electrode is "<<Dp<<", which is too high. The maximum is "<<sini[2*nch + 10] <<", the original diffusion constant"<<endl<<flush;
	bool dnmin = Dn <= sini[2*nch + 11] / 5;
	if (dnmin)
		cerr<<"Error in State::validState. The cell has degraded too much and the diffusion constant of the negative electrode is "<<Dn<<", which is too low. The minimum is "<<sini[2*nch + 11] / 5 <<", 1/5 of the original diffusion constant"<<endl<<flush;
		// errors will happen if the diffusion constant becomes 0 or negative. otherwise no direct problems are expected
		// on a longer term, you will get discretisation errors (see above, too bad diffusion -> too large surface concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool dnmax = Dn > sini[2*nch + 11]*1.001;
	if (dnmax)
		cerr<<"Error in State::validState. The diffusion constant of the negative electrode is "<<Dn<<", which is too high. The maximum is "<<sini[2*nch + 11] <<", the original effective diffusion constant"<<endl<<flush;

	// specific resistance of the electrodes (one value for both electrodes combined)
	bool rmin = r <= 0;
	if(rmin)
		cerr<<"Error in State::validState. The specific resistance is "<<r<<", which is too low, it must be strictly positive"<<endl<<flush;
	bool rmax = r > 1000*sini[2*nch + 12];
	if(rmax)
		cerr<<"Error in State::validState. The cell has degraded too much and the specific resistance is "<<r<<", which is too high. The maximum is "<<1000*sini[2*nch + 12] <<", 1000 times the original specific resistance"<<endl<<flush;

	// thickness of the plated litium layer
	bool delpl = delta_pl < 0;
	if (delpl)
		cerr<<"Error in State::validState. The thickness of the plated lithium is "<<delta_pl<<", which is too low. Only strictly positive values are allowed"<<endl<<flush;
		// 0 is allowed
		// a negative value might lead to a further decrease in thickness, so it will keep getting more and more negative

	// throw an error if one of the states was invalid
	if (Tmin || Tmax || del || li || tpmin || tpmax || tnmin || tnmax || epmin || epmax || enmin || enmax || apmin || apmax || anmin || anmax
			|| csmin || csmax || dpmin || dpmax || dnmin || dnmax || rmin || rmax || delpl)
		throw 15;

}

State::~State() {

}

} /* namespace std */
