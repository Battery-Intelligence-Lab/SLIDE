/*
 * Cell.cpp
 *
 * Implements the functions for the parent class of the Cells
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Cell.hpp"
#include "ReadCSVfiles.h"
#include "Interpolation.h"
#include <cassert>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <math.h>

using namespace std;

double Cell::getNominalCap(){
	/*
	 * returns nominal capacity of the cell [Ah]
	 */
	return nomCapacity;
}

double Cell::getVmax(){
	/*
	 * Returns the maximum voltage of this cell [V] above which operation is not allowed
	 */
	return Vmax;
}

double Cell::getVmin(){
	/*
	 * Returns the minimum voltage of this cell [V] below which operation is not allowed
	 */
	return Vmin;
}

double Cell::getT(){
	/*
	 * returns the uniform battery temperature in [K]
	 */
	return s.getT();
}

double Cell::getTenv(){
	/*
	 * Returns the environmental temperature [K]
	 */
	return T_env;
}

void Cell::getTemperatures(double* Tenv, double* Tref){
	/*
	 * Function to get the environmental and reference temperatures
	 *
	 * OUT
	 * Tenv 	environmental temperature [K]
	 * Tref 	reference temperature [K] at which cell parameters are measured
	 */
	*Tenv = T_env;
	*Tref = T_ref;
}

double Cell::getR(){
	/*
	 * Return the total cell DC resistance [Ohm]
	 * The total cell resistance is the sum of the resistance of the electrodes and the SEI layer
	 * 		the resistance of the electrodes increases due to loss of active material (LAM)
	 * 		the resistance of the SEI layer increases as the layer becomes thicker
	 * it is assumed the plated lithium does not add resistance
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getR starting"<<endl;

	return (Rsei*s.getDelta() + s.getR(elec_surf));

	if(verbose >= printCellFunctions)
		cout<<"Cell::getR starting"<<endl;

}

double Cell::getAnodeSurface(){
	/*
	 * Return the pure surface area of the anode
	 * This is the product of the effective surface area (an) with the electrode volume
	 */
	return s.getAn()*s.getThickn()*elec_surf;
}

double Cell::getI(){
	/*
	 * Returns the cell current [A]
	 * 	positive for discharging
	 * 	negative for charging
	 */
	return Icell;
}

void Cell::getStates(int nin, double states[]){
	/*
	 * Returns all states which describe the status of the cell in the array.
	 * see State::getStates
	 * 		zp[nch] zn[nch] T delta LLI thickp thickn ep en ap an CS Dp Dn r delta_pl
	 *
	 * IN
	 * nin 		length of the array provided
	 *
	 * OUT
	 * states 	array of length ns (defined on top of State.hpp)
	 *
	 * THROWS
	 * 100 		the array provided has the wrong length
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getStates(int, double[]) starting"<<endl;

	if (nin != ns){
		cerr<<"ERROR in Cell::getStates. The length of the array provided is "<<nin<<" instead of "<<ns<<endl<<flush;
		throw 100;
	}

	s.getStates(nin, states);

	if(verbose >= printCellFunctions)
		cout<<"Cell::getStates(int, double[]) terminating"<<endl;
}

void Cell::getStates(State& si, double* I){
	/*
	 * Copies all the states which describe the status of the cell to the State-object
	 *
	 * OUT
	 * si 	a reference to a state-boject in which the states will be written
	 * 		the & ensures this is a call by reference
	 * 		i.e. you pass a pointer to the memory location such that changes in this function to si are
	 * 			reflected in the State object which you used to call this function.
	 * 		else it would be a call by value, so changes in this function to si won't be reflected in the object you used to call this function.
	 * 		see https://stackoverflow.com/questions/21215409/does-c-pass-objects-by-value-or-reference
	 * I 	the actual cell current [A]
	 * 			positive for discharging
	 * 			negative for charging
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getStates(State, double*) starting"<<endl;

	s.getStates(si);
	*I = Icell;

	if(verbose >= printCellFunctions)
		cout<<"Cell::getStates(State, double*) starting"<<endl;
}

void Cell::getCSurf(double* cps, double* cns){
	/*
	 * Calculates the surface concentration at each particle.
	 * Uses the matrices from the Model struct with the spatial discretisation for the solid diffusion PDE.
	 *
	 * OUT
	 * cps 	surface li-concentration at the positive particle [mol m-3]
	 * cns 	surface li-concentration at the negative particle [mol m-3]
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getCSurf starting"<<endl;

	// get the transformed concentrations at the inner nodes from the state-object of the cell
	double zp[nch], zn[nch];
	s.getZ(nch, zp, zn);

	// Calculate the diffusion constant at the battery temperature using an Arrhenius relation
	double Dpt = s.getDp()*exp(Dp_T/Rg*(1/T_ref - 1/s.getT())); 	// diffusion constant of the positive particle [m s-1]
	double Dnt = s.getDn()*exp(Dn_T/Rg*(1/T_ref - 1/s.getT()));		// diffusion constant of the negative particle [m s-1]

	// Calculate the molar flux on the surfaces
	double jp = -Icell/(s.getAp()*elec_surf*s.getThickp())/(n*F);	// molar flux on the positive particle [mol m-2 s-1]
	double jn = Icell/(s.getAn()*elec_surf*s.getThickn())/(n*F);	// molar flux on the negative particle [mol m-2 s-1]

	// Calculate the surface concentration at the positive particle
	// 	cp_surf = M.Cp[0][:] * zp[:] + M.Dp*jp/Dpt
	double cp_surf = 0;
	for (int j=0;j<nch;j++)
		cp_surf += M.Cp[0][j]*zp[j];
	cp_surf += M.Dp[0]*jp/Dpt;

	// Calculate the surface concentration at the negative particle
	// 	cn_surf = M.Cn[0][:] * zn[:] + M.Dn*jn/Dnt
	double cn_surf = 0;
	for (int j=0;j<nch;j++)
		cn_surf += M.Cn[0][j]*zn[j];
	cn_surf += M.Dn[0]*jn/Dnt;

	// Make the output parameters
	*cns = cn_surf;
	*cps = cp_surf;

	if(verbose >= printCellFunctions)
		cout<<"Cell::getCSurf terminating"<<endl;
}

void Cell::getC(int nin, double cp[], double cn[]){
	/*
	 * Calculates the lithium concentration at each positive chebyshev node (0 <= x <= 1), including the centre and surface nodes.
	 * Uses the matrices from the state space matrices.
	 * See the equations in the explanatory documents.
	 *
	 * IN
	 * nin 	length of the arrays provided
	 *
	 * OUT
	 * cp 	li-concentration at each chebyshev node in the positive electrode [mol m-3], length of the array should be nch+2
	 * 			cp[0]			concentration at the surface of the sphere
	 * 			cp[1 to nch]	concentration at the inner nodes
	 * 			cp[nch + 1]		concentration at the centre of the sphere
	 * cn 	li-concentration at each chebyshev node in the negative electrode [mol m-3], length of the array should be nch+2
	 * 			cn[0]			concentration at the surface of the sphere
	 * 			cn[1 to nch]	concentration at the inner nodes
	 * 			cn[nch + 1]		concentration at the centre of the sphere
	 *
	 * THROWS
	 * 100 	the arrays provided have the wrong length
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getC starting"<<endl;

	if (nin != nch+2){												// nch is the number of inner nodes, +2 for the node at the centre and surface
		cerr<<"ERROR in Cell::getC. The length of the array provided is "<<nin<<" instead of "<<nch+2<<endl<<flush;
		throw 100;
	}

	// get the transformed concentrations at the inner nodes from the state-object of the cell
	double zp[nch], zn[nch];
	s.getZ(nch, zp, zn);

	// Calculate the diffusion constant at the battery temperature using an Arrhenius relation
	double Dpt = s.getDp()*exp(Dp_T/Rg*(1/T_ref - 1/s.getT())); 	// diffusion constant of the positive particle [m s-1]
	double Dnt = s.getDn()*exp(Dn_T/Rg*(1/T_ref - 1/s.getT()));		// diffusion constant of the negative particle [m s-1]

	// Calculate the molar flux on the surfaces
	double jp = -Icell/(s.getAp()*elec_surf*s.getThickp())/(n*F);	// molar flux on the positive particle [mol m-2 s-1]
	double jn = Icell/(s.getAn()*elec_surf*s.getThickn())/(n*F);	// molar flux on the negative particle [mol m-2 s-1]

	// Calculate concentration at the surface and inner nodes using the matrices from the spatial discretisation of the solid diffusion PDE
	// 	cp = M.Cp[:][:] * zp[:] + M.Dp*jp/Dpt
	// 	cn = M.Cn[:][:] * zn[:] + M.Dn*jn/Dnt
	double cpt, cnt;
	for(int i=0;i<nch+1;i++){										// loop to calculate at each surface + inner node
		cpt = 0.0;
		cnt = 0.0;
		for(int j=0;j<nch;j++){
			cpt += M.Cp[i][j]*zp[j];
			cnt += M.Cn[i][j]*zn[j];
		}
		cp[i] = cpt + M.Dp[i]*jp / Dpt;
		cn[i] = cnt + M.Dn[i]*jn / Dnt;
	}

	// Calculate the concentration at centre node using the boundary condition (the concentration gradient at the centre has to be 0 due to symmetry)
	// cp_centre = -1/2 (M.Cc[:]*cp +jp*Rp/Dpt)
	// cn_centre = -1/2 (M.Cc[:]*cn +jn*Rn/Dnt)
	double DM = 2.0;												// we need a constant of 2 in the equations
	cpt = 0.0;
	cnt = 0.0;
	for(int i=0;i<nch+1;i++){
		cpt += M.Cc[i]*cp[i];
		cnt += M.Cc[i]*cn[i];
	}
	cp[nch+1] = (-1.0/DM)*(cpt + jp*Rp/Dpt);
	cn[nch+1] = (-1.0/DM)*(cnt + jn*Rn/Dnt);

	if(verbose >= printCellFunctions)
		cout<<"Cell::getC terminating"<<endl;
}

bool Cell::getVoltage(bool print, double* V, double* OCVp, double* OCVn, double* etap, double* etan, double* Rdrop, double* Temp){
	/*
	 * Function to calculate the cell voltage and give detailed feedback about the cell
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if an error will be critical or not
	 *
	 * OUT
	 * V 		battery voltage [V]
	 * OCVp		cathode potential [V]
	 * OCVn		anode potential [V]
	 * etap 	overpotential at positive electrode [V], < 0 on discharge
	 * etan 	overpotential at negative electrode [V],  > 0 on discharge
	 * Rdrop 	resistive voltage drop [V], < 0 on discharge, > 0 on charge
	 * Temp 	cell temperature [K]
	 * bool 	indicates if the voltage is in the allowed range Vmin <= V <= Vmax
	 *
	 * THROWS
	 * 101		invalid surface concentration detected
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getVoltage starting"<<endl;

	// Get the surface concentrations
	double cps, cns;
	getCSurf(&cps, &cns);

	// check if the surface concentration is within the allowed range
	// 	0 < cp < Cmaxpos
	// 	0 < cn < Cmaxneg
	// don't allow 0 or Cmax because in that case, i0 will be 0, and the overpotentials will have 1/0 = inf or nan
	if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg){
		if (print){		// print error message unless you want to suppress the output
			cerr<<"ERROR in Cell::getVoltage: concentration out of bounds. the positive lithium fraction is "<<cps/Cmaxpos<<" and the negative lithium fraction is "<<cns/Cmaxneg;
			cerr<<" they should both be between 0 and 1"<<endl<<flush;
		}
		*V = nan("double");			// set the voltage to nan (Not A Number)
		throw 101;
		return false;				// the voltage is not within the limits
	}
	else{
		// Calculate the li-fraction (instead of the li-concentration)
		double zp_surf = (cps / Cmaxpos);
		double zn_surf = (cns / Cmaxneg);

		// Calculate the electrode potentials
		double OCV_p;				// cathode potential [V]
		double OCV_n;				// anode potential [V]
		double dOCV;				// entropic coefficient of the total cell voltage [V/K]
		bool bound = true;			// in linear interpolation, throw an error if you are out of the allowed range
		try{
			dOCV = linInt(print, bound, dOCV_tot_x, dOCV_tot_y, dOCV_tot_n, zp_surf);
			OCV_n = linInt(print, bound, OCV_neg_x, OCV_neg_y, OCV_neg_n, zn_surf);
			OCV_p = linInt(print, bound, OCV_pos_x, OCV_pos_y, OCV_pos_n, zp_surf);
		}
		catch(int e){
			if(print)
				cout<<"error in Cell::getVoltage when getting the electrode potentials "<<e<<". Throwing it up"<<endl<<flush;
			throw e;
		}

		// Calculate the rate constants at the cell's temperature using an Arrhenius relation
		double kpt = kp*exp(kp_T/Rg*(1/T_ref - 1/s.getT()));
		double knt = kn*exp(kn_T/Rg*(1/T_ref - 1/s.getT()));

		// Calculate the overpotential using the Bulter-Volmer equation
		// 		if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
		// 		and asinh(x) = ln(x + sqrt(1+x^2)
		double i_app = Icell/elec_surf;										// current density on the electrodes [I m-2]
		double i0p = kpt*n*F*sqrt(C_elec)*sqrt(cps)*sqrt(Cmaxpos-cps); 		// exchange current density of the positive electrode
		double i0n = knt*n*F*sqrt(C_elec)*sqrt(cns)*sqrt(Cmaxneg-cns);		// exchange current density of the negative electrode
		double xp = -0.5*i_app/(s.getAp()*s.getThickp()) / i0p;				// x for the cathode
		double xn = 0.5*i_app/(s.getAn()*s.getThickn()) / i0n;				// x for the anode
		double etapi = (2*Rg*s.getT())/(n*F) * log(xp + sqrt(1+xp*xp));		// cathode overpotential [V], < 0 on discharge
		double etani = (2*Rg*s.getT())/(n*F) * log(xn + sqrt(1+xn*xn));		// anode overpotential [V],  > 0 on discharge

		// Calculate the cell voltage
		// the cell OCV at the reference temperature is OCV_p - OCV_n
		// this OCV is adapted to the actual cell temperature using the entropic coefficient dOCV * (T - Tref)
		// then the overpotentials and the resistive voltage drop are added
		*V = (OCV_p-OCV_n + (s.getT()-T_ref)*dOCV) + (etapi - etani)  - getR()*Icell;

		// make the output variables
		*OCVp = OCV_p;
		*OCVn = OCV_n;
		*etap = etapi;
		*etan = etani;
		*Rdrop = - getR()*Icell;
		*Temp = s.getT();

		if(verbose >= printCellFunctions)
			cout<<"Cell::getVoltage terminating with V = "<<*V<<" and valid is "<< (Vmin <= *V && *V <= Vmax)<<endl;

		// Return whether this voltage is within the allowed limits or not
		return Vmin <= *V && *V <= Vmax;
	}
}

void Cell::getDaiStress(int nin, double* sigma_p, double* sigma_n, double sigma_r_p[], double sigma_r_n[],double sigma_t_p[], double sigma_t_n[], double sigma_h_p[], double sigma_h_n[]){
	/*
	 * Calculates the radial and tangential stress for each positive Chebyshev node according to the formula by
	 * Dai, Cai, White, Journal of Power sources 247, 2014
	 *
	 * It takes quite long to calculate the stress, so only call this function when needed.
	 *
	 * IN
	 * nin 			length of the arrays provided
	 *
	 * OUT
	 * sigma_p		maximum hydrostatic stress in the positive particle, can be both positive and negative [Pa]
	 * sigma_n		maximum hydrostatic stress in the negative particle, can be both positive and negative [Pa]
	 * sigma_r_p	array with the radial stress at each positive chebyshev node in the positive electrode, length nch+2, [Pa]
	 * sigma_r_n	array with the radial stress at each positive chebyshev node in the negative electrode, length nch+2, [Pa]
	 * sigma_t_p	array with the tangential stress at each positive chebyshev node in the positive electrode, length nch+2, [Pa]
	 * sigma_t_n	array with the tangential stress at each positive chebyshev node in the negative electrode, length nch+2, [Pa]
	 * sigma_h_p	array with the hydrostatic stress at each positive chebyshev node in the positive electrode, length nch+2, [Pa]
	 * sigma_h_n	array with the hydrostatic stress at each positive chebyshev node in the negative electrode, length nch+2, [Pa]
	 * 				[0]			stress at the surface of the sphere
	 * 				[1 to nch]	stress at the inner nodes
	 * 				[nch + 1]	stress at the centre of the sphere
	 *
	 * THROWS
	 * 100 	the arrays provided have the wrong length
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getDaiStress starting"<<endl;

	// check the arrays have the correct length
	if (nin != nch+2){
		cerr<<"ERROR in Cell::getDaiStress. The length of the array provided is "<<nin<<" instead of "<<nch+2<<endl<<flush;
		throw 100;
	}

	// Get the locations of the Chebyshev nodes
	double xp[nch+2]; 										// location (x-value) of the positive Chebyshev nodes
	xp[0] = 1;
	for(int i=0;i<nch;i++)
		xp[i+1] = M.xch[i];
	xp[nch+1] = 0;
	double xtot[2*nch+3];									// location (x-value) of the positive and negative Chebyshev nodes [-surface .. centre .. +surface]
	xtot[nch+1] = xp[nch+1];								// centre node
	for(int i=0;i<nch+1;i++){
		xtot[i] = -xp[i];									// negative nodes
		xtot[nch+2+i] = xp[nch-i];							// positive nodes
	}

	// get concentrations at each Chebyshev node
	// Due to symmetry, the concentration at the negative point is the same as the concentration of the positive point: c(-x) = c(x)
	double cp[nch+2], cn[nch+2]; 							// positive and negative nodes, [+surface .. inner .. centre]
	getC(nch+2,cp,cn);
	double CP[2*nch+3], CN[2*nch+3];						// concentrations at all nodes, [-surface .. inner .. centre .. inner .. +surface]
	CP[nch+1] = cp[nch+1];									// cathode centre node
	CN[nch+1] = cn[nch+1];									// anode centre node
	for(int i=0;i<nch+1;i++){
		CP[i] = cp[i];										// cathode negative points
		CN[i] = cn[i];										// anode negative points
		CP[nch+2+i] = cp[nch-i]; 							// cathode positive points
		CN[nch+2+i] = cn[nch-i]; 							// anode positive points
	}

	// The formula's to calculate the stress have integrals.
	// Integrals of Chebyshev points can be calculated using the Q-matrix in the state space struct (M)
	// The integral from the negative surface to node i is given by row i of the product Q*f
	// 		with Q the Chebyshev integration matrix
	// 			 f the value of the function you want to integrate, evaluated at every node
	// All integrals have to start from the negative surface (because you must cover the entire Chebyshev domain)
	// 	so if you need the integral of a function f from the centre (x = 0) to a positive point in the sphere (x = i)
	// 		int(f, x = 0 .. i) = int(f, x=-1 .. i) - int(f, x=-1 .. 0)
	// E.g. to get the integral of the positive li-concentration from the centre until the 4th positive Chebyshev node:
	// 		F = Q * CP 				[-surface .. centre .. +surface]
	// 		int(c(x), x = 0 .. i(4)) = int(c(x), x=-1 .. i(4)) - int(c(x), x=-1 .. 0)
	// 							  	 = F[nch + 1 + 4] 		   - F[nch+1]
	// 		(remember that the centre node is at [nch+1])

	// Calculate the matrix-vector product of the required functions as given in the paper by Dai et al. (concentration * radius^2)
	// we need to remember the transformation of variables from x (-1<x<1) to r (-R<r<R)
	// 		int( c * r^2 dr) = int(c * (x*R)^2 * d(R*x)) = int(c x^2 R^3 dx)
	// so the matrix-vector product we need is F = Q * (concentration * x^2 * R^3)
	double Fp[2*nch+3], Fn[2*nch+3]; 						// arrays with the product for the positive and negative electrode
	for(int i=0;i<2*nch+3;i++){								// loop for each row (one row = one node)
		Fp[i] = 0; 											// calculate the matrix-vector product for row i as you would do it by hand:
		Fn[i] = 0;											// 	F(i) = sum(Q(i,j)*C(j)*x(j)^2*R^3, j=0..2*nch+2)
		for(int j=0;j<2*nch+3;j++){							// loop through the columns to calculate the sum
			Fp[i] += M.Q[i][j]*(CP[j]*xtot[j]*xtot[j]);		// 		Q(i,j)*C(j)*x(j)^2
			Fn[i] += M.Q[i][j]*(CN[j]*xtot[j]*xtot[j]);
		}
		Fp[i] = Fp[i]*Rp*Rp*Rp;								// *R^3 (R is constant so it can be out of the sum)
		Fn[i] = Fn[i]*Rn*Rn*Rn;
	}

	// Calculate the integral from the centre to the positive surface, which is a constant present in all equations
	double ap = Fp[2*nch+2]-Fp[nch+1];						// int( cp*r^2, r=0..Rp )
	double an = Fn[2*nch+2]-Fn[nch+1];						// int( cn*r^2, r=0..Rn )

	// Calculate the equations for all nodes
	double rp;												// radius of positive node i in the positive particle
	double rn;												// radius of positive node i in the negative particle
	double bp;												// integral from the centre to positive node i, int(cp*zp^2, zp=0..rp(i) )
	double bn;												// integral from the centre to positive node i, int(cn*zn^2, zn=0..rn(i) )
	double srp[nch+2];										// radial stress in the positive particle at the positive nodes [centre .. +surface]
	double srn[nch+2];										// radial stress in the negative particle at the positive nodes [centre .. +surface]
	double stp[nch+2];										// tangential stress in the positive particle at the positive nodes [centre .. +surface]
	double stn[nch+2];										// tangential stress in the negative particle at the positive nodes [centre .. +surface]
	for(int i=0;i<nch+2;i++){								// loop for the positive nodes
		rp = Rp*xtot[nch+1+i]; 								// r(i) = R * x(i)
		rn = Rn*xtot[nch+1+i];
		bp = Fp[nch+1+i]-Fp[nch+1];							// int( cp*z^2, z=0..r(i) ) = F[nch+1+i] - F[nch+1]
		bn = Fn[nch+1+i]-Fn[nch+1];

		// Implement the equations from Dai et al.
		if (i == 0){ 										// centre node -> special formula (31 & 33) in Dai, Cai, White
			srp[i] = 2*omegap*Ep/(9*(1-nup))*( 3/pow(Rp,3)*ap - CP[nch+1]);
			srn[i] = 2*omegan*En/(9*(1-nun))*( 3/pow(Rn,3)*an - CN[nch+1]);

			stp[i] = 2*omegap*Ep/(9*(1-nup))*( 3/pow(Rp,3)*ap - CP[nch+1]);
			stn[i] = 2*omegan*En/(9*(1-nun))*( 3/pow(Rn,3)*an - CN[nch+1]);
		}
		else { 												// other nodes -> equation 13 in Dai, Cai, White
			srp[i] = 2*omegap*Ep/(3*(1-nup))*( 1/pow(Rp,3)*ap - 1/pow(rp,3)*bp );//ap = int (c x^2, x=0..R), bp = int (c x^2 , x=0..r)
			srn[i] = 2*omegan*En/(3*(1-nun))*( 1/pow(Rn,3)*an - 1/pow(rn,3)*bn );

			stp[i] = omegap*Ep/(3*(1-nup))*( 2/pow(Rp,3)*ap + 1/pow(rp,3)*bp - cp[i] );
			stn[i] = omegan*En/(3*(1-nun))*( 2/pow(Rn,3)*an + 1/pow(rn,3)*bn - cn[i] );
		}
	}

	// Flip all arrays to get the opposite order (now it is [centre .. +surface] and we want [+surface .. centre]
	// and store in the output arrays
	for(int i=0;i<nch+2;i++){								// loop for the positive nodes
		sigma_r_p[i] = srp[nch+2-1-i];
		sigma_r_n[i] = srn[nch+2-1-i];
		sigma_t_p[i] = stp[nch+2-1-i];
		sigma_t_n[i] = stn[nch+2-1-i];
	}

	// Make the hydrostatic stress sh = (sr + 2sp)/3
	int sp = 0;												// node with the maximum hydrostatic stress in the positive particle
	int sn = 0;												// node with the maximum hydrostatic stress in the negative
	for(int i=0;i<nch+2;i++){								// loop for all nodes
		sigma_h_p[i] = (sigma_r_p[i] + 2*sigma_t_p[i]) / 3; // calculate hydrostatic stress
		sigma_h_n[i] = (sigma_r_n[i] + 2*sigma_t_n[i]) / 3;

		// find the maximum (in absolute value) of the stresses
		if (abs(sigma_h_p[i]) > abs(sigma_h_p[sp]))
			sp = i;
		if (abs(sigma_h_n[i]) > abs(sigma_h_n[sp]))
			sn = i;
	}
	*sigma_p = sigma_h_p[sp];
	*sigma_n = sigma_h_n[sn];

	if(verbose >= printCellFunctions)
		cout<<"Cell::getDaiStress terminating"<<endl;
}

void Cell::updateDaiStress(){
	/*
	 * Function which will update the values stored in the stress variables relating with Dai's stress model
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::updateDaiStress starting"<<endl;

	// Make variables to store the stress
	double sigma_r_p[nch+2],sigma_r_n[nch+2], sigma_t_p[nch+2],sigma_t_n[nch+2], sigma_h_p[nch+2],sigma_h_n[nch+2];
	double sigma_p, sigma_n;

	// Get the stress
	try{
		getDaiStress(nch+2, &sigma_p, &sigma_n, sigma_r_p, sigma_r_n, sigma_t_p, sigma_t_n, sigma_h_p, sigma_h_n);
	}
	catch(int e){
		cout<<"Error in Cell::getDaiStress when calling updateDaiStress: "<<e<<". throwing it on"<<endl<<flush;
		throw e;
	}

	// Update the stored values
	s_dai_p = sigma_p;
	s_dai_n = sigma_n;

	// indicate that the values in the class variables are updated
	s_dai_update = true;

	if(verbose >= printCellFunctions)
		cout<<"Cell::updateDaiStress terminating"<<endl;
}

void Cell::getLaresgoitiStress(bool print, double* sigma_n){
	/*
	 * Calculate the stress according to Laresgoiti's stress model
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 *
	 * OUT
	 * sigma_n	stress in the negative particle [MPa]
	 *
	 * THROWS
	 * 101 		the surface concentration is out of bounds
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::getLaresgoitiStress starting"<<endl;


	// Arrays with the stress from Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015 (figure 5)
	double xx[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};				// li-fraction in the graphite
	double yy[11] = {0.0, 5.5, 8.5, 9.5, 10.0, 10.0, 10.5, 13.0, 16.5, 21.0, 23.5};			// stress [MPa]

	// Get the surface concentration
	double cps, cns, zn_surf;
	getCSurf(&cps, &cns);			// get the surface lithium concentration [mol m-3]
	zn_surf = (cns / Cmaxneg);		// lithium fraction on negative surface [0 1]

	// check if the surface concentration is within the allowed range
	// 	0 < cp < Cmaxpos
	// 	0 < cn < Cmaxneg
	if (cps <0 || cns <0 || cps> Cmaxpos || cns > Cmaxneg){
		if(print){
			cerr<<"ERROR in Cell::getLaresgoitiStress: concentration out of bounds. the positive lithium fraction is "<<cps/Cmaxpos<<" and the negative lithium fraction is "<<cns/Cmaxneg;
			cerr<<"they should both be between 0 and 1"<<endl<<flush;
		}
		throw 101;
	}

	// Interpolate linearly to get the stress
	double s;
	try{
		s = linInt(print, true, xx,yy,11, zn_surf); 			// throw an error if you are out of the allowed range
	}
	catch (int e){
		if(print)
			cout<<"Error in Cell::getLaresgoitiStress when interpolating in the arrays "<<e<<". Throwing it on"<<endl<<flush;
		throw e;
	}

	// Make the output variable
	*sigma_n = s;

	if(verbose >= printCellFunctions)
		cout<<"Cell::getLaresgoitiStress terminating"<<endl;
}

void Cell::updateLaresgoitiStress(bool print){
	/*
	 * Function which will update the values stored in the stress variables relating with Laresgoiti's stress model
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::updateLaresgoitiStress starting"<<endl;

	double s;
	getLaresgoitiStress(print, &s);

	// Update the stored value
	s_lares_n = s;
	s_lares_update = true;			// indicate that the values in the class variables are updated

	if(verbose >= printCellFunctions)
		cout<<"Cell::updateLaresgoitiStress terminating"<<endl;
}

void Cell::setTenv(double Tenv){
	/*
	 * Sets the environmental temperature
	 *
	 * IN
	 * Tenv 	environmental temperature, 273 <= T <= 333 [K]
	 *
	 * THROWS
	 * 102 		illegal value of T
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::setTenv starting"<<endl;

	// check if the environmental temperature is in the allowed limits
	bool T = (Tenv < TMIN) || (Tenv > TMAX); 			// the temperature limits are defined in State.hpp
	if(T){
		cerr<<"ERROR in Cell::setTenv, illegal value of environmental temperature "<<Tenv<<"K. The value has to be between The value has to be between "<< TMIN<< "and "<<TMAX<<endl<<flush;
		throw 102;
	}

	// update the temperature
	T_env = Tenv;

	if(verbose >= printCellFunctions)
		cout<<"Cell::setTenv terminating"<<endl;
}

void Cell::setVlimits(double VMAX, double VMIN){
	/*
	 * sets the voltage limits of this cell to the given values.
	 *
	 * IN
	 * VMAX		maximum voltage of this cell, [V]
	 * VMIN 	minimum voltage of this cell, [V]
	 *
	 * THROWS
	 * 103 		illegal value of the voltage limits
	 * 			Both values have to be between the maximum and minimum of the OCV curve of the cell
	 * 			We are not checking the full OCV curve because that would take too long.
	 * 			Instead, we check that the values are positive and that VMAX is below the maximum of the cathode OCV.
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::setVlimits starting with upper limit "<<VMAX<<" and lower limit "<<VMIN<<endl;

	bool vmax = VMAX < 0 || VMAX > OCV_pos_y[0];
	if (vmax)
		cerr<<"ERROR in Cell::setVlimits. The value of the maximum voltage is "<<VMAX<<"V but it has to be positive and lower than the maximum value of the OCV curve of the cathode, which is "<<OCV_pos_y[OCV_pos_n-1]<<endl<<flush;
	bool vmin = VMIN < 0;
	if(vmin)
		cerr<<"ERROR in Cell::setVlimits. The value of the minimum voltage is "<<VMIN<<"V but it has to be positive"<<endl<<flush;
	if (vmax || vmin)
		throw 103;

	Vmax = VMAX;
	Vmin = VMIN;

	if(verbose >= printCellFunctions)
		cout<<"Cell::setVlimits terminating with upper limit "<<Vmax<<" and lower limit "<<Vmin<<endl;
}

void Cell::setT(double Ti ){
	/*
	 * Set the temperature of the battery
	 *
	 * IN
	 * Ti		uniform cell temperature, 273 <= T <= 333 [K]
	 *
	 * THROWS
	 * 102 		illegal value of T
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::setT starting"<<endl;

	// check that the temperature is in the allowed limits
	bool T = (Ti < TMIN) || (Ti > TMAX); 				// the temperature limits are defined in State.hpp
	if(T){
		cerr<<"ERROR in Cell::setT, illegal value of cell temperature "<<T<<"K. The value has to be between "<< TMIN<< "and "<<TMAX<<endl<<flush;
		throw 102;
	}

	s.setT(Ti);

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	s_dai_update = false;
	s_lares_update = false;

	if(verbose >= printCellFunctions)
		cout<<"Cell::setT terminating"<<endl;
}

void Cell::setStates(int nin, double states[]){
	/*
	 * sets all battery state variables
	 *
	 * IN
	 * nin 		length of the array
	 * states	array with the battery states
	 *
	 * THROWS
	 * 100 		the array provided has the wrong length
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::setStates(int, double[]) starting"<<endl;

	if (nin != ns){
		cerr<<"ERROR in Cell::setStates. The length of the array provided is "<<nin<<" instead of "<<ns<<endl<<flush;
		throw 100;
	}

	s.getStates(n, states);

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	s_dai_update = false;
	s_lares_update = false;

	if(verbose >= printCellFunctions)
		cout<<"Cell::setStates(int, double[]) terminating"<<endl;
}
void Cell::setStates(State si, double I){
	/*
	 * Set all battery state variables to the ones in the State-object provided.
	 * Set the cell current to the value provided (without ramping the current since we are changing the state any way)
	 *
	 * IN
	 * si 	new state of the cell
	 * I 	new current of the cell [A], positive for discharging, negative for charging
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::setStates(State, double) starting"<<endl;

	s.setStates(si);
	Icell = I;

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	s_dai_update = false;
	s_lares_update = false;

	if(verbose >= printCellFunctions)
		cout<<"Cell::setStates(State, double) terminating"<<endl;
}

void Cell::setC(double cp0, double cn0){
	/*
	 * Function to set the value of the transformed concentrations to the values
	 * corresponding with a uniform (user-specified) lithium concentration.
	 * The cell's current is set to 0 because the concentration is fully uniform (which means the gradient at the surface is 0, so the current must be 0)
	 *
	 * IN
	 * cp0	lithium fraction in the positive electrode 0 <= cp0 <= 1
	 * cn0	lithium fraction in the negative electrode 0 <= cn0 <= 1
	 *
	 * THROWS
	 * 104 	illegal input lithium fractions
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::setC starting"<<endl;

	bool pp = (cp0 < 0) || (cp0 > 1);
	if (pp)
		cerr<<"ERROR in Cell::setC, illegal input li-fraction for the positive electrode : "<<cp0<<". The value has to be between 0 and 1"<<endl<<flush;
	bool nn = (cn0 < 0) || (cn0 > 1);
	if (nn)
		cerr<<"ERROR in Cell::setC, illegal input li-fraction for the negative electrode : "<<cn0<<". The value has to be between 0 and 1"<<endl<<flush;
	if (pp || nn)
		throw 104;

	// Calculate the corresponding li-concentrations in [mol m-3]
	double cp = cp0*Cmaxpos;
	double cn = cn0*Cmaxneg;

	// Do the first transformation, to u(i) = radius(i) * concentration = x(i) * R * concentration(i)
	double uneg[nch], upos[nch];
	for(int i=0;i<nch;i++){
		uneg[i] = cn * M.xch[i] * Rn;
		upos[i] = cp * M.xch[i] * Rp;
	}

	// The second transformation is to the eigenspace: z = V * u with V the inverse of the matrix with the eigenvectors.
	// As explained, we know that there is one eigenvector corresponding to a uniform concentration
	// So we need to calculate only this one nonzero value for the (twice) transformed concentration
	// The location of the uniform eigenvector (which has a 0 eigenvalue) is written in M.Input[3]
	int ind = M.Input[3];
	double znu = 0;						// (twice) transformed negative uniform concentration
	double zpu = 0;						// (twice) transformed positive uniform concentration
	for (int i=0;i<nch;i++){			// loop to calculate the row of V * u corresponding to the uniform concentration
		znu += M.Vn[ind][i]*uneg[i];
		zpu += M.Vp[ind][i]*upos[i];
	}

	// Make the full arrays for the (twice) transformed concentration
	double zp[nch], zn[nch];
	for (int i=0;i<nch;i++){			// set all values to 0
		zp[i] = 0;
		zn[i] = 0;
	}
	zp[ind] = zpu;						// set the non-zero value
	zn[ind] = znu;
	s.setZ(nch, zp, zn);

	// Set the cell current to 0 to reflect the boundary condition for a fully uniform concentration
	Icell = 0;

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	s_dai_update = false;
	s_lares_update = false;

	if(verbose >= printCellFunctions)
		cout<<"Cell::setC terminating"<<endl;
}

void Cell::setI(bool print, bool check, double I){
	/*
	 * Function to set the cell current to the specified value.
	 * The current is slowly ramped from the actual cell current to the specified value.
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * check 	if true, a check is done to see if the battery state at the end is valid or not (an error is thrown if not)
	 * 			if false, the current is just set without checking if the state is valid or not
	 * I 		value to which the current should be set [A]
	 * 			> 0 for discharge
	 * 			< 0 for charge
	 *
	 * THROWS
	 * 105 		check == true and the specified current could not be set without violating the cell's limits.
	 * 			the original battery state and current are restored
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::setI starting with current "<<I<<endl;

	// check if the specified value is different from the actual cell current
	if(abs(Icell - I)<pow(10,-10))
		return;									// the values are the same -> we don't need to do anything

	// Store the old current and state to restore it if needed
	double Iold;
	State sold;
	getStates(sold, &Iold);

	// settings
	bool blockDegradation = true; 				// don't account for degradation while the current is ramping
	bool reached = false; 						// boolean indicating if we have reached the set current
	int sign;									// sign whether we should increase or decrease the current
	if (I > Icell)
		sign = 1;
	else
		sign = -1;

	// loop to ramp the current
	while(!reached){
		// increase the current
		Icell += sign*dIcell;

		// check if you have reached the specified current
		if(sign*Icell > sign * I){				// increase I: Icell > I, decrease I: Icell < I
			Icell = I;
			reached = true;
		}

		// take one small time step
		ETI(print, dt_I, blockDegradation);
	}

	// Check the cell's conditions are still valid if we want to check the final state
	double v, ocvp, ocvn, etap, etan, rdrop, tem;
	bool valid;
	if(check){
		// check the state
		try{
			s.validState();						// throws an error if the state is illegal
		}
		catch(int e){
			if(print)
				cout<<"Cell::setI illegal state after setting the current to "<<Icell<<", error: "<<e<<". Throwing an error"<<endl<<flush;
			setStates(sold, Iold);				// restore the original battery state and current
			throw e;
		}

		// check the voltage
		try{
			valid = getVoltage(print, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem); // throws an error if the surface concentration is out of bounds
		}
		catch(int e){
			valid = false;						// in that case, the voltage is illegal
		}
		if(!valid){
			if(print)
				cerr<<"Cell::setI Illegal voltage after trying to set the current to "<<Icell<<", the voltage is: "<<v<<"V. Throwing an error"<<endl<<flush;
			setStates(sold, Iold);				// restore the original battery state and current
			throw 105;
		}
	}

	// the stress values stored in the class variables for stress are no longer valid because the cell current has changed
	s_dai_update = false;
	s_lares_update = false;

	if(verbose >= printCellFunctions){
		if(check)
			cout<<"Cell::setI terminating with current "<<I<<" and voltage "<<v<<endl;
		else
			cout<<"Cell::setI terminating with current "<<I<<" without checking the voltage"<<endl;
	}
}

void Cell::SEI(double OCVnt, double etan, double* isei, double* den){
	/*
	 * Function to calculate the degradation effects of growth of the SEI layer
	 *
	 * IN
	 * OCVnt 	the OCV of the negative electrode at the battery temperature [V]
	 * etan 	the overpotential at the negative electrode [V]
	 *
	 * OUT
	 * isei 	current density for the SEI side-reaction [A m-2]
	 * den 		decrease in the active volume fraction as a result of SEI growth [sec-1]
	 *
	 * THROWS
	 * 106 		illegal value in id or por
	 * 107		too many degradation models
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::SEI starting"<<endl;

	// variables
	double kseit, Dseit;				// temperature-dependent SEI parameters, using Arrhenius' relation
	double isei1, isei2, isei3;			// temporary parameters to calculate the SEI growth
	double is = 0;						// SEI side reaction current density of all models combined

	// check that the number of degradation models is not too long
	if(deg_id.SEI_n > deg_id.len){
		cerr<<"ERROR in Cell::SEI the user wants to use more than "<<deg_id.len<<" degradation models. throwing an error"<<endl<<flush;
		cout<<"The maximum length of the array with the degradation identifiers is "<<deg_id.len<<". If you want to use more models,"
				" you have to increase the value of 'len' in the struct DEG_ID, defined in Cell.hpp"<<endl<<flush;
		throw 107;
	}

	// Loop for each model to use
	for(int i=0;i<deg_id.SEI_n;i++){

		// Use a switch to calculate the magnitude of the SEI growth according to this degradation model
		switch (deg_id.SEI_id[i]) {
		case 0 :						// no SEI growth
			is += 0;
			break;
		case 1 :						// Kinetic model according to Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004
			kseit = seiparam.sei1k*exp(seiparam.sei1k_T/Rg*(1/T_ref - 1/s.getT())); 		// Arrhenius relation for the rate parameter at the cell temperature
			is += nsei*F*kseit*exp(-nsei*F/(Rg*s.getT()) * alphasei * (OCVnt + etan - OCVsei + Rsei*s.getDelta()*Icell)); // Add the effect of this model
				// eta_sei = OCVneg + etaneg - OCVsei + Rsei*I
				// isei = nFk exp(-nF/RT alpha eta_sei)
				// on charge, I < 0 and etan < 0.
				// so higher charging current -> more negative term in exponential -> larger isei
			break;
		case 2 :						// kinetics and diffusion according to Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
			kseit = seiparam.sei2k*exp(seiparam.sei2k_T/Rg*(1/T_ref - 1/s.getT())); 		// Arrhenius relation for the rate parameter at the cell temperature
			Dseit = seiparam.sei2D*exp(seiparam.sei2D_T/Rg*(1/T_ref - 1/s.getT())); 		// Arrhenius relation for the diffusion constant at the cell temperature

			// derivation of the formula:
				// start equations (use the symbols from Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
				// but with opposite sign for j
					// j = nFk c exp(..)
					// j/nF = - D/delta (c - c0)
				// j/nF = D/delta (- j/(nFk exp(..)) + c0)
				// j * ( 1/(nFk exp(..)) + delta/DnF ) = c0
				// j = c0 / ( 1/(nFk exp(..)) + delta/DnF )
			isei2 = nsei*F*kseit* exp(-nsei*F/(Rg*s.getT()) * alphasei * (OCVnt + etan - OCVsei + Rsei*s.getDelta()*Icell));
			isei3 = s.getDelta()/(nsei * F * Dseit);
			is += c_elec0 / (1/isei2 + isei3);												// Add the effects of this model
			break;
		case 3 :						// model from Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
			kseit = seiparam.sei3k*exp(seiparam.sei3k_T/Rg*(1/T_ref - 1/s.getT())); 		// Arrhenius relation for the rate parameter at the cell temperature
			Dseit = seiparam.sei3D*exp(seiparam.sei3D_T/Rg*(1/T_ref - 1/s.getT())); 		// Arrhenius relation for the diffusion constant at the cell temperature

			// Use equation [22] from the paper
			isei1 = 0.134461 * exp(-nsei*F*(etan+Rsei*s.getDelta()*Icell)/(Rg*s.getT())); 	// the parameter a_L_K is set to 0.134461 but this constant can be lumped into the rate- and diffusion constants
			isei2 = nsei*F*kseit*exp(-nsei*F/(Rg*s.getT()) * alphasei * (OCVnt -OCVsei));
			isei3 = s.getDelta()/(nsei * F * Dseit);
			is += isei1 / (1 / isei2 + isei3);												// Add the effects of this model
			break;
		default :						// unknown degradation model
			cerr<<"ERROR in Cell::SEI, unknown SEI degradation model with identifier "<<deg_id.SEI_id[i]<<". Only values 0 to 3 are allowed. Throw an error"<<endl<<flush;
			throw 106;
			break;
		}
	} // end loop for all the models you want to use

	// Make the output for the SEI side reaction current density
	*isei = is;

	// Calculate how much we decrease the volume fraction due to SEI growth
	if(deg_id.SEI_porosity == 0) 		// don't decrease volume fraction
		*den = 0;
	else if (deg_id.SEI_porosity == 1){	// decrease volume fraction according to Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
		double jn = Icell/elec_surf/(s.getAn()*n*F*s.getThickn()); // molar flux on the negative particle
		*den = - seiparam.sei_porosity*(jn * Vmain + *isei * Vsei);
			// note: they use J = volumetric current [A/m3] -> they multiply with 'an' but we already have density [A/m2]
			// - because they use the porosity while we use the volume fraction
	}
	else{						// unknown degradation model
		cerr<<"ERROR in Cell::SEI, unknown value for decreasing the volume fraction "<<deg_id.SEI_porosity<<". Only values 0 or 1 are allowed. Throw an error"<<endl<<flush;
		throw 106;
	}

	if(verbose >= printCellFunctions)
		cout<<"Cell::SEI terminating"<<endl;
}

void Cell::CS(double OCVnt, double etan, double* isei_multiplyer, double* dCS, double* dDn){
	/*
	 * function to calculate the degradation effect of surface cracking due to fatigue
	 *
	 * IN
	 * OCVnt 			the OCV of the negative electrode at the battery temperature [V]
	 * etan 			the overpotential at the negative electrode [V]
	 *
	 * OUT
	 * isei_multiplyer 	Extra SEI side reaction due to crack growth as a fraction of the original SEI side reaction current [-].
	 * 						i.e. total SEI growth = (1+isei_multiplyer)*isei
	 * dCS 				increase in surface area due to fatigue in this time step [m2 sec-1]
	 * dDn 				decrease in the diffusion constant of the graphite due to surface cracks [m s-1 s-1]
	 *
	 * THROWS
	 * 106 				illegal value in id or d
	 * 107				too many degradation models
	 * 108 				the stress values are not up to date
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::CS starting"<<endl;

	// output parameters
	double ism = 0; 						// isei multiplier from all models to be considered
	double dcs = 0;							// increase in surface area from all models to be considered

	if(deg_id.CS_n > deg_id.len){
		cerr<<"ERROR in Cell::CS the user wants to use more than "<<deg_id.len<<" degradation models. throwing an error"<<endl<<flush;
		cout<<"The maximum length of the array with the degradation identifiers is "<<deg_id.len<<". If you want to use more models,"
				" you have to increase the value of 'len' in the struct DEG_ID, defined in Cell.hpp"<<endl<<flush;
		throw 107;
	}

	double ASn = getAnodeSurface(); 		// active surface area of the anode [m2]
											// this active area is used to translate absolute values (such as currents) to relative values (such as current densities)

	// Loop for each model we want to use
	for(int i=0;i<deg_id.CS_n;i++){

		// a switch to calculate the effect according to model i
		switch (deg_id.CS_id[i]){
		case 0 : 							// no surface cracks
			ism += 0;
			dcs += 0;
			break;
		case 1 : 							// Laresgoiti's stress and crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
											// this model calculates crack growth due to temporal variations in the li-concentration
			// check the calculated stress values are up to date
			if(!s_lares_update){
				cerr<<"ERROR in Cell::CS. The stress values for Laresgoiti's stress model are not updated. Throwing an error"<<endl<<flush;
																			// if you see this error, you have to call Cell::updateLaresgoitiStress(), which calculates the stress and stores the values, before you call this function
				throw 108;
			}

			// Implement the equation from the paper
			// capacity loss is m-power of the so-called stress amplitude (sigma_max - sigma_min)/2
			// sigma_max and sigma_min are the max and min stresses 'of the cyclic signal' i.e. within one charge/discharge
			// assume m = 1, then (max - min) = (max - t1) + (t1-t2) + (t2-t3) + ... + (tn - min)
			// so stress amplitude can be substituted by (the stress in the previous time step) - (the stress in this time step)
			dcs += csparam.CS1alpha*abs(s_lares_n - s_lares_n_prev)/2; 		// equations (22)+ (27) from the paper
			ism += s.getCS()/ASn;											// increase SEI growth proportionally the crack surface
																			// current density on particle = I /(elec_surf * thick * a)
																			// 		isei also acts on this scale since it is an extra boundary condition ( itot = (jn + isei) =  (surface gradient)/nF )
																			// 		crack growth -> increase isei_on_particle -> (jn + isei*(initial+crack_surface)/initial)*nF
																			// 			such that if the crack surface area is the same as the initial electrode surface area, we double isei
			break;
		case 2 :							// Laresgoiti's crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
											// but with stress model from Dai, Cai, White, Journal of Power sources 247, 2014
												// instead of Laresgoiti's stress from figure 5
											// this model calculates crack growth due to spatial (Dai) and temporal (Laresgoiti) variations in the li-concentration
			// ensure the stress values are up to date
			if(!s_dai_update){
				cerr<<"ERROR in Cell::CS. The stress values for Dai's stress model are not updated. Throwing an error"<<endl<<flush;
																			// if you see this error, you have to call Cell::updateDaiStress(), which calculates the stress and stores the values before you call this function
				throw 108;
			}

			// Add the effects of this model
			dcs += csparam.CS2alpha*abs(s_dai_n - s_dai_n_prev)/2; 			// equations (22)+ (27) from the paper but with Dai's stress
			ism += s.getCS()/ASn;											// increase SEI growth proportionally the crack surface
			break;
		case 3 :							// model by Deshpande & Bernardi,Journal of the Electrochemical Society 164 (2), 2017
											// this model is adapted to calculate crack growth due to spatial variations in the li-concentration
			// get concentrations
			double cp[nch+2], cn[nch+2];
			getC(nch+2,cp,cn);

			// Add the effects of this model
			dcs += csparam.CS3alpha * pow((cn[0]-cn[nch+1])/Cmaxneg,2); 	// equations (8) + (21)
																			// Note that eqn (8) refers to the change with respect to time while here we use the spatial variation
																			// This makes the model capture more of the spatial variation in stress
																			// Laresgoiti's model already accounted for temporal variation, so simply use that if you are interested in temporal rather than spatial variation
			ism += s.getCS()/ASn; 											// increase SEI growth proportionally the crack surface
			break;
		case 4 : 							// model from Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
			// equation (1a): CS = Amax(1 - exp(-m * Ah)) with m a fitting parameter and Ah the charge throughput up to now
			// 	dCS/dAh = m Amax exp(-m Ah) = m Amax - m Amax + m Amax exp(-m Ah) = m Amax - m CS = m(Amax - CS)
			// 	dCS/dt = dCS/dAh * dAh/dt = dCS/dAh * abs(I) = m(Amax - CS)*abs(I)
			// 		where we use the absolute value of I because 'Ah' is the total charge throughput, i.e. int ( abs(I) dt )

			double Amax;													// 'maximum crack surface area', a fitting parameters
			Amax = max(csparam.CS4Amax, s.getCS()); 						// avoid negative crack growth if the crack surface becomes slightly larger than Amax
																			// this is possible due to discrete time steps: CS(t) is just smaller, but CS (t+1) = CS(t) + dCS*dt is just larger

			// Add the effects of this model
			dcs += csparam.CS4alpha*(Amax - s.getCS()) * abs(Icell); 		// see above, with m = csparam.CS4
			ism += s.getCS()/ASn;											// increase SEI growth proportionally the crack surface
			break;
		case 5 :							// model from Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
			double etasei;													// overpotential for the crack-side-reaction = overpotential for the SEI reaction
			double  kcr;													// rate constant for the side reaction

			etasei = (OCVnt + etan - OCVsei + Rsei*s.getDelta()*Icell); 	// overpotential [V], equation (6)

			// get surface concentration
			double cps, cns;
			getCSurf(&cps, &cns);											// get the surface lithium concentration

			// Calculate the rate constant, equation (11) with an Arrhenius relation for the temperature (which wasn't considered by Ekstrom)
			if (Icell > 0)
				kcr = 0;
			else if (cns / Cmaxneg < 0.3)
				kcr = 2*csparam.CS5k * exp(csparam.CS5k_T/Rg*(1/T_ref - 1/s.getT()));
			else if (cns / Cmaxneg < 0.7)
				kcr = 0;
			else
				kcr = csparam.CS5k * exp(csparam.CS5k_T/Rg*(1/T_ref - 1/s.getT()));

			// Add the effects of this model
			dcs += nsei * F * kcr * exp(-alphasei * nsei * F/(Rg * s.getT()) * etasei); // equation (9)
			ism += s.getCS()/ASn;											// increase SEI growth proportionally the crack surface
			break;
		default :							// unknown degradation model
			cerr<<"ERROR in Cell::CS, unknown crack growth model with identifier "<<deg_id.CS_id[i]<<". Only values 0 to 5 are allowed. Throw an error"<<endl<<flush;
			throw 106;
			break;
		}
	}

	// Make the output variables
	*isei_multiplyer = ism;
	*dCS = dcs;

	// Decrease the negative diffusion constant if needed
	if(deg_id.CS_diffusion == 0)			// don't decrease the negative diffusion constant
		*dDn = 0;
	else if (deg_id.CS_diffusion == 1){		// decrease it according to Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
		// equation (2) D(t) = D0 (1 - CS)^gamma
		// 	but this can become negative is CS is larger than 1, we assume there should be a term /Amax in eqn (2): D(t) = D0 (1 - (CS/Amax))^gamma, which becomes 0 if CS = Amax
		// so dD/dt = - gamma D0 (1-CS/Amax)^(gamma-1) 1/Amax dCS/dt

		double Dnmax;														// cap the decrease rate at a maximum value of 2e-7 (2e-5% = kill the battery in  about 1000h)
		double Amax;														// 'maximum crack surface area', a fitting parameters
		Amax = max(csparam.CS4Amax, s.getCS()); 							// avoid increasing diffusion coefficient if the crack surface becomes larger than Amax
																			// this is possible if the user chooses a different CS growth model, which does give larger crack surfaces (i.e. not CS4 which is Barai's crack growth model)
		Dnmax = csparam.CS_diffusion * pow(1-s.getCS()/Amax,csparam.CS_diffusion-1) / Amax *dcs;
		Dnmax = min(Dnmax, 2*pow(10,-7));
		*dDn = -Dnmax*s.getDn();
	}
	else{								// unknown degradation model
		cerr<<"ERROR in Cell::CS, unknown value for decreasing the diffusion constant "<<deg_id.CS_diffusion<<". Only values 0 or 1 are allowed. Throw an error"<<endl<<flush;
		throw 106;
	}

	if(verbose >= printCellFunctions)
		cout<<"Cell::CS terminating"<<endl;
}

void Cell::LAM(bool print, double zp_surf, double etap,
		double* dthickp, double* dthickn, double* dap, double* dan, double* dep, double* den){
	/*
	 * Function to calculate the effect of loss of active material (LAM).
	 *
	 * LAM is simulated by decreasing the amount of active material.
	 * This will increase the current density on the particle for the same overall cell current (because there is less material to 'spread' it over).
	 * This will mean that for the same overall cell current,
	 * 		there will be a larger change in lithium concentration,
	 * 			so there is a larger change in open circuit voltage,
	 * 			so a smaller capacity before a voltage limit is reached
	 * 			so the capacity decreases,
	 * 		and there will be a larger resistive voltage drop, so the 'effective' resistance increases
	 *
	 * There are a couple of variables describing the amount of active material:
	 * 	elec_surf 	the (geometric) surface of the electrode, i.e. the product of the height and length of the electrode [m2]
	 * 	thick		the thickness of the electrode material [m]
	 * 	R 			radius of the particle of the single particle model [m]
	 * 	e			the volume fraction of active material	[-]
	 * 	a 			the effective surface, i.e. the surface per unit of electrode volume [m2 m-3]
	 * 				a = 3*e/R
	 *
	 * 	The current density is calculated as:
	 * 	i = I / (a * thick * elec_surf) = I / (3* e / R * thick * elec_surf)
	 *
	 * 	A decrease in any of these geometric parameters will have exactly the same effect on the battery:
	 * 	you can double i by halving a, or by halving thick, or by halving elec_surf.
	 * 	and you can halve a by halving e or doubling R.
	 * 	So it doesn't really matter which of the geometric parameters you decrease to account for LAM.
	 * 	However, this model doesn't allow to change elec_surf and R because these values are needed on multiple locations in the code.
	 * 	E.g. the value of R is needed to calculate values of the matrices used for the diffusion state space model (because the spatial discretisation depends on R).
	 * 	So if you were to change R, you have to re-calculate the matrices, which would needlessly complicate the model.
	 * 	Therefore, degradation can only decrease the values of 'thick', 'a' and 'e'
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * zp_surf 		li-fraction at the surface of the positive particle [-]
	 * etap 		overpotential at the positive electrode [V]
	 *
	 * OUT
	 * dthickp 		change in electrode thickness of the positive electrode [m s-1]
	 * dthickn		change in electrode thickness of the negative electrode [m s-1]
	 * dap			change in effective electrode surface of the positive electrode [m2 m-3 s-1]
	 * dan			change in effective electrode surface of the negative electrode [m2 m-3 s-1]
	 * dep			change in volume fraction of active material in the positive electrode [s-1]
	 * den			change in volume fraction of active material in the negative electrode [s-1]
	 *
	 * THROWS
	 * 106 			illegal value in id
	 * 107			too many degradation models
	 * 108 			the stress values are not updated
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::LAM starting"<<endl;

	// output parameters
	double dthickpp = 0;
	double dthicknn = 0;
	double dapp = 0;
	double dann = 0;
	double depp = 0;
	double denn = 0;

	if(deg_id.LAM_n > deg_id.len){
		cerr<<"ERROR in Cell::LAM the user wants to use more than "<<deg_id.len<<" degradation models. throwing an error"<<endl<<flush;
		cout<<"The maximum length of the array with the degradation identifiers is "<<deg_id.len<<". If you want to use more models,"
				" you have to increase the value of 'len' in the struct DEG_ID, defined in Cell.hpp"<<endl<<flush;
		throw 107;
	}

	// loop for each model to use
	for(int i=0;i<deg_id.LAM_n;i++){

		// calculate the effect of this model
		switch (deg_id.LAM_id[i]){
		case 0 : 							// no LAM
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
			depp += 0;
			denn += 0;
			break;
		case 1 : 							// Stress model from Dai, Cai, White, Journal of Power sources 247, 2014
											// LAM equation similar to CS equation from Laresgoiti
											// (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
			// ensure the stress values are up to date
			if(!s_dai_update){
				cerr<<"ERROR in Cell::LAM. The stress values for Dai's stress model are not updated. Throwing an error"<<endl<<flush;
																			// if you see this error, you have to call Cell::updateDaiStress(), which calculates the stress and stores the values before calling this function
				throw 108;
			}

			// Laresgoiti's equation to link stress to LAM
			dthickpp += - lamparam.lam1p * abs(s_dai_p - s_dai_p_prev)/2; 	// with Dai's stress model (values stored in s_dai_p)
			dthicknn += - lamparam.lam1n * abs(s_dai_n - s_dai_n_prev)/2;
			// assume the other effects are 0
			dapp += 0;
			dann += 0;
			depp += 0;
			denn += 0;
			break;
		case 2:								// Model by Delacourt & Safari, Journal of the Electrochemical Society 159 (8), 2012
			// Get the molar flux on each particle
			double i_app, jp, jn;
			i_app = Icell/elec_surf;										// current density on the electrode [A m-2]
			jp = -i_app/(s.getAp()*n*F*s.getThickp()); 						// molar flux on the positive particle [mol m-2 s-1]
			jn = i_app/(s.getAn()*n*F*s.getThickn()); 						// molar flux on the negative particle [mol m-2 s-1]

			// Use Arrhenius relations to update the fitting parameters for the cell temperature
			double ap, bp, an, bn;
			ap = lamparam.lam2ap * exp(lamparam.lam2t / Rg * (1/T_ref - 1/s.getT()));
			an = lamparam.lam2an * exp(lamparam.lam2t / Rg * (1/T_ref - 1/s.getT()));
			bp = lamparam.lam2bp * exp(lamparam.lam2t / Rg * (1/T_ref - 1/s.getT()));
			bn = lamparam.lam2bn * exp(lamparam.lam2t / Rg * (1/T_ref - 1/s.getT()));

			// Add the effects of this model
			depp += ap * abs(jp) + bp *sqrt(abs(jp));						// equation (5) from the paper
			denn += an * abs(jn) + bn *sqrt(abs(jn));
			// assume the other effects are 0
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
			break;
		case 3 :								// Model by Kindermann, Keil, Frank, Jossen, Journal of the Electrochemical Society 164 (12), 2017
			double OCVpt; 													// cathode potential
			double etap_LAM;												// overpotential for the NMC dissolution reaction
			double kt;														// temperature dependent rate constant
			double idiss;													// current density of the NMC dissolution reaction

			try{
				OCVpt = linInt(print, true, OCV_pos_x, OCV_pos_y, OCV_pos_n, zp_surf); 	// get OCV of positive electrode, throw error if out of bounds
				// this should be updated for the cell's temperature using the entropic coefficient of the cathode
				// but I couldn't find any data on this, so I have ignored the effect
			}
			catch(int e){
				cout<<"Error in Cell::LAM when calculating the cathode potential for LAM: "<<e<<". Throwing it on"<<endl<<flush;
				throw e;
			}

			etap_LAM = OCVpt + etap - OCVnmc; 								// equation (9) from the paper
			kt = lamparam.lam3k * exp(lamparam.lam3k_T/Rg * (1/T_ref - 1/s.getT())); // Arrhenius law
			idiss = - kt * exp(n*F/Rg/s.getT()*etap_LAM) / (n*F); 			// equation (8) from the paper
			idiss = max(idiss, -5*pow(10,-6)); 								// cap the effect at 5e-6 to avoid a very fast drop of capacity (which could  cause an error)
																			// a value of 5e-6 gives dap = -3.5. The initial value is about 17000, so the cell is dead in 5,000 seconds
																			// so this cap is quite high

			// Add the effects of this model
			depp += idiss;
			// assume the other effects are 0
			denn += 0;
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
			break;
		case 4 :								// Model by Narayanrao, Joglekar, Inguva, Journal of the Electrochemical Society 160 (1), 2012
			// Add the effects of this model
			dapp += - lamparam.lam4p*s.getAp(); 						// equation (7) from the paper
			dann += - lamparam.lam4n*s.getAn();
			// assume the other effects are 0
			dthickpp += 0;
			dthicknn += 0;
			depp += 0;
			denn += 0;
			break;
		default:								// unknown degradation model
			cerr<<"ERROR in Cell::LAM, unknown LAM degradation model with identifier "<<deg_id.LAM_id[i]<<". Only values 0 to 4 are allowed. Throw an error"<<endl<<flush;
			throw 106;
			break;
		}
	}

	// Make the output variables
	*dthickp = dthickpp;
	*dthickn = dthicknn;
	*dap = dapp;
	*dan = dann;
	*dep = depp;
	*den = denn;

	if(verbose >= printCellFunctions)
		cout<<"Cell::LAM terminating"<<endl;
}

void Cell::LiPlating(double OCVnt, double etan, double* ipl){
	/*
	 * Function to simulate the effect of lithium plating
	 *
	 * IN
	 * OCVnt 	the OCV of the negative electrode at the current battery temperature [V]
	 * etan 	the overpotential at the negative electrode [V]
	 *
	 * OUT
	 * ipl 		current density for the plating side-reaction [A m-2]
	 *
	 * THROWS
	 * 106 		illegal value in id
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::LiPlating starting"<<endl;

	// Arrhenius relation for temperature-dependent plating parameters
	double kplt = plparam.pl1k*exp(plparam.pl1k_T/Rg*(1/T_ref - 1/s.getT()));	// Rate constant

	if(deg_id.pl_id == 0)								// no plating
		*ipl = 0;
	else if(deg_id.pl_id == 1)							// Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
		*ipl = npl*F*kplt*exp(-n*F/(Rg*s.getT()) * alphapl * (OCVnt + etan - OCVpl + Rsei*s.getDelta()*Icell));
	else{
		cerr<<"ERROR in Cell::LiPlating, illegal degradation model identifier "<<deg_id.pl_id<<", only values 0 and 1 are allowed. Throwing an error"<<endl<<flush;
		throw 106;
	}

	if(verbose >= printCellFunctions)
		cout<<"Cell::LiPlating starting"<<endl;
}

// state space model
void Cell::dState(bool print, bool blockDegradation, int electr, int nin, double dstates[]){
	/*
	 * function calculating the time derivatives of the battery states.
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * blockDegradation 	if true, battery degradation is ignored (i.e. the time derivatives of those states are 0)
	 * electr 	integer describing which electrodes should be considered in the case you want to do half-cell cycling.
	 * 			half-cell cycling is only allowed if blockDegradation is true (i.e. you ignore degradation)
	 * 			this is because some of the degradation mechanisms could give NaN or Inf due to a divide by 0
	 * 				1 					only the positive electrode is considered [half-cell cycling with only the positive electrode]
	 * 				2 					only the negative electrode is considered [half-cell cycling with only the negative electrode]
	 * 				any other value		both electrodes are considered [normal operation]
	 * nin	length of the output array
	 *
	 * OUT
	 * dstates	change in the states
	 * 		dzp			time derivative of the transformed concentration at the positive inner nodes of the positive electrode (dzp/dt)
	 * 		dzn			time derivative of the transformed concentration at the positive inner nodes of the negative electrode (dzn/dt)
	 * 		dT			time derivative of the battery temperature [K s-1] (dT/dt)
	 * 		ddelta 		time derivative of the SEI thickness [m s-1] (ddelta/dt)
	 * 		dLLI 		time derivative of the lost lithium inventory [C s-1] (dLLI/dt)
	 * 		dthickp 	time derivative of the thickness of the positive electrode [m s-1] (dthickp/dt), <0 (dthickp/dt)
	 * 		dthickn		time derivative of the thickness of the negative electrode [m s-1] (dthickn/dt), <0 (dthickn/dt)
	 * 		dep			time derivative of the volume fraction in the positive electrode [s-1] (dep/dt)
	 * 		den			time derivative of the volume fraction in the negative electrode [s-1] (den/dt)
	 * 		dap			time derivative of the effective surface area of the positive electrode [m2 m-3 s-1] (dap/dt)
	 * 		dan			time derivative of the effective surface area of the negative electrode [m2 m-3 s-1] (dan/dt)
	 * 		dCS 		time derivative of the crack surface [m2 s-1], dCS/st > 0 (dCS/dt)
	 * 		dDp 		time derivative of the diffusion constant at reference temperature of the positive electrode [m s-1 s-1] (dDp/dt)
	 * 		dDn			time derivative of the diffusion constant at reference temperature of the negative electrode [m s-1 s-1] (dDn/dt)
	 * 		dR			time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
	 * 		ddelta_pl 	time derivative of the thickness of the plated lithium layer [m s-1] (ddelta_pl/dt)
	 *
	 * THROWS
	 * 100 	the array provided has the wrong length
	 * 101 	the surface li-concentration is out of bounds
	 * 109	illegal (combination of) input parameters
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::dState starting"<<endl;

	// check the input parameters
	if (nin != ns){
		cerr<<"ERROR in Cell::dState. The length of the array provided is "<<nin<<" instead of "<<ns<<endl<<flush;
		throw 100;
	}
	if ( (electr == 1 || electr == 2) && !blockDegradation){
		cerr<<"ERROR in Cell::dState. you are cycling with only one electrode "<<electr<<" but you are also accounting for degradation."
				" That is not allowed, half-cell cycling is only allowed if no degradation is considered. Either change the value of 'electr' to 3"
				" such that you cycle with the full cell, or set 'blockDegradation' to true to ignore degradation"<<endl<<flush;
		throw 109;
	}

	// get surface concentrations
		double zp[nch], zn[nch];
		double cps, cns;
		s.getZ(nch, zp, zn);										// get the transformed concentrations at the positive inner nodes
		getCSurf(&cps, &cns);										// get the surface lithium concentration
		double zp_surf = (cps / Cmaxpos);							// lithium fraction (0 to 1)
		double zn_surf = (cns / Cmaxneg);

		// check if the surface concentration is within the allowed range
		// 	0 < cp < Cmaxpos
		// 	0 < cn < Cmaxneg
		if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg){
			if(print){
				cerr<<"ERROR in Cell::dState: concentration out of bounds. the positive lithium fraction is "<<cps/Cmaxpos<<" and the negative lithium fraction is "<<cns/Cmaxneg;
				cerr<<"they should both be between 0 and 1"<<endl<<flush;
			}
			throw 101;
		}

	// current density
		double i_app = Icell/elec_surf;								// current density on the electrode surfaces [A m-2]
		double jp = -i_app/(s.getAp()*n*F*s.getThickp()); 			// molar flux on the positive particle [mol m-2 s-1]
		double jn = i_app/(s.getAn()*n*F*s.getThickn());			// molar flux on the negative particle [mol m-2 s-1]
		if (electr == 1) 											// only consider positive electrode, ignore the negative electrode
			jn = 0;
		if (electr == 2)											// only consider negative electrode, ignore the positive electrode
			jp = 0;

	// Arrhenius relation for temperature-dependent parameters
		double Dpt = s.getDp()*exp(Dp_T/Rg*(1/T_ref - 1/s.getT())); // Diffusion constant at the positive electrode at the cell's temperature [m s-1]
		double Dnt = s.getDn()*exp(Dn_T/Rg*(1/T_ref - 1/s.getT()));	// Diffusion constant at the negative electrode at the cell's temperature [m s-1]
		double kpt = kp*exp(kp_T/Rg*(1/T_ref - 1/s.getT()));		// Rate constant at the positive electrode at the cell's temperature [m s-1]
		double knt = kn*exp(kn_T/Rg*(1/T_ref - 1/s.getT()));		// Rate constant at the negative electrode at the cell's temperature [m s-1]
		if (electr == 1){ 											// only consider positive electrode, ignore the negative electrode
			Dnt = 0;
			knt = 0;
		}
		if (electr == 2){											// only consider negative electrode, ignore the positive electrode
			Dpt = 0;
			kpt = 0;
		}

	// Calculate the effect of the main li-reaction on the (transformed) concentration
		double ctep, cten, dzp[nch], dzn[nch];
		for (int j=0;j<nch;j++){									// loop for each row of the matrix-vector product A * z
			ctep = M.Ap[j]*zp[j]; 									// A is diagonal, so the array M.A has only the diagonal elements
			cten = M.An[j]*zn[j];
			dzp[j] = (Dpt*ctep+M.Bp[j]*jp);							// dz/dt = D * A * z + B * j
			dzn[j] = (Dnt*cten+M.Bn[j]*jn);
		}

	// Calculate the overpotential using the Bulter-Volmer equation
		// if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
		// and asinh(x) = ln(x + sqrt(1+x^2)
		double i0p = kpt*n*F*sqrt(C_elec)*sqrt(cps)*sqrt(Cmaxpos-cps); 	// exchange current density of the positive electrode
		double i0n = knt*n*F*sqrt(C_elec)*sqrt(cns)*sqrt(Cmaxneg-cns);	// exchange current density of the negative electrode
		double xp = -0.5*i_app/(s.getAp()*s.getThickp()) / i0p;			// x for the cathode
		double xn = 0.5*i_app/(s.getAn()*s.getThickn()) / i0n;			// x for the anode
		double etap = (2*Rg*s.getT())/(n*F) * log(xp + sqrt(1+xp*xp));	// cathode overpotential [V], < 0 on discharge
		double etan = (2*Rg*s.getT())/(n*F) * log(xn + sqrt(1+xn*xn));	// anode overpotential [V],  > 0 on discharge
		if (electr == 1) 												// only consider positive electrode, ignore the negative electrode
			etan = 0;
		if (electr == 2)												// only consider negative electrode, ignore the positive electrode
			etap = 0;

	// Calculate the entropic coefficient
		double dOCV;
		bool bound = true;												// in linear interpolation, throw an error if you are outside of the allowed range of the data
		try{
			dOCV = linInt(print, bound, dOCV_tot_x, dOCV_tot_y, dOCV_tot_n, zp_surf);	// entropic coefficient of the entire cell OCV [V K-1]
		}
		catch(int e){
			if(print)
				cout<<"Error in Cell::dState when calculating the entropic coefficient "<<e<<". Throwing it on"<<endl<<flush;
			throw e;
		}

	// temperature
		// Calculate the thermal sources/sinks/transfers per unit of volume of the battery
		// The battery volume is given by the product of the cell thickness and the electrode surface
		double Qrev = -i_app/L*s.getT()*dOCV;						// reversible heat due to entropy changes [W m-3]
		double Qrea = i_app/L * (etan-etap);						// reaction heat due to the kinetics [W m-3]
		double Qohm = pow(Icell,2)*getR() / (L*elec_surf);			// Ohmic heat due to electrode resistance [W m-3]
		double Qc = -Qch*SAV*(s.getT() - T_env);					// cooling with the environment [W m-3]

	// If we ignore degradation in this time step, we have calculated everything we need
		if (blockDegradation){
			for (int j=0;j<nch;j++){
				dstates[j] = dzp[j];								// dzp		diffusion
				dstates[nch+j] = dzn[j];							// dzn
			}
			dstates[2*nch + 0] = 1/(rho*Cp)*(Qrev+Qrea+Qohm+Qc);	// dT		cell temperature
			dstates[2*nch + 1] = 0;									// ddelta	SEI thickness
			dstates[2*nch + 2] = 0;									// dLLI		lost lithium
			dstates[2*nch + 3] = 0;									// dthickp 	electrode thickness
			dstates[2*nch + 4] = 0;									// dthickn
			dstates[2*nch + 5] = 0;									// dep		volume fraction of active material
			dstates[2*nch + 6] = 0;									// den
			dstates[2*nch + 7] = 0;									// dap		effective surface are, a = 3 e/R
			dstates[2*nch + 8] = 0;									// dan
			dstates[2*nch + 9] = 0;									// dCS		surface area of the cracks
			dstates[2*nch + 10] = 0;								// dDp 		diffusion constant
			dstates[2*nch + 11] = 0;								// dDn
			dstates[2*nch + 12]= 0;									// dR 		electrode resistance
			dstates[2*nch + 13]= 0;									// ddelta_pl thickness of the plated lithium

			// check that we have updated all the defined states
			assert(ns == 2*nch + 13 + 1); 	// check that ns has the correct value.
			// if this assertion fails, the user has changed something in the code at some point, without accounting for this change somewhere else.

			if(verbose >= printCellFunctions)
				cout<<"Cell::dState terminating without degradation"<<endl;

			return;													// stop calculating
		}

	// calculate the anode potential (needed for various degradation models)
		double OCV_n, dOCVn;										// anode potential at reference temperature and entropic coefficient
		try{
			dOCVn = linInt(print, bound, dOCV_neg_x, dOCV_neg_y, dOCV_neg_n, zn_surf);// entropic coefficient of the anode potential [V K-1]
			OCV_n = linInt(print, bound, OCV_neg_x, OCV_neg_y, OCV_neg_n, zn_surf);	// anode potential [V]
		}
		catch(int e){
			if(print)
				cout<<"Error in Cell::dState when calculating the anode potential "<<e<<". Throwing it on"<<endl<<flush;
			throw e;
		}
		double OCVnt = OCV_n + (s.getT()-T_ref)*dOCVn;				// anode potential at the cell's temperature [V]

	// SEI growth
		double isei;												// current density of the SEI growth side reaction [A m-2]
		double den_sei;												// decrease in volume fraction due to SEI growth [s-1]
		double dznsei[nch];											// additional diffusion in the anode due to isei
		try{
			SEI(OCVnt, etan, &isei, &den_sei);
		}
		catch(int e){
			if(print)
				cout<<"Error in Cell::dState when calculating the effect of SEI growth: "<<e<<". Throwing it on"<<endl<<flush;
			throw e;
		}

		// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF)
		for (int j=0;j<nch;j++)
			dznsei[j] = (M.Bn[j]*isei/(nsei*F));

	// crack growth leading to additional exposed surface area
		double isei_multiplyer;										// relative increase in isei due to additional SEI growth on the extra exposed surface area [-]
		double dCS;													// increase in crack surface area [m2 s-1]
		double dDn;													// change in negative diffusion constant [m s-1 s-1]
		double dznsei_CS[nch];										// additional diffusion in the anode due to extra SEI growth on the crack surface
		try{
			CS(OCVnt, etan, &isei_multiplyer, &dCS, &dDn);
		}
		catch(int e){
			if(print)
				cout<<"Error in Cell::dState when calculating the effect of crack growth: "<<e<<". Throwing it on"<<endl<<flush;
			throw e;
		}

		// crack surface leads to extra SEI growth because the exposed surface area increases.
		// (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF + isei_CS/nF)
		double isei_CS = isei*isei_multiplyer;						// extra SEI side reaction current density due to the crack surface area [A m-2]
		for (int j=0;j<nch;j++)
			dznsei_CS[j] = (M.Bn[j]*isei_CS/(nsei*F));

	// loss of active material LAM
		double dthickp, dthickn, dap, dan, dep, den;				// change in geometric parameters describing the amount of active material
		try{
			LAM(print, zp_surf, etap, &dthickp, &dthickn, &dap, &dan, &dep, &den);
		}
		catch(int e){
			if(print)
				cout<<"Error in Cell::dState when calculating the LAM: "<<e<<". Throwing it on"<<endl<<flush;
			throw e;
		}

	// lithium plating
		double ipl;													// current density of the plating side reaction [A m-2]
		double dzn_pl[nch];											// additional diffusion in the anode due to ipl
		try{
			LiPlating(OCVnt, etan, &ipl);
		}
		catch(int e){
			if(print)
				cout<<"Error in Cell::dState when calculating the lithium plating: "<<e<<". Throwing it on"<<endl<<flush;
			throw e;
		}

		// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + ipl/nF)
		for (int j=0;j<nch;j++)
			dzn_pl[j] = (M.Bn[j]*ipl/(npl*F));

	// time derivatives
		for (int j=0;j<nch;j++){
			dstates[j] = dzp[j];																	// dzp 		diffusion
			dstates[nch+j] = (dzn[j] + dznsei[j] + dznsei_CS[j] + dzn_pl[j]);						// dzn		jtot = jn + isei/nF + isei_CS/nF + ipl/nF
		}
		dstates[2*nch + 0] = 1/(rho*Cp)*(Qrev+Qrea+Qohm+Qc);										// dT 		cell temperature
		dstates[2*nch + 1] = isei / (nsei*F*rhosei);												// ddelta	thickness of the SEI layer
			// delta uses only isei (and not isei + isei_CS) since crack growth increases the area, not the thickness
		dstates[2*nch + 2] = (isei + isei_CS + ipl)*elec_surf*s.getThickn()*s.getAn();				// dLLI 	loss of lithium
			// i_sei = density => * active surface area = * (surf*thick*specific_surf_neg)
		dstates[2*nch + 3] = dthickp;																// dthickp 	electrode thickness
		dstates[2*nch + 4] = dthickn;																// dthickn
		dstates[2*nch + 5] = dep;																	// dep		volume fraction of active material
		dstates[2*nch + 6] = den + den_sei;															// den
		dstates[2*nch + 7] = dap + 3/Rp*dstates[2*nch + 5];											// dap		effective surface area, a = 3 e/R -> da/dt = da/dt + 3/R de/dt
		dstates[2*nch + 8] = dan + 3/Rn*dstates[2*nch + 6];											// dan
		dstates[2*nch + 9] = dCS;																	// dCS 		surface area of the cracks
		dstates[2*nch + 10] = 0;																	// dDp 		diffusion constant
		dstates[2*nch + 11] = dDn;																	// dDn
		dstates[2*nch + 12]= 0;																		// dR 		specific electrode resistance
		dstates[2*nch + 13]= ipl / (npl*F*rhopl);													// ddelta_pl thickness of the plated lithium

		// check that we have updated all the defined states
		assert(ns == 2*nch + 13 + 1); 	// check that ns has the correct value.
		// if this assertion fails, the user has changed something in the code at some point, without accounting for this change somewhere else.

		if(verbose >= printCellFunctions)
			cout<<"Cell::dState terminating with degradation"<<endl;
		return;
}

void Cell::ETI(bool print, double dti, bool blockDegradation){
	/*
	 * Performs forward Euler time integration over one time step of dti seconds
	 * s(t+1) = s(t) + ds/dt * dti
	 *
	 * IN
	 * print 			boolean indicating if we want to print error messages or not
	 * 						if true, error messages are printed
	 * 					if false no error messages are printed (but the error will still be thrown)
	 * 					we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * dti 				time step over which time integradation should be done [s]
	 * 					it should be small enough to ensure numerical stability and accuracy
	 * 					1-5 seconds in usually ok (depending on the magnitude of the current, the temperature, etc.)
	 * blockDegradation if true, degradation is not accounted for in this time step
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::ETI starting"<<endl;

	// Update the stress values stored in the attributes with the stress of the previous time step
	s_dai_p_prev = s_dai_p;							// Dai's stress in the positive particle in the previous time step
	s_dai_n_prev = s_dai_n;							// Dai's stress in the negative particle in the previous time step
	s_lares_n_prev = s_lares_n;						// Laresgoiti's stress in the negative particle in the previous time step

	// Calculate the stress values stored in the attributes for the stress in this time step
	if(s_dai) 										// only a few degradation models actually need the stress according to Dai, so only calculate it if needed
		updateDaiStress();
	if(s_lares)										// only a few degradation models need the stress according to Laresgoiti
		updateLaresgoitiStress(print);

	// Calculate the time derivatives
	double states[ns], dstates[ns];					// arrays with the state and dstate/dt
	s.getStates(ns, states);						// get states
	dState(print, blockDegradation, 0, ns, dstates);	// calculate time derivatives, electr = 0 to account for both electrodes (i.e. cycle the full cell)

	// forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
	for (int s=0;s<ns;s++)							// loop for all states
		states[s] = states[s] + dti * dstates[s];
	s.setStates(ns, states);						// store new states, checks if they are illegal (throws an error in that case)

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	s_dai_update = false;
	s_lares_update = false;

	if(verbose >= printCellFunctions)
		cout<<"Cell::ETI terminating"<<endl;
}

void Cell::ETI_electr(bool print, double I, double dti, bool blockDegradation, bool pos){
	/*
	 * Performs forward Euler time integration over one time step of dti seconds
	 * Considers only one electrode, i.e. you can do half-cell cycling with only one electrode
	 * This allows to calculate the half-cell OCV curves (versus a reference electrode of metallic lithium)
	 *
	 * USE THIS FUNCTION WITH CARE:
	 * It doesn't check the input current
	 * It might get the cell to an illegal situation
	 * You can't undo the effects of this function, unless you had stored the states before calling this function
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * I				current applied during this time step [A]
	 * 						> 0 for discharge
	 * 						< 0 for charge
	 * dti				time step [s]
	 * blockDegradation if true, degradation is not accounted for in this time step. Must be true
	 * pos 				if true, the positive electrode is considered (half-cell cycling with the positive electrode)
	 * 					if false, the negative electrode is considered (half-cell cycling with the negative electrode)
	 *
	 * THROWS
	 * 109 				illegal input parameters
	 */

	if(verbose >= printCellFunctions)
		cout<<"Cell::ETI_electr starting"<<endl;

	if(!blockDegradation){
		cerr<<"ERROR in Cell::ETI_electr, you are cycling only one electrode but want to account for degradation. This is not allowed"<<endl<<flush;
		 // half-cell cycling is only allowed if blockDegradation is true (i.e. you ignore degradation)
		 // this is because some of the degradation mechanisms could give NaN or Inf due to a divide by 0
		// So you can only call this function with 'true' for blockDegradation
		throw 109;
	}

	// Set the specified current
	Icell = I;												// don't call setI because that function ramps the current on both electrodes
															// so here 'cheat it' and directly set the current
															// this means you avoid the checks done in setI, so you don't know if the current is feasible or not

	// Calculate the time derivatives
	double states[ns], dstates[ns];							// arrays with the state and dstate/dt
	s.getStates(ns, states);								// get states
	if (pos)
		dState(print, blockDegradation, 1, ns, dstates);	// calculate time derivatives of the positive electrode
	else
		dState(print, blockDegradation, 2, ns, dstates);	// calculate time derivatives of the negative electrode

	// forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
	for (int s=0;s<ns;s++)									// loop for all states
		states[s] = states[s] + dti * dstates[s];
	s.setStates(ns, states);								// store new states

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	s_dai_update = false;
	s_lares_update = false;

	if(verbose >= printCellFunctions)
		cout<<"Cell::ETI_electr terminating"<<endl;
}
