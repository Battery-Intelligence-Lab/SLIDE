/*
 * Cell.cpp
 *
 * Implements the functions for the parent class of the Cells
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>

#include "cell.hpp"
#include "read_CSVfiles.h"
#include "interpolation.h"
#include "util.hpp"
#include "constants.hpp"
#include "param/cell_param.hpp"

void Cell::getStates(slide::State &si, double *I)
{
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

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getStates(State, double*) starting\n";

	try
	{
		validState(); // throw 15 if the state is illegal
	}
	catch (int e)
	{
		std::cout << "Error in State::setStates(double states[]), the suggested state is illegal: " << e << ". throwing it on.\n";
		throw e;
	}
	//	validState(); // throw 15 if the state is illegal
	si = s;
	*I = Icell;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getStates(State, double*) terminating.\n";
}

void Cell::getCSurf(double *cps, double *cns)
{
	/*
	 * Calculates the surface concentration at each particle.
	 * Uses the matrices from the Model struct with the spatial discretisation for the solid diffusion PDE.
	 *
	 * OUT
	 * cps 	surface li-concentration at the positive particle [mol m-3]
	 * cns 	surface li-concentration at the negative particle [mol m-3]
	 */

	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getCSurf starting\n";

	// get the transformed concentrations at the inner nodes from the state-object of the cell

	// Calculate the diffusion constant at the battery temperature using an Arrhenius relation
	const double Dpt = s.get_Dp() * std::exp(Dp_T / Rg * (1 / T_ref - 1 / s.get_T())); // diffusion constant of the positive particle [m s-1]
	const double Dnt = s.get_Dn() * std::exp(Dn_T / Rg * (1 / T_ref - 1 / s.get_T())); // diffusion constant of the negative particle [m s-1]

	// Calculate the molar flux on the surfaces
	const double jp = -Icell / (s.get_ap() * elec_surf * s.get_thickp()) / (n * F); // molar flux on the positive particle [mol m-2 s-1]
	const double jn = Icell / (s.get_an() * elec_surf * s.get_thickn()) / (n * F);	// molar flux on the negative particle [mol m-2 s-1]

	// Calculate the surface concentration at the positive particle
	// 	cp_surf = M.Cp[0][:] * zp[:] + M.Dp*jp/Dpt
	double cp_surf = 0;
	for (int j = 0; j < settings::nch; j++)
		cp_surf += M.Cp[0][j] * s.get_zp(j);
	cp_surf += M.Dp[0] * jp / Dpt;

	// Calculate the surface concentration at the negative particle
	// 	cn_surf = M.Cn[0][:] * zn[:] + M.Dn*jn/Dnt
	double cn_surf = 0;
	for (int j = 0; j < settings::nch; j++)
		cn_surf += M.Cn[0][j] * s.get_zn(j);
	cn_surf += M.Dn[0] * jn / Dnt;

	// Make the output parameters
	*cns = cn_surf;
	*cps = cp_surf;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getCSurf terminating\n";
}

void Cell::getC(double cp[], double cn[])
{
	/*
	 * Calculates the lithium concentration at each positive Chebyshev node (0 <= x <= 1), including the centre and surface nodes.
	 * Uses the matrices from the state space matrices.
	 * See the equations in the explanatory documents.
	 *
	 * OUT
	 * cp 	li-concentration at each Chebyshev node in the positive electrode [mol m-3], length of the array should be nch+2
	 * 			cp[0]			concentration at the surface of the sphere
	 * 			cp[1 to nch]	concentration at the inner nodes
	 * 			cp[nch + 1]		concentration at the centre of the sphere
	 * cn 	li-concentration at each Chebyshev node in the negative electrode [mol m-3], length of the array should be nch+2
	 * 			cn[0]			concentration at the surface of the sphere
	 * 			cn[1 to nch]	concentration at the inner nodes
	 * 			cn[nch + 1]		concentration at the centre of the sphere
	 *
	 * THROWS
	 * 100 	the arrays provided have the wrong length
	 */

	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getC starting.\n";

	// get the transformed concentrations at the inner nodes from the state-object of the cell
	//const auto &zp{s.get_zp()}, &zn{s.get_zn()};

	// Calculate the diffusion constant at the battery temperature using an Arrhenius relation
	const double Dpt = s.get_Dp() * std::exp(Dp_T / Rg * (1 / T_ref - 1 / s.get_T())); // diffusion constant of the positive particle [m s-1]
	const double Dnt = s.get_Dn() * std::exp(Dn_T / Rg * (1 / T_ref - 1 / s.get_T())); // diffusion constant of the negative particle [m s-1]

	// Calculate the molar flux on the surfaces
	const double jp = -Icell / (s.get_ap() * elec_surf * s.get_thickp()) / (n * F); // molar flux on the positive particle [mol m-2 s-1]
	const double jn = Icell / (s.get_an() * elec_surf * s.get_thickn()) / (n * F);	// molar flux on the negative particle [mol m-2 s-1]

	// Calculate concentration at the surface and inner nodes using the matrices from the spatial discretisation of the solid diffusion PDE
	// 	cp = M.Cp[:][:] * zp[:] + M.Dp*jp/Dpt
	// 	cn = M.Cn[:][:] * zn[:] + M.Dn*jn/Dnt
	for (int i = 0; i < settings::nch + 1; i++) // Problem here!!!!!!! #CHECK
	{											// loop to calculate at each surface + inner node
		double cpt{0}, cnt{0};
		for (int j = 0; j < settings::nch; j++)
		{
			cpt += M.Cp[i][j] * s.get_zp(j);
			cnt += M.Cn[i][j] * s.get_zn(j);
		}
		cp[i] = cpt + M.Dp[i] * jp / Dpt;
		cn[i] = cnt + M.Dn[i] * jn / Dnt;
	}

	// Calculate the concentration at centre node using the boundary condition (the concentration gradient at the centre has to be 0 due to symmetry)
	// cp_centre = -1/2 (M.Cc[:]*cp +jp*Rp/Dpt)
	// cn_centre = -1/2 (M.Cc[:]*cn +jn*Rn/Dnt)
	constexpr double DM = 2.0; // we need a constant of 2 in the equations
	double cpt{0}, cnt{0};
	for (int i = 0; i < settings::nch + 1; i++)
	{
		cpt += M.Cc[i] * cp[i];
		cnt += M.Cc[i] * cn[i];
	}
	cp[settings::nch + 1] = (-1.0 / DM) * (cpt + jp * Rp / Dpt);
	cn[settings::nch + 1] = (-1.0 / DM) * (cnt + jn * Rn / Dnt);

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getC terminating.\n";
}

bool Cell::getVoltage(bool print, double *V, double *OCVp, double *OCVn, double *etap, double *etan, double *Rdrop, double *Temp)
{
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

	// #HOTFUNC
	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getVoltage starting\n";

	// Get the surface concentrations
	double cps, cns;
	getCSurf(&cps, &cns);

	// check if the surface concentration is within the allowed range
	// 	0 < cp < Cmaxpos
	// 	0 < cn < Cmaxneg
	// don't allow 0 or Cmax because in that case, i0 will be 0, and the overpotentials will have 1/0 = inf or nan
	if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg)
	{
		if (print)
		{ // print error message unless you want to suppress the output
			std::cerr << "ERROR in Cell::getVoltage: concentration out of bounds. the positive lithium fraction is " << cps / Cmaxpos
					  << " and the negative lithium fraction is " << cns / Cmaxneg << " they should both be between 0 and 1.\n";
		}
		*V = nan("double"); // set the voltage to nan (Not A Number)
		throw 101;
		return false; // the voltage is not within the limits
	}
	else
	{
		// Calculate the li-fraction (instead of the li-concentration)
		const double zp_surf = (cps / Cmaxpos);
		const double zn_surf = (cns / Cmaxneg);

		// Calculate the electrode potentials
		double OCV_p;	   // cathode potential [V]
		double OCV_n;	   // anode potential [V]
		double dOCV;	   // entropic coefficient of the total cell voltage [V/K]
		bool bound = true; // in linear interpolation, throw an error if you are out of the allowed range
		try
		{
			dOCV = OCV_curves.linInt_dOCV_tot(zp_surf, print, bound); // Question: Why do we use zp_surf?
			OCV_n = OCV_curves.linInt_OCV_neg(zn_surf, print, bound);
			OCV_p = OCV_curves.linInt_OCV_pos(zp_surf, print, bound);
		}
		catch (int e)
		{
			if (print)
				std::cout << "error in Cell::getVoltage when getting the electrode potentials " << e << ". Throwing it up.\n";
			throw e;
		}

		// Calculate the rate constants at the cell's temperature using an Arrhenius relation
		const double kpt = kp * std::exp(kp_T / Rg * (1 / T_ref - 1 / s.get_T()));
		const double knt = kn * std::exp(kn_T / Rg * (1 / T_ref - 1 / s.get_T()));

		// Calculate the overpotential using the Bulter-Volmer equation
		// 		if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
		// 		and asinh(x) = ln(x + sqrt(1+x^2)
		const double i_app = Icell / elec_surf;											   // current density on the electrodes [I m-2]
		const double i0p = kpt * n * F * sqrt(C_elec * cps * (Cmaxpos - cps));			   // exchange current density of the positive electrode
		const double i0n = knt * n * F * sqrt(C_elec * cns * (Cmaxneg - cns));			   // exchange current density of the negative electrode
		const double xp = -0.5 * i_app / (s.get_ap() * s.get_thickp()) / i0p;			   // x for the cathode
		const double xn = 0.5 * i_app / (s.get_an() * s.get_thickn()) / i0n;			   // x for the anode
		const double etapi = (2 * Rg * s.get_T()) / (n * F) * log(xp + sqrt(1 + xp * xp)); // cathode overpotential [V], < 0 on discharge
		const double etani = (2 * Rg * s.get_T()) / (n * F) * log(xn + sqrt(1 + xn * xn)); // anode overpotential [V],  > 0 on discharge

		// Calculate the cell voltage
		// the cell OCV at the reference temperature is OCV_p - OCV_n
		// this OCV is adapted to the actual cell temperature using the entropic coefficient dOCV * (T - Tref)
		// then the overpotentials and the resistive voltage drop are added
		*V = (OCV_p - OCV_n + (s.get_T() - T_ref) * dOCV) + (etapi - etani) - getR() * Icell;

		// make the output variables
		*OCVp = OCV_p;
		*OCVn = OCV_n;
		*etap = etapi;
		*etan = etani;
		*Rdrop = -getR() * Icell;
		*Temp = s.get_T();

		if constexpr (settings::verbose >= printLevel::printCellFunctions)
			std::cout << "Cell::getVoltage terminating with V = " << *V << " and valid is " << (Vmin <= *V && *V <= Vmax) << '\n';

		// Return whether this voltage is within the allowed limits or not
		return Vmin <= *V && *V <= Vmax;
	}
}

void Cell::getDaiStress(double *sigma_p, double *sigma_n, sigma_type &sigma_r_p, sigma_type &sigma_r_n,
						sigma_type &sigma_t_p, sigma_type &sigma_t_n, sigma_type &sigma_h_p, sigma_type &sigma_h_n)
{
	/*
	 * Calculates the radial and tangential stress for each positive Chebyshev node according to the formula by
	 * Dai, Cai, White, Journal of Power sources 247, 2014
	 *
	 * It takes quite long to calculate the stress, so only call this function when needed.
	 *
	 * OUT
	 * sigma_p		maximum hydrostatic stress in the positive particle, can be both positive and negative [Pa]
	 * sigma_n		maximum hydrostatic stress in the negative particle, can be both positive and negative [Pa]
	 * sigma_r_p	array with the radial stress at each positive Chebyshev node in the positive electrode, length nch+2, [Pa]
	 * sigma_r_n	array with the radial stress at each positive Chebyshev node in the negative electrode, length nch+2, [Pa]
	 * sigma_t_p	array with the tangential stress at each positive Chebyshev node in the positive electrode, length nch+2, [Pa]
	 * sigma_t_n	array with the tangential stress at each positive Chebyshev node in the negative electrode, length nch+2, [Pa]
	 * sigma_h_p	array with the hydrostatic stress at each positive Chebyshev node in the positive electrode, length nch+2, [Pa]
	 * sigma_h_n	array with the hydrostatic stress at each positive Chebyshev node in the negative electrode, length nch+2, [Pa]
	 * 				[0]			stress at the surface of the sphere
	 * 				[1 to nch]	stress at the inner nodes
	 * 				[nch + 1]	stress at the centre of the sphere
	 *
	 * THROWS
	 * 100 	the arrays provided have the wrong length
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getDaiStress starting\n";

	// Get the locations of the Chebyshev nodes
	std::array<double, settings::nch + 2> xp; // location (x-value) of the positive Chebyshev nodes
	xp[0] = 1;
	std::copy(M.xch.begin(), M.xch.end(), xp.begin() + 1);
	xp[settings::nch + 1] = 0;

	double xtot[2 * settings::nch + 3];				 // location (x-value) of the positive and negative Chebyshev nodes [-surface .. centre .. +surface]
	xtot[settings::nch + 1] = xp[settings::nch + 1]; // centre node
	for (int i = 0; i < settings::nch + 1; i++)
	{
		xtot[i] = -xp[i];									 // negative nodes
		xtot[settings::nch + 2 + i] = xp[settings::nch - i]; // positive nodes
	}

	// get concentrations at each Chebyshev node
	// Due to symmetry, the concentration at the negative point is the same as the concentration of the positive point: c(-x) = c(x)
	double cp[settings::nch + 2], cn[settings::nch + 2]; // positive and negative nodes, [+surface .. inner .. centre]
	getC(cp, cn);
	double CP[2 * settings::nch + 3], CN[2 * settings::nch + 3]; // concentrations at all nodes, [-surface .. inner .. centre .. inner .. +surface]
	CP[settings::nch + 1] = cp[settings::nch + 1];				 // cathode centre node
	CN[settings::nch + 1] = cn[settings::nch + 1];				 // anode centre node
	for (int i = 0; i < settings::nch + 1; i++)
	{
		CP[i] = cp[i];									   // cathode negative points
		CN[i] = cn[i];									   // anode negative points
		CP[settings::nch + 2 + i] = cp[settings::nch - i]; // cathode positive points
		CN[settings::nch + 2 + i] = cn[settings::nch - i]; // anode positive points
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
	double Fp[2 * settings::nch + 3], Fn[2 * settings::nch + 3]; // arrays with the product for the positive and negative electrode
	for (int i = 0; i < 2 * settings::nch + 3; i++)
	{			   // loop for each row (one row = one node)
		Fp[i] = 0; // calculate the matrix-vector product for row i as you would do it by hand:
		Fn[i] = 0; // 	F(i) = sum(Q(i,j)*C(j)*x(j)^2*R^3, j=0..2*nch+2)
		for (int j = 0; j < 2 * settings::nch + 3; j++)
		{													  // loop through the columns to calculate the sum
			Fp[i] += M.Q[i][j] * (CP[j] * xtot[j] * xtot[j]); // 		Q(i,j)*C(j)*x(j)^2
			Fn[i] += M.Q[i][j] * (CN[j] * xtot[j] * xtot[j]);
		}
		Fp[i] = Fp[i] * Rp * Rp * Rp; // *R^3 (R is constant so it can be out of the sum)
		Fn[i] = Fn[i] * Rn * Rn * Rn;
	}

	// Calculate the integral from the centre to the positive surface, which is a constant present in all equations
	double ap = Fp[2 * settings::nch + 2] - Fp[settings::nch + 1]; // int( cp*r^2, r=0..Rp )
	double an = Fn[2 * settings::nch + 2] - Fn[settings::nch + 1]; // int( cn*r^2, r=0..Rn )

	// Calculate the equations for all nodes
	double srp[settings::nch + 2]; // radial stress in the positive particle at the positive nodes [centre .. +surface]
	double srn[settings::nch + 2]; // radial stress in the negative particle at the positive nodes [centre .. +surface]
	double stp[settings::nch + 2]; // tangential stress in the positive particle at the positive nodes [centre .. +surface]
	double stn[settings::nch + 2]; // tangential stress in the negative particle at the positive nodes [centre .. +surface]
	for (int i = 0; i < settings::nch + 2; i++)
	{ // loop for the positive nodes

		const double rp = Rp * xtot[settings::nch + 1 + i];					 // r(i) = R * x(i) radius of positive node i in the positive particle
		const double rn = Rn * xtot[settings::nch + 1 + i];					 // radius of positive node i in the negative particle
		const double bp = Fp[settings::nch + 1 + i] - Fp[settings::nch + 1]; // integral from the centre to positive node i int(cp*zp^2, zp=0..rp(i)) = F[nch+1+i] - F[nch+1]
		const double bn = Fn[settings::nch + 1 + i] - Fn[settings::nch + 1]; // integral from the centre to positive node i int(cn*zn^2, zn=0..rn(i))

		// Implement the equations from Dai et al.
		if (i == 0)
		{ // centre node -> special formula (31 & 33) in Dai, Cai, White
			srp[i] = 2 * sparam.omegap * sparam.Ep / (9 * (1 - sparam.nup)) * (3 / pow(Rp, 3.0) * ap - CP[settings::nch + 1]);
			srn[i] = 2 * sparam.omegan * sparam.En / (9 * (1 - sparam.nun)) * (3 / pow(Rn, 3.0) * an - CN[settings::nch + 1]);

			stp[i] = 2 * sparam.omegap * sparam.Ep / (9 * (1 - sparam.nup)) * (3 / pow(Rp, 3.0) * ap - CP[settings::nch + 1]);
			stn[i] = 2 * sparam.omegan * sparam.En / (9 * (1 - sparam.nun)) * (3 / pow(Rn, 3.0) * an - CN[settings::nch + 1]);
		}
		else
		{																													   // other nodes -> equation 13 in Dai, Cai, White
			srp[i] = 2 * sparam.omegap * sparam.Ep / (3 * (1 - sparam.nup)) * (1 / pow(Rp, 3.0) * ap - 1 / pow(rp, 3.0) * bp); //ap = int (c x^2, x=0..R), bp = int (c x^2 , x=0..r)
			srn[i] = 2 * sparam.omegan * sparam.En / (3 * (1 - sparam.nun)) * (1 / pow(Rn, 3.0) * an - 1 / pow(rn, 3.0) * bn);

			stp[i] = sparam.omegap * sparam.Ep / (3 * (1 - sparam.nup)) * (2 / pow(Rp, 3.0) * ap + 1 / pow(rp, 3.0) * bp - cp[i]);
			stn[i] = sparam.omegan * sparam.En / (3 * (1 - sparam.nun)) * (2 / pow(Rn, 3.0) * an + 1 / pow(rn, 3.0) * bn - cn[i]);
		}
	}

	// Flip all arrays to get the opposite order (now it is [centre .. +surface] and we want [+surface .. centre]
	// and store in the output arrays
	for (int i = 0; i < settings::nch + 2; i++)
	{ // loop for the positive nodes
		sigma_r_p[i] = srp[settings::nch + 2 - 1 - i];
		sigma_r_n[i] = srn[settings::nch + 2 - 1 - i];
		sigma_t_p[i] = stp[settings::nch + 2 - 1 - i];
		sigma_t_n[i] = stn[settings::nch + 2 - 1 - i];
	}

	// Make the hydrostatic stress sh = (sr + 2sp)/3
	int sp = 0; // node with the maximum hydrostatic stress in the positive particle
	int sn = 0; // node with the maximum hydrostatic stress in the negative
	for (int i = 0; i < settings::nch + 2; i++)
	{														  // loop for all nodes
		sigma_h_p[i] = (sigma_r_p[i] + 2 * sigma_t_p[i]) / 3; // calculate hydrostatic stress
		sigma_h_n[i] = (sigma_r_n[i] + 2 * sigma_t_n[i]) / 3;

		// find the maximum (in absolute value) of the stresses
		if (std::abs(sigma_h_p[i]) > std::abs(sigma_h_p[sp]))
			sp = i;
		if (std::abs(sigma_h_n[i]) > std::abs(sigma_h_n[sp]))
			sn = i;
	}
	*sigma_p = sigma_h_p[sp];
	*sigma_n = sigma_h_n[sn];

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getDaiStress terminating.\n";
}

void Cell::updateDaiStress()
{
	/*
	 * Function which will update the values stored in the stress variables relating with Dai's stress model
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::updateDaiStress starting.\n";

	// Make variables to store the stress
	std::array<double, settings::nch + 2> sigma_r_p, sigma_r_n, sigma_t_p, sigma_t_n, sigma_h_p, sigma_h_n;

	// Get the stress
	try
	{
		getDaiStress(&sparam.s_dai_p, &sparam.s_dai_n, sigma_r_p, sigma_r_n, sigma_t_p, sigma_t_n, sigma_h_p, sigma_h_n);
	}
	catch (int e)
	{
		std::cout << "Error in Cell::getDaiStress when calling updateDaiStress: " << e << ". throwing it on.\n";
		throw e;
	}

	// indicate that the values in the class variables are updated
	sparam.s_dai_update = true;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::updateDaiStress terminating.\n";
}

void Cell::getLaresgoitiStress(bool print, double *sigma_n)
{
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

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getLaresgoitiStress starting.\n";

	// Arrays with the stress from Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015 (figure 5)
	const std::array<double, 11> xx{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};		   // li-fraction in the graphite
	const std::array<double, 11> yy{0.0, 5.5, 8.5, 9.5, 10.0, 10.0, 10.5, 13.0, 16.5, 21.0, 23.5}; // stress [MPa]

	const bool is_xx_fixed = true; // Change this if xx is not fixed time step.

	// Get the surface concentration
	double cps, cns, zn_surf;
	getCSurf(&cps, &cns);	   // get the surface lithium concentration [mol m-3]
	zn_surf = (cns / Cmaxneg); // lithium fraction on negative surface [0 1]

	// check if the surface concentration is within the allowed range
	// 	0 < cp < Cmaxpos
	// 	0 < cn < Cmaxneg
	if (cps < 0 || cns < 0 || cps > Cmaxpos || cns > Cmaxneg)
	{
		if (print)
		{
			std::cerr << "ERROR in Cell::getLaresgoitiStress: concentration out of bounds. the positive lithium fraction is " << cps / Cmaxpos
					  << " and the negative lithium fraction is " << cns / Cmaxneg << "they should both be between 0 and 1.\n";
		}
		throw 101;
	}

	// Interpolate linearly to get the stress
	double s;
	try
	{
		s = linInt(print, true, xx, yy, 11, zn_surf, is_xx_fixed); // throw an error if you are out of the allowed range
	}
	catch (int e)
	{
		if (print)
			std::cout << "Error in Cell::getLaresgoitiStress when interpolating in the arrays "
					  << e << ". Throwing it on.\n";
		throw e;
	}

	// Make the output variable
	*sigma_n = s;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::getLaresgoitiStress terminating.\n";
}

void Cell::updateLaresgoitiStress(bool print)
{
	/*
	 * Function which will update the values stored in the stress variables relating with Laresgoiti's stress model
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::updateLaresgoitiStress starting.\n";

	double s;
	getLaresgoitiStress(print, &s);

	// Update the stored value
	sparam.s_lares_n = s;
	sparam.s_lares_update = true; // indicate that the values in the class variables are updated

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::updateLaresgoitiStress terminating.\n";
}

void Cell::setTenv(double Tenv)
{
	/*
	 * Sets the environmental temperature
	 *
	 * IN
	 * Tenv 	environmental temperature, 273 <= T <= 333 [K]
	 *
	 * THROWS
	 * 102 		illegal value of T
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setTenv starting.\n";

	// check if the environmental temperature is in the allowed limits
	bool T = (Tenv < settings::Tmin_K) || (Tenv > settings::Tmax_K); // the temperature limits are defined in State.hpp
	if (T)
	{
		std::cerr << "ERROR in Cell::setTenv, illegal value of environmental temperature " << Tenv
				  << "K. The value has to be between The value has to be between " << settings::Tmin_K << "and " << settings::Tmax_K << ".\n";
		throw 102;
	}

	// update the temperature
	T_env = Tenv;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setTenv terminating.\n";
}

void Cell::setVlimits(double VMAX, double VMIN)
{
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

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setVlimits starting with upper limit " << VMAX << " and lower limit " << VMIN << '\n';

	bool vmax = VMAX < 0 || VMAX > OCV_curves.OCV_pos_y[0];
	if (vmax)
		std::cerr << "ERROR in Cell::setVlimits. The value of the maximum voltage is " << VMAX
				  << "V but it has to be positive and lower than the maximum value of the OCV curve of the cathode, which is "
				  << OCV_curves.OCV_pos_y.back() << ".\n";
	bool vmin = VMIN < 0;
	if (vmin)
		std::cerr << "ERROR in Cell::setVlimits. The value of the minimum voltage is "
				  << VMIN << "V but it has to be positive.\n";

	if (vmax || vmin)
		throw 103;

	Vmax = VMAX;
	Vmin = VMIN;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setVlimits terminating with upper limit " << Vmax << " and lower limit " << Vmin << '\n';
}

void Cell::setT(double Ti)
{
	/*
	 * Set the temperature of the battery
	 *
	 * IN
	 * Ti		uniform cell temperature, 273 <= T <= 333 [K]
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setT starting.\n";

	s.setT(Ti); // checks if T is in limits.

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	sparam.s_dai_update = false;
	sparam.s_lares_update = false;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setT terminating\n";
}

void Cell::setStates(slide::states_type &&states)
{
	/*
	 * sets all battery state variables
	 *
	 * IN
	 * states	array with the battery states
	 *
	 * THROWS
	 * 100 		the array provided has the wrong length
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setStates(int, double[]) starting.\n";

	try
	{
		validState(); // throw 15 if the state is illegal
	}
	catch (int e)
	{
		std::cout << "Error in State::setStates(double states[]), the suggested state is illegal: " << e << ". throwing it on.\n";
		throw e;
	}

	s.setStates(std::forward<slide::states_type>(states));

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	sparam.s_dai_update = false;
	sparam.s_lares_update = false;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setStates(int, double[]) terminating.\n";
}
void Cell::setStates(const slide::State &si, double I)
{
	/*
	 * Set all battery state variables to the ones in the State-object provided.
	 * Set the cell current to the value provided (without ramping the current since we are changing the state any way)
	 *
	 * IN
	 * si 	new state of the cell
	 * I 	new current of the cell [A], positive for discharging, negative for charging
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setStates(State, double) starting.\n";

	s = si;
	Icell = I;

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	sparam.s_dai_update = false;
	sparam.s_lares_update = false;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setStates(State, double) terminating.\n";
}

void Cell::setC(double cp0, double cn0)
{
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
	// #NOTHOTFUNCTION
	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setC starting.\n";

	bool pp = (cp0 < 0) || (cp0 > 1);
	if (pp)
		std::cerr << "ERROR in Cell::setC, illegal input li-fraction for the positive electrode : " << cp0
				  << ". The value has to be between 0 and 1.\n";

	bool nn = (cn0 < 0) || (cn0 > 1);
	if (nn)
		std::cerr << "ERROR in Cell::setC, illegal input li-fraction for the negative electrode : " << cn0
				  << ". The value has to be between 0 and 1.\n";
	if (pp || nn)
		throw 104;

	// Calculate the corresponding li-concentrations in [mol m-3]
	const double cp = cp0 * Cmaxpos;
	const double cn = cn0 * Cmaxneg;

	// Do the first transformation, to u(i) = radius(i) * concentration = x(i) * R * concentration(i)
	double uneg[settings::nch], upos[settings::nch];
	for (int i = 0; i < settings::nch; i++)
	{
		uneg[i] = cn * M.xch[i] * Rn;
		upos[i] = cp * M.xch[i] * Rp;
	}

	// The second transformation is to the eigenspace: z = V * u with V the inverse of the matrix with the eigenvectors.
	// As explained, we know that there is one eigenvector corresponding to a uniform concentration
	// So we need to calculate only this one nonzero value for the (twice) transformed concentration
	// The location of the uniform eigenvector (which has a 0 eigenvalue) is written in M.Input[3]
	const int ind = M.Input[3];
	double znu = 0; // (twice) transformed negative uniform concentration
	double zpu = 0; // (twice) transformed positive uniform concentration

	for (int i = 0; i < settings::nch; i++)
	{ // loop to calculate the row of V * u corresponding to the uniform concentration
		znu += M.Vn[ind][i] * uneg[i];
		zpu += M.Vp[ind][i] * upos[i];
	}

	// Make the full arrays for the (twice) transformed concentration
	slide::z_type zp{}, zn{}; // set all values to 0
	zp[ind] = zpu;			  // set the non-zero value
	zn[ind] = znu;
	s.setZ(zp, zn);

	// Set the cell current to 0 to reflect the boundary condition for a fully uniform concentration
	Icell = 0;

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	sparam.s_dai_update = false;
	sparam.s_lares_update = false;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setC terminating.\n";
}

void Cell::setI(bool print, bool check, double I)
{
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
	// #HOTFUNC -> ~500k calls over 20M.
	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::setI starting with current " << I << '\n';

	// check if the specified value is different from the actual cell current
	if (std::abs(Icell - I) < 1e-10)
		return; // the values are the same -> we don't need to do anything

	// Store the old current and state to restore it if needed
	double Iold;
	slide::State sold;
	getStates(sold, &Iold);

	// settings
	bool blockDegradation = true; // don't account for degradation while the current is ramping
	bool reached = false;		  // boolean indicating if we have reached the set current
	int sign;					  // sign whether we should increase or decrease the current
	if (I > Icell)
		sign = 1;
	else
		sign = -1;

	// loop to ramp the current
	while (!reached)
	{
		// increase the current
		Icell += sign * dIcell;

		// check if you have reached the specified current
		if (sign * Icell > sign * I)
		{ // increase I: Icell > I, decrease I: Icell < I
			Icell = I;
			reached = true;
		}

		// take one small time step
		ETI(print, dt_I, blockDegradation);
	}

	// Check the cell's conditions are still valid if we want to check the final state
	double v, ocvp, ocvn, etap, etan, rdrop, tem;
	bool valid;
	if (check)
	{
		// check the state
		try
		{
			validState(); // throws an error if the state is illegal
		}
		catch (int e)
		{
			if (print)
				std::cout << "Cell::setI illegal state after setting the current to " << Icell << ", error: " << e << ". Throwing an error.\n";
			setStates(sold, Iold); // restore the original battery state and current
			throw e;
		}

		// check the voltage
		try
		{
			valid = getVoltage(print, &v, &ocvp, &ocvn, &etap, &etan, &rdrop, &tem); // throws an error if the surface concentration is out of bounds
		}
		catch (int e)
		{
			valid = false; // in that case, the voltage is illegal
		}
		if (!valid)
		{
			if (print)
				std::cerr << "Cell::setI Illegal voltage after trying to set the current to " << Icell << ", the voltage is: " << v << "V. Throwing an error.\n";
			setStates(sold, Iold); // restore the original battery state and current
			throw 105;
		}
	}

	// the stress values stored in the class variables for stress are no longer valid because the cell current has changed
	sparam.s_dai_update = false;
	sparam.s_lares_update = false;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
	{
		if (check)
			std::cout << "Cell::setI terminating with current " << I << " and voltage " << v << '\n';
		else
			std::cout << "Cell::setI terminating with current " << I << " without checking the voltage.\n";
	}
}

void Cell::SEI(double OCVnt, double etan, double *isei, double *den)
{
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

	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::SEI starting\n";

	// variables
	double kseit, Dseit;		// temperature-dependent SEI parameters, using Arrhenius' relation
	double isei1, isei2, isei3; // temporary parameters to calculate the SEI growth
	double is = 0;				// SEI side reaction current density of all models combined

	// check that the number of degradation models is not too long
	if (deg_id.SEI_n > deg_id.len)
	{
		std::cerr << "ERROR in Cell::SEI the user wants to use more than " << deg_id.len << " degradation models. throwing an error.\n";
		std::cout << "The maximum length of the array with the degradation identifiers is " << deg_id.len
				  << ". If you want to use more models, you have to increase the value of 'len' in the struct DEG_ID, defined in Cell.hpp.\n";
		throw 107;
	}

	// Loop for each model to use
	for (int i = 0; i < deg_id.SEI_n; i++)
	{

		// Use a switch to calculate the magnitude of the SEI growth according to this degradation model
		switch (deg_id.SEI_id[i])
		{
		case 0: // no SEI growth
			is += 0;
			break;
		case 1:																																	 // Kinetic model according to Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004
			kseit = seiparam.sei1k * std::exp(seiparam.sei1k_T / Rg * (1 / T_ref - 1 / s.get_T()));												 // Arrhenius relation for the rate parameter at the cell temperature
			is += nsei * F * kseit * std::exp(-nsei * F / (Rg * s.get_T()) * alphasei * (OCVnt + etan - OCVsei + Rsei * s.get_delta() * Icell)); // Add the effect of this model
																																				 // eta_sei = OCVneg + etaneg - OCVsei + Rsei*I
																																				 // isei = nFk std::exp(-nF/RT alpha eta_sei)
																																				 // on charge, I < 0 and etan < 0.
																																				 // so higher charging current -> more negative term in exponential -> larger isei
			break;
		case 2:																						// kinetics and diffusion according to Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
			kseit = seiparam.sei2k * std::exp(seiparam.sei2k_T / Rg * (1 / T_ref - 1 / s.get_T())); // Arrhenius relation for the rate parameter at the cell temperature
			Dseit = seiparam.sei2D * std::exp(seiparam.sei2D_T / Rg * (1 / T_ref - 1 / s.get_T())); // Arrhenius relation for the diffusion constant at the cell temperature

			// derivation of the formula:
			// start equations (use the symbols from Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
			// but with opposite sign for j
			// j = nFk c exp(..)
			// j/nF = - D/delta (c - c0)
			// j/nF = D/delta (- j/(nFk exp(..)) + c0)
			// j * ( 1/(nFk exp(..)) + delta/DnF ) = c0
			// j = c0 / ( 1/(nFk exp(..)) + delta/DnF )
			isei2 = nsei * F * kseit * std::exp(-nsei * F / (Rg * s.get_T()) * alphasei * (OCVnt + etan - OCVsei + Rsei * s.get_delta() * Icell));
			isei3 = s.get_delta() / (nsei * F * Dseit);
			is += c_elec0 / (1 / isei2 + isei3); // Add the effects of this model
			break;
		case 3:																						// model from Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
			kseit = seiparam.sei3k * std::exp(seiparam.sei3k_T / Rg * (1 / T_ref - 1 / s.get_T())); // Arrhenius relation for the rate parameter at the cell temperature
			Dseit = seiparam.sei3D * std::exp(seiparam.sei3D_T / Rg * (1 / T_ref - 1 / s.get_T())); // Arrhenius relation for the diffusion constant at the cell temperature

			// Use equation [22] from the paper
			isei1 = 0.134461 * std::exp(-nsei * F * (etan + Rsei * s.get_delta() * Icell) / (Rg * s.get_T())); // the parameter a_L_K is set to 0.134461 but this constant can be lumped into the rate- and diffusion constants
			isei2 = nsei * F * kseit * std::exp(-nsei * F / (Rg * s.get_T()) * alphasei * (OCVnt - OCVsei));
			isei3 = s.get_delta() / (nsei * F * Dseit);
			is += isei1 / (1 / isei2 + isei3); // Add the effects of this model
			break;
		default: // unknown degradation model
			std::cerr << "ERROR in Cell::SEI, unknown SEI degradation model with identifier " << deg_id.SEI_id[i] << ". Only values 0 to 3 are allowed. Throw an error"
					  << '\n'
					  << std::flush;
			throw 106;
			break;
		}
	} // end loop for all the models you want to use

	// Make the output for the SEI side reaction current density
	*isei = is;

	// Calculate how much we decrease the volume fraction due to SEI growth
	if (deg_id.SEI_porosity == 0) // don't decrease volume fraction
		*den = 0;
	else if (deg_id.SEI_porosity == 1)
	{																		   // decrease volume fraction according to Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
		double jn = Icell / elec_surf / (s.get_an() * n * F * s.get_thickn()); // molar flux on the negative particle
		*den = -seiparam.sei_porosity * (jn * Vmain + *isei * Vsei);
		// note: they use J = volumetric current [A/m3] -> they multiply with 'an' but we already have density [A/m2]
		// - because they use the porosity while we use the volume fraction
	}
	else
	{ // unknown degradation model
		std::cerr << "ERROR in Cell::SEI, unknown value for decreasing the volume fraction " << deg_id.SEI_porosity << ". Only values 0 or 1 are allowed. Throw an error.\n";
		throw 106;
	}

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::SEI terminating\n";
}

void Cell::CS(double OCVnt, double etan, double *isei_multiplyer, double *dCS, double *dDn)
{
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

	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::CS starting.\n";

	// output parameters
	double ism = 0; // isei multiplier from all models to be considered
	double dcs = 0; // increase in surface area from all models to be considered

	if (deg_id.CS_n > deg_id.len)
	{
		std::cerr << "ERROR in Cell::CS the user wants to use more than " << deg_id.len << " degradation models. throwing an error.\n";
		std::cout << "The maximum length of the array with the degradation identifiers is " << deg_id.len << ". If you want to use more models,"
				  << " you have to increase the value of 'len' in the struct DEG_ID, defined in Cell.hpp.\n";
		throw 107;
	}

	double ASn = getAnodeSurface(); // active surface area of the anode [m2]
									// this active area is used to translate absolute values (such as currents) to relative values (such as current densities)

	// Loop for each model we want to use
	for (int i = 0; i < deg_id.CS_n; i++)
	{

		// a switch to calculate the effect according to model i
		switch (deg_id.CS_id[i])
		{
		case 0: // no surface cracks
			ism += 0;
			dcs += 0;
			break;
		case 1: // Laresgoiti's stress and crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
				// this model calculates crack growth due to temporal variations in the li-concentration
			// check the calculated stress values are up to date
			if (!sparam.s_lares_update)
			{
				std::cerr << "ERROR in Cell::CS. The stress values for Laresgoiti's stress model are not updated. Throwing an error.\n";
				// if you see this error, you have to call Cell::updateLaresgoitiStress(), which calculates the stress and stores the values, before you call this function
				throw 108;
			}

			// Implement the equation from the paper
			// capacity loss is m-power of the so-called stress amplitude (sigma_max - sigma_min)/2
			// sigma_max and sigma_min are the max and min stresses 'of the cyclic signal' i.e. within one charge/discharge
			// assume m = 1, then (max - min) = (max - t1) + (t1-t2) + (t2-t3) + ... + (tn - min)
			// so stress amplitude can be substituted by (the stress in the previous time step) - (the stress in this time step)
			dcs += csparam.CS1alpha * std::abs(sparam.s_lares_n - sparam.s_lares_n_prev) / 2; // equations (22)+ (27) from the paper
			ism += s.get_CS() / ASn;														  // increase SEI growth proportionally the crack surface
																							  // current density on particle = I /(elec_surf * thick * a)
																							  // 		isei also acts on this scale since it is an extra boundary condition ( itot = (jn + isei) =  (surface gradient)/nF )
																							  // 		crack growth -> increase isei_on_particle -> (jn + isei*(initial+crack_surface)/initial)*nF
																							  // 			such that if the crack surface area is the same as the initial electrode surface area, we double isei
			break;
		case 2: // Laresgoiti's crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
				// but with stress model from Dai, Cai, White, Journal of Power sources 247, 2014
				// instead of Laresgoiti's stress from figure 5
				// this model calculates crack growth due to spatial (Dai) and temporal (Laresgoiti) variations in the li-concentration
			// ensure the stress values are up to date
			if (!sparam.s_dai_update)
			{
				std::cerr << "ERROR in Cell::CS. The stress values for Dai's stress model are not updated. Throwing an error.\n";
				// if you see this error, you have to call Cell::updateDaiStress(), which calculates the stress and stores the values before you call this function
				throw 108;
			}

			// Add the effects of this model
			dcs += csparam.CS2alpha * std::abs(sparam.s_dai_n - sparam.s_dai_n_prev) / 2; // equations (22)+ (27) from the paper but with Dai's stress
			ism += s.get_CS() / ASn;													  // increase SEI growth proportionally the crack surface
			break;
		case 3: // model by Deshpande & Bernardi,Journal of the Electrochemical Society 164 (2), 2017
				// this model is adapted to calculate crack growth due to spatial variations in the li-concentration
			// get concentrations
			double cp[settings::nch + 2], cn[settings::nch + 2];
			getC(cp, cn);

			// Add the effects of this model
			dcs += csparam.CS3alpha * pow((cn[0] - cn[settings::nch + 1]) / Cmaxneg, 2.0); // equations (8) + (21)
																						   // Note that eqn (8) refers to the change with respect to time while here we use the spatial variation
																						   // This makes the model capture more of the spatial variation in stress
																						   // Laresgoiti's model already accounted for temporal variation, so simply use that if you are interested in temporal rather than spatial variation
			ism += s.get_CS() / ASn;													   // increase SEI growth proportionally the crack surface
			break;
		case 4:
		{
			// model from Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
			// equation (1a): CS = Amax(1 - exp(-m * Ah)) with m a fitting parameter and Ah the charge throughput up to now
			// 	dCS/dAh = m Amax exp(-m Ah) = m Amax - m Amax + m Amax exp(-m Ah) = m Amax - m CS = m(Amax - CS)
			// 	dCS/dt = dCS/dAh * dAh/dt = dCS/dAh * abs(I) = m(Amax - CS)*abs(I)
			// 		where we use the absolute value of I because 'Ah' is the total charge throughput, i.e. int ( abs(I) dt )

			const double Amax = std::max(csparam.CS4Amax, s.get_CS());
			// 'maximum crack surface area', a fitting parameters
			// avoid negative crack growth if the crack surface becomes slightly larger than Amax
			// this is possible due to discrete time steps: CS(t) is just smaller, but CS (t+1) = CS(t) + dCS*dt is just larger

			// Add the effects of this model
			dcs += csparam.CS4alpha * (Amax - s.get_CS()) * std::abs(Icell); // see above, with m = csparam.CS4
			ism += s.get_CS() / ASn;										 // increase SEI growth proportionally the crack surface
		}
		break;
		case 5:
		{
			// model from Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
			// overpotential for the crack-side-reaction = overpotential for the SEI reaction
			const double etasei = (OCVnt + etan - OCVsei + Rsei * s.get_delta() * Icell); // overpotential [V], equation (6)

			// get surface concentration
			double cps, cns;
			getCSurf(&cps, &cns); // get the surface lithium concentration

			double kcr; // rate constant for the side reaction
			// Calculate the rate constant, equation (11) with an Arrhenius relation for the temperature (which wasn't considered by Ekstrom)
			if (Icell > 0)
				kcr = 0;
			else if (cns / Cmaxneg < 0.3)
				kcr = 2 * csparam.CS5k * std::exp(csparam.CS5k_T / Rg * (1 / T_ref - 1 / s.get_T()));
			else if (cns / Cmaxneg < 0.7)
				kcr = 0;
			else
				kcr = csparam.CS5k * std::exp(csparam.CS5k_T / Rg * (1 / T_ref - 1 / s.get_T()));

			// Add the effects of this model
			dcs += nsei * F * kcr * std::exp(-alphasei * nsei * F / (Rg * s.get_T()) * etasei); // equation (9)
			ism += s.get_CS() / ASn;															// increase SEI growth proportionally the crack surface
		}
		break;
		default: // unknown degradation model
			std::cerr << "ERROR in Cell::CS, unknown crack growth model with identifier " << deg_id.CS_id[i] << ". Only values 0 to 5 are allowed. Throw an error.\n";
			throw 106;
			break;
		}
	}

	// Make the output variables
	*isei_multiplyer = ism;
	*dCS = dcs;

	// Decrease the negative diffusion constant if needed
	if (deg_id.CS_diffusion == 0) // don't decrease the negative diffusion constant
		*dDn = 0;
	else if (deg_id.CS_diffusion == 1)
	{ // decrease it according to Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
		// equation (2) D(t) = D0 (1 - CS)^gamma
		// 	but this can become negative is CS is larger than 1, we assume there should be a term /Amax in eqn (2): D(t) = D0 (1 - (CS/Amax))^gamma, which becomes 0 if CS = Amax
		// so dD/dt = - gamma D0 (1-CS/Amax)^(gamma-1) 1/Amax dCS/dt

		// 'maximum crack surface area', a fitting parameters
		// avoid increasing diffusion coefficient if the crack surface becomes larger than Amax
		// this is possible if the user chooses a different CS growth model, which does give larger crack surfaces (i.e. not CS4 which is Barai's crack growth model)
		const double Amax = std::max(csparam.CS4Amax, s.get_CS());

		// cap the decrease rate at a maximum value of 2e-7 (2e-5% = kill the battery in  about 1000h)
		const double Dnmax = std::min(2e-7, csparam.CS_diffusion * pow(1 - s.get_CS() / Amax, csparam.CS_diffusion - 1) / Amax * dcs);
		*dDn = -Dnmax * s.get_Dn();
	}
	else
	{ // unknown degradation model
		std::cerr << "ERROR in Cell::CS, unknown value for decreasing the diffusion constant " << deg_id.CS_diffusion << ". Only values 0 or 1 are allowed. Throw an error.\n";
		throw 106;
	}

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::CS terminating\n";
}

void Cell::LAM(bool print, double zp_surf, double etap,
			   double *dthickp, double *dthickn, double *dap, double *dan, double *dep, double *den)
{
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

	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::LAM starting.\n";

	// output parameters
	double dthickpp = 0;
	double dthicknn = 0;
	double dapp = 0;
	double dann = 0;
	double depp = 0;
	double denn = 0;

	if (deg_id.LAM_n > deg_id.len)
	{
		std::cerr << "ERROR in Cell::LAM the user wants to use more than " << deg_id.len << " degradation models. throwing an error.\n";
		std::cout << "The maximum length of the array with the degradation identifiers is " << deg_id.len << ". If you want to use more models,"
				  << " you have to increase the value of 'len' in the struct DEG_ID, defined in Cell.hpp.\n";
		throw 107;
	}

	// loop for each model to use
	for (int i = 0; i < deg_id.LAM_n; i++)
	{

		// calculate the effect of this model
		switch (deg_id.LAM_id[i])
		{
		case 0: // no LAM
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
			depp += 0;
			denn += 0;
			break;
		case 1: // Stress model from Dai, Cai, White, Journal of Power sources 247, 2014
				// LAM equation similar to CS equation from Laresgoiti
				// (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
			// ensure the stress values are up to date
			if (!sparam.s_dai_update)
			{
				std::cerr << "ERROR in Cell::LAM. The stress values for Dai's stress model are not updated. Throwing an error.\n";
				// if you see this error, you have to call Cell::updateDaiStress(), which calculates the stress and stores the values before calling this function
				throw 108;
			}

			// Laresgoiti's equation to link stress to LAM
			dthickpp += -lamparam.lam1p * std::abs(sparam.s_dai_p - sparam.s_dai_p_prev) / 2; // with Dai's stress model (values stored in s_dai_p)
			dthicknn += -lamparam.lam1n * std::abs(sparam.s_dai_n - sparam.s_dai_n_prev) / 2;
			// assume the other effects are 0
			dapp += 0;
			dann += 0;
			depp += 0;
			denn += 0;
			break;
		case 2: // Model by Delacourt & Safari, Journal of the Electrochemical Society 159 (8), 2012
		{		// Get the molar flux on each particle

			const double i_app = Icell / elec_surf;							  // current density on the electrode [A m-2]
			const double jp = -i_app / (s.get_ap() * n * F * s.get_thickp()); // molar flux on the positive particle [mol m-2 s-1]
			const double jn = i_app / (s.get_an() * n * F * s.get_thickn());  // molar flux on the negative particle [mol m-2 s-1]

			// Use Arrhenius relations to update the fitting parameters for the cell temperature

			const double ap = lamparam.lam2ap * std::exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / s.get_T()));
			const double an = lamparam.lam2an * std::exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / s.get_T()));
			const double bp = lamparam.lam2bp * std::exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / s.get_T()));
			const double bn = lamparam.lam2bn * std::exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / s.get_T()));

			// Add the effects of this model
			depp += ap * std::abs(jp) + bp * sqrt(std::abs(jp)); // equation (5) from the paper
			denn += an * std::abs(jn) + bn * sqrt(std::abs(jn));
			// assume the other effects are 0
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
		}
		break;
		case 3: // Model by Kindermann, Keil, Frank, Jossen, Journal of the Electrochemical Society 164 (12), 2017
		{
			double OCVpt; // cathode potential
			try
			{
				OCVpt = OCV_curves.linInt_OCV_pos(zp_surf, print);
				// get OCV of positive electrode, throw error if out of bounds
				// this should be updated for the cell's temperature using the entropic coefficient of the cathode
				// but I couldn't find any data on this, so I have ignored the effect
			}
			catch (int e)
			{
				std::cout << "Error in Cell::LAM when calculating the cathode potential for LAM: " << e << ". Throwing it on.\n";
				throw e;
			}

			// overpotential for the NMC dissolution reaction
			const double etap_LAM = OCVpt + etap - OCVnmc; // equation (9) from the paper

			// temperature dependent rate constant
			const double kt = lamparam.lam3k * std::exp(lamparam.lam3k_T / Rg * (1 / T_ref - 1 / s.get_T())); // Arrhenius law

			// current density of the NMC dissolution reaction
			const double idiss = std::max(-5e-6, -kt * std::exp(n * F / Rg / s.get_T() * etap_LAM) / (n * F)); // equation (8) from the paper
			// cap the effect at 5e-6 to avoid a very fast drop of capacity (which could  cause an error)
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
		}
		break;
		case 4: // Model by Narayanrao, Joglekar, Inguva, Journal of the Electrochemical Society 160 (1), 2012
			// Add the effects of this model
			dapp += -lamparam.lam4p * s.get_ap(); // equation (7) from the paper
			dann += -lamparam.lam4n * s.get_an();
			// assume the other effects are 0
			dthickpp += 0;
			dthicknn += 0;
			depp += 0;
			denn += 0;
			break;
		default: // unknown degradation model
			std::cerr << "ERROR in Cell::LAM, unknown LAM degradation model with identifier " << deg_id.LAM_id[i] << ". Only values 0 to 4 are allowed. Throw an error.\n";
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

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::LAM terminating\n";
}

void Cell::LiPlating(double OCVnt, double etan, double *ipl)
{
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
	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::LiPlating starting\n";

	// Arrhenius relation for temperature-dependent plating parameters
	const double kplt = plparam.pl1k * std::exp(plparam.pl1k_T / Rg * (1 / T_ref - 1 / s.get_T())); // Rate constant

	if (deg_id.pl_id == 0) // no plating
		*ipl = 0;
	else if (deg_id.pl_id == 1) // Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
		*ipl = npl * F * kplt * std::exp(-n * F / (Rg * s.get_T()) * alphapl * (OCVnt + etan - OCVpl + Rsei * s.get_delta() * Icell));
	else
	{
		std::cerr << "ERROR in Cell::LiPlating, illegal degradation model identifier " << deg_id.pl_id << ", only values 0 and 1 are allowed. Throwing an error.\n";
		throw 106;
	}

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::LiPlating starting\n";
}

// state space model
slide::states_type Cell::dState(bool print, bool blockDegradation, int electr)
{
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

	using namespace PhyConst;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::dState starting\n";

	if ((electr == 1 || electr == 2) && !blockDegradation)
	{
		std::cerr << "ERROR in Cell::dState. you are cycling with only one electrode " << electr
				  << " but you are also accounting for degradation."
					 " That is not allowed, half-cell cycling is only allowed if no degradation is considered. Either change the value of 'electr' to 3"
					 " such that you cycle with the full cell, or set 'blockDegradation' to true to ignore degradation.\n";
		throw 109;
	}

	// get surface concentrations
	double cps, cns;

	getCSurf(&cps, &cns); // get the surface lithium concentration

	// check if the surface concentration is within the allowed range
	// 	0 < cp < Cmaxpos
	// 	0 < cn < Cmaxneg
	if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg)
	{
		if (print)
		{
			std::cerr << "ERROR in Cell::dState: concentration out of bounds. the positive lithium fraction is " << cps / Cmaxpos << " and the negative lithium fraction is " << cns / Cmaxneg;
			std::cerr << "they should both be between 0 and 1.\n";
		}
		throw 101;
	}

	const double zp_surf = (cps / Cmaxpos); // lithium fraction (0 to 1)
	const double zn_surf = (cns / Cmaxneg);

	// (electr == 1)? only consider positive electrode, ignore the negative electrode : only consider negative electrode, ignore the positive electrode

	// current density
	const double i_app = Icell / elec_surf;												  // current density on the electrode surfaces [A m-2]
	const double jp = (electr == 2) ? 0 : -i_app / (s.get_ap() * n * F * s.get_thickp()); // molar flux on the positive particle [mol m-2 s-1]
	const double jn = (electr == 1) ? 0 : i_app / (s.get_an() * n * F * s.get_thickn());  // molar flux on the negative particle [mol m-2 s-1]

	// Arrhenius relation for temperature-dependent parameters
	const double Dpt = (electr == 2) ? 0 : s.get_Dp() * std::exp(Dp_T / Rg * (1 / T_ref - 1 / s.get_T())); // Diffusion constant at the positive electrode at the cell's temperature [m s-1]
	const double Dnt = (electr == 1) ? 0 : s.get_Dn() * std::exp(Dn_T / Rg * (1 / T_ref - 1 / s.get_T())); // Diffusion constant at the negative electrode at the cell's temperature [m s-1]
	const double kpt = (electr == 2) ? 0 : kp * std::exp(kp_T / Rg * (1 / T_ref - 1 / s.get_T()));		   // Rate constant at the positive electrode at the cell's temperature [m s-1]
	const double knt = (electr == 1) ? 0 : kn * std::exp(kn_T / Rg * (1 / T_ref - 1 / s.get_T()));		   // Rate constant at the negative electrode at the cell's temperature [m s-1]

	// Calculate the effect of the main li-reaction on the (transformed) concentration
	slide::z_type dzp, dzn;
	for (int j = 0; j < settings::nch; j++)
	{											   // loop for each row of the matrix-vector product A * z
		const double ctep = M.Ap[j] * s.get_zp(j); // A is diagonal, so the array M.A has only the diagonal elements
		const double cten = M.An[j] * s.get_zn(j);
		dzp[j] = (Dpt * ctep + M.Bp[j] * jp); // dz/dt = D * A * z + B * j
		dzn[j] = (Dnt * cten + M.Bn[j] * jn);
	}

	// Calculate the overpotential using the Bulter-Volmer equation
	// if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
	// and asinh(x) = ln(x + sqrt(1+x^2)
	const double i0p = kpt * n * F * sqrt(C_elec * cps * (Cmaxpos - cps));								  // exchange current density of the positive electrode
	const double i0n = knt * n * F * sqrt(C_elec * cns * (Cmaxneg - cns));								  // exchange current density of the negative electrode
	const double xp = -0.5 * i_app / (s.get_ap() * s.get_thickp()) / i0p;								  // x for the cathode
	const double xn = 0.5 * i_app / (s.get_an() * s.get_thickn()) / i0n;								  // x for the anode
	const double etap = (electr == 2) ? 0 : (2 * Rg * s.get_T()) / (n * F) * log(xp + sqrt(1 + xp * xp)); // cathode overpotential [V], < 0 on discharge
	const double etan = (electr == 1) ? 0 : (2 * Rg * s.get_T()) / (n * F) * log(xn + sqrt(1 + xn * xn)); // anode overpotential [V],  > 0 on discharge

	// Calculate the entropic coefficient
	double dOCV;
	bool bound = true; // in linear interpolation, throw an error if you are outside of the allowed range of the data
	try
	{
		dOCV = OCV_curves.linInt_dOCV_tot(zp_surf, print, bound); // entropic coefficient of the entire cell OCV [V K-1]
	}
	catch (int e)
	{
		if (print)
			std::cout << "Error in Cell::dState when calculating the entropic coefficient " << e << ". Throwing it on.\n";
		throw e;
	}

	// temperature
	// Calculate the thermal sources/sinks/transfers per unit of volume of the battery
	// The battery volume is given by the product of the cell thickness and the electrode surface
	const double Qrev = -i_app / L * s.get_T() * dOCV;			  // reversible heat due to entropy changes [W m-3]
	const double Qrea = i_app / L * (etan - etap);				  // reaction heat due to the kinetics [W m-3]
	const double Qohm = Icell * Icell * getR() / (L * elec_surf); // Ohmic heat due to electrode resistance [W m-3]
	const double Qc = -Qch * SAV * (s.get_T() - T_env);			  // cooling with the environment [W m-3]

	// If we ignore degradation in this time step, we have calculated everything we need
	if (blockDegradation)
	{
		slide::states_type dstates{};												 // Initialize as zero.
		std::copy(dzp.begin(), dzp.end(), dstates.begin());							 // first nch dstates are d_zp,
		std::copy(dzn.begin(), dzn.end(), dstates.begin() + settings::nch);			 // first nch dstates are d_zn,
		dstates[2 * settings::nch + 0] = 1 / (rho * Cp) * (Qrev + Qrea + Qohm + Qc); // dT		cell temperature

		// Others are zero : ddelta	SEI thickness, dLLI	lost lithium, dthickp/dthickn 	electrode thickness,
		// dep/den volume fraction of active material, dap/dan effective surface are, a = 3 e/R, dCS surface area of the cracks,
		// dDp/dDn diffusion constant, dR electrode resistance, ddelta_pl thickness of the plated lithium

		if constexpr (settings::verbose >= printLevel::printCellFunctions)
			std::cout << "Cell::dState terminating without degradation.\n";

		return dstates; // stop calculating
	}

	// calculate the anode potential (needed for various degradation models)
	double OCV_n, dOCVn; // anode potential at reference temperature and entropic coefficient
	try
	{
		dOCVn = OCV_curves.linInt_dOCV_neg(zn_surf, print, bound); // entropic coefficient of the anode potential [V K-1]
		OCV_n = OCV_curves.linInt_OCV_neg(zn_surf, print, bound);  // anode potential [V]
	}
	catch (int e)
	{
		if (print)
			std::cout << "Error in Cell::dState when calculating the anode potential " << e << ". Throwing it on.\n";

		throw e;
	}
	const double OCVnt = OCV_n + (s.get_T() - T_ref) * dOCVn; // anode potential at the cell's temperature [V]

	// SEI growth
	double isei;				  // current density of the SEI growth side reaction [A m-2]
	double den_sei;				  // decrease in volume fraction due to SEI growth [s-1]
	double dznsei[settings::nch]; // additional diffusion in the anode due to isei
	try
	{
		SEI(OCVnt, etan, &isei, &den_sei);
	}
	catch (int e)
	{
		if (print)
			std::cout << "Error in Cell::dState when calculating the effect of SEI growth: " << e << ". Throwing it on.\n";
		throw e;
	}

	// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF)
	for (int j = 0; j < settings::nch; j++)
		dznsei[j] = (M.Bn[j] * isei / (nsei * F));

	// crack growth leading to additional exposed surface area
	double isei_multiplyer;			 // relative increase in isei due to additional SEI growth on the extra exposed surface area [-]
	double dCS;						 // increase in crack surface area [m2 s-1]
	double dDn;						 // change in negative diffusion constant [m s-1 s-1]
	double dznsei_CS[settings::nch]; // additional diffusion in the anode due to extra SEI growth on the crack surface
	try
	{
		CS(OCVnt, etan, &isei_multiplyer, &dCS, &dDn);
	}
	catch (int e)
	{
		if (print)
			std::cout << "Error in Cell::dState when calculating the effect of crack growth: " << e << ". Throwing it on.\n";
		throw e;
	}

	// crack surface leads to extra SEI growth because the exposed surface area increases.
	// (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF + isei_CS/nF)
	const double isei_CS = isei * isei_multiplyer; // extra SEI side reaction current density due to the crack surface area [A m-2]
	for (int j = 0; j < settings::nch; j++)
		dznsei_CS[j] = (M.Bn[j] * isei_CS / (nsei * F));

	// loss of active material LAM
	double dthickp, dthickn, dap, dan, dep, den; // change in geometric parameters describing the amount of active material
	try
	{
		LAM(print, zp_surf, etap, &dthickp, &dthickn, &dap, &dan, &dep, &den);
	}
	catch (int e)
	{
		if (print)
			std::cout << "Error in Cell::dState when calculating the LAM: " << e << ". Throwing it on.\n";
		throw e;
	}

	// lithium plating
	double ipl;					  // current density of the plating side reaction [A m-2]
	double dzn_pl[settings::nch]; // additional diffusion in the anode due to ipl
	try
	{
		LiPlating(OCVnt, etan, &ipl);
	}
	catch (int e)
	{
		if (print)
			std::cout << "Error in Cell::dState when calculating the lithium plating: " << e << ". Throwing it on.\n";
		throw e;
	}

	// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + ipl/nF)
	for (int j = 0; j < settings::nch; j++)
		dzn_pl[j] = (M.Bn[j] * ipl / (npl * F));

	// time derivatives
	slide::states_type dstates{}; // Initialize as zero.
	for (int j = 0; j < settings::nch; j++)
	{
		dstates[j] = dzp[j];														  // dzp 		diffusion
		dstates[settings::nch + j] = (dzn[j] + dznsei[j] + dznsei_CS[j] + dzn_pl[j]); // dzn		jtot = jn + isei/nF + isei_CS/nF + ipl/nF
	}
	dstates[2 * settings::nch + 0] = 1 / (rho * Cp) * (Qrev + Qrea + Qohm + Qc);					   // dT 		cell temperature
	dstates[2 * settings::nch + 1] = isei / (nsei * F * rhosei);									   // ddelta	thickness of the SEI layer
																									   // delta uses only isei (and not isei + isei_CS) since crack growth increases the area, not the thickness
	dstates[2 * settings::nch + 2] = (isei + isei_CS + ipl) * elec_surf * s.get_thickn() * s.get_an(); // dLLI 	loss of lithium
																									   // i_sei = density => * active surface area = * (surf*thick*specific_surf_neg)
	dstates[2 * settings::nch + 3] = dthickp;														   // dthickp 	electrode thickness
	dstates[2 * settings::nch + 4] = dthickn;														   // dthickn
	dstates[2 * settings::nch + 5] = dep;															   // dep		volume fraction of active material
	dstates[2 * settings::nch + 6] = den + den_sei;													   // den
	dstates[2 * settings::nch + 7] = dap + 3 / Rp * dep;											   // dap		effective surface area, a = 3 e/R -> da/dt = da/dt + 3/R de/dt
	dstates[2 * settings::nch + 8] = dan + 3 / Rn * (den + den_sei);								   // dan
	dstates[2 * settings::nch + 9] = dCS;															   // dCS 		surface area of the cracks
	dstates[2 * settings::nch + 10] = 0;															   // dDp 		diffusion constant
	dstates[2 * settings::nch + 11] = dDn;															   // dDn
	dstates[2 * settings::nch + 12] = 0;															   // dR 		specific electrode resistance
	dstates[2 * settings::nch + 13] = ipl / (npl * F * rhopl);										   // ddelta_pl thickness of the plated lithium

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::dState terminating with degradation.\n";
	return dstates;
}

void validState(slide::State &s, slide::State &s_ini)
{
	/* Moved from State to Cell. 
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
	if (s.get_T() < settings::Tmin_K) // check the temperature is above 0 degrees
	{
		std::cerr << "Error in State::setT. The temperature " << s.get_T() << "K is too low. The minimum value is " << settings::Tmin_K << std::endl;
		throw 15;
	}
	else if (s.get_T() > settings::Tmax_K) // check the temperature is below 60 degrees
	{
		std::cerr << "Error in State::setT. The temperature " << s.get_T() << "K is too high. The maximum value is " << settings::Tmax_K << std::endl;
		throw 15;
	}

	// thickness of the SEI layer
	const bool del = s.get_delta() <= 0;
	if (del)
	{
		std::cerr << "Error in State::validState. The SEI thickness delta is " << s.get_delta() << ", which is too low. Only strictly positive values are allowed"
				  << '\n';
		throw 15;
	}
	// a value of 0 gives problems in some equations, which have a term 1/delta, which would give nan or inf if delta = 0
	// a negative value might lead to a further decrease in SEI thickness, so it will keep getting more and more negative

	// lost lithium
	bool li = s.get_LLI() < 0;
	if (li)
	{
		std::cerr << "Error in State::validState. The lost lithium LLI is " << s.get_LLI() << ", which is too low. Only non-negative values are allowed"
				  << '\n';
		throw 15;
	}
	// a value of 0 is allowed (it won't lead to errors in the code)

	// thickness of the electrodes
	bool tpmin = s.get_thickp() <= s_ini.get_thickp() / 5;
	if (tpmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the thickness of the positive electrode is " << s.get_thickp() << ", which is too low. The minimum is " << s_ini.get_thickp() / 5 << ", 1/5 of the original thickness"
				  << '\n';
		throw 15;
	}
	// errors will happen if the thickness becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool tpmax = s.get_thickp() > s_ini.get_thickp() * 1.001; // leave a 0.1% margin for numerical errors
	if (tpmax)
	{
		std::cerr << "Error in State::validState. The thickness of the positive electrode is " << s.get_thickp()
				  << ", which is too high. The maximum is " << s_ini.get_thickp() << ", the original thickness\n";
		throw 15;
	}
	bool tnmin = s.get_thickn() <= s_ini.get_thickn() / 5;
	if (tnmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the thickness of the negative electrode is " << s.get_thickn()
				  << ", which is too low. The minimum is " << s_ini.get_thickn() / 5 << ", 1/5 of the original thickness\n";
		throw 15;
	}

	// errors will happen if the thickness becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool tnmax = s.get_thickn() > s_ini.get_thickn() * 1.001;
	if (tnmax)
	{
		std::cerr << "Error in State::validState. The thickness of the negative electrode is " << s.get_thickn()
				  << ", which is too high. The maximum is " << s_ini.get_thickn() << ", the original thickness"
				  << '\n';

		throw 15;
	}

	// volume fraction of the active material in the electrodes
	bool epmin = s.get_ep() <= s_ini.get_ep() / 5;
	if (epmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the volume fraction of the positive electrode is "
				  << s.get_ep() << ", which is too low. The minimum is " << s_ini.get_ep() / 5 << ", 1/5 of the original volume fraction.\n";
		throw 15;
	}

	// errors will happen if the volume fraction becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool epmax = s.get_ep() > s_ini.get_ep() * 1.001;
	if (epmax)
	{
		std::cerr << "Error in State::validState. The volume fraction of the positive electrode is " << s.get_ep()
				  << ", which is too high. The maximum is " << s_ini.get_ep() << ", the original volume fraction.\n";
		throw 15;
	}

	bool enmin = s.get_en() <= s_ini.get_en() / 5;
	if (enmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the volume fraction of the negative electrode is " << s.get_en()
				  << ", which is too low. The minimum is " << s_ini.get_en() / 5 << ", 1/5 of the original volume fraction.\n";
		throw 15;
	}

	// errors will happen if the volume fraction becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool enmax = s.get_en() > s_ini.get_en() * 1.001;
	if (enmax)
	{
		std::cerr << "Error in State::validState. The volume fraction of the negative electrode is " << s.get_en()
				  << ", which is too high. The maximum is " << s_ini.get_en() << ", the original volume fraction.\n";
		throw 15;
	}

	// effective surface area of the porous electrodes
	bool apmin = s.get_ap() <= s_ini.get_ap() / 5;
	if (apmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the effective surface of the positive electrode is " << s.get_ap()
				  << ", which is too low. The minimum is " << s_ini.get_ap() / 5 << ", 1/5 of the original effective surface.\n";
		throw 15;
	}
	// errors will happen if the effective surface becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool apmax = s.get_ap() > s_ini.get_ap() * 1.001;
	if (apmax)
	{
		std::cerr << "Error in State::validState. The effective surface of the positive electrode is " << s.get_ap()
				  << ", which is too high. The maximum is " << s_ini.get_ap() << ", the original effective surface.\n";
		throw 15;
	}

	bool anmin = s.get_an() <= s_ini.get_an() / 5;
	if (anmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the effective surface of the negative electrode is " << s.get_an()
				  << ", which is too low. The minimum is " << s_ini.get_an() / 5 << ", 1/5 of the original effective surface.\n";
		throw 15;
	}

	// errors will happen if the effective surface becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too little active material -> too large concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool anmax = s.get_an() > s_ini.get_an() * 1.001;
	if (anmax)
	{
		std::cerr << "Error in State::validState. The effective surface of the negative electrode is " << s.get_an()
				  << ", which is too high. The maximum is " << s_ini.get_an() << ", the original effective surface.\n";
		throw 15;
	}

	// surface area of the cracks growing at the particle surface
	bool csmin = s.get_CS() <= 0;
	if (csmin)
	{
		std::cerr << "Error in State::validState. The crack surface area is " << s.get_CS()
				  << ", which is too low. It must be strictly positive.\n";
		throw 15;
	}

	// don't allow 0 because it might give problems in some of the equations, which might grow CS proportional to the existing CS (so a value of 0 gives 0 increase)
	// a negative value might lead to a further decrease in CS, so it will keep getting more and more negative
	bool csmax = s.get_CS() > 1e4 * s_ini.get_CS();
	if (csmax)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the crack surface area is " << s.get_CS()
				  << ", which is too high. The maximum is " << 1e4 * s_ini.get_CS() << ", 10,000 times the original crack surface area.\n";
		throw 15;
	}
	// normally, the initial value is 1% of the total real electrode surface area, so 10,000*initial value = 100 * total electrode surface area
	// but in theory no maximum value will give errors in the code

	// diffusion constant at reference temperature for the electrodes
	bool dpmin = s.get_Dp() <= s_ini.get_Dp() / 5;
	if (dpmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the diffusion constant of the positive electrode is " << s.get_Dp()
				  << ", which is too low. The minimum is " << s_ini.get_Dp() / 5 << ", 1/5 of the original diffusion constant.\n";
		throw 15;
	}

	// errors will happen if the diffusion constant becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too bad diffusion -> too large surface concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool dpmax = s.get_Dp() > s_ini.get_Dp() * 1.001;
	if (dpmax)
	{
		std::cerr << "Error in State::validState. The diffusion constant of the positive electrode is " << s.get_Dp()
				  << ", which is too high. The maximum is " << s_ini.get_Dp() << ", the original diffusion constant.\n";
		throw 15;
	}

	bool dnmin = s.get_Dn() <= s_ini.get_Dn() / 5;
	if (dnmin)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the diffusion constant of the negative electrode is " << s.get_Dn()
				  << ", which is too low. The minimum is " << s_ini.get_Dn() / 5 << ", 1/5 of the original diffusion constant.\n";
		throw 15;
	}
	// errors will happen if the diffusion constant becomes 0 or negative. otherwise no direct problems are expected
	// on a longer term, you will get discretisation errors (see above, too bad diffusion -> too large surface concentration swings in one time step -> numerical errors / concentration out of bound, etc
	bool dnmax = s.get_Dn() > s_ini.get_Dn() * 1.001;
	if (dnmax)
	{
		std::cerr << "Error in State::validState. The diffusion constant of the negative electrode is " << s.get_Dn() << ", which is too high. The maximum is "
				  << s_ini.get_Dn() << ", the original effective diffusion constant.\n";
		throw 15;
	}

	// specific resistance of the electrodes (one value for both electrodes combined)
	bool rmin = s.get_r() <= 0;
	if (rmin)
	{
		std::cerr << "Error in State::validState. The specific resistance is " << s.get_r() << ", which is too low, it must be strictly positive.\n";
		throw 15;
	}

	bool rmax = s.get_r() > 1000 * s_ini.get_r();
	if (rmax)
	{
		std::cerr << "Error in State::validState. The cell has degraded too much and the specific resistance is " << s.get_r()
				  << ", which is too high. The maximum is " << 1000 * s_ini.get_r() << ", 1000 times the original specific resistance.\n";
		throw 15;
	}

	// thickness of the plated litium layer
	bool delpl = s.get_delta_pl() < 0;
	if (delpl)
	{
		std::cerr << "Error in State::validState. The thickness of the plated lithium is " << s.get_delta_pl()
				  << ", which is too low. Only strictly positive values are allowed.\n";
		throw 15;

		// 0 is allowed
		// a negative value might lead to a further decrease in thickness, so it will keep getting more and more negative
	}

	// throw an error if one of the states was invalid
}

void Cell::ETI(bool print, double dti, bool blockDegradation)
{
	/*
	 * Performs forward Euler time integration over one time step of dti seconds
	 * s(t+1) = s(t) + ds/dt * dti
	 *
	 * IN
	 * print 			boolean indicating if we want to print error messages or not
	 * 					if true, error messages are printed
	 * 					if false no error messages are printed (but the error will still be thrown)
	 * 					we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * dti 				time step over which time integradation should be done [s]
	 * 					it should be small enough to ensure numerical stability and accuracy
	 * 					1-5 seconds in usually ok (depending on the magnitude of the current, the temperature, etc.)
	 * blockDegradation if true, degradation is not accounted for in this time step
	 */

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::ETI starting.\n";

	// Update the stress values stored in the attributes with the stress of the previous time step
	sparam.s_dai_p_prev = sparam.s_dai_p;	  // Dai's stress in the positive particle in the previous time step
	sparam.s_dai_n_prev = sparam.s_dai_n;	  // Dai's stress in the negative particle in the previous time step
	sparam.s_lares_n_prev = sparam.s_lares_n; // Laresgoiti's stress in the negative particle in the previous time step

	// Calculate the stress values stored in the attributes for the stress in this time step
	if (sparam.s_dai) // only a few degradation models actually need the stress according to Dai, so only calculate it if needed
		updateDaiStress();
	if (sparam.s_lares) // only a few degradation models need the stress according to Laresgoiti
		updateLaresgoitiStress(print);

	// Calculate the time derivatives
	slide::states_type states = s.getStates_arr();

	// calculate time derivatives, electr = 0 to account for both electrodes (i.e. cycle the full cell)
	const slide::states_type dstates = dState(print, blockDegradation, 0); // arrays with the state and dstate/dt

	// forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
	for (int i = 0; i < settings::ns; i++) // loop for all states
		states[i] += dti * dstates[i];

	setStates(std::move(states)); // store new states, checks if they are illegal (throws an error in that case)

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	sparam.s_dai_update = false;
	sparam.s_lares_update = false;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::ETI terminating.\n";
}

void Cell::ETI_electr(bool print, double I, double dti, bool blockDegradation, bool pos)
{
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

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::ETI_electr starting\n";

	if (!blockDegradation)
	{
		std::cerr << "ERROR in Cell::ETI_electr, you are cycling only one electrode but want to account for degradation. This is not allowed.\n";
		// half-cell cycling is only allowed if blockDegradation is true (i.e. you ignore degradation)
		// this is because some of the degradation mechanisms could give NaN or Inf due to a divide by 0
		// So you can only call this function with 'true' for blockDegradation
		throw 109;
	}

	// Set the specified current
	Icell = I; // don't call setI because that function ramps the current on both electrodes
			   // so here 'cheat it' and directly set the current
			   // this means you avoid the checks done in setI, so you don't know if the current is feasible or not

	// Calculate the time derivatives
	slide::states_type states = s.getStates_arr();

	// calculate time derivatives of the positive or negative electrode
	const slide::states_type dstates = pos ? dState(print, blockDegradation, 1) : dState(print, blockDegradation, 2); // arrays with the state and dstate/dt

	// forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
	for (int i = 0; i < settings::ns; i++) // loop for all states
		states[i] += dti * dstates[i];

	s.setStates(std::move(states)); // store new states

	// the stress values stored in the class variables for stress are no longer valid because the state has changed
	sparam.s_dai_update = false;
	sparam.s_lares_update = false;

	if constexpr (settings::verbose >= printLevel::printCellFunctions)
		std::cout << "Cell::ETI_electr terminating.\n";
}

void Cell::checkModelparam()
{
	// check if the inputs to the Matlab code are the same as the ones here in the C++ code
	// input:
	// 		M.Input[0] has to be the same as nch (defined in State.hpp)
	// 		M.Input[1] has to be the same as Rp (defined earlier in this constructor)
	// 		M.Input[2] has to be the same as Rn (defined earlier in this constructor)
	// 		M.input[3] has to give the location of the 0 eigenvalue
	bool Mnch = (M.Input[0] - settings::nch) / M.Input[0] > 1e-10; // allow a relative difference of e-10 due to numerical errors
	if (Mnch)
		std::cerr << "ERROR in Cell_KokamNMC::Cell_KokamNMC: the value of nch used in the Matlab script " << M.Input[0]
				  << " is not the same as the value of nch used in the c++ code " << settings::nch << ".\n";

	bool Mrp = (M.Input[1] - Rp) / M.Input[1] > 1e-10; // allow a relative difference of e-10 due to numerical errors
	if (Mrp)
		std::cerr << "ERROR in Cell_KokamNMC::Cell_KokamNMC: the value of Rp used in the Matlab script " << M.Input[1]
				  << " is not the same as the value of Rp used in the c++ code " << Rp << ".\n";
	bool Mrn = (M.Input[2] - Rn) / M.Input[2] > 1e-10; // allow a relative difference of e-10 due to numerical errors
	if (Mrn)
		std::cerr << "ERROR in Cell_KokamNMC::Cell_KokamNMC: the value of Rn used in the Matlab script " << M.Input[2]
				  << " is not the same as the value of Rn used in the c++ code " << Rn << ".\n";
	int a = M.Input[3];
	bool Meig = std::abs(M.An[a]) > 1e-10 || std::abs(M.Ap[a]) > 1e-10; // allow a relative difference of e-10 due to numerical errors
	if (Meig)
		std::cerr << "ERROR in Cell_KokamNMC::Cell_KokamNMC: the row of the 0-eigenvalue is " << M.Input[3]
				  << " but that row has a positive eigenvalue of " << M.Ap[a] << " and negative eigenvalue of " << M.An[a] << ". They are not 0.\n";
	if (Mnch || Mrp || Mrn || Meig)
	{
		std::cout << "The Matlab script modelSetup.m produces matrices used by the C++ code for the discretisation of the solid diffusion equation."
					 " Matlab needs some input parameters, such as the number of nodes and the radius of the particles. These parameters are specified on top of the Matlab scripts."
					 " These values are also defined in the C++ code. Of course, both values have to be the same."
					 " It turned out this was not the case, so you either have to change the values in the Matlab script or the ones in the C++ code."
					 " We are throwing an error.\n";
		throw 110;
	}
}