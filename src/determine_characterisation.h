/*
 * determineRateParameters.h
 *
 * Header file for the functions used to find the parameters which will match measured CCCV cycles.
 * This are the diffusion constants, rate constants and DC resistance.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <string>

#include "cell_fit.hpp"
#include "slide_aux.hpp"
#include "model.h"
#include "cycler.hpp"

bool CCCV_fit(Cell_Fit c1, double Crate, double Ccut, double Tref, double Dp, double Dn, double kp, double kn, double R, const struct OCVparam &ocvfit, const struct slide::Model &M,
			  slide::vec_XYdata &Vsim, slide::vec_XYdata &Tsim);

void CCCV(double Crate, double Ccut, double Tref, double Dp, double Dn, double kp, double kn, double R, const struct OCVparam &ocvfit,
		  const struct slide::Model &M, slide::vec_XYdata &Vsim, slide::vec_XYdata &Tsim);

void fitDiffusionAndRate(int hierarchy, int ir, double R, slide::fixed_data<double> Dp_space, slide::fixed_data<double> Dn_space,
						 slide::fixed_data<double> kp_space, slide::fixed_data<double> kn_space,
						 std::vector<slide::vec_XYdata> &Vdata_all, double weights[],
						 double Crates[], double Ccuts[], double Tref, const struct OCVparam &ocvfit,
						 double *err, std::array<double, 5> &par, std::array<int, 5> &parindex);

void hierarchicalCharacterisationFit(int hmax, slide::fixed_data<double> r_space, slide::fixed_data<double> Dp_space,
									 slide::fixed_data<double> Dn_space, slide::fixed_data<double> kp_space,
									 slide::fixed_data<double> kn_space, std::vector<slide::vec_XYdata> &Vdata_all,
									 double weights[], double Crates[], double Ccuts[], double Tref,
									 const struct OCVparam &ocvfit, double *err, std::array<double, 5> &par);

void writeCharacterisationParam(int h, const std::array<double, 5> &par, double err);

void estimateCharacterisation();

// Define a struct with the parameters of the OCV curve
// These parameters are calculated by the functions in determineOCV.cpp
struct OCVparam
{
	double elec_surf; // electrode surface
	double ep;		  // volume fraction of active material in the cathode
	double en;		  // volume fraction of active material in the anode
	double thickp;	  // thickness of the cathode
	double thickn;	  // thickness of the anode

	std::string namepos; // name of the CSV file with the cathode OCV curve
	std::string nameneg; // name of the CSV file with the anode OCV curve
	int np;				 // number of points in the cathode OCV curve
	int nn;				 // number of points in the anode OCV curve

	double lifracpini; // lithium fraction in the cathode at 50% soC
	double lifracnini; // lithium fraction in the anode at 50% SoC
	double cmaxp;	   // maximum lithium concentration in the cathode [mol m-3]
	double cmaxn;	   // maximum lithium concentration in the anode [mol m-3]
	double cap;		   // the capacity of the cell [Ah]
	double Vmax;	   // maximum voltage of the cell [V]
	double Vmin;	   // minimum voltage of the cell [V]
};