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

#ifndef SRC_DETERMINECHARACTERISATION_H_
#define SRC_DETERMINECHARACTERISATION_H_


#include <string>

using namespace std;

void CCCV(double Crate, double Ccut, double Tref, double Dp, double Dn, double kp, double kn, double R, struct OCVparam ocvfit,const struct Model& M,
		int nin, double Vsim[][2], double Tsim[][2], int* nsim);
void fitDiffusionAndRate(int hierarchy, int ir, double R, int nDstep, double Dpstep, double Dnstep, bool logDstep, double Dpmin, double Dnmin,
		int nkstep, double kpstep, double knstep, bool logkstep, double kpmin, double knmin,
		int nCCCV, double weights[], string names[], int lengths[],
		double Crates[], double Ccuts[], double Tref, struct OCVparam ocvfit,
		double* err, double par[], double parindex[]);
void oneLevelCharacterisationFit(int hierarchy, int nrstep, double rstep, double rmin,
		int nDstep, double Dpstep, double Dnstep, bool logDstep, double Dpmin, double Dnmin,
		int nkstep, double kpstep, double knstep, bool logkstep, double kpmin, double knmin,
		int nCCCV, double weights[], string names[], int lengths[],
		double Crates[], double Ccuts[], double Tref, struct OCVparam ocvfit,
		double* err, double par[], int parindex[]);
void hierarchicalCharacterisationFit(int hmax, int nrstep, double Rstep, double Rmin,
		int nDstep, double Dpstep, double Dnstep, bool logDstep, double Dpmin, double Dnmin,
		int nkstep, double kpstep, double knstep, bool logkstep, double kpmin, double knmin,
		int nCCCV, double weights[], string names[], int lengths[],
		double Crates[], double Ccuts[], double Tref, struct OCVparam ocvfit,
		double* err, double par[]);
void estimateCharacterisation();

// Define a struct with the parameters of the OCV curve
// These parameters are calculated by the functions in determineOCV.cpp
struct OCVparam{
	double elec_surf;		// electrode surface
	double ep;				// volume fraction of active material in the cathode
	double en;				// volume fraction of active material in the anode
	double thickp;			// thickness of the cathode
	double thickn;			// thickness of the anode

	string namepos;			// name of the CSV file with the cathode OCV curve
	string nameneg;			// name of the CSV file with the anode OCV curve
	int np;					// number of points in the cathode OCV curve
	int nn;					// number of points in the anode OCV curve

	double lifracpini;		// lithium fraction in the cathode at 50% soC
	double lifracnini;		// lithium fraction in the anode at 50% SoC
	double cmaxp;			// maximum lithium concentration in the cathode [mol m-3]
	double cmaxn;			// maximum lithium concentration in the anode [mol m-3]
	double cap;				// the capacity of the cell [Ah]
	double Vmax;			// maximum voltage of the cell [V]
	double Vmin; 			// minimum voltage of the cell [V]
};

#endif /* SRC_DETERMINECHARACTERISATION_H_ */
