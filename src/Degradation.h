/*
 * Degradation.h
 *
 * Header file for the degradation simulations.
 * The exact degradation procedures are defined in the Cycler.
 * In the functions defined here, the functions from the Cycler are called various times with slightly different parameters (e.g. different temperatures, C rates, etc).
 * As such, we can simulate the effect the different parameters have on battery degradation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */


#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <assert.h>
#include <memory>
#include <thread>

#include "Cell.hpp"

#ifndef SRC_DEGRADATION_H_
#define SRC_DEGRADATION_H_

using namespace std;

string print_DEG_ID(DEG_ID degid);

// Auxiliary functions for multi-threaded simulations
void Calendar_one(const struct Model& M, const struct DEG_ID& degid, int cellType, int verbose,														// simulate one calendar ageing experiment
		double V, double Ti, int Time, int mode, int timeCycleData, int timeCheck, struct checkUpProcedure proc, string name);
void Cycle_one(const struct Model& M, const struct DEG_ID& degid, int cellType, int verbose, double Vma, double Vmi,								// simulate one cycle ageing experiment
		double Ccha, bool CVcha, double Icutcha, double Cdis, bool CVdis, double Icutdis, double Ti, int timeCycleData, int nrCycles, int nrCap, struct checkUpProcedure proc, string name);
void Profile_one(const struct Model& M, const struct DEG_ID& degid, int cellType, int verbose, string profName, int n, int limit,					// simulate one drive cycle ageing experiment
		double Vma, double Vmi, double Ti, int timeCycleData, int nrProfiles, int nrCap, struct checkUpProcedure proc, string name);

// Degradation experiments
void CycleAgeing(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose);										// simulate a range of cycle ageing experiments (different temperatures, SoC windows, currents)
void CalendarAgeig(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose);										// simulate a range of calendar ageing experiments (different temperatures, SoC levels)
void ProfileAgeing(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose);										// simulate a range of drive cycle experiments (different cycles, different temperatures, etc.)

#endif /* SRC_DEGRADATION_H_ */
