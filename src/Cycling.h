/*
 * Cycling.h
 *
 * Header file for the functions which simulate some cycling regimes for cell
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

#ifndef SRC_CYCLING_H_
#define SRC_CYCLING_H_

using namespace std;

void CCCV(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose);				// function to load a cell with a few CCCV cycles
void FollowCurrent(const struct Model& M, string pref, const struct DEG_ID& degid, int cellType, int verbose);		// function to let a cell follow a predefined current profile

#endif /* SRC_CYCLING_H_ */
