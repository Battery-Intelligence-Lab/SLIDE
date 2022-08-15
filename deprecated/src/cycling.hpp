/*
 * Cycling.hpp
 *
 * Header file for the functions which simulate some cycling regimes for cell
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <memory>

#include "model.h"

void CCCV(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, const int verbose);    //!< function to load a cell with a few CCCV cycles
void FollowCurrent(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, int verbose); //!< function to let a cell follow a predefined current profile
