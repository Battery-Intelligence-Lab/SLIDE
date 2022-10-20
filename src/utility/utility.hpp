/*
 * utility.hpp
 *
 * All utility functions. So the code is less verbose.
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "timing.hpp"
#include "parallelisation.hpp"
#include "util.hpp"
#include "interpolation.hpp"
#include "free_functions.hpp"
#include "predicate_functions.hpp"
#include "slide_aux.hpp"

#include "io/read_CSVfiles.hpp"
