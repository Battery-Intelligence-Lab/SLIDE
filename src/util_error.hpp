/*
 * Util_error.hpp
 *
 * Utility functions for error. So the code is less verbose. 
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

#include "cell.hpp"

class Cell;

namespace slide::util::error
{
    void checkInputParam_CalAge(Cell &c, double V, double Ti, int Time, int timeCheck, int mode);
    void checkInputParam_CC_V_CV_I(Cell &c, double Crate, double Vset, double Ccut);
    void checkInputParam_CycAge(Cell &c, double Vma, double Vmi, double Ccha, double Ccutcha,
                                double Cdis, double Ccutdis, double Ti, int nrCycles, int nrCap);

}