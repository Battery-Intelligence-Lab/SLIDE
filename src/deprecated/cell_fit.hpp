/*
 * CellFit.hpp
 *
 * One of the child classes of Cell.
 *
 * This cell is used to fit the cycling parameters of a cell to data.
 * It is done by the functions in determineCharacterisation.cpp
 * This cell should never be used by other functions (e.g. it should never be used for degradation simulations)
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "Cell.hpp"
#include "settings/settings.hpp"

class Cell_Fit : public Cell
{
public:
  Cell_Fit(const struct slide::Model_SPM &, int verbosei); // standard constructor

  void setOCVcurve(const std::string &namepos, const std::string &nameneg);                                         // sets the OCV curve of the cell to the given value
  void setInitialConcentration(double cmaxp, double cmaxn, double lifracp, double lifracn);                         // sets the initial concentration
  void setGeometricParameters(double capnom, double elec_surf, double ep, double en, double thickp, double thickn); // sets the geometric parameters related to the amount of active material
  void setRamping(double Istep, double tstep);                                                                      // sets the ramping parameters

  void setCharacterisationParam(double Dp, double Dn, double kp, double kn, double Rdc); // sets the parameters related to the characterisation of the cell
};
