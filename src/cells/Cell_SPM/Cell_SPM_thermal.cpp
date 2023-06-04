/*
 * Cell_SPM_thermal.cpp
 *
 * Implements the functions for the parent class of the Cells
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "Cell_SPM.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <utility>

namespace slide {
double Cell_SPM::getThermalSurface()
{
  /*
   * Assume heat transfer is on one side of the pouch (one of the faces).
   * We assume that on the other side is another cell, with which this one will also exchange heat
   *
   * So even though the total surface area = 2*Acell (+small edges on top, bottom and sides), Atherm is just Acell
   */
  return geo.Acell;
}

double Cell_SPM::thermalModel_cell()
{
  /*
   * Calculate the thermal model for this cell on its own.
   * We do not account for neighbouring cells, and this cell is cooled convectively by the environment
   */

  //!< total heat generation since last time the temperature was updated
  double Etot = Therm_Qgen;

  //!< cooling with the environment
  double Qc = Qch * getThermalSurface() * (T_env - T()) * Therm_time; //!< cooling with the environment [W m-3]
  Etot += Qc;

  //!< Calculate the new temperature
  //!< rho * cp * dT/dt = Qtot / V
  //!< 		where 	Qtot is total power in W
  //!< 				V is the cell's volume L * elec_surf
  //!< so integrated over time this is
  //!< rho * cp * (Tnew - Told) = Etot / V
  double Tnew = T() + Etot / (rho * Cp * geo.L * geo.elec_surf);

  //!< Check the new temperature is valid, and if so, set it
  if (Tnew > Tmax() || Tnew < Tmin() || std::isnan(Tnew)) {
    if constexpr (settings::printBool::printCrit) {
      std::cerr << "ERROR in Cell_SPM::thermalModel, the new temperature of " << Tnew
                << " is outside the allowed range from " << Tmin() << " to " << Tmax()
                << ". The time since the last time this function was called is " << Therm_time << '\n';

      std::cout << "Total thermal energy " << Etot << ". internal heat generation " << Therm_Qgen << '\n'
                << "giving change in temperature: " << Etot / (rho * Cp * geo.L * geo.elec_surf) << '\n';
    }
    throw 9;
  }

  //!< setting the temperature is done by the parent module. else some cells will update their T before others, and we get inconsistencies
  //!< 		e.g. exchange from cell 2 to this cell will not be same as from this cell to cell 2 since T of this cell would have changed.

  //!< Reset the cumulative thermal variables
  Therm_Qgentot += Therm_Qgen;
  Therm_Qgen = 0;
  Therm_time = 0;

  return Tnew;
}

double Cell_SPM::thermalModel_coupled(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim)
{
  /*
   * Calculate the thermal model for this cell and update its cell temperature.
   * The heat exchange in [W] with element i is given by
   * 		Qc = Kneighbours[i]*A*(Tneighbours[i] - getT());
   * We assume all temperatures were constant for the past period, such that the thermal energy in [J] is given by
   * 		Ec = Qc * Therm_time
   * This is added up with the internal heath generated, and the cell's temperature is returned
   *
   * IN
   * Nneigtbours 		the number of neighbouring cells / cooling systems etc.
   * Tneighbours  	array with the temperature of the neighbouring  elements [K]
   * Kneighbours		array with the heat transfer constants of the neighbouring elements, k or h
   * Aneighb 			array with the surface area of the neigbouring cell
   * 						the heat transfer happens over the smaller one of Aneighb[i] and the area of this cell
   * tim 				the time since the last thermal model calculation [for verification]
   *
   * throws
   * 8 	invalid time keeping
   * 9 	invalid module temperature
   */

  //!< check the time since the last checkup is not too large.
  //!< if the parent has not called the thermal model for a while, the equation becomes unstable
  //!< 		cause E = time * kA dT, so even a small dT will cause a huge E, and therefore a very large temperature swint
  //
  if (Therm_time > 15 * 60) //!< #TODO magic number.
  {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Warning in Cell_SPM::thermalModel, the time since this function was called last is very large, "
                << Therm_time << " which might lead to excessive temperature variations" << '\n';
  }

  //!< then check whether our internal time keeping matches up with the external one
  if (std::abs(Therm_time - tim) > 1) {
    if constexpr (settings::printBool::printCrit) {
      std::cerr << "ERROR in Cell_SPM::thermalModel, according to the cell's internal timing, " << Therm_time
                << "s have passed since the last thermal model solution. The external time provided was "
                << tim << "s, which is more than 1s difference. Throwing an error.\n";
    }
    throw 8;
  }

  //!< calculate the total thermal balance
  double Etot = Therm_Qgen;
  double Atherm;
  for (int i = 0; i < Nneighbours; i++) {
    Atherm = std::min(Aneighb[i], getThermalSurface());
    Etot += Kneighbours[i] * Atherm * (Tneighbours[i] - T()) * Therm_time;
  }

  //!< Calculate the new temperature
  //!< rho * cp * dT/dt = Qtot / V
  //!< 		where 	Qtot is total power in W
  //!< 				V is the cell's volume L * elec_surf
  //!< so integrated over time this is
  //!< rho * cp * (Tnew - Told) = Etot / V
  const double Tnew = T() + Etot / (rho * Cp * geo.L * geo.elec_surf);

  //!< Check the new temperature is valid, and if so, set it
  if (Tnew > Tmax() || Tnew < Tmin() || std::isnan(Tnew)) {
    if constexpr (settings::printBool::printCrit) {
      std::cerr << "ERROR in Cell_SPM::thermalModel, the new temperature of " << Tnew
                << " is outside the allowed range from " << Tmin() << " to " << Tmax()
                << ". The time since the last time this function was called is " << Therm_time << '\n';

      std::cout << "Total thermal energy " << Etot << ". internal heat generation "
                << Therm_Qgen << " and external contributions as below:\n";

      for (int i = 0; i < Nneighbours; i++)
        std::cout << Aneighb[i] << ", " << Kneighbours[i] << "," << Tneighbours[i] << ","
                  << Kneighbours[i] * Atherm * (Tneighbours[i] - T()) * Therm_time << '\n';

      std::cout << "giving change in temperature: " << Etot / (rho * Cp * geo.L * geo.elec_surf) << '\n';
    }
    throw 9;
  }

  //!< setting the temperature is done by the parent module. else some cells will update their T before others, and we get inconsistencies
  //!< 		e.g. exchange from cell 2 to this cell will not be same as from this cell to cell 2 since T of this cell would have changed.

  //!< Reset the cumulative thermal variables
  Therm_Qgentot += Therm_Qgen;
  Therm_Qgen = 0;
  Therm_time = 0;

  return Tnew;
}

//!< total heat generation [J] since start of cell's life
double Cell_SPM::thermal_getTotalHeat() { return Therm_Qgentot; }

/**
 * Calculate the thermal model for this cell and update its cell temperature.
 * The heat exchange in [W] with element i is given by
 * 		Qc = Kneighbours[i]*A*(Tneighbours[i] - getT());
 * We assume all temperatures were constant for the past period, such that the thermal energy in [J] is given by
 * 		Ec = Qc * Therm_time
 * This is added up with the internal heath generated, and the cell's temperature is returned
 *
 * IN
 * Nneigtbours 		the number of neighbouring cells / cooling systems etc.
 * Tneighbours  	array with the temperature of the neighbouring  elements [K]
 * Kneighbours		array with the heat transfer constants of the neighbouring elements, k or h
 * Aneighb 			array with the surface area of the neigbouring cell
 * 						the heat transfer happens over the smaller one of Aneighb[i] and the area of this cell
 * tim 				the time since the last thermal model calculation [for verification]
 *
 * throws
 * 8 	invalid time keeping
 * 9 	invalid module temperature
 **/
double Cell_SPM::thermalModel(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim)
{
  double Tnew;
  try {
    if constexpr (settings::T_MODEL == 0) //!< #TODO thermal model implementation should be outside.
      Tnew = T();
    else if constexpr (settings::T_MODEL == 1)
      Tnew = thermalModel_cell();
    else if constexpr (settings::T_MODEL == 2)
      Tnew = thermalModel_coupled(Nneighbours, Tneighbours, Kneighbours, Aneighb, tim);
  } catch (int e) {

    //!< indicate we have ran the thermal model
    Therm_Qgentot += Therm_Qgen;
    Therm_Qgen = 0;
    Therm_time = 0;
    throw e;
  }
  //!< setting the temperature is done by the parent module. else some cells will update their T before others,
  //!< and we get inconsistencies e.g. exchange from cell 2 to this cell will not be same as from this cell to
  //!< cell 2 since T of this cell would have changed.
  //!< Reset the cumulative thermal variables
  Therm_Qgentot += Therm_Qgen;
  Therm_Qgen = 0;
  Therm_time = 0;

  return Tnew;
}

} // namespace slide