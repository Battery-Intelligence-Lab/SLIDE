/*
 * Geometry_SPM.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

namespace slide::param {
//!< Define a structure with the fitting parameters of the SEI growth models (SEI)
struct Geometry_SPM
{
  double L{ 1.6850e-4 };          //!< thickness of one layer of the cell [m]
  double Acell{ 0.1 * 0.2 };      //!< geometric surface area of the cell [m2] / width * height of the pouch.
  double elec_surf{ Acell * 31 }; //!< geometric surface area of the electrodes (electrode height * electrode width*layers) [m2] Doubly coated so multiply by two when putting into PyBAMM
  double SAV{ 252.9915 };         //!< surface area to volume-ratio of the cell [m2/m3]

  double Rp{ 8.5e-6 }, Rn{ 1.25e-5 }; //!< radius of the positive/neg sphere of the Single Particle model [m]
  //!< do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied MATLAB scripts.
  //!< See the word document '2 overview of the code', section 'MATLAB setup before running the C++ code'

  //!< other geometric parameters are part of State because they can change over the battery'
};
} // namespace slide::param