/*
 * state.cpp
 *
 * Implements a class State which defines the state-variables of a cell for the state-space model formulation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "state.hpp"
#include "../utility/utility.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <array>

namespace slide {

void State::initialise(z_type &zpi, z_type &zni, double Ti, double deltai, double LLIi,
                       double thickpi, double thickni, double epi, double eni, double api, double ani,
                       double CSi, double Dpi, double Dni, double ri, double deltalii)
{
  /*
   * Initialise the state variables to the given values.
   * Use this function only once, immediately after calling the constructor.
   * Calling it again at a later point will generate an error (11)
   *
   * IN
   * nin 		length of the arrays with the (transformed) concentration
   * zpi 		transformed concentration at the positive inner Chebyshev nodes of the positive particle
   * zni		transformed concentration at the positive inner Chebyshev nodes of the negative particle
   * Ti 		cell temperature [K]
   * deltai 	thickness of the SEI layer [m]
   * LLIi 	lost lithium inventory [As]
   * thickpi 	thickness of the cathode [m]
   * thickni 	thickness of the anode [m]
   * epi 		volume fraction of active material in the cathode [-]
   * eni 		volume fraction of active material in the anode [-]
   * api 		effective surface area of the porous cathode [m2 m-3]
   * ani 		effective surface area of the porous anode [m2 m-3]
   * CSi 		surface area of the cracks at the surface of the negative particle [m2]
   * Dpi 		diffusion constant of the cathode at reference temperature [m s-1]
   * Dni 		diffusion constant of the anode at reference temperature [m s-1]
   * ri 		specific resistance of both electrodes combined [Ohm m2]
   * deltalii thickness of the plated lithium layer [m]
   *
   * Note on ri:
   * this is the resistance times the electrode surface, averaged between both electrodes.
   * The total cell resistance is (see Cell::getR() ): r /( (thickp*ap*elec_surf + thickn*an*elec_surf)/2 )
   * 	with r the resistance times the average electrode surface
   * 		 thicki the thickness of electrode i
   * 		 ai the specific surface area of electrode i
   * 		 elec_surf the geometric surface area of the electrode (height of the electrode * width of the electrode)
   * so if the measured DC resistance of the cell is R, the value of r can be calculated using:
   * 		 ri = R * ( (thickp*ap*elec_surf + thickn*an*elec_surf)/2 )
   *
   * THROWS
   * 10		the arrays have the wrong length
   * 11 		this function is called when the initial states have already been initialised
   * 12 		the state suggested is illegal
   */
  //!< Set the state variables
  std::copy(zpi.begin(), zpi.end(), this->begin());
  std::copy(zni.begin(), zni.end(), this->begin() + zpi.size());
  get_T() = Ti;
  get_delta() = deltai;
  get_LLI() = LLIi;
  get_thickp() = thickpi;
  get_thickn() = thickni;
  get_ep() = epi;
  get_en() = eni;
  get_ap() = api;
  get_an() = ani;
  get_CS() = CSi;
  get_Dp() = Dpi;
  get_Dn() = Dni;
  get_r() = ri;
  get_delta_pl() = deltalii;
}

void State::setT(double Ti)
{
  /*
   * Sets the temperature to the given value
   *
   * IN
   * Ti		temperature [K], must be between 0 and 60 degrees, so 273 and 273+60
   *
   * THROWS
   * 13 		illegal value of T
   */

  //!< the temperature limits are defined in State.hpp
  if (Ti < settings::Tmin_K) //!< check the temperature is above 0 degrees
  {
    std::cerr << "Error in State::setT. The temperature " << Ti << "K is too low. The minimum value is " << settings::Tmin_K << std::endl;
    throw 13;
  } else if (Ti > settings::Tmax_K) //!< check the temperature is below 60 degrees
  {
    std::cerr << "Error in State::setT. The temperature " << Ti << "K is too high. The maximum value is " << settings::Tmax_K << std::endl;
    throw 13;
  }

  get_T() = Ti;
}

void State::overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani)
{
  /*
   * Function to overwrite the geometric parameters of the state.
   * It also overwrites the initial states, so use it with extreme caution.
   * It should only be called when you are parametrising a cell (determineCharacterisation.cpp), never while cycling a cell.
   *
   * IN
   * thickpi 	thickness of the cathode [m]
   * thickni 	thickness of the anode [m]
   * epi 		volume fraction of active material in the cathode [-]
   * eni 		volume fraction of active material in the anode [-]
   * api 		effective surface area of the porous cathode [m2 m-3]
   * ani 		effective surface area of the porous anode [m2 m-3]
   */

  //!< set the states
  this->get_thickp() = thickpi;
  this->get_thickn() = thickni;
  this->get_ep() = epi;
  this->get_en() = eni;
  this->get_ap() = api;
  this->get_an() = ani;
}

void State::overwriteCharacterisationStates(double Dpi, double Dni, double ri)
{
  /*
   * Function to overwrite the parameters related to the characterisation of the cell.
   * The states and initial states are overwritten so use this function with caution.
   * It should only be called when you are parametrising a cell (determineCharacterisation.cpp), never while cycling a cell.
   *
   * IN
   * Dpi	diffusion constant of the cathode at rate temperature [m s-1]
   * Dni 	diffusion constant of the anode at rate temperature [m s-1]
   * r 	specific resistance of the combined electrodes [Ohm m2]
   */

  //!< Set the states
  this->get_Dp() = Dpi;
  this->get_Dn() = Dni;
  this->get_r() = ri;
}

} // namespace slide
