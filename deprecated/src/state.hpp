/*
 * State.hpp
 *
 * Defines a class State which defines the state-variables of a cell for the state-space model formulation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <array>
#include <utility>

#include "settings/settings.hpp"

namespace slide {
template <size_t Nstates>
class State : public std::array<double, Nstates>
{
};

using z_type = std::array<double, settings::nch>;
using states_type = std::array<double, settings::ns>;
using settings::nch;

class State_SPM : public std::array<double, settings::ns>
{

public:
  //!< State() = default; //!<Default constructor which DOESN'T initialise the states. All states are set to 0
  //!< State(const slide::states_type &s) : x{s} {}
  void initialise(z_type &zpi, z_type &zni, double Ti, double deltai, double LLIi,
                  double thickpi, double thickni, double epi, double eni, double api, double ani,
                  double CSi, double Dpi, double Dni, double Ri, double deltalii); //!< initialises all state variables to the given values

  //!< inline const auto &getStates_arr() const { return x; } //!<get states in array format.

  inline const auto &get_zp(int i) { return this->at(i); }       //!< z_type, get the transformed concentrations
  inline const auto &get_zn(int i) { return this->at(nch + i); } //!< z_type, get the transformed concentrations

  inline auto &get_T() { return this->at(2 * nch + 0); }         //!< get the temperature [K]
  inline auto &get_delta() { return this->at(2 * nch + 1); }     //!< get the SEI thickness
  inline auto &get_LLI() { return this->at(2 * nch + 2); }       //!< get the lost lithium
  inline auto &get_thickp() { return this->at(2 * nch + 3); }    //!< get the thickness of the cathode
  inline auto &get_thickn() { return this->at(2 * nch + 4); }    //!< get the thickness of the anode
  inline auto &get_ep() { return this->at(2 * nch + 5); }        //!< get the volume fraction of active material of the cathode
  inline auto &get_en() { return this->at(2 * nch + 6); }        //!< get the volume fraction of active material of the anode
  inline auto &get_ap() { return this->at(2 * nch + 7); }        //!< get the effective surface of the cathode
  inline auto &get_an() { return this->at(2 * nch + 8); }        //!< get the effective surface of the anode
  inline auto &get_CS() { return this->at(2 * nch + 9); }        //!< get the crack surface area
  inline auto &get_Dp() { return this->at(2 * nch + 10); }       //!< get the diffusion constant of the cathode at reference temperature [m/s]
  inline auto &get_Dn() { return this->at(2 * nch + 11); }       //!< get the diffusion constant of the anode at reference temperature [m/s]
  inline auto &get_r() { return this->at(2 * nch + 12); }        //!< get the specific resistance (resistance times real surface area of the combined electrodes) [Ohm m2]
  inline auto &get_delta_pl() { return this->at(2 * nch + 13); } //!< get the thickness of the plated lithium layer

  double getR(double elec_surf) //!< get the total DC resistance of the electrodes
  {
    /*
     * IN
     * elec_surf geometric surface area of the electrodes [m2]
     */
    const double surfp = get_ap() * get_thickp() * elec_surf; //!< real surface area of the positive electrode
    const double surfn = get_an() * get_thickn() * elec_surf; //!< real surface area of the negative electrode
    const double surf = (surfp + surfn) / 2;                  //!< mean surface area
    const double r_elec = get_r() / surf;                     //!< resistance of the electrodes [Ohm]
    return r_elec;
  }

  void setT(double Ti);                           //!< set the temperature
  void setZ(const z_type &zpi, const z_type &zni) //!< set the transformed concentrations
  {

    std::copy(zpi.begin(), zpi.end(), this->begin());
    std::copy(zni.begin(), zni.end(), this->begin() + nch);
  }                                                                                                              //!< set the transformed concentration
  void overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani); //!< overwrite the states related to the geometry of a cell
  void overwriteCharacterisationStates(double Dpi, double Dni, double ri);                                       //!< overwrite the states related to the characterisation of a cell

  //!< private:
  //!<  battery states
  //!<  zp[nch] zn[nch] T delta LLI thickp thickn ep en ap an CS Dp Dn R delta_liPlating
  //!<  slide::states_type x{}; //!<Array to hold all states.
  //!<  						//	slide::State *sini_ptr{nullptr}; //!<array ptr with the initial battery states
};
} // namespace slide