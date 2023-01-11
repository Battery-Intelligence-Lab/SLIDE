/*
 * State_SPM.hpp
 *
 * Implements a class State which defines the state-variables of a cell for the state-space model formulation
 *
 *  Created on: 12 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../../settings/settings.hpp"
#include "../../types/State.hpp"

#include <cstdlib>
#include <array>
#include <span>

namespace slide {
class State_SPM : public State<29> //!< #TODO how can we make this so it takes 29=N_states from enum?
{
public:
  constexpr static auto nch = settings::nch;

  enum Index : size_t //!< Index variables for:
  {
    i_I,
    i_V,
    i_T,                 //!< cell temperature [K]
    i_delta,             //!< thickness of the SEI layer [m]
    i_LLI,               //!< lost lithium inventory [As]
    i_thickp,            //!< thickness of the cathode [m]
    i_thickn,            //!< thickness of the anode [m]
    i_ep,                //!< volume fraction of active material in the cathode [-]
    i_en,                //!< volume fraction of active material in the anode [-]
    i_ap,                //!< effective surface area of the porous cathode [m2 m-3]
    i_an,                //!< effective surface area of the porous anode [m2 m-3]
    i_CS,                //!< surface area of the cracks at the surface of the negative particle [m2]
    i_Dp,                //!< diffusion constant of the cathode at reference temperature [m s-1]
    i_Dn,                //!< diffusion constant of the anode at reference temperature [m s-1]
    i_delta_pl,          //!< thickness of the plated lithium layer [m]
    i_zp,                //!< transformed concentration at the positive inner Chebyshev nodes of the positive particle
    i_zn = i_zp + nch,   //!< transformed concentration at the positive inner Chebyshev nodes of the negative particle
    i_rDCp = i_zn + nch, //!< (ONLY CATHODE) specific resistance of both electrodes combined [Ohm m2]
    i_rDCn,              //!< (ONLY ANODE) specific resistance of both electrodes combined [Ohm m2]
    i_rDCcc,             //!< (ONLY SEPARATOR) specific resistance of both electrodes combined [Ohm m2]
    i_SOC,
    N_states,
    N_save = i_zp //!< Save until i_zp
  };

  /* (OLD NOTE) Note on ri:
   * this is the resistance times the electrode surface, averaged between both electrodes.
   * The total cell resistance is (see Cell::getR() ): r /( (thickp*ap*elec_surf + thickn*an*elec_surf)/2 )
   * 	with r the resistance times the average electrode surface
   * 		 thicki the thickness of electrode i
   * 		 ai the specific surface area of electrode i
   * 		 elec_surf the geometric surface area of the electrode (height of the electrode * width of the electrode)
   * so if the measured DC resistance of the cell is R, the value of r can be calculated using:
   * 		 ri = R * ( (thickp*ap*elec_surf + thickn*an*elec_surf)/2 )
   */

  using z_type = std::array<value_type, nch>;
  using states_type = std::array<value_type, N_states>;

  //!< State() = default; //!< Default constructor which DOESN'T initialise the states. All states are set to 0
  //!< State(const slide::states_type &s) : x{s} {}

  inline auto &zp(size_t i) { return (*this)[i_zp + i]; } //!< z_type, transformed li concentration at the positive inner nodes of the positive particle
  inline auto &zn(size_t i) { return (*this)[i_zn + i]; } //!< z_type, transformed li concentration at the positive inner nodes of the negative particle
  inline auto &z(size_t i) { return (*this)[i_zp + i]; }  //!< Both z_p and z_n;

  inline auto zp() { return std::span<double>(&zp(0), &zn(0)); }
  inline auto zn() { return std::span<double>(&zn(0), &zn(0) + nch); }
  inline auto z() { return std::span<double>(&zp(0), &zn(0) + nch); }


  inline auto &T() { return (*this)[i_T]; }               //!< cell temperature [K]
  inline auto &delta() { return (*this)[i_delta]; }       //!< thickness of the SEI layer [m]
  inline auto &LLI() { return (*this)[i_LLI]; }           //!< lost lithium [As]
  inline auto &thickp() { return (*this)[i_thickp]; }     //!< thickness of the cathode [m]
  inline auto &thickn() { return (*this)[i_thickn]; }     //!< thickness of the anode [m]
  inline auto &ep() { return (*this)[i_ep]; }             //!< volume fraction of active material in the cathode [-]
  inline auto &en() { return (*this)[i_en]; }             //!< volume fraction of active material in the anode [-]
  inline auto &ap() { return (*this)[i_ap]; }             //!< effective surface area of the cathode [m2 m-3]
  inline auto &an() { return (*this)[i_an]; }             //!< effective surface area of the anode [m2 m-3]
  inline auto &CS() { return (*this)[i_CS]; }             //!< surface area of the cracks at the surface of the negative particle [m2]
  inline auto &Dp() { return (*this)[i_Dp]; }             //!< diffusion constant at reference temperature of the cathode [m s-1]
  inline auto &Dn() { return (*this)[i_Dn]; }             //!< diffusion constant at reference temperature of the anode [m s-1]
  inline auto &rDCp() { return (*this)[i_rDCp]; }         //!< specific resistance (resistance times real surface area of the combined electrodes) [Ohm m2]
  inline auto &rDCn() { return (*this)[i_rDCn]; }         //!< specific resistance (resistance times real surface area of the combined electrodes) [Ohm m2]
  inline auto &rDCcc() { return (*this)[i_rDCcc]; }       //!< specific resistance (resistance times real surface area of the combined electrodes) [Ohm m2]
  inline auto &delta_pl() { return (*this)[i_delta_pl]; } //!< thickness of the plated lithium layer [m]
  inline auto &SOC() { return (*this)[i_SOC]; }           //!< thickness of the plated lithium layer [m]
  inline auto &I() { return (*this)[i_I]; }               //!< current [A]
  inline auto &V() { return (*this)[i_V]; }               //!< voltage [V]

  //!< void setT(double Ti);                                                                                          //!< set the temperature
  void overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani); //!< overwrite the states related to the geometry of a cell
  void overwriteCharacterisationStates(double Dpi, double Dni, double ri);                                       //!< overwrite the states related to the characterisation of a cell

  std::span<double> viewGeometricStates() { return std::span<double>(&delta(), &delta_pl() + 1); } //!< #Check and fix. Why this also checks SOC?

  // Const methods:

  inline auto I() const { return (*this)[i_I]; } //!< current [A]
};

inline void State_SPM::overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani)
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
  this->thickp() = thickpi;
  this->thickn() = thickni;
  this->ep() = epi;
  this->en() = eni;
  this->ap() = api;
  this->an() = ani;
}

inline void State_SPM::overwriteCharacterisationStates(double Dpi, double Dni, double ri)
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
  this->Dp() = Dpi;
  this->Dn() = Dni;
  this->rDCp() = ri; //!< Considering both anode and cathode r is same.
  this->rDCn() = ri;
  //!< #TODO -> r was //!< the specific resistance (resistance times real surface area of the combined electrodes) [Ohm m2]
}
} // namespace slide