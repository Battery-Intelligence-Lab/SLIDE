/**
 * @file State_SPM.hpp
 * @brief  Implements a class State which defines the state-variables of a cell for the state-space model formulation
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 12 Apr 2022
 */

#pragma once

#include "Pair.hpp"
#include "settings.hpp"
#include "State.hpp"

#include <cstdlib>
#include <array>
#include <span>

namespace slide {
class State_SPM : public State<19 + 2 * settings::nch> //!< #TODO how can we make this so it takes 29=N_states from enum?
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
  inline auto &z(size_t i, Domain domain) { return (*this)[i_zp + nch * domain + i]; }

  template <size_t domain>
  inline auto &z(size_t i)
  {
    if constexpr (domain == 0)
      return zp(i);
    else
      return zn(i);
  }

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

  // Additional definitions for domain-based access:
  inline auto &thick(size_t domain) { return (*this)[i_thickp + domain]; }
  inline auto &e(size_t domain) { return (*this)[i_ep + domain]; }
  inline auto &a(size_t domain) { return (*this)[i_ap + domain]; }
  inline auto &D(size_t domain) { return (*this)[i_Dp + domain]; }
  inline auto &rDC(size_t domain) { return (*this)[i_rDCp + domain]; }


  //!< void setT(double Ti);                                                                                          //!< set the temperature
  void overwriteGeometricStates(DPair thick_, DPair e_, DPair a_); //!< overwrite the states related to the geometry of a cell
  void overwriteCharacterisationStates(DPair D_, double ri);       //!< overwrite the states related to the characterisation of a cell

  std::span<double> viewGeometricStates() { return std::span<double>(&delta(), &delta_pl() + 1); } //!< #Check and fix. Why this also checks SOC?

  // Const methods:

  inline auto I() const { return (*this)[i_I]; } //!< current [A]
};

inline void State_SPM::overwriteGeometricStates(DPair thick_, DPair e_, DPair a_)
{
  /*
   * Function to overwrite the geometric parameters of the state.
   * It also overwrites the initial states, so use it with extreme caution.
   * It should only be called when you are parametrising a cell (determineCharacterisation.cpp), never while cycling a cell.
   *
   * IN
   * thick_ 	thickness of the cathode/anode [m]
   * e_ 		volume fraction of active material in the cathode/anode [-]
   * a_ 		effective surface area of the porous cathode/anode [m2 m-3]
   */

  //!< set the states
  for (auto dom : { pos, neg }) {
    this->thick(dom) = thick_[dom];
    this->e(dom) = e_[dom];
    this->a(dom) = a_[dom];
  }
}

inline void State_SPM::overwriteCharacterisationStates(DPair D_, double ri)
{
  /*
   * Function to overwrite the parameters related to the characterisation of the cell.
   * The states and initial states are overwritten so use this function with caution.
   * It should only be called when you are parametrising a cell (determineCharacterisation.cpp), never while cycling a cell.
   *
   * IN
   * D_	diffusion constant of the cathode/anode at rate temperature [m s-1]
   * r 	specific resistance of the combined electrodes [Ohm m2]
   */

  //!< Set the states
  for (auto dom : { pos, neg }) {
    this->D(dom) = D_[dom];
    this->rDC(dom) = ri; //!< Considering both anode and cathode r is same.
  }

  //!< #TODO -> r was //!< the specific resistance (resistance times real surface area of the combined electrodes) [Ohm m2]
}
} // namespace slide