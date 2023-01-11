/*
 * State_ECM.hpp
 *
 *  Created on: 14 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../../types/State.hpp"

namespace slide {
struct State_ECM : public State<4>
{
  enum Index : size_t //!< Index variables for:
  {
    i_I,
    i_Ir,
    i_SOC,
    i_T,      //!< cell temperature [K]
    N_states, // Do not use N_states for total states, use .size()
  };

  // Const methods:
  inline auto I() const { return (*this)[i_I]; } //!< Current, [A], + for discharge, - for charge

  // Non-const methods:
  inline auto &I() { return (*this)[i_I]; }     //!< Current, [A], + for discharge, - for charge
  inline auto &Ir() { return (*this)[i_Ir]; }   //!< Current through the parallel resistance, [I]
  inline auto &SOC() { return (*this)[i_SOC]; } //!< state of charge [0-1]
  inline auto &T() { return (*this)[i_T]; }     //!< temperature, [K]
};
} // namespace slide