/*
 * State_ECM.hpp
 *
 *  Created on: 14 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../../types/State.hpp"

namespace slide {

template <size_t N_RC = 1>
struct State_ECM : public State<3 + N_RC>
{
  enum Index : size_t //!< Index variables for:
  {
    i_T, //!< cell temperature [K]
    i_SOC,
    i_I,
    i_Ir,
    N_states, // #TODO Do not use N_states for total states, use .size()
  };

  // Const methods:
  inline auto I() const { return (*this)[i_I]; } //!< Current, [A], + for discharge, - for charge

  // Non-const methods:
  inline auto &I() { return (*this)[i_I]; }                      //!< Current, [A], + for discharge, - for charge
  inline auto &Ir(size_t i = 0) { return (*this)[i_I + i + 1]; } //!< Current through the parallel resistance, [I]
  inline auto &SOC() { return (*this)[i_SOC]; }                  //!< state of charge [0-1]
  inline auto &T() { return (*this)[i_T]; }                      //!< temperature, [K]
};

using State_Bucket = State_ECM<0>;

} // namespace slide