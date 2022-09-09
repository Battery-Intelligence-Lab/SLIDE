/*
 * State_Bucket.hpp
 *
 *  Created on: 12 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../../types/State.hpp"

namespace slide {
struct State_Bucket : public State<3>
{
  enum Index : size_t //!< Index variables for:
  {
    i_I,
    i_SOC,
    i_T,      //!< cell temperature [K]
    N_states, // Do not use N_states for total states, use .size()
  };

  inline auto &I() { return (*this)[i_I]; }     //!< Current, [I]
  inline auto &SOC() { return (*this)[i_SOC]; } //!< state of charge [0-1]
  inline auto &T() { return (*this)[i_T]; }     //!< temperature, [K]

  constexpr static auto description(size_t i)
  {
    if (i >= N_states)
      return State::description(i);
    else {
      if (i == i_I)
        return "Current [A]";
      else if (i == i_SOC)
        return "SOC [-]";
      else if (i == i_T)
        return "Temperature [K]";
      else
        return "";
    }
  }
};
} // namespace slide