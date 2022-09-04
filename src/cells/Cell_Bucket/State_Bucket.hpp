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
    i_I = State::N_states,
    i_SOC,
    i_T, //!< cell temperature [K]
    N_states,
  };

  constexpr static std::array<const char *, N_states> description{
    "Current [A]",
    "SOC [-]",
    "Temperature [K]"
  };

  inline auto &I() { return (*this)[i_I]; }     //!< Current, [I]
  inline auto &SOC() { return (*this)[i_SOC]; } //!< state of charge [0-1]
  inline auto &T() { return (*this)[i_T]; }     //!< temperature, [K]
};
} // namespace slide