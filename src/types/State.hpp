/*
 * State.hpp
 *
 * Generic State class to hold time, Ah, Wh.
 *
 *  Created on: 02 Sep 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once


namespace slide {
struct State
{
  enum Index : size_t //!< Index variables for:
  {
    i_time,
    i_Ah,
    i_Wh,
    N_states,
  };

  inline auto &I() { return (*this)[i_I]; }     //!< Current, [I]
  inline auto &SOC() { return (*this)[i_SOC]; } //!< state of charge [0-1]
  inline auto &T() { return (*this)[i_T]; }     //!< temperature, [K]
};
} // namespace slide