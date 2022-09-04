/*
 * State.hpp
 *
 * Generic State class to hold time, Ah, Wh.
 *
 *  Created on: 02 Sep 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <array>

#include "../settings/settings.hpp"

namespace slide {

template <size_t N, size_t N_cumulative = settings::data::N_CumulativeData>
struct State : public std::array<double, N + N_cumulative>
{
  enum Index : size_t //!< Index variables for:
  {
    i_time,
    i_Ah,
    i_Wh,
    N_states = N_cumulative,
  };

  inline auto &time() { return (*this)[i_time]; } //!< Current, [I]
  inline auto &Ah() { return (*this)[i_Ah]; }     //!< state of charge [0-1]
  inline auto &Wh() { return (*this)[i_Wh]; }     //!< temperature, [K]
};

} // namespace slide