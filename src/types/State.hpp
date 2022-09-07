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

  inline auto &time() { return (N_states != 3) ? 0 : (*this)[i_time]; } //!< Current, [I]
  inline auto &Ah() { return (N_states != 3) ? 0 : (*this)[i_Ah]; }     //!< state of charge [0-1]
  inline auto &Wh() { return (N_states != 3) ? 0 : (*this)[i_Wh]; }     //!< temperature, [K]

  constexpr static auto description(size_t i)
  {
    if (N_states != 3)
      return "";
    else {
      if (i == i_time)
        return "time [s]";
      else if (i == i_Ah)
        return "Current throughput [Ah]";
      else if (i == i_Wh)
        return "Energy throughput [Wh]";
      else
        return "";
    }
  }
};

} // namespace slide