/*
 * State.hpp
 *
 * Generic State class to hold time, Ah, Wh.
 *
 *  Created on: 02 Sep 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../settings/settings.hpp"

#include <array>

namespace slide {

template <size_t N, size_t N_cumulative = settings::data::N_CumulativeData>
struct State : public std::array<double, N + N_cumulative>
{
  enum Index : size_t //!< Index variables for:
  {
    i_time = N,
    i_Ah,
    i_Wh,
    N_states,
  };

  inline auto &time() { return (N_cumulative != 3) ? (*this)[0] : (*this)[i_time]; } //!< time [s]
  inline auto &Ah() { return (N_cumulative != 3) ? (*this)[0] : (*this)[i_Ah]; }     //!< Current throughput [Ah]
  inline auto &Wh() { return (N_cumulative != 3) ? (*this)[0] : (*this)[i_Wh]; }     //!< Energy throughput [Wh]

  constexpr static auto description(size_t i)
  {
    if constexpr (N_cumulative != 3)
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