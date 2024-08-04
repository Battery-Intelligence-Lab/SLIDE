/**
 * @file State.hpp
 * @brief Generic State class to hold time, Ah, Wh.
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 02 Sep 2022
 */

#pragma once

#include "../settings/settings.hpp"
#include "../utility/array_util.hpp"

#include <array>

namespace slide {

template <size_t N>
struct State : public std::array<double, N + 3>
{
  enum Index : size_t //!< Index variables for:
  {
    i_time = N,
    i_Ah,
    i_Wh,
    N_states,
  };

  constexpr auto &time() { return (*this)[i_time]; } //!< time [s]
  constexpr auto &Ah() { return (*this)[i_Ah]; }     //!< Current throughput [Ah]
  constexpr auto &Wh() { return (*this)[i_Wh]; }     //!< Energy throughput [Wh]

  constexpr auto time() const { return (*this)[i_time]; } //!< time [s]
  constexpr auto Ah() const { return (*this)[i_Ah]; }     //!< Current throughput [Ah]
  constexpr auto Wh() const { return (*this)[i_Wh]; }     //!< Energy throughput [Wh]

  constexpr static auto description(size_t i)
  {
    if (i == i_time)
      return "time [s]";
    else if (i == i_Ah)
      return "Current throughput [Ah]";
    else if (i == i_Wh)
      return "Energy throughput [Wh]";
    else
      return "";
  }

  constexpr auto reset() { this->fill(0); }
};

using ThroughputData = State<0>; // time, Ah, Wh

} // namespace slide