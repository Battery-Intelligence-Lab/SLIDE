/*
 * State_Bucket.hpp
 *
 *  Created on: 12 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <array>
#include <cstdlib> //!< for size_t

namespace slide {
struct State_Bucket : public std::array<double, 3>
{
  enum Index : size_t //!< Index variables for:
  {
    i_I,
    i_SOC,
    i_T, //!< cell temperature [K]
    N_states,
  };

  inline auto &I() { return (*this)[i_I]; }     //!< Current, [I]
  inline auto &SOC() { return (*this)[i_SOC]; } //!< state of charge [0-1]
  inline auto &T() { return (*this)[i_T]; }     //!< temperature, [K]
};
} // namespace slide