/*
 * State_ECM.hpp
 *
 *  Created on: 14 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <array>

namespace slide {
struct State_ECM : public std::array<double, 4>
{
  inline auto &I() { return this->at(0); }   //!< Current, [A], + for discharge, - for charge
  inline auto &Ir() { return this->at(1); }  //!< Current through the parallel resistance, [I]
  inline auto &SOC() { return this->at(2); } //!< state of charge [0-1]
  inline auto &T() { return this->at(3); }   //!< temperature, [K]
};
} // namespace slide