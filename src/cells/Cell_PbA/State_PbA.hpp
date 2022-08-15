/*
 * State_PbA.hpp
 *
 *  Created on: 12 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <array>

namespace slide {
struct State_PbA : public std::array<double, 4>
{
  inline auto &I() { return this->at(0); }       //!< < Current, [I]
  inline auto &SOC() { return this->at(1); }     //!< < state of charge [0-1]
  inline auto &T() { return this->at(2); }       //!< < temperature, [K]
  inline auto &Delta_W() { return this->at(3); } //!< < An effective layer thickness, Delta W
  inline auto &rho_c() { return this->at(4); }
  inline auto &rho_d() { return this->at(5); }
  inline auto &Znom() { return this->at(6); }
  inline auto &Zw() { return this->at(7); }
  inline auto &Delta_tSOC() { return this->at(8); } //!< < time since the last full charge
  inline auto &SOC_min() { return this->at(9); }    //!< < the lowest SOC since the last full charge
  inline auto &n_bad() { return this->at(10); }     //!< < Number of bad charges
  inline auto &f_stratification() { return this->at(11); }
};
} // namespace slide