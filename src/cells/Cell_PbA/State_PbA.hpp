/*
 * State_PbA.hpp
 *
 *  Created on: 12 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../../types/State.hpp"

namespace slide {

struct State_PbA : public State<5>
{
  enum Index : size_t //!< Index variables for:
  {
    i_I,
    i_SOC,
    i_T, //!< cell temperature [K]
    i_Delta_W,
    i_c_H2SO4,
    i_rho_c,
    i_rho_d,
    i_Znom,
    i_Zw,
    i_Delta_tSOC,
    i_SOC_min,
    i_n_bad,
    i_f_stratification,
    N_states, // Do not use N_states for total states, use .size()
  };

  auto &I() { return (*this)[i_I]; }             //!< Current, [I]
  auto &SOC() { return (*this)[i_SOC]; }         //!< state of charge [0-1]
  auto &T() { return (*this)[i_T]; }             //!< temperature, [K]
  auto &Delta_W() { return (*this)[i_Delta_W]; } //!< < An effective layer thickness, Delta W
  auto &c_H2SO4() { return (*this)[i_c_H2SO4]; }
  auto &rho_c() { return (*this)[i_T]; }
  auto &rho_d() { return (*this)[i_T]; }
  auto &Znom() { return (*this)[i_T]; }
  auto &Zw() { return (*this)[i_T]; }
  auto &Delta_tSOC() { return (*this)[i_T]; } //!< < time since the last full charge
  auto &SOC_min() { return (*this)[i_T]; }    //!< < the lowest SOC since the last full charge
  auto &n_bad() { return (*this)[i_T]; }      //!< < Number of bad charges
  auto &f_stratification() { return (*this)[i_T]; }
};


} // namespace slide