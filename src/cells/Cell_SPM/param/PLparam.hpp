/*
 * PLparam.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

namespace slide::param {
//!< Define a structure with the fitting parameters of the li-plating models (PL)
struct PLparam
{
  double pl1k{ 4.5e-10 };       //!< rate constant of the li-plating side reaction at reference temperature in the 1st model
  double pl1k_T{ -2.014008e5 }; //!< activation energy of pl1k

  [[nodiscard]] constexpr auto begin() noexcept { return &pl1k; }
  [[nodiscard]] constexpr auto end() noexcept { return &pl1k_T + 1; }

  auto &operator*=(double a)
  {
    for (auto &item : *this)
      item *= a;

    return *this;
  }
};
} // namespace slide::param