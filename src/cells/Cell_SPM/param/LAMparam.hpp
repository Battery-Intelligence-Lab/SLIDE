/*
 * LAMparam.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

namespace slide::param {
//!< Define a structure with the fitting parameters of the models for loss of active material (LAM)
struct LAMparam //!< #TODO if begin and end is ok for alignment.
{
  double lam1p; //!< fitting parameter for the positive electrode for the 1st LAM model
  double lam1n; //!< fitting parameter for the negative electrode for the 1st LAM model

  double lam2ap; //!< fitting parameter 1 at reference temperature for the positive electrode for the 2nd LAM model
  double lam2bp; //!< fitting parameter 2 at reference temperature for the positive electrode for the 2nd LAM model
  double lam2an; //!< fitting parameter 1 at reference temperature for the negative electrode for the 2nd LAM model
  double lam2bn; //!< fitting parameter 2 at reference temperature for the negative electrode for the 2nd LAM model
  double lam2t;  //!< activation energy for all the parameters of the 2nd LAM model

  double lam3k;   //!< rate constant at reference temperature for the cathode dissolution side reaction
  double lam3k_T; //!< activation energy for lam3k

  double lam4p; //!< fitting parameter for the positive electrode for the 4th LAM model
  double lam4n; //!< fitting parameter for the negative electrode for the 4th LAM model

  [[nodiscard]] constexpr auto begin() noexcept { return &lam1p; }
  [[nodiscard]] constexpr auto end() noexcept { return &lam4n + 1; }

  auto &operator*=(double a)
  {
    for (auto &item : *this)
      item *= a;

    return *this;
  }
};
} // namespace slide::param