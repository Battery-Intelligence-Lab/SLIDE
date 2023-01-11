/*
 * CSparam.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */
#pragma once

namespace slide::param {
//!< Define a structure with the fitting parameters of the surface crack growth models (CS)
struct CSparam
{
  using value_type = double;
  value_type CS1alpha{ 4.25e-5 };   //!< fitting parameter of the 1st surface crack growth model
  value_type CS2alpha{ 6.3e-7 };    //!< fitting parameter of the 2nd surface crack growth model
  value_type CS3alpha{ 2.31e-16 };  //!< fitting parameter of the 3rd surface crack growth model
  value_type CS4alpha{ 4.3306e-8 }; //!< fitting parameter of the 4th surface crack growth model

  value_type CS4Amax; //!< maximum crack growth surface for the 4th surface crack growth model 5 * getAn() * getThickn() * elec_surf

  value_type CS5k{ 1e-18 };     //!< rate parameter of the 5th surface crack growth model at reference temperature
  value_type CS5k_T{ -127040 }; //!< activation energy of CS5k

  value_type CS_diffusion{ 2 }; //!< fitting parameter to decrease the diffusion constant due to surface cracks

  [[nodiscard]] constexpr auto begin() noexcept { return &CS1alpha; }
  [[nodiscard]] constexpr auto end() noexcept { return &CS_diffusion + 1; }

  auto &operator*=(double a)
  {
    for (auto &item : *this)
      item *= a;

    return *this;
  }
};
} // namespace slide::param