/*
 * SEIparam.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

namespace slide::param {
struct SEIparam
{
  using value_type = double;
  value_type sei_porosity{ 3 * 7.5e-7 }; //!< proportionality constant between the SEI growth and the decrease in volume fraction of active material
                                         //!< (because the SEI layer blocks the pores at the surface)

  value_type sei1k{ 0.075e-14 }; //!< rate parameter of the SEI side reaction in the 1st SEI model
  value_type sei1k_T{ 130e3 };   //!< activation energy of sei1k

  value_type sei2k{ 2.75e-11 };  //!< rate parameter of the SEI side reaction in the 2nd SEI model
  value_type sei2k_T{ 130e3 };   //!< activation energy of sei2k
  value_type sei2D{ 1.125e-14 }; //!< diffusion constant of the SEI layer in the 2nd SEI model
  value_type sei2D_T{ 20e3 };    //!< activation energy of sei2D

  value_type sei3k{ 1.1458e-15 }; //!< rate parameter of the SEI side reaction in the 3rd SEI model
  value_type sei3k_T{ 65e3 };     //!< activation energy of sei3k
  value_type sei3D{ 0.25e-15 };   //!< diffusion constant of the SEI layer in the 3rd SEI model
  value_type sei3D_T{ 200e3 };    //!< activation energy of sei3D

  value_type sei4k{ 3.75e-15 };         //!< rate parameter of the SEI side reaction in the 3rd SEI model
  value_type sei4k_T{ 130000.0 / 1.5 }; //!< activation energy of sei3k
  value_type sei4D{ 0.5e-16 / 15.0 };   //!< diffusion constant of the SEI layer in the 3rd SEI model
  value_type sei4D_T{ 80000 };          //!< activation energy of sei3D

  auto begin() noexcept { return &sei_porosity; }
  auto end() noexcept { return &sei4D_T + 1; }

  auto &operator*=(double a)
  {
    for (auto &item : *this)
      item *= a;

    return *this;
  }
};


} // namespace slide::param
