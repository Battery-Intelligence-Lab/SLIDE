/*
 * SEIparam_default.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

#include "SEIparam.hpp"

namespace slide::param::def {
//!< fitting parameters of the models
constexpr SEIparam SEIparam_Kokam{
  3 * 7.5e-7, //!< sei_porosity

  0.075e-14, //!< sei1k
  130e3,     //!< sei1k_T

  2.75e-11,  //!< sei2k
  130e3,     //!< sei2k_T
  1.125e-14, //!< sei2D
  20e3,      //!< sei2D_T

  1.1458e-15, //!< sei3k
  65e3,       //!< sei3k_T
  0.25e-15,   //!< sei3D
  200e3,      //!< sei3D_T

  1.5 * 0.5e-14 / 2.0, //!< sei4k -> ageing fit.
  130000.0 / 1.5,      //!< sei4k_T -> ageing fit.
  0.5e-16 / 15.0,      //!< sei4D -> ageing fit.
  200000.0 / 2.5       //!< sei4D_T -> ageing fit.
};

constexpr SEIparam SEIparam_Kokam_old{
  7.5e-7,     //!< sei_porosity
  0.075e-14,  //!< sei1k
  130e3,      //!< sei1k_T
  2.75e-11,   //!< sei2k
  130e3,      //!< sei2k_T
  1.125e-14,  //!< sei2D
  20e3,       //!< sei2D_T
  1.1458e-15, //!< sei3k
  65e3,       //!< sei3k_T
  0.25e-15,   //!< sei3D
  200e3       //!< sei3D_T
};

constexpr SEIparam SEIparam_LGCChemNMC{
  7.5e-7,    //!< sei_porosity
  0.075e-14, //!< sei1k
  130e3,     //!< sei1k_T
  2.75e-11,  //!< sei2k
  130e3,     //!< sei2k_T
  2.5e-15,   //!< sei2D
  200e3,     //!< sei2D_T
  1e-11,     //!< sei3k
  0.0,       //!< sei3k_T
  1.05e-16,  //!< sei3D
  20e3       //!< sei3D_T
};

constexpr SEIparam SEIparam_User{
  7.5e-7,     //!< sei_porosity
  0.075e-14,  //!< sei1k
  130e3,      //!< sei1k_T
  2.75e-11,   //!< sei2k
  130e3,      //!< sei2k_T
  1.125e-14,  //!< sei2D
  20e3,       //!< sei2D_T
  1.1458e-15, //!< sei3k
  65e3,       //!< sei3k_T
  0.25e-15,   //!< sei3D
  200e3       //!< sei3D_T
};
} // namespace slide::param::def