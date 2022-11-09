/*
 * LAMparam_default.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

#include "LAMparam.hpp"

namespace slide::param::def {
constexpr LAMparam LAMparam_Kokam{
  3.4985e-9 / 10.0, //!< lam1p
  5.88e-13 * 2e4,   //!< lam1n
  -1.675e-5,        //!< lam2ap
  0.0,              //!< lam2bp
  -1.675e-5,        //!< lam2an
  0.0,              //!< am2bn
  54611.0,          //!< lam2t
  12.5e-6,          //!< lam3k
  27305.0,          //!< lam3k_T
  8.3333e-9,        //!< lam4p
  8.3333e-9         //!< lam4n
};

constexpr LAMparam LAMparam_LGCChemNMC{
  2.6031e-9,   //!< lam1p
  1.0417e-12,  //!< lam1n
  3.015e-11,   //!< lam2ap
  -1.72125e-6, //!< lam2bp
  3.015e-11,   //!< lam2an
  -1.72125e-6, //!< lam2bn
  54611.0,     //!< lam2t
  1.21e-6,     //!< lam3k
  27305.0,     //!< lam3k_T
  7.5e-9,      //!< lam4p
  7.5e-9       //!< lam4n
};

constexpr LAMparam LAMparam_User{
  3.4985e-9, //!< lam1p
  5.88e-13,  //!< lam1n
  -1.675e-5, //!< lam2ap
  0.0,       //!< lam2bp
  -1.675e-5, //!< lam2an
  0.0,       //!< am2bn
  54611.0,   //!< lam2t
  12.5e-6,   //!< lam3k
  27305.0,   //!< lam3k_T
  8.3333e-9, //!< lam4p
  8.3333e-9  //!< lam4n
};

} // namespace slide::param::def