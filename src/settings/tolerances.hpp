/*
 * tolerances.hpp
 *
 * Author : Jorn Reniers, Volkan Kumtepeli
 *
 * Defines constants to be used in the program.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

namespace slide {
namespace tol {
  //!< Tolerance values for different functions.
  constexpr double findCVcurrent = 0.001;     //!< Relative tolerance for finding CV phase current.
  constexpr double findCVcurrent_bin = 0.001; //!< Relative tolerance for finding CV phase current (False position method)
} // namespace tol

} // namespace slide

namespace slide::settings //!< #TODO if it is in slide::settings it gives 600+ errors. Why?
{
constexpr double MODULE_P_V_ABSTOL = 0.01;   //!< tolerance on the difference between parallel-connected cells
                                             //!< i.e. the difference in voltage between two parallel connected cell has to be below this value
constexpr double MODULE_P_V_RELTOL = 0.001;  //!< tolerance on the relative difference between parallel-connected cells
constexpr double MODULE_P_I_ABSTOL = 0.005;  //!< absolute tolerance on the difference between the set current and real current in a parallel modules
                                             //!< i.e. in Module_p::setCurrent(), Inew - I < Itol
constexpr double MODULE_P_I_RELTOL = 0.0005; //!< relative tolerance on the difference between the set current and real current in a parallel modules
                                             //!< i.e. in Module_p::setCurrent(), abs(Inew - I) < Itol*Inew
} // namespace slide::settings
