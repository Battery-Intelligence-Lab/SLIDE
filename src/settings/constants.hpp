/*
 * constants.hpp
 *
 * Author : Volkan Kumtepeli
 *
 * Defines constants to be used in the program.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

namespace slide {
constexpr double TIME_INF = 20 * 365 * 24 * 3600; // Very long time for no time limit.
}

namespace slide::PhyConst {
constexpr double Kelvin = 273.0;
constexpr double F = 96487;  //!< Faraday's constant
constexpr double Rg = 8.314; //!< ideal gas constant
constexpr double pi = 3.141592653589793;
} // namespace slide::PhyConst