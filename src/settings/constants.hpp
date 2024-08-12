/**
 * @file constants.hpp
 * @brief Defines constants to be used in the program.
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 11 Jun 2020
 */

#pragma once

#include <numbers>

namespace slide {
inline constexpr double TIME_INF = 20 * 365 * 24 * 3600; // Very long time for no time limit.
inline constexpr double M_PI = std::numbers::pi;
} // namespace slide

namespace slide::PhyConst {
inline constexpr double Kelvin = 273.0;
inline constexpr double F = 96487;     //!< Faraday's constant
inline constexpr double Rg = 8.314;    //!< ideal gas constant
inline constexpr double C_elec = 1000; //!< Li- concentration in electrolyte [mol m-3] standard concentration of 1 molar
inline constexpr double n = 1;         //!< number of electrons involved in the main reaction [-] #TODO if really constant?

} // namespace slide::PhyConst