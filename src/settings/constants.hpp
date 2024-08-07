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
inline constexpr double F = 96487;  //!< Faraday's constant
inline constexpr double Rg = 8.314; //!< ideal gas constant
} // namespace slide::PhyConst