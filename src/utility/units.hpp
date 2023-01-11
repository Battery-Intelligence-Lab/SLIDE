/*
 * units.hpp
 *
 * Definition of different units
 *  Created on: 30 Jun 2022
 *   Author(s): Volkan Kumtepeli
 */

#pragma once

#include "../settings/constants.hpp"

namespace slide {

// Current:
constexpr auto operator""_kA(long double d) { return 1e3 * d; }  //!< kiloAmperes
constexpr auto operator""_A(long double d) { return d; }         //!< Amperes
constexpr auto operator""_mA(long double d) { return d * 1e-3; } //!< miliAmperes
constexpr auto operator""_uA(long double d) { return d * 1e-6; } //!< microAmperes

// Voltage:
constexpr auto operator""_kV(long double d) { return 1e3 * d; }  //!< kiloVolts
constexpr auto operator""_V(long double d) { return d; }         //!< Volts
constexpr auto operator""_mV(long double d) { return d * 1e-3; } //!< miliVolts
constexpr auto operator""_uV(long double d) { return d * 1e-6; } //!< microVolts


// Temperature:
constexpr inline double C_to_Kelvin(auto Celsius) { return static_cast<double>(PhyConst::Kelvin + Celsius); } // #TODO check if it better to have auto or double
constexpr inline double K_to_Celsius(auto Kelvin) { return static_cast<double>(Kelvin - PhyConst::Kelvin); }

constexpr double operator""_degC(long double d) { return C_to_Kelvin(d); }      //!< degrees Celsius
constexpr double operator""_K(long double d) { return static_cast<double>(d); } //!< Kelvins

constexpr double operator""_degC(size_t d) { return C_to_Kelvin(d); }      //!< degrees Celsius
constexpr double operator""_K(size_t d) { return static_cast<double>(d); } //!< Kelvins


} // namespace slide
