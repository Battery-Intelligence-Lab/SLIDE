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
constexpr inline auto C_to_Kelvin(double Celsius) { return PhyConst::Kelvin + Celsius; }
constexpr inline auto K_to_Celsius(double Kelvin) { return Kelvin - PhyConst::Kelvin; }

constexpr auto operator""_degC(long double d) { return C_to_Kelvin(d); } //!< degrees Celsius
constexpr auto operator""_K(long double d) { return d; }                 //!< Kelvins


} // namespace slide
