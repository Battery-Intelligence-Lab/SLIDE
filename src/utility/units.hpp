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
constexpr double operator""_kA(long double d) { return 1e3 * d; }  //!< kiloAmperes
constexpr double operator""_A(long double d) { return d; }         //!< Amperes
constexpr double operator""_mA(long double d) { return d * 1e-3; } //!< miliAmperes
constexpr double operator""_uA(long double d) { return d * 1e-6; } //!< microAmperes

// Voltage:
constexpr double operator""_kV(long double d) { return 1e3 * d; }  //!< kiloVolts
constexpr double operator""_V(long double d) { return d; }         //!< Volts
constexpr double operator""_mV(long double d) { return d * 1e-3; } //!< miliVolts
constexpr double operator""_uV(long double d) { return d * 1e-6; } //!< microVolts


// Temperature:
constexpr inline double C_to_Kelvin(double Celsius) { return PhyConst::Kelvin + Celsius; }
constexpr inline double K_to_Celsius(double Kelvin) { return Kelvin - PhyConst::Kelvin; }

constexpr double operator""_degC(long double d) { return C_to_Kelvin(d); } //!< degrees Celsius
constexpr double operator""_K(long double d) { return d; }                 //!< Kelvins


} // namespace slide
