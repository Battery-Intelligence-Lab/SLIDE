/**
 * @file Dual.hpp
 * @brief Dependency-free first-order dual number for forward sensitivities.
 */

#pragma once

#include "Numeric.hpp"

#include <cmath>

namespace slide::core {

struct Dual
{
  real_t value{};
  real_t derivative{};

  constexpr Dual() = default;
  constexpr Dual(real_t primal) : value{ primal } {}
  constexpr Dual(real_t primal, real_t tangent)
    : value{ primal }, derivative{ tangent }
  {}

  constexpr explicit operator real_t() const { return value; }

  constexpr Dual &operator+=(Dual other)
  {
    value += other.value;
    derivative += other.derivative;
    return *this;
  }
  constexpr Dual &operator-=(Dual other)
  {
    value -= other.value;
    derivative -= other.derivative;
    return *this;
  }
  constexpr Dual &operator*=(Dual other)
  {
    derivative = derivative * other.value + value * other.derivative;
    value *= other.value;
    return *this;
  }
  constexpr Dual &operator/=(Dual other)
  {
    derivative = (derivative * other.value - value * other.derivative)
                 / (other.value * other.value);
    value /= other.value;
    return *this;
  }
};

constexpr Dual operator+(Dual left, Dual right) { return left += right; }
constexpr Dual operator-(Dual left, Dual right) { return left -= right; }
constexpr Dual operator*(Dual left, Dual right) { return left *= right; }
constexpr Dual operator/(Dual left, Dual right) { return left /= right; }
constexpr Dual operator-(Dual value)
{
  return { -value.value, -value.derivative };
}

inline Dual exp(Dual value)
{
  const real_t result = std::exp(value.value);
  return { result, result * value.derivative };
}

inline Dual expm1(Dual value)
{
  return { std::expm1(value.value), std::exp(value.value) * value.derivative };
}

inline Dual sqrt(Dual value)
{
  const real_t result = std::sqrt(value.value);
  return { result, value.derivative / (2.0 * result) };
}

/** Differentiate a variable base raised to a scalar exponent. */
inline Dual pow(Dual value, real_t exponent)
{
  if (exponent == 0.0)
    return { 1.0, 0.0 };
  if (exponent == 1.0)
    return value;
  const real_t result = std::pow(value.value, exponent);
  return { result,
           exponent * std::pow(value.value, exponent - 1.0)
             * value.derivative };
}

inline Dual asinh(Dual value)
{
  return { std::asinh(value.value),
           value.derivative / std::sqrt(1.0 + value.value * value.value) };
}

inline Dual abs(Dual value)
{
  return value.value < 0.0 ? -value : value;
}

} // namespace slide::core
