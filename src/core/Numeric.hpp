/**
 * @file Numeric.hpp
 * @brief Scalar-generic numeric helpers that remain valid under fast-math builds.
 */

#pragma once

#include "StateArena.hpp"

#include <bit>
#include <cstdint>

namespace slide::core {

template <class Real>
constexpr real_t primal_value(const Real &value)
{
  return static_cast<real_t>(value);
}

/**
 * IEEE-754 finite check that clang cannot discard under `-ffast-math`/`-Ofast`.
 * `std::isfinite` is not a validation primitive in those builds because the compiler is
 * permitted to assume NaN and infinity never occur and fold the check to true.
 */
inline bool is_finite(const real_t &value) noexcept
{
  constexpr std::uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
  // Both the opaque reference boundary and volatile integer barrier are intentional. Clang
  // can otherwise attach `nofpclass` to a by-value argument and fold even an integer bit-cast
  // classification to true under -ffinite-math-only (part of -Ofast).
  volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  return (bits & exponent_mask) != exponent_mask;
}

inline bool is_finite(const volatile real_t &value) noexcept
{
  const real_t copied = value;
  return is_finite(copied);
}

/**
 * IEEE-754 strictly-positive finite check with the same fast-math barriers as
 * is_finite(). This is required for derived coefficients: under an FTZ
 * environment, a mathematically positive reciprocal can materialise as +0,
 * while a floating `value > 0` check can be removed by finite-math reasoning.
 */
inline bool is_strictly_positive_finite(const real_t &value) noexcept
{
  constexpr std::uint64_t sign_mask = UINT64_C(0x8000000000000000);
  constexpr std::uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
  constexpr std::uint64_t magnitude_mask = UINT64_C(0x7fffffffffffffff);
  volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  return (bits & sign_mask) == 0 && (bits & exponent_mask) != exponent_mask
         && (bits & magnitude_mask) != 0;
}

template <class Real>
inline bool is_finite_primal(const Real &value) noexcept
{
  return is_finite(primal_value(value));
}

} // namespace slide::core
