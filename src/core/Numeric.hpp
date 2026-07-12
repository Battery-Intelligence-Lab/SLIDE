/**
 * @file Numeric.hpp
 * @brief Scalar-generic numeric helpers that remain valid under fast-math builds.
 */

#pragma once

#include "StateArena.hpp"

#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>

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

/**
 * Multiply two finite, nonnegative scalars without ever evaluating an
 * overflowing product or an overflowing overflow guard.  `maximum / right`
 * is safe only when `right > 1`; when either operand is at most one, their
 * product cannot exceed the other finite operand.  Equality at the rounded
 * overflow boundary is rejected conservatively.
 */
inline bool try_multiply_nonnegative(const real_t &left,
                                     const real_t &right,
                                     real_t &product) noexcept
{
  if (!(is_finite(left) && left >= 0.0
        && is_finite(right) && right >= 0.0))
    return false;

  constexpr real_t maximum = std::numeric_limits<real_t>::max();
  if (left > 1.0 && right > 1.0 && left >= maximum / right)
    return false;

  const real_t candidate = left * right;
  product = candidate;
  return true;
}

template <class Real>
inline bool is_finite_primal(const Real &value) noexcept
{
  if constexpr (std::is_same_v<std::remove_cvref_t<Real>, real_t>) {
    // Preserve the opaque reference boundary for the native scalar.  Routing a double
    // through primal_value() creates a by-value SSA edge on which ThinLTO may attach
    // `nofpclass` under finite-math, defeating the integer classifier in is_finite().
    return is_finite(value);
  } else {
    const real_t primal = primal_value(value);
    return is_finite(primal);
  }
}

} // namespace slide::core
