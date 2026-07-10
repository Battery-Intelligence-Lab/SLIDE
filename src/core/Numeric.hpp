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
inline bool is_finite(real_t value) noexcept
{
  constexpr std::uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
  // The volatile integer barrier is intentional. Clang otherwise recognises even an integer
  // bit_cast classification and folds it to true under -ffinite-math-only (part of -Ofast).
  volatile std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  return (bits & exponent_mask) != exponent_mask;
}

template <class Real>
inline bool is_finite_primal(const Real &value) noexcept
{
  return is_finite(primal_value(value));
}

} // namespace slide::core
