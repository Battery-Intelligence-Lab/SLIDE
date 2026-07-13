/**
 * @file PackSolverInternal.hpp
 * @brief Internal numeric helpers shared only with focused solver tests.
 * @details Owns M0.3 diagnostic-bound saturation; hot, allocation-free, and
 *          deliberately excluded from the public PackSolver surface.
 * @surface internal
 */

#pragma once

#include "PackSolver.hpp"
#include "Numeric.hpp"

#include <algorithm>
#include <span>

namespace slide::core::detail {

[[nodiscard]] real_t conservativePackRoundoffBound(real_t accumulation_ratio,
                                                   real_t current_scale,
                                                   real_t operation_scale) noexcept;

// Shared by the direct sparse solve and by the two matrix-free strategies. `addCompensatedFinite`
// keeps its strict-FP pragmas with its definition: the Kahan evaluation order is the contract,
// and the Release configuration's fast-math would otherwise reassociate it away.
inline bool knownSolveMode(PackSolveMode mode)
{
  switch (mode) {
  case PackSolveMode::sparse_newton:
  case PackSolveMode::ladder:
  case PackSolveMode::relaxation:
    return true;
  }
  return false;
}

inline bool addFinite(real_t &target, real_t increment)
{
  if (!is_finite(increment))
    return false;
  const real_t updated = target + increment;
  if (!is_finite(updated))
    return false;
  target = updated;
  return true;
}

#if defined(_MSC_VER) && !defined(__clang__)
#pragma float_control(precise, on, push)
#endif
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-fast-math")))
#endif
inline bool addCompensatedFinite(real_t &sum, real_t &compensation, real_t value)
{
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
  if (!is_finite(value))
    return false;
  // Volatile stages preserve Kahan's evaluation order under the Release
  // configuration's finite-math/reassociation flags.
  volatile real_t adjusted = value - compensation;
  volatile real_t updated = sum + adjusted;
  volatile real_t next_compensation = (updated - sum) - adjusted;
  if (!is_finite(adjusted) || !is_finite(updated)
      || !is_finite(next_compensation))
    return false;
  sum = updated;
  compensation = next_compensation;
  return true;
}
#if defined(_MSC_VER) && !defined(__clang__)
#pragma float_control(pop)
#endif

inline bool finiteCandidate(std::span<const real_t> current,
                     std::span<const real_t>
                       node_voltage,
                     real_t terminal_voltage)
{
  if (!is_finite(terminal_voltage))
    return false;
  const auto finite = [](const real_t &value) { return is_finite(value); };
  return std::all_of(current.begin(), current.end(), finite)
         && std::all_of(node_voltage.begin(), node_voltage.end(), finite);
}

} // namespace slide::core::detail
