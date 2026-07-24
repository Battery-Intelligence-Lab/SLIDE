/**
 * @file PackSolverInternal.hpp
 * @brief Internal numeric helpers shared only with focused solver tests.
 * @details Owns PLAN.md §3.4 affine branch/current algebra, checked hot
 *          arithmetic, and M0.3 diagnostic-bound saturation; allocation-free
 *          and deliberately excluded from the public PackSolver surface.
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

[[nodiscard]] inline real_t branchDrop(
  const CompiledElectricalBranch &branch,
  std::span<const real_t>
    node_voltage) noexcept
{
  return node_voltage[branch.node_positive]
         - node_voltage[branch.node_negative];
}

struct BranchAffine
{
  real_t resistance{};
  real_t source{};
};

[[nodiscard]] inline BranchAffine branchAffine(
  const CompiledElectricalBranch &branch,
  std::span<const real_t>
    cell_ocv,
  std::span<const real_t>
    cell_resistance) noexcept
{
  if (branch.kind == ElectricalBranchKind::cell)
    return { .resistance = cell_resistance[branch.cell],
             .source = cell_ocv[branch.cell] };
  return { .resistance = branch.resistance, .source = real_t{} };
}

[[nodiscard]] inline real_t branchCurrentNumerator(
  const real_t &drop,
  const real_t &source) noexcept
{
  return drop - source;
}

[[nodiscard]] inline real_t branchCurrentOut(
  const real_t &numerator,
  const real_t &resistance) noexcept
{
  return numerator / resistance;
}

[[nodiscard]] inline bool cellCurrentFromDrop(
  const real_t &drop,
  const real_t &ocv,
  const real_t &resistance,
  real_t &current) noexcept
{
  const real_t numerator = ocv - drop;
  if (!is_finite(drop) || !is_finite(numerator))
    return false;
  current = numerator / resistance;
  return is_finite(current);
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

// Keep these strict-FP pragmas with the compensated finite sum: the Kahan
// evaluation order is the contract, and fast math would reassociate it.
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
