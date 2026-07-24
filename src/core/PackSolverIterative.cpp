/**
 * @file PackSolverIterative.cpp
 * @brief The two matrix-free pack strategies: the resistor-free ladder and PI relaxation.
 *
 * Owns: `PackSolver::solveLadder` and `PackSolver::solveRelaxation`. Implements PLAN.md §3.4.
 * Hot: both run inside the pack step. Neither touches `SolverWorkspace::Impl` -- they never
 * assemble or factorise a matrix, which is exactly what separates them from the sparse Newton
 * solve in `PackSolver.cpp`; the shared checked arithmetic lives in `PackSolverInternal.hpp`.
 * @surface internal
 */

#include "PackSolver.hpp"
#include "PackSolverInternal.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

namespace slide::core {

using detail::addCompensatedFinite;
using detail::addFinite;
using detail::branchAffine;
using detail::branchCurrentNumerator;
using detail::branchCurrentOut;
using detail::branchDrop;
using detail::cellCurrentFromDrop;

slide::Status PackSolver::solveLadder(real_t applied_current)
{
  const auto &netlist = topology_.electrical;
  // configure() validates the ladder and solveImpl() dispatches here only for
  // a ladder-compatible mode.
  assert(netlist.series_parallel_ladder && netlist.ladder_offsets.size() >= 2
         && netlist.ladder_nodes.size() == netlist.ladder_offsets.size());
  const std::size_t layers = netlist.ladder_offsets.size() - 1;
  assert(layer_voltage_.size() == layers);
  for (std::size_t layer = 0; layer < layers; ++layer) {
    real_t conductance_sum{};
    real_t source_sum{};
    for (std::uint32_t i = netlist.ladder_offsets[layer];
         i < netlist.ladder_offsets[layer + 1];
         ++i) {
      const auto cell = netlist.ladder_cells[i];
      const real_t conductance = 1.0 / resistance_[cell];
      const real_t source = ocv_[cell] * conductance;
      if (!is_strictly_positive_finite(conductance)
          || !is_finite(source)
          || !addFinite(conductance_sum, conductance)
          || !addFinite(source_sum, source))
        return slide::Status::Invalid_states;
    }
    // Validated offsets give every layer at least one cell; accepted
    // conductances are positive and every addition above remained finite.
    assert(is_finite(conductance_sum) && conductance_sum > 0.0);
    const real_t numerator = source_sum - applied_current;
    if (!is_finite(numerator))
      return slide::Status::Invalid_states;
    layer_voltage_[layer] = numerator / conductance_sum;
    if (!is_finite(layer_voltage_[layer]))
      return slide::Status::Invalid_states;
    for (std::uint32_t i = netlist.ladder_offsets[layer];
         i < netlist.ladder_offsets[layer + 1];
         ++i) {
      const auto cell = netlist.ladder_cells[i];
      if (!cellCurrentFromDrop(layer_voltage_[layer],
                               ocv_[cell],
                               resistance_[cell],
                               candidate_current_[cell]))
        return slide::Status::Invalid_states;
    }
  }
  std::fill(candidate_node_voltage_.begin(), candidate_node_voltage_.end(), 0.0);
  real_t voltage{};
  for (std::size_t reverse = layers; reverse > 0; --reverse) {
    if (!addFinite(voltage, layer_voltage_[reverse - 1]))
      return slide::Status::Invalid_states;
    candidate_node_voltage_[netlist.ladder_nodes[reverse - 1]] = voltage;
  }
  candidate_terminal_voltage_ = voltage;
  assert(is_finite(candidate_terminal_voltage_));
  return slide::Status::Success;
}

slide::Status PackSolver::solveRelaxation(real_t applied_current)
{
  const auto &netlist = topology_.electrical;
  // configure() starts with zeros and only a fully finite candidate is ever
  // published for a warm start.
  assert(std::all_of(candidate_node_voltage_.begin(),
                     candidate_node_voltage_.end(),
                     [](const real_t &value) { return is_finite(value); }));
  std::fill(relaxation_.diagonal.begin(), relaxation_.diagonal.end(), 0.0);
  std::fill(relaxation_.diagonal_compensation.begin(),
            relaxation_.diagonal_compensation.end(),
            0.0);
  std::fill(relaxation_.rhs.begin(), relaxation_.rhs.end(), 0.0);
  std::fill(relaxation_.rhs_compensation.begin(),
            relaxation_.rhs_compensation.end(),
            0.0);
  auto stamp = [&](const CompiledElectricalBranch &branch, real_t resistance, real_t source) {
    const real_t conductance = 1.0 / resistance;
    if (!is_strictly_positive_finite(conductance)
        || !is_finite(source))
      return false;
    const auto p = branch.node_positive;
    const auto n = branch.node_negative;
    const real_t positive_source = candidate_node_voltage_[n] + source;
    const real_t negative_source = candidate_node_voltage_[p] - source;
    if (!is_finite(positive_source) || !is_finite(negative_source))
      return false;
    const real_t positive_rhs = conductance * positive_source;
    const real_t negative_rhs = conductance * negative_source;
    return is_finite(positive_rhs) && is_finite(negative_rhs)
           && addCompensatedFinite(relaxation_.diagonal[p],
                                   relaxation_.diagonal_compensation[p],
                                   conductance)
           && addCompensatedFinite(relaxation_.diagonal[n],
                                   relaxation_.diagonal_compensation[n],
                                   conductance)
           && addCompensatedFinite(relaxation_.rhs[p],
                                   relaxation_.rhs_compensation[p],
                                   positive_rhs)
           && addCompensatedFinite(relaxation_.rhs[n],
                                   relaxation_.rhs_compensation[n],
                                   negative_rhs);
  };
  for (const auto &branch : netlist.branches) {
    const auto affine = branchAffine(branch, ocv_, resistance_);
    // Linearization and netlist validation established this immediately
    // before entering the numeric kernel.
    assert(is_finite(affine.resistance) && affine.resistance > 0.0);
    if (!stamp(branch, affine.resistance, affine.source))
      return slide::Status::Invalid_states;
  }
  if (!addCompensatedFinite(relaxation_.rhs[netlist.terminal_positive],
                            relaxation_.rhs_compensation[netlist.terminal_positive],
                            -applied_current)
      || !addCompensatedFinite(
        relaxation_.rhs[netlist.terminal_negative],
        relaxation_.rhs_compensation[netlist.terminal_negative],
        applied_current))
    return slide::Status::Invalid_states;

  for (std::uint32_t node = 0; node < netlist.node_count; ++node) {
    if (node == netlist.terminal_negative) {
      relaxation_.target[node] = real_t{};
      continue;
    }
    if (!(is_strictly_positive_finite(relaxation_.diagonal[node])
          && is_finite(relaxation_.rhs[node])))
      return slide::Status::Numerical_failure;
    relaxation_.target[node] = relaxation_.rhs[node]
                               / relaxation_.diagonal[node];
    if (!is_finite(relaxation_.target[node]))
      return slide::Status::Invalid_states;
  }
  for (std::uint32_t node = 0; node < netlist.node_count; ++node)
    if (node != netlist.terminal_negative) {
      const real_t difference = relaxation_.target[node]
                                - candidate_node_voltage_[node];
      const real_t correction = relaxation_alpha_ * difference;
      if (!is_finite(difference) || !is_finite(correction)
          || !addFinite(candidate_node_voltage_[node], correction))
        return slide::Status::Invalid_states;
    }

  for (const auto &branch : netlist.branches)
    if (branch.kind == ElectricalBranchKind::cell) {
      const real_t voltage = branchDrop(branch, candidate_node_voltage_);
      if (!cellCurrentFromDrop(voltage,
                               ocv_[branch.cell],
                               resistance_[branch.cell],
                               candidate_current_[branch.cell]))
        return slide::Status::Invalid_states;
    }
  candidate_terminal_voltage_ = candidate_node_voltage_[netlist.terminal_positive]
                                - candidate_node_voltage_[netlist.terminal_negative];
  // terminal_negative is the skipped reference target and remains zero.
  assert(candidate_node_voltage_[netlist.terminal_negative] == 0.0);
  assert(is_finite(candidate_terminal_voltage_));

  std::fill(relaxation_.residual.begin(), relaxation_.residual.end(), 0.0);
  std::fill(relaxation_.residual_compensation.begin(),
            relaxation_.residual_compensation.end(),
            0.0);
  real_t roundoff_operation_scale = std::abs(applied_current);
  real_t roundoff_current_scale = std::abs(applied_current);
  for (const auto &branch : netlist.branches) {
    const real_t voltage = branchDrop(branch, candidate_node_voltage_);
    if (!is_finite(voltage))
      return slide::Status::Invalid_states;
    const auto affine = branchAffine(branch, ocv_, resistance_);
    const real_t numerator =
      branchCurrentNumerator(voltage, affine.source);
    // Cell drops used the opposite subtraction above, and resistor sources
    // are positive zero; no node or source changed in between.
    assert(is_finite(numerator));
    const real_t branch_current =
      branchCurrentOut(numerator, affine.resistance);
    real_t operation_scale{};
    if (!addFinite(operation_scale,
                   std::abs(candidate_node_voltage_[branch.node_positive]))
        || !addFinite(operation_scale,
                      std::abs(candidate_node_voltage_[branch.node_negative]))
        || !addFinite(operation_scale, std::abs(affine.source)))
      return slide::Status::Invalid_states;
    operation_scale /= affine.resistance;
    const real_t current_magnitude = std::abs(branch_current);
    if (!is_finite(branch_current) || !is_finite(operation_scale)
        || !is_finite(current_magnitude)
        || !addFinite(roundoff_operation_scale, operation_scale)
        || !addFinite(roundoff_current_scale, current_magnitude)
        || !addCompensatedFinite(
          relaxation_.residual[branch.node_positive],
          relaxation_.residual_compensation[branch.node_positive],
          branch_current)
        || !addCompensatedFinite(
          relaxation_.residual[branch.node_negative],
          relaxation_.residual_compensation[branch.node_negative],
          -branch_current))
      return slide::Status::Invalid_states;
  }
  if (!addCompensatedFinite(
        relaxation_.residual[netlist.terminal_positive],
        relaxation_.residual_compensation[netlist.terminal_positive],
        applied_current)
      || !addCompensatedFinite(
        relaxation_.residual[netlist.terminal_negative],
        relaxation_.residual_compensation[netlist.terminal_negative],
        -applied_current))
    return slide::Status::Invalid_states;
  real_t drift{};
  for (std::uint32_t node = 0; node < netlist.node_count; ++node) {
    if (node == netlist.terminal_negative)
      continue;
    const real_t residual = relaxation_.residual[node];
    const real_t magnitude = std::abs(residual);
    assert(is_finite(magnitude)); // abs preserves finiteness
    drift = std::max(drift, magnitude);
  }
  // This is a conservative diagnostic bound, not solver state. If its bound
  // arithmetic is itself unrepresentable, saturation remains conservative.
  diagnostics_.constraint_bound = std::max(
    diagnostics_.constraint_bound,
    detail::conservativePackRoundoffBound(
      static_cast<real_t>(netlist.branches.size() + 1)
        * std::numeric_limits<real_t>::epsilon(),
      roundoff_current_scale,
      roundoff_operation_scale));
  diagnostics_.constraint_drift = drift;
  return slide::Status::Success;
}


} // namespace slide::core
