/**
 * @file PackSolver.cpp
 * @brief Affine Thevenin pack solving with persistent sparse and ladder workspaces.
 */

#include "PackSolver.hpp"

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace slide::core {

slide::Status TheveninBatchView::linearize(std::span<const real_t> current,
                                           std::span<real_t>
                                             ocv,
                                           std::span<real_t>
                                             resistance) const
{
  if (!valid() || static_cast<int>(current.size()) != lanes_
      || ocv.size() != current.size() || resistance.size() != current.size())
    return slide::Status::Invalid_parameters;
  return linearize_(object_, current, ocv, resistance);
}

slide::Status PackTheveninSystem::configure(
  std::span<const CompiledCell> cells,
  std::span<const std::string>
    batch_archetypes,
  std::span<const TheveninBatchView>
    batches)
{
  if (batches.size() != batch_archetypes.size() || cells.empty())
    return slide::Status::Invalid_parameters;
  std::vector<int> required_lanes(batches.size());
  for (const auto &cell : cells) {
    if (cell.location.batch >= batches.size())
      return slide::Status::Invalid_parameters;
    auto &lanes = required_lanes[cell.location.batch];
    lanes = std::max(lanes, static_cast<int>(cell.location.lane + 1));
  }
  std::vector<BatchScratch> scratch(batches.size());
  for (std::size_t batch = 0; batch < batches.size(); ++batch) {
    if (!batches[batch].valid() || batches[batch].lanes() != required_lanes[batch])
      return slide::Status::Invalid_parameters;
    scratch[batch].view = batches[batch];
    scratch[batch].current.resize(static_cast<std::size_t>(required_lanes[batch]));
    scratch[batch].ocv.resize(static_cast<std::size_t>(required_lanes[batch]));
    scratch[batch].resistance.resize(static_cast<std::size_t>(required_lanes[batch]));
  }
  cells_.assign(cells.begin(), cells.end());
  batches_ = std::move(scratch);
  return slide::Status::Success;
}

slide::Status PackTheveninSystem::linearize(std::span<const real_t> cell_current,
                                            std::span<real_t>
                                              cell_ocv,
                                            std::span<real_t>
                                              cell_resistance)
{
  if (cell_current.size() != cells_.size() || cell_ocv.size() != cells_.size()
      || cell_resistance.size() != cells_.size())
    return slide::Status::Invalid_parameters;
  for (std::size_t cell = 0; cell < cells_.size(); ++cell) {
    const auto location = cells_[cell].location;
    batches_[location.batch].current[location.lane] = cell_current[cell];
  }
  for (auto &batch : batches_) {
    const auto status = batch.view.linearize(batch.current, batch.ocv, batch.resistance);
    if (status != slide::Status::Success)
      return status;
  }
  for (std::size_t cell = 0; cell < cells_.size(); ++cell) {
    const auto location = cells_[cell].location;
    const auto ocv = batches_[location.batch].ocv[location.lane];
    const auto resistance = batches_[location.batch].resistance[location.lane];
    if (!is_finite(ocv) || !is_finite(resistance) || !(resistance > 0.0))
      return slide::Status::Invalid_states;
    cell_ocv[cell] = ocv;
    cell_resistance[cell] = resistance;
  }
  return slide::Status::Success;
}

struct SolverWorkspace::Impl
{
  using Matrix = Eigen::SparseMatrix<real_t>;
  Matrix matrix{};
  Eigen::SparseLU<Matrix> factorization{};
  Eigen::VectorXd rhs{};
  Eigen::VectorXd unknown_voltage{};
  std::vector<int> node_to_unknown{};
};

SolverWorkspace::SolverWorkspace() : impl_{ std::make_unique<Impl>() } {}
SolverWorkspace::~SolverWorkspace() = default;
SolverWorkspace::SolverWorkspace(SolverWorkspace &&) noexcept = default;
SolverWorkspace &SolverWorkspace::operator=(SolverWorkspace &&) noexcept = default;

slide::Status SolverWorkspace::configure(const CompiledElectricalNetlist &netlist,
                                         std::size_t cell_count)
{
  if (!netlist.connected || netlist.node_count < 2 || cell_count == 0)
    return slide::Status::Invalid_parameters;
  auto candidate = std::make_unique<Impl>();
  candidate->node_to_unknown.assign(netlist.node_count, -1);
  int unknowns{};
  for (std::uint32_t node = 0; node < netlist.node_count; ++node)
    if (node != netlist.terminal_negative)
      candidate->node_to_unknown[node] = unknowns++;
  std::vector<Eigen::Triplet<real_t>> triplets;
  triplets.reserve(netlist.nodal_sparsity.size() * 2);
  for (const auto &[node_a, node_b] : netlist.nodal_sparsity) {
    const int a = candidate->node_to_unknown[node_a];
    const int b = candidate->node_to_unknown[node_b];
    if (a >= 0 && b >= 0) {
      triplets.emplace_back(a, b, 1.0);
      if (a != b)
        triplets.emplace_back(b, a, 1.0);
    }
  }
  candidate->matrix.resize(unknowns, unknowns);
  candidate->matrix.setFromTriplets(triplets.begin(), triplets.end());
  candidate->matrix.makeCompressed();
  candidate->factorization.analyzePattern(candidate->matrix);
  candidate->rhs.resize(unknowns);
  candidate->unknown_voltage.resize(unknowns);
  impl_ = std::move(candidate);
  factorized_resistance_.assign(cell_count, 0.0);
  valid_ = false;
  age_ = 0;
  numeric_factorizations_ = 0;
  symbolic_factorizations_ = 1;
  return slide::Status::Success;
}

slide::Status PackSolver::configure(const CompiledPackTopology &topology,
                                    std::span<const TheveninBatchView>
                                      batches)
{
  PackTheveninSystem thevenin;
  auto status = thevenin.configure(topology.cells, topology.batch_archetypes, batches);
  if (status != slide::Status::Success)
    return status;
  SolverWorkspace workspace;
  status = workspace.configure(topology.electrical, topology.cells.size());
  if (status != slide::Status::Success)
    return status;

  topology_ = topology;
  thevenin_ = std::move(thevenin);
  workspace_ = std::move(workspace);
  const auto cells = topology.cells.size();
  const auto nodes = topology.electrical.node_count;
  solution_.cell_current.assign(cells, 0.0);
  solution_.node_voltage.assign(nodes, 0.0);
  current_guess_.assign(cells, 0.0);
  candidate_current_.assign(cells, 0.0);
  ocv_.assign(cells, 0.0);
  resistance_.assign(cells, 0.0);
  candidate_node_voltage_.assign(nodes, 0.0);
  rollback_cell_current_.assign(cells, 0.0);
  rollback_node_voltage_.assign(nodes, 0.0);
  layer_voltage_.assign(topology.electrical.ladder_offsets.empty()
                          ? 0
                          : topology.electrical.ladder_offsets.size() - 1,
                        0.0);
  configured_ = true;
  has_solution_ = false;
  diagnostics_ = {};
  return slide::Status::Success;
}

slide::Status PackSolver::solve(real_t applied_current, PackSolveMode mode,
                                real_t current_tolerance, int max_iterations)
{
  return solveImpl(applied_current, mode, current_tolerance, max_iterations, true);
}

slide::Status PackSolver::solveImpl(real_t applied_current, PackSolveMode mode,
                                    real_t current_tolerance, int max_iterations,
                                    bool allow_source_stepping)
{
  if (!configured_ || !is_finite(applied_current) || !(current_tolerance > 0.0)
      || !is_finite(current_tolerance) || max_iterations <= 0)
    return slide::Status::Invalid_parameters;
  if (mode == PackSolveMode::ladder && !topology_.electrical.series_parallel_ladder)
    return slide::Status::Invalid_parameters;
  if (has_solution_)
    std::copy(solution_.cell_current.begin(), solution_.cell_current.end(), current_guess_.begin());
  else
    std::fill(current_guess_.begin(), current_guess_.end(), 0.0);
  if (has_solution_)
    std::copy(solution_.node_voltage.begin(), solution_.node_voltage.end(), candidate_node_voltage_.begin());
  else
    std::fill(candidate_node_voltage_.begin(), candidate_node_voltage_.end(), 0.0);

  diagnostics_.iterations = 0;
  diagnostics_.jacobian_refreshes = 0;
  diagnostics_.source_steps = 0;
  real_t previous_residual = std::numeric_limits<real_t>::max();
  int consecutive_divergence{};
  const int factorization_at_entry = workspace_.numericFactorizations();
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    auto status = thevenin_.linearize(current_guess_, ocv_, resistance_);
    if (status != slide::Status::Success)
      return status;
    status = mode == PackSolveMode::sparse_newton ? solveSparse(applied_current,
                                                                previous_residual,
                                                                iteration,
                                                                consecutive_divergence)
                                                  : solveLadder(applied_current);
    if (status != slide::Status::Success)
      return status;
    ++diagnostics_.iterations;
    if (mode == PackSolveMode::sparse_newton) {
      if (is_finite(previous_residual) && previous_residual > 0.0
          && residual_norm_ >= previous_residual)
        ++consecutive_divergence;
      else
        consecutive_divergence = 0;
      previous_residual = residual_norm_;
      diagnostics_.residual_norm = residual_norm_;
    }
    real_t max_change{};
    for (std::size_t cell = 0; cell < current_guess_.size(); ++cell)
      max_change = std::max(max_change,
                            std::abs(candidate_current_[cell] - current_guess_[cell]));
    std::copy(candidate_current_.begin(), candidate_current_.end(), current_guess_.begin());
    if (max_change <= current_tolerance) {
      solution_.cell_current = candidate_current_;
      solution_.node_voltage = candidate_node_voltage_;
      solution_.terminal_voltage = candidate_terminal_voltage_;
      has_solution_ = true;
      if (workspace_.numericFactorizations() == factorization_at_entry)
        ++workspace_.age_;
      diagnostics_.numeric_factorizations = workspace_.numericFactorizations();
      diagnostics_.symbolic_factorizations = workspace_.symbolicFactorizations();
      return slide::Status::Success;
    }
  }
  workspace_.invalidate();
  if (allow_source_stepping && applied_current != 0.0) {
    const bool rollback_has_solution = has_solution_;
    std::copy(solution_.cell_current.begin(), solution_.cell_current.end(), rollback_cell_current_.begin());
    std::copy(solution_.node_voltage.begin(), solution_.node_voltage.end(), rollback_node_voltage_.begin());
    const real_t rollback_terminal_voltage = solution_.terminal_voltage;
    has_solution_ = false;

    constexpr int source_steps = 8;
    auto status = solveImpl(0.0, mode, current_tolerance, max_iterations, false);
    for (int step = 1; status == slide::Status::Success && step <= source_steps; ++step)
      status = solveImpl(applied_current * static_cast<real_t>(step)
                           / static_cast<real_t>(source_steps),
                         mode,
                         current_tolerance,
                         max_iterations,
                         false);
    if (status == slide::Status::Success) {
      diagnostics_.source_steps = source_steps;
      return status;
    }

    std::copy(rollback_cell_current_.begin(), rollback_cell_current_.end(), solution_.cell_current.begin());
    std::copy(rollback_node_voltage_.begin(), rollback_node_voltage_.end(), solution_.node_voltage.begin());
    solution_.terminal_voltage = rollback_terminal_voltage;
    has_solution_ = rollback_has_solution;
    workspace_.invalidate();
  }
  return slide::Status::Numerical_failure;
}

slide::Status PackSolver::solveSparse(real_t applied_current,
                                      real_t previous_residual,
                                      int iteration,
                                      int consecutive_divergence)
{
  auto &impl = *workspace_.impl_;
  impl.rhs.setZero();
  auto addResidual = [&](std::uint32_t node, real_t value) {
    const int unknown = impl.node_to_unknown[node];
    if (unknown >= 0)
      impl.rhs[unknown] += value;
  };
  for (const auto &branch : topology_.electrical.branches) {
    const real_t voltage = candidate_node_voltage_[branch.node_positive]
                           - candidate_node_voltage_[branch.node_negative];
    const real_t branch_current = branch.kind == ElectricalBranchKind::cell
                                    ? (voltage - ocv_[branch.cell])
                                        / resistance_[branch.cell]
                                    : voltage / branch.resistance;
    addResidual(branch.node_positive, branch_current);
    addResidual(branch.node_negative, -branch_current);
  }
  addResidual(topology_.electrical.terminal_positive, applied_current);
  addResidual(topology_.electrical.terminal_negative, -applied_current);
  residual_norm_ = impl.rhs.lpNorm<Eigen::Infinity>();
  if (!is_finite(residual_norm_))
    return slide::Status::Invalid_states;

  const real_t contraction = is_finite(previous_residual) && previous_residual > 0.0
                               ? residual_norm_ / previous_residual
                               : 0.0;
  const bool refresh_requested = contraction > 0.5 || iteration > 4
                                 || consecutive_divergence >= 2;
  const bool same_resistance = workspace_.valid_
                               && std::equal(resistance_.begin(), resistance_.end(), workspace_.factorized_resistance_.begin());
  if (!workspace_.valid_ || (refresh_requested && !same_resistance)) {
    std::fill(impl.matrix.valuePtr(),
              impl.matrix.valuePtr() + impl.matrix.nonZeros(),
              0.0);
    auto stamp = [&](std::uint32_t positive, std::uint32_t negative, real_t conductance) {
      const int p = impl.node_to_unknown[positive];
      const int n = impl.node_to_unknown[negative];
      if (p >= 0)
        impl.matrix.coeffRef(p, p) += conductance;
      if (n >= 0)
        impl.matrix.coeffRef(n, n) += conductance;
      if (p >= 0 && n >= 0) {
        impl.matrix.coeffRef(p, n) -= conductance;
        impl.matrix.coeffRef(n, p) -= conductance;
      }
    };
    for (const auto &branch : topology_.electrical.branches) {
      const real_t resistance = branch.kind == ElectricalBranchKind::cell
                                  ? resistance_[branch.cell]
                                  : branch.resistance;
      if (!(is_finite(resistance) && resistance > 0.0))
        return slide::Status::Invalid_states;
      stamp(branch.node_positive, branch.node_negative, 1.0 / resistance);
    }
    impl.factorization.factorize(impl.matrix);
    if (impl.factorization.info() != Eigen::Success) {
      workspace_.invalidate();
      return slide::Status::Numerical_failure;
    }
    workspace_.factorized_resistance_ = resistance_;
    workspace_.valid_ = true;
    workspace_.age_ = 0;
    ++workspace_.numeric_factorizations_;
    if (refresh_requested)
      ++diagnostics_.jacobian_refreshes;
  }

  impl.rhs *= -1.0;
  impl.unknown_voltage = impl.factorization.solve(impl.rhs);
  if (impl.factorization.info() != Eigen::Success || !impl.unknown_voltage.allFinite())
    return slide::Status::Numerical_failure;

  real_t damping = 1.0;
  for (std::uint32_t node = 0; node < topology_.electrical.node_count; ++node) {
    const int unknown = impl.node_to_unknown[node];
    if (unknown >= 0)
      candidate_node_voltage_[node] += impl.unknown_voltage[unknown];
  }
  for (const auto &branch : topology_.electrical.branches)
    if (branch.kind == ElectricalBranchKind::cell) {
      const real_t voltage = candidate_node_voltage_[branch.node_positive]
                             - candidate_node_voltage_[branch.node_negative];
      candidate_current_[branch.cell] = (ocv_[branch.cell] - voltage)
                                        / resistance_[branch.cell];
      const real_t change = std::abs(candidate_current_[branch.cell]
                                     - current_guess_[branch.cell]);
      const real_t limit = 2.0 * std::max({ real_t{ 1.0 }, std::abs(applied_current), std::abs(current_guess_[branch.cell]) });
      if (change > limit)
        damping = std::min(damping, limit / change);
    }
  if (damping < 1.0) {
    for (std::uint32_t node = 0; node < topology_.electrical.node_count; ++node) {
      const int unknown = impl.node_to_unknown[node];
      if (unknown >= 0)
        candidate_node_voltage_[node] -= (1.0 - damping) * impl.unknown_voltage[unknown];
    }
    for (const auto &branch : topology_.electrical.branches)
      if (branch.kind == ElectricalBranchKind::cell) {
        const real_t voltage = candidate_node_voltage_[branch.node_positive]
                               - candidate_node_voltage_[branch.node_negative];
        candidate_current_[branch.cell] = (ocv_[branch.cell] - voltage)
                                          / resistance_[branch.cell];
      }
  }
  candidate_terminal_voltage_ = candidate_node_voltage_[topology_.electrical.terminal_positive]
                                - candidate_node_voltage_[topology_.electrical.terminal_negative];
  return slide::Status::Success;
}

slide::Status PackSolver::solveLadder(real_t applied_current)
{
  const auto &netlist = topology_.electrical;
  if (!netlist.series_parallel_ladder || netlist.ladder_offsets.size() < 2
      || netlist.ladder_nodes.size() != netlist.ladder_offsets.size())
    return slide::Status::Invalid_parameters;
  const std::size_t layers = netlist.ladder_offsets.size() - 1;
  if (layer_voltage_.size() != layers)
    return slide::Status::Invalid_parameters;
  for (std::size_t layer = 0; layer < layers; ++layer) {
    real_t conductance_sum{};
    real_t source_sum{};
    for (std::uint32_t i = netlist.ladder_offsets[layer];
         i < netlist.ladder_offsets[layer + 1];
         ++i) {
      const auto cell = netlist.ladder_cells[i];
      const real_t conductance = 1.0 / resistance_[cell];
      conductance_sum += conductance;
      source_sum += ocv_[cell] * conductance;
    }
    if (!(is_finite(conductance_sum) && conductance_sum > 0.0))
      return slide::Status::Invalid_states;
    layer_voltage_[layer] = (source_sum - applied_current) / conductance_sum;
    for (std::uint32_t i = netlist.ladder_offsets[layer];
         i < netlist.ladder_offsets[layer + 1];
         ++i) {
      const auto cell = netlist.ladder_cells[i];
      candidate_current_[cell] = (ocv_[cell] - layer_voltage_[layer]) / resistance_[cell];
    }
  }
  std::fill(candidate_node_voltage_.begin(), candidate_node_voltage_.end(), 0.0);
  real_t voltage{};
  for (std::size_t reverse = layers; reverse > 0; --reverse) {
    voltage += layer_voltage_[reverse - 1];
    candidate_node_voltage_[netlist.ladder_nodes[reverse - 1]] = voltage;
  }
  candidate_terminal_voltage_ = voltage;
  return slide::Status::Success;
}

} // namespace slide::core
