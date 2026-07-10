/**
 * @file PackSolver.cpp
 * @brief Affine Thevenin pack solving with persistent sparse and ladder workspaces.
 */

#include "PackSolver.hpp"
#include "PackTopologyInternal.hpp"

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace slide::core {
namespace {

  bool knownSolveMode(PackSolveMode mode)
  {
    switch (mode) {
    case PackSolveMode::sparse_newton:
    case PackSolveMode::ladder:
    case PackSolveMode::relaxation:
      return true;
    }
    return false;
  }

  bool addFinite(real_t &target, real_t increment)
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
  bool addCompensatedFinite(real_t &sum, real_t &compensation, real_t value)
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

  bool finiteCandidate(std::span<const real_t> current,
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

} // namespace

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
  for (std::size_t batch = 0; batch < batch_archetypes.size(); ++batch) {
    if (batch_archetypes[batch].empty()
        || std::find(batch_archetypes.begin(),
                     batch_archetypes.begin() + static_cast<std::ptrdiff_t>(batch),
                     batch_archetypes[batch])
             != batch_archetypes.begin() + static_cast<std::ptrdiff_t>(batch))
      return slide::Status::Invalid_parameters;
  }
  std::vector<int> required_lanes(batches.size());
  for (const auto &cell : cells) {
    if (cell.location.batch >= batches.size()
        || cell.location.lane >= static_cast<std::uint32_t>(
             std::numeric_limits<int>::max())
        || cell.archetype != batch_archetypes[cell.location.batch])
      return slide::Status::Invalid_parameters;
    auto &lanes = required_lanes[cell.location.batch];
    lanes = std::max(lanes, static_cast<int>(cell.location.lane + 1));
  }
  std::vector<std::vector<unsigned char>> occupied(batches.size());
  for (std::size_t batch = 0; batch < batches.size(); ++batch)
    occupied[batch].resize(static_cast<std::size_t>(required_lanes[batch]));
  for (const auto &cell : cells) {
    auto &slot = occupied[cell.location.batch][cell.location.lane];
    if (slot != 0)
      return slide::Status::Invalid_parameters;
    slot = 1;
  }
  for (const auto &batch : occupied)
    if (std::any_of(batch.begin(), batch.end(), [](unsigned char used) {
          return used == 0;
        }))
      return slide::Status::Invalid_parameters;
  std::vector<BatchScratch> scratch(batches.size());
  for (std::size_t batch = 0; batch < batches.size(); ++batch) {
    if (!batches[batch].valid() || batches[batch].lanes() != required_lanes[batch])
      return slide::Status::Invalid_parameters;
    for (std::size_t prior = 0; prior < batch; ++prior)
      if (batches[batch].identity() == batches[prior].identity())
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
                                              cell_resistance,
                                            BatchExecutor &executor)
{
  if (cell_current.size() != cells_.size() || cell_ocv.size() != cells_.size()
      || cell_resistance.size() != cells_.size())
    return slide::Status::Invalid_parameters;
  for (std::size_t cell = 0; cell < cells_.size(); ++cell) {
    const auto location = cells_[cell].location;
    batches_[location.batch].current[location.lane] = cell_current[cell];
  }
  const auto status = executor.parallelFor(
    batches_.size(),
    [&](std::size_t index) {
      auto &batch = batches_[index];
      return batch.view.linearize(batch.current, batch.ocv, batch.resistance);
    });
  if (status != slide::Status::Success)
    return status;
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
  const auto validation = detail::validateElectricalNetlist(netlist, cell_count);
  if (validation != slide::Status::Success)
    return validation;
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
                                      batches,
                                    unsigned workers)
{
  PackTheveninSystem thevenin;
  auto status = thevenin.configure(
    topology.cells, topology.batch_archetypes, batches);
  if (status != slide::Status::Success)
    return status;
  SolverWorkspace workspace;
  status = workspace.configure(topology.electrical, topology.cells.size());
  if (status != slide::Status::Success)
    return status;
  BatchExecutor batch_executor;
  status = batch_executor.configure(batches.size(), workers);
  if (status != slide::Status::Success)
    return status;

  topology_ = topology;
  thevenin_ = std::move(thevenin);
  batch_executor_ = std::move(batch_executor);
  workspace_ = std::move(workspace);
  const auto cells = topology.cells.size();
  const auto nodes = topology.electrical.node_count;
  solution_.cell_current.assign(cells, 0.0);
  solution_.node_voltage.assign(nodes, 0.0);
  solution_.terminal_voltage = 0.0;
  current_guess_.assign(cells, 0.0);
  candidate_current_.assign(cells, 0.0);
  ocv_.assign(cells, 0.0);
  resistance_.assign(cells, 0.0);
  candidate_node_voltage_.assign(nodes, 0.0);
  rollback_cell_current_.assign(cells, 0.0);
  rollback_node_voltage_.assign(nodes, 0.0);
  relaxation_diagonal_.assign(nodes, 0.0);
  relaxation_rhs_.assign(nodes, 0.0);
  relaxation_target_.assign(nodes, 0.0);
  relaxation_compensation_.assign(nodes, 0.0);
  layer_voltage_.assign(topology.electrical.ladder_offsets.empty()
                          ? 0
                          : topology.electrical.ladder_offsets.size() - 1,
                        0.0);
  candidate_terminal_voltage_ = 0.0;
  residual_norm_ = 0.0;
  configured_ = true;
  relaxation_alpha_ = topology.electrical.series_parallel_ladder ? 1.0 : 2.0 / 3.0;
  has_solution_ = false;
  diagnostics_ = {};
  return slide::Status::Success;
}

slide::Status PackSolver::solve(real_t applied_current, PackSolveMode mode,
                                real_t current_tolerance, int max_iterations)
{
  return solveImpl(applied_current, mode, current_tolerance, max_iterations, true);
}

slide::Status PackSolver::setRelaxationGain(real_t alpha)
{
  if (!configured_ || !is_finite(alpha) || !(alpha > 0.0 && alpha <= 1.0))
    return slide::Status::Invalid_parameters;
  relaxation_alpha_ = alpha;
  return slide::Status::Success;
}

slide::Status PackSolver::solveImpl(real_t applied_current, PackSolveMode mode,
                                    real_t current_tolerance, int max_iterations,
                                    bool allow_source_stepping)
{
  if (!configured_ || !is_finite(applied_current) || !(current_tolerance > 0.0)
      || !is_finite(current_tolerance) || max_iterations <= 0)
    return slide::Status::Invalid_parameters;
  if (!knownSolveMode(mode))
    return slide::Status::Invalid_parameters;
  if (mode == PackSolveMode::ladder && !topology_.electrical.series_parallel_ladder)
    return slide::Status::Invalid_parameters;
  if (mode == PackSolveMode::relaxation && !topology_.electrical.index1_candidate)
    return slide::Status::Invalid_parameters;
  if (has_solution_)
    std::copy(solution_.cell_current.begin(), solution_.cell_current.end(), current_guess_.begin());
  else
    std::fill(current_guess_.begin(), current_guess_.end(), 0.0);
  if (has_solution_)
    std::copy(solution_.node_voltage.begin(), solution_.node_voltage.end(), candidate_node_voltage_.begin());
  else
    std::fill(candidate_node_voltage_.begin(), candidate_node_voltage_.end(), 0.0);

  diagnostics_ = {};
  diagnostics_.relaxation_gain = mode == PackSolveMode::relaxation
                                   ? relaxation_alpha_
                                   : 0.0;
  diagnostics_.constraint_bound = mode == PackSolveMode::relaxation
                                    ? current_tolerance
                                    : 0.0;
  residual_norm_ = 0.0;
  real_t previous_residual = std::numeric_limits<real_t>::max();
  int consecutive_divergence{};
  const int factorization_at_entry = workspace_.numericFactorizations();
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    auto status = thevenin_.linearize(
      current_guess_, ocv_, resistance_, batch_executor_);
    if (status != slide::Status::Success) {
      workspace_.invalidate();
      return status;
    }
    status = mode == PackSolveMode::sparse_newton
               ? solveSparse(applied_current, previous_residual, iteration, consecutive_divergence)
             : mode == PackSolveMode::ladder
               ? solveLadder(applied_current)
               : solveRelaxation(applied_current);
    if (status != slide::Status::Success) {
      workspace_.invalidate();
      return status;
    }
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
    if (!finiteCandidate(candidate_current_, candidate_node_voltage_, candidate_terminal_voltage_)) {
      workspace_.invalidate();
      return slide::Status::Invalid_states;
    }
    for (std::size_t cell = 0; cell < current_guess_.size(); ++cell) {
      const real_t difference = candidate_current_[cell] - current_guess_[cell];
      if (!is_finite(difference)) {
        workspace_.invalidate();
        return slide::Status::Invalid_states;
      }
      const real_t change = std::abs(difference);
      if (!is_finite(change)) {
        workspace_.invalidate();
        return slide::Status::Invalid_states;
      }
      max_change = std::max(max_change, change);
    }
    if (mode != PackSolveMode::sparse_newton)
      diagnostics_.residual_norm = max_change;
    std::copy(candidate_current_.begin(), candidate_current_.end(), current_guess_.begin());
    const bool constraint_converged = mode != PackSolveMode::relaxation
                                      || (is_finite(diagnostics_.constraint_drift)
                                          && diagnostics_.constraint_drift
                                               <= current_tolerance);
    if (max_change <= current_tolerance && constraint_converged) {
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
    if (!is_finite(value))
      return false;
    const int unknown = impl.node_to_unknown[node];
    if (unknown >= 0)
      return addFinite(impl.rhs[unknown], value);
    return true;
  };
  for (const auto &branch : topology_.electrical.branches) {
    const real_t voltage = candidate_node_voltage_[branch.node_positive]
                           - candidate_node_voltage_[branch.node_negative];
    if (!is_finite(voltage))
      return slide::Status::Invalid_states;
    real_t branch_current{};
    if (branch.kind == ElectricalBranchKind::cell) {
      const real_t numerator = voltage - ocv_[branch.cell];
      if (!is_finite(numerator))
        return slide::Status::Invalid_states;
      branch_current = numerator / resistance_[branch.cell];
    } else {
      branch_current = voltage / branch.resistance;
    }
    if (!is_finite(branch_current)
        || !addResidual(branch.node_positive, branch_current)
        || !addResidual(branch.node_negative, -branch_current))
      return slide::Status::Invalid_states;
  }
  if (!addResidual(topology_.electrical.terminal_positive, applied_current)
      || !addResidual(topology_.electrical.terminal_negative, -applied_current))
    return slide::Status::Invalid_states;
  residual_norm_ = impl.rhs.lpNorm<Eigen::Infinity>();
  if (!is_finite(residual_norm_))
    return slide::Status::Invalid_states;

  real_t contraction{};
  if (is_finite(previous_residual) && previous_residual > 0.0) {
    contraction = residual_norm_ / previous_residual;
    if (!is_finite(contraction))
      return slide::Status::Invalid_states;
  }
  const bool refresh_requested = contraction > 0.5 || iteration > 4
                                 || consecutive_divergence >= 2;
  const bool same_resistance = workspace_.valid_
                               && std::equal(resistance_.begin(), resistance_.end(), workspace_.factorized_resistance_.begin());
  if (!workspace_.valid_ || (refresh_requested && !same_resistance)) {
    std::fill(impl.matrix.valuePtr(),
              impl.matrix.valuePtr() + impl.matrix.nonZeros(),
              0.0);
    auto addCoefficient = [&](int row, int column, real_t value) {
      real_t &coefficient = impl.matrix.coeffRef(row, column);
      return addFinite(coefficient, value);
    };
    auto stamp = [&](std::uint32_t positive, std::uint32_t negative, real_t conductance) {
      if (!is_finite(conductance))
        return false;
      const int p = impl.node_to_unknown[positive];
      const int n = impl.node_to_unknown[negative];
      if (p >= 0 && !addCoefficient(p, p, conductance))
        return false;
      if (n >= 0 && !addCoefficient(n, n, conductance))
        return false;
      if (p >= 0 && n >= 0) {
        if (!addCoefficient(p, n, -conductance)
            || !addCoefficient(n, p, -conductance))
          return false;
      }
      return true;
    };
    for (const auto &branch : topology_.electrical.branches) {
      const real_t resistance = branch.kind == ElectricalBranchKind::cell
                                  ? resistance_[branch.cell]
                                  : branch.resistance;
      if (!(is_finite(resistance) && resistance > 0.0))
        return slide::Status::Invalid_states;
      const real_t conductance = 1.0 / resistance;
      if (!is_finite(conductance)
          || !stamp(branch.node_positive, branch.node_negative, conductance))
        return slide::Status::Invalid_states;
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
  for (Eigen::Index i = 0; i < impl.rhs.size(); ++i)
    if (!is_finite(impl.rhs[i]))
      return slide::Status::Invalid_states;
  impl.unknown_voltage = impl.factorization.solve(impl.rhs);
  if (impl.factorization.info() != Eigen::Success)
    return slide::Status::Numerical_failure;
  for (Eigen::Index i = 0; i < impl.unknown_voltage.size(); ++i)
    if (!is_finite(impl.unknown_voltage[i]))
      return slide::Status::Invalid_states;

  real_t damping = 1.0;
  for (std::uint32_t node = 0; node < topology_.electrical.node_count; ++node) {
    const int unknown = impl.node_to_unknown[node];
    if (unknown >= 0
        && !addFinite(candidate_node_voltage_[node], impl.unknown_voltage[unknown]))
      return slide::Status::Invalid_states;
  }
  for (const auto &branch : topology_.electrical.branches)
    if (branch.kind == ElectricalBranchKind::cell) {
      const real_t voltage = candidate_node_voltage_[branch.node_positive]
                             - candidate_node_voltage_[branch.node_negative];
      const real_t numerator = ocv_[branch.cell] - voltage;
      if (!is_finite(voltage) || !is_finite(numerator))
        return slide::Status::Invalid_states;
      candidate_current_[branch.cell] = numerator / resistance_[branch.cell];
      const real_t difference = candidate_current_[branch.cell]
                                - current_guess_[branch.cell];
      if (!is_finite(candidate_current_[branch.cell]) || !is_finite(difference))
        return slide::Status::Invalid_states;
      const real_t change = std::abs(difference);
      const real_t magnitude = std::max({ real_t{ 1.0 },
                                          std::abs(applied_current),
                                          std::abs(current_guess_[branch.cell]) });
      const real_t limit = magnitude > std::numeric_limits<real_t>::max() / 2.0
                             ? std::numeric_limits<real_t>::max()
                             : 2.0 * magnitude;
      if (!is_finite(change) || !is_finite(magnitude) || !is_finite(limit))
        return slide::Status::Invalid_states;
      if (change > limit)
        damping = std::min(damping, limit / change);
    }
  if (!is_finite(damping))
    return slide::Status::Invalid_states;
  if (damping < 1.0) {
    for (std::uint32_t node = 0; node < topology_.electrical.node_count; ++node) {
      const int unknown = impl.node_to_unknown[node];
      if (unknown >= 0) {
        const real_t correction = -(1.0 - damping) * impl.unknown_voltage[unknown];
        if (!is_finite(correction)
            || !addFinite(candidate_node_voltage_[node], correction))
          return slide::Status::Invalid_states;
      }
    }
    for (const auto &branch : topology_.electrical.branches)
      if (branch.kind == ElectricalBranchKind::cell) {
        const real_t voltage = candidate_node_voltage_[branch.node_positive]
                               - candidate_node_voltage_[branch.node_negative];
        const real_t numerator = ocv_[branch.cell] - voltage;
        if (!is_finite(voltage) || !is_finite(numerator))
          return slide::Status::Invalid_states;
        candidate_current_[branch.cell] = numerator / resistance_[branch.cell];
        if (!is_finite(candidate_current_[branch.cell]))
          return slide::Status::Invalid_states;
      }
  }
  candidate_terminal_voltage_ = candidate_node_voltage_[topology_.electrical.terminal_positive]
                                - candidate_node_voltage_[topology_.electrical.terminal_negative];
  if (!is_finite(candidate_terminal_voltage_))
    return slide::Status::Invalid_states;
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
      const real_t source = ocv_[cell] * conductance;
      if (!is_finite(conductance) || !is_finite(source)
          || !addFinite(conductance_sum, conductance)
          || !addFinite(source_sum, source))
        return slide::Status::Invalid_states;
    }
    if (!(is_finite(conductance_sum) && conductance_sum > 0.0))
      return slide::Status::Invalid_states;
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
      const real_t current_numerator = ocv_[cell] - layer_voltage_[layer];
      if (!is_finite(current_numerator))
        return slide::Status::Invalid_states;
      candidate_current_[cell] = current_numerator / resistance_[cell];
      if (!is_finite(candidate_current_[cell]))
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
  if (!is_finite(candidate_terminal_voltage_))
    return slide::Status::Invalid_states;
  return slide::Status::Success;
}

slide::Status PackSolver::solveRelaxation(real_t applied_current)
{
  const auto &netlist = topology_.electrical;
  const auto finite = [](const real_t &value) { return is_finite(value); };
  if (!std::all_of(candidate_node_voltage_.begin(),
                   candidate_node_voltage_.end(),
                   finite))
    return slide::Status::Invalid_states;
  std::fill(relaxation_diagonal_.begin(), relaxation_diagonal_.end(), 0.0);
  std::fill(relaxation_rhs_.begin(), relaxation_rhs_.end(), 0.0);
  std::fill(relaxation_target_.begin(), relaxation_target_.end(), 0.0);
  std::fill(relaxation_compensation_.begin(),
            relaxation_compensation_.end(),
            0.0);
  auto stamp = [&](const CompiledElectricalBranch &branch, real_t resistance, real_t source) {
    const real_t conductance = 1.0 / resistance;
    if (!is_finite(conductance) || !is_finite(source))
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
           && addCompensatedFinite(relaxation_diagonal_[p],
                                   relaxation_target_[p],
                                   conductance)
           && addCompensatedFinite(relaxation_diagonal_[n],
                                   relaxation_target_[n],
                                   conductance)
           && addCompensatedFinite(relaxation_rhs_[p],
                                   relaxation_compensation_[p],
                                   positive_rhs)
           && addCompensatedFinite(relaxation_rhs_[n],
                                   relaxation_compensation_[n],
                                   negative_rhs);
  };
  for (const auto &branch : netlist.branches) {
    const bool cell = branch.kind == ElectricalBranchKind::cell;
    const real_t resistance = cell ? resistance_[branch.cell] : branch.resistance;
    if (!(is_finite(resistance) && resistance > 0.0))
      return slide::Status::Invalid_states;
    if (!stamp(branch, resistance, cell ? ocv_[branch.cell] : 0.0))
      return slide::Status::Invalid_states;
  }
  if (!addCompensatedFinite(relaxation_rhs_[netlist.terminal_positive],
                            relaxation_compensation_[netlist.terminal_positive],
                            -applied_current)
      || !addCompensatedFinite(
        relaxation_rhs_[netlist.terminal_negative],
        relaxation_compensation_[netlist.terminal_negative],
        applied_current))
    return slide::Status::Invalid_states;

  for (std::uint32_t node = 0; node < netlist.node_count; ++node) {
    if (node == netlist.terminal_negative) {
      relaxation_target_[node] = 0.0;
      continue;
    }
    if (!(is_finite(relaxation_diagonal_[node])
          && relaxation_diagonal_[node] > 0.0 && is_finite(relaxation_rhs_[node])))
      return slide::Status::Numerical_failure;
    relaxation_target_[node] = relaxation_rhs_[node]
                               / relaxation_diagonal_[node];
    if (!is_finite(relaxation_target_[node]))
      return slide::Status::Invalid_states;
  }
  for (std::uint32_t node = 0; node < netlist.node_count; ++node)
    if (node != netlist.terminal_negative) {
      const real_t difference = relaxation_target_[node]
                                - candidate_node_voltage_[node];
      const real_t correction = relaxation_alpha_ * difference;
      if (!is_finite(difference) || !is_finite(correction)
          || !addFinite(candidate_node_voltage_[node], correction))
        return slide::Status::Invalid_states;
    }

  for (const auto &branch : netlist.branches)
    if (branch.kind == ElectricalBranchKind::cell) {
      const real_t voltage = candidate_node_voltage_[branch.node_positive]
                             - candidate_node_voltage_[branch.node_negative];
      const real_t numerator = ocv_[branch.cell] - voltage;
      if (!is_finite(voltage) || !is_finite(numerator))
        return slide::Status::Invalid_states;
      candidate_current_[branch.cell] = numerator / resistance_[branch.cell];
      if (!is_finite(candidate_current_[branch.cell]))
        return slide::Status::Invalid_states;
    }
  candidate_terminal_voltage_ = candidate_node_voltage_[netlist.terminal_positive]
                                - candidate_node_voltage_[netlist.terminal_negative];
  if (!is_finite(candidate_terminal_voltage_))
    return slide::Status::Invalid_states;

  std::fill(relaxation_target_.begin(), relaxation_target_.end(), 0.0);
  std::fill(relaxation_compensation_.begin(),
            relaxation_compensation_.end(),
            0.0);
  real_t roundoff_operation_scale = std::abs(applied_current);
  real_t roundoff_current_scale = std::abs(applied_current);
  for (const auto &branch : netlist.branches) {
    const real_t voltage = candidate_node_voltage_[branch.node_positive]
                           - candidate_node_voltage_[branch.node_negative];
    if (!is_finite(voltage))
      return slide::Status::Invalid_states;
    real_t branch_current{};
    if (branch.kind == ElectricalBranchKind::cell) {
      const real_t numerator = voltage - ocv_[branch.cell];
      if (!is_finite(numerator))
        return slide::Status::Invalid_states;
      branch_current = numerator / resistance_[branch.cell];
    } else {
      branch_current = voltage / branch.resistance;
    }
    const real_t resistance = branch.kind == ElectricalBranchKind::cell
                                ? resistance_[branch.cell]
                                : branch.resistance;
    const real_t source = branch.kind == ElectricalBranchKind::cell
                            ? ocv_[branch.cell]
                            : 0.0;
    real_t operation_scale{};
    if (!addFinite(operation_scale,
                   std::abs(candidate_node_voltage_[branch.node_positive]))
        || !addFinite(operation_scale,
                      std::abs(candidate_node_voltage_[branch.node_negative]))
        || !addFinite(operation_scale, std::abs(source)))
      return slide::Status::Invalid_states;
    operation_scale /= resistance;
    const real_t current_magnitude = std::abs(branch_current);
    if (!is_finite(branch_current) || !is_finite(operation_scale)
        || !is_finite(current_magnitude)
        || !addFinite(roundoff_operation_scale, operation_scale)
        || !addFinite(roundoff_current_scale, current_magnitude)
        || !addCompensatedFinite(
          relaxation_target_[branch.node_positive],
          relaxation_compensation_[branch.node_positive],
          branch_current)
        || !addCompensatedFinite(
          relaxation_target_[branch.node_negative],
          relaxation_compensation_[branch.node_negative],
          -branch_current))
      return slide::Status::Invalid_states;
  }
  if (!addCompensatedFinite(
        relaxation_target_[netlist.terminal_positive],
        relaxation_compensation_[netlist.terminal_positive],
        applied_current)
      || !addCompensatedFinite(
        relaxation_target_[netlist.terminal_negative],
        relaxation_compensation_[netlist.terminal_negative],
        -applied_current))
    return slide::Status::Invalid_states;
  real_t drift{};
  for (std::uint32_t node = 0; node < netlist.node_count; ++node) {
    if (node == netlist.terminal_negative)
      continue;
    const real_t residual = relaxation_target_[node];
    const real_t magnitude = std::abs(residual);
    if (!is_finite(magnitude))
      return slide::Status::Invalid_states;
    drift = std::max(drift, magnitude);
  }
  const real_t epsilon = std::numeric_limits<real_t>::epsilon();
  const real_t accumulation_ratio = static_cast<real_t>(netlist.branches.size() + 1)
                                    * epsilon;
  if (!(is_finite(accumulation_ratio) && accumulation_ratio < 1.0))
    return slide::Status::Invalid_states;
  const real_t accumulation_estimate = accumulation_ratio
                                       / (1.0 - accumulation_ratio)
                                       * roundoff_current_scale;
  const real_t operation_estimate = 8.0 * epsilon * roundoff_operation_scale;
  const real_t roundoff_estimate = accumulation_estimate + operation_estimate;
  if (!is_finite(accumulation_estimate) || !is_finite(operation_estimate)
      || !is_finite(roundoff_estimate))
    return slide::Status::Invalid_states;
  diagnostics_.constraint_bound = std::max(diagnostics_.constraint_bound,
                                           roundoff_estimate);
  diagnostics_.constraint_drift = drift;
  return slide::Status::Success;
}

} // namespace slide::core
