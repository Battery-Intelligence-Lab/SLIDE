/**
 * @file PackSolver.cpp
 * @brief Affine Thevenin pack solving with persistent sparse and ladder workspaces.
 */

#include "PackSolver.hpp"
#include "PackSolverInternal.hpp"
#include "PackTopologyInternal.hpp"

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>

namespace slide::core {
using detail::addCompensatedFinite;
using detail::addFinite;
using detail::finiteCandidate;
using detail::knownSolveMode;


namespace detail {

  real_t conservativePackRoundoffBound(real_t accumulation_ratio,
                                       real_t current_scale,
                                       real_t operation_scale) noexcept
  {
    assert(is_finite(current_scale) && current_scale >= 0.0);
    assert(is_finite(operation_scale) && operation_scale >= 0.0);
    const real_t epsilon = std::numeric_limits<real_t>::epsilon();
    if (!(is_finite(accumulation_ratio) && accumulation_ratio >= 0.0
          && accumulation_ratio < 1.0))
      return std::numeric_limits<real_t>::max();
    const real_t accumulation = accumulation_ratio
                                / (1.0 - accumulation_ratio) * current_scale;
    const real_t operations = 8.0 * epsilon * operation_scale;
    const real_t total = accumulation + operations;
    return is_finite(accumulation) && is_finite(operations) && is_finite(total)
             ? total
             : std::numeric_limits<real_t>::max();
  }

} // namespace detail

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
  try {
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

    const auto cells = topology.cells.size();
    const auto nodes = static_cast<std::size_t>(topology.electrical.node_count);
    const auto layers = topology.electrical.ladder_offsets.empty()
                          ? std::size_t{}
                          : topology.electrical.ladder_offsets.size() - 1;
    CompiledPackTopology candidate_topology = topology;
    PackSolution solution{ .cell_current = std::vector<real_t>(cells, 0.0),
                           .node_voltage = std::vector<real_t>(nodes, 0.0),
                           .terminal_voltage = 0.0 };
    std::vector<real_t> current_guess(cells, 0.0);
    std::vector<real_t> candidate_current(cells, 0.0);
    std::vector<real_t> ocv(cells, 0.0);
    std::vector<real_t> resistance(cells, 0.0);
    std::vector<real_t> candidate_node_voltage(nodes, 0.0);
    std::vector<real_t> layer_voltage(layers, 0.0);
    std::vector<real_t> rollback_cell_current(cells, 0.0);
    std::vector<real_t> rollback_node_voltage(nodes, 0.0);
    std::vector<real_t> relaxation_diagonal(nodes, 0.0);
    std::vector<real_t> relaxation_rhs(nodes, 0.0);
    std::vector<real_t> relaxation_target(nodes, 0.0);
    std::vector<real_t> relaxation_compensation(nodes, 0.0);

    static_assert(std::is_nothrow_move_assignable_v<CompiledPackTopology>);
    static_assert(std::is_nothrow_move_assignable_v<PackTheveninSystem>);
    static_assert(std::is_nothrow_move_assignable_v<BatchExecutor>);
    static_assert(std::is_nothrow_move_assignable_v<SolverWorkspace>);
    static_assert(std::is_nothrow_move_assignable_v<PackSolution>);
    static_assert(std::is_nothrow_move_assignable_v<std::vector<real_t>>);

    // Every fallible operation is complete. Publish the candidate state only
    // through no-throw moves so failed reconfiguration preserves the old solver.
    topology_ = std::move(candidate_topology);
    thevenin_ = std::move(thevenin);
    batch_executor_ = std::move(batch_executor);
    workspace_ = std::move(workspace);
    solution_ = std::move(solution);
    current_guess_ = std::move(current_guess);
    candidate_current_ = std::move(candidate_current);
    ocv_ = std::move(ocv);
    resistance_ = std::move(resistance);
    candidate_node_voltage_ = std::move(candidate_node_voltage);
    layer_voltage_ = std::move(layer_voltage);
    rollback_cell_current_ = std::move(rollback_cell_current);
    rollback_node_voltage_ = std::move(rollback_node_voltage);
    relaxation_diagonal_ = std::move(relaxation_diagonal);
    relaxation_rhs_ = std::move(relaxation_rhs);
    relaxation_target_ = std::move(relaxation_target);
    relaxation_compensation_ = std::move(relaxation_compensation);
    candidate_terminal_voltage_ = 0.0;
    residual_norm_ = 0.0;
    configured_ = true;
    relaxation_alpha_ = topology.electrical.series_parallel_ladder
                          ? 1.0
                          : 2.0 / 3.0;
    has_solution_ = false;
    diagnostics_ = {};
    return slide::Status::Success;
  } catch (const std::bad_alloc &) {
    return slide::Status::Numerical_failure;
  } catch (const std::length_error &) {
    return slide::Status::Invalid_parameters;
  }
}

slide::Status PackSolver::solveImpl(real_t applied_current,
                                    PackSolveMode mode,
                                    real_t current_tolerance,
                                    int max_iterations,
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
    // The selected private kernel validates every candidate assignment before
    // success; topology validation proves one branch per cell and full node
    // coverage. Keep the shared postcondition executable in Debug.
    assert(finiteCandidate(candidate_current_, candidate_node_voltage_, candidate_terminal_voltage_));
    for (std::size_t cell = 0; cell < current_guess_.size(); ++cell) {
      const real_t difference = candidate_current_[cell] - current_guess_[cell];
      if (!is_finite(difference)) {
        workspace_.invalidate();
        return slide::Status::Invalid_states;
      }
      const real_t change = std::abs(difference);
      assert(is_finite(change)); // abs preserves finiteness
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
      status = solveImpl(applied_current
                           * (static_cast<real_t>(step)
                              / static_cast<real_t>(source_steps)),
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
    // First use starts from zero; every successful kernel validates all branch
    // drops before publication, so a warm-start drop is finite here.
    assert(is_finite(voltage));
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
  // Infinity norm is the maximum absolute value of the finite assembled RHS.
  assert(is_finite(residual_norm_));

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
      // Linearization and netlist validation established this immediately
      // before entering the numeric kernel.
      assert(is_finite(resistance) && resistance > 0.0);
      const real_t conductance = 1.0 / resistance;
      if (!is_strictly_positive_finite(conductance)
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
    assert(is_finite(impl.rhs[i])); // finite negation is exact
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
      assert(is_finite(change) && is_finite(magnitude) && is_finite(limit));
      if (change > limit)
        damping = std::min(damping, limit / change);
    }
  assert(is_finite(damping) && damping >= 0.0 && damping <= 1.0);
  if (damping < 1.0) {
    bool damped_candidate_valid = true;
    for (std::uint32_t node = 0; node < topology_.electrical.node_count; ++node) {
      const int unknown = impl.node_to_unknown[node];
      if (unknown >= 0) {
        const real_t correction = -(1.0 - damping) * impl.unknown_voltage[unknown];
        if (!is_finite(correction)
            || !addFinite(candidate_node_voltage_[node], correction)) {
          damped_candidate_valid = false;
          break;
        }
      }
    }
    for (const auto &branch : topology_.electrical.branches)
      if (damped_candidate_valid
          && branch.kind == ElectricalBranchKind::cell) {
        const real_t voltage = candidate_node_voltage_[branch.node_positive]
                               - candidate_node_voltage_[branch.node_negative];
        const real_t numerator = ocv_[branch.cell] - voltage;
        if (is_finite(voltage) && is_finite(numerator))
          candidate_current_[branch.cell] = numerator / resistance_[branch.cell];
        if (!is_finite(voltage) || !is_finite(numerator)
            || !is_finite(candidate_current_[branch.cell])) {
          damped_candidate_valid = false;
          break;
        }
      }
    if (!damped_candidate_valid)
      return slide::Status::Invalid_states;
  }
  for (const auto &branch : topology_.electrical.branches)
    if (branch.kind == ElectricalBranchKind::resistor) {
      const real_t voltage = candidate_node_voltage_[branch.node_positive]
                             - candidate_node_voltage_[branch.node_negative];
      const real_t current = voltage / branch.resistance;
      if (!is_finite(voltage) || !is_finite(current))
        return slide::Status::Invalid_states;
    }
  candidate_terminal_voltage_ = candidate_node_voltage_[topology_.electrical.terminal_positive]
                                - candidate_node_voltage_[topology_.electrical.terminal_negative];
  // SolverWorkspace excludes terminal_negative from the unknowns, so it stays
  // exactly zero; terminal_positive was checked on every correction.
  assert(candidate_node_voltage_[topology_.electrical.terminal_negative]
         == 0.0);
  assert(is_finite(candidate_terminal_voltage_));
  return slide::Status::Success;
}

} // namespace slide::core
