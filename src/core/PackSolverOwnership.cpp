/**
 * @file PackSolverOwnership.cpp
 * @brief No-throw ownership transfer for configured pack solvers.
 */

#include "PackSolver.hpp"

#include <type_traits>
#include <utility>

namespace slide::core {
namespace {

  template <class... Types>
  inline constexpr bool nothrow_movable =
    (std::is_nothrow_move_constructible_v<Types> && ...)
    && (std::is_nothrow_move_assignable_v<Types> && ...);

  static_assert(std::is_nothrow_default_constructible_v<PackSolution>
                && std::is_nothrow_default_constructible_v<
                  PackSolveDiagnostics>);

} // namespace

PackSolver::PackSolver(PackSolver &&other) noexcept
  : topology_{ std::move(other.topology_) },
    thevenin_{ std::move(other.thevenin_) },
    batch_executor_{ std::move(other.batch_executor_) },
    workspace_{ std::move(other.workspace_) },
    solution_{ std::exchange(other.solution_, PackSolution{}) },
    diagnostics_{
      std::exchange(other.diagnostics_, PackSolveDiagnostics{})
    },
    current_guess_{ std::move(other.current_guess_) },
    candidate_current_{ std::move(other.candidate_current_) },
    ocv_{ std::move(other.ocv_) },
    resistance_{ std::move(other.resistance_) },
    candidate_node_voltage_{
      std::move(other.candidate_node_voltage_)
    },
    layer_voltage_{ std::move(other.layer_voltage_) },
    rollback_cell_current_{ std::move(other.rollback_cell_current_) },
    rollback_node_voltage_{ std::move(other.rollback_node_voltage_) },
    relaxation_{ std::move(other.relaxation_) },
    candidate_terminal_voltage_{
      std::exchange(other.candidate_terminal_voltage_, real_t{})
    },
    residual_norm_{ std::exchange(other.residual_norm_, real_t{}) },
    relaxation_alpha_{
      std::exchange(other.relaxation_alpha_, real_t{ 2.0 / 3.0 })
    },
    configured_{ std::exchange(other.configured_, false) },
    has_solution_{ std::exchange(other.has_solution_, false) }
{
  static_assert(nothrow_movable<
                CompiledPackTopology,
                PackTheveninSystem,
                BatchExecutor,
                SolverWorkspace,
                PackSolution,
                PackSolveDiagnostics,
                std::vector<real_t>,
                RelaxationScratch,
                real_t,
                bool>);
}

PackSolver &PackSolver::operator=(PackSolver &&other) noexcept
{
  if (this != &other) {
    topology_ = std::move(other.topology_);
    thevenin_ = std::move(other.thevenin_);
    batch_executor_ = std::move(other.batch_executor_);
    workspace_ = std::move(other.workspace_);
    solution_ = std::exchange(other.solution_, PackSolution{});
    diagnostics_ =
      std::exchange(other.diagnostics_, PackSolveDiagnostics{});
    current_guess_ = std::move(other.current_guess_);
    candidate_current_ = std::move(other.candidate_current_);
    ocv_ = std::move(other.ocv_);
    resistance_ = std::move(other.resistance_);
    candidate_node_voltage_ =
      std::move(other.candidate_node_voltage_);
    layer_voltage_ = std::move(other.layer_voltage_);
    rollback_cell_current_ = std::move(other.rollback_cell_current_);
    rollback_node_voltage_ =
      std::move(other.rollback_node_voltage_);
    relaxation_ = std::move(other.relaxation_);
    candidate_terminal_voltage_ =
      std::exchange(other.candidate_terminal_voltage_, real_t{});
    residual_norm_ = std::exchange(other.residual_norm_, real_t{});
    relaxation_alpha_ =
      std::exchange(other.relaxation_alpha_, real_t{ 2.0 / 3.0 });
    configured_ = std::exchange(other.configured_, false);
    has_solution_ = std::exchange(other.has_solution_, false);
  }
  return *this;
}

} // namespace slide::core
