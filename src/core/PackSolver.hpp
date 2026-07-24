/**
 * @file PackSolver.hpp
 * @brief Batch Thevenin interface, persistent sparse workspace, and Mode A/B pack solvers.
 * @surface api
 */

#pragma once

#include "PackTopology.hpp"
#include "ThreadPool.hpp"

#include <memory>
#include <span>
#include <vector>

namespace slide::core {

class TheveninBatchView
{
public:
  using LinearizeFn = slide::Status (*)(void *, std::span<const real_t>,
                                        std::span<real_t>, std::span<real_t>);

  TheveninBatchView() = default;
  TheveninBatchView(void *object, int lanes, LinearizeFn linearize)
    : object_{ object }, n_lanes_{ lanes }, linearize_{ linearize }
  {}

  template <class Batch>
  static TheveninBatchView bind(Batch &batch, int lanes)
  {
    return { &batch, lanes, [](void *object, std::span<const real_t> current, std::span<real_t> ocv, std::span<real_t> resistance) {
              return static_cast<Batch *>(object)->linearizeThevenin(current, ocv, resistance);
            } };
  }

  bool valid() const { return object_ != nullptr && n_lanes_ > 0 && linearize_ != nullptr; }
  const void *identity() const { return object_; }
  int n_lanes() const { return n_lanes_; }
  [[nodiscard]] slide::Status linearize(std::span<const real_t> current,
                                        std::span<real_t> ocv,
                                        std::span<real_t> resistance) const;

private:
  void *object_{};
  int n_lanes_{};
  LinearizeFn linearize_{};
};

class PackTheveninSystem
{
public:
  [[nodiscard]] slide::Status configure(std::span<const CompiledCell> cells,
                                        std::span<const std::string> batch_archetypes,
                                        std::span<const TheveninBatchView> batches);
  [[nodiscard]] slide::Status linearize(std::span<const real_t> cell_current,
                                        std::span<real_t> cell_ocv,
                                        std::span<real_t> cell_resistance,
                                        BatchExecutor &executor);

private:
  struct BatchScratch
  {
    TheveninBatchView view{};
    std::vector<real_t> current{};
    std::vector<real_t> ocv{};
    std::vector<real_t> resistance{};
  };
  std::vector<CompiledCell> cells_{};
  std::vector<BatchScratch> batches_{};
};

enum class PackSolveMode : unsigned char { sparse_newton,
                                           ladder,
                                           relaxation };

struct PackSolveDiagnostics
{
  int iterations{};
  int numeric_factorizations{};
  int symbolic_factorizations{};
  int jacobian_refreshes{};
  int source_steps{};
  real_t residual_norm{};
  real_t constraint_drift{};
  /** max(user KCL tolerance, floating-point roundoff estimate); not a general contraction proof */
  real_t constraint_bound{};
  real_t relaxation_gain{};
};

struct PackSolution
{
  std::vector<real_t> cell_current{}; //!< positive discharge current [A]
  std::vector<real_t> node_voltage{};
  real_t terminal_voltage{};
};

class SolverWorkspace
{
public:
  SolverWorkspace();
  ~SolverWorkspace();
  SolverWorkspace(SolverWorkspace &&) noexcept;
  SolverWorkspace &operator=(SolverWorkspace &&) noexcept;
  SolverWorkspace(const SolverWorkspace &) = delete;
  SolverWorkspace &operator=(const SolverWorkspace &) = delete;

  [[nodiscard]] slide::Status configure(const CompiledElectricalNetlist &netlist,
                                        std::size_t cell_count);
  void invalidate() { valid_ = false; }
  bool valid() const { return valid_; }
  int age() const { return age_; }
  int numericFactorizations() const { return numeric_factorizations_; }
  int symbolicFactorizations() const { return symbolic_factorizations_; }

private:
  struct Impl;
  std::unique_ptr<Impl> impl_{};
  std::vector<real_t> factorized_resistance_{};
  bool valid_{};
  int age_{};
  int numeric_factorizations_{};
  int symbolic_factorizations_{};

  friend class PackSolver;
};

class PackSolver
{
public:
  [[nodiscard]] slide::Status configure(const CompiledPackTopology &topology,
                                        std::span<const TheveninBatchView> batches,
                                        unsigned workers = 0);
  [[nodiscard]] slide::Status solve(real_t applied_current,
                                    PackSolveMode mode = PackSolveMode::sparse_newton,
                                    real_t current_tolerance = 1e-10,
                                    int max_iterations = 8);
  [[nodiscard]] slide::Status setRelaxationGain(real_t alpha);
  void invalidate()
  {
    workspace_.invalidate();
    has_solution_ = false;
  }

  const PackSolution &solution() const { return solution_; }
  const PackSolveDiagnostics &diagnostics() const { return diagnostics_; }
  const SolverWorkspace &workspace() const { return workspace_; }
  unsigned batchWorkerCount() const { return batch_executor_.workerCount(); }

private:
  [[nodiscard]] slide::Status solveImpl(real_t applied_current,
                                        PackSolveMode mode,
                                        real_t current_tolerance,
                                        int max_iterations,
                                        bool allow_source_stepping);
  [[nodiscard]] slide::Status solveSparse(real_t applied_current,
                                          real_t previous_residual,
                                          int iteration,
                                          int consecutive_divergence);
  [[nodiscard]] slide::Status solveLadder(real_t applied_current);
  [[nodiscard]] slide::Status solveRelaxation(real_t applied_current);

  CompiledPackTopology topology_{};
  PackTheveninSystem thevenin_{};
  BatchExecutor batch_executor_{};
  SolverWorkspace workspace_{};
  PackSolution solution_{};
  PackSolveDiagnostics diagnostics_{};
  std::vector<real_t> current_guess_{};
  std::vector<real_t> candidate_current_{};
  std::vector<real_t> ocv_{};
  std::vector<real_t> resistance_{};
  std::vector<real_t> candidate_node_voltage_{};
  std::vector<real_t> layer_voltage_{};
  std::vector<real_t> rollback_cell_current_{};
  std::vector<real_t> rollback_node_voltage_{};
  struct RelaxationScratch
  {
    std::vector<real_t> diagonal{};
    std::vector<real_t> diagonal_compensation{};
    std::vector<real_t> rhs{};
    std::vector<real_t> rhs_compensation{};
    std::vector<real_t> target{};
    std::vector<real_t> residual{};
    std::vector<real_t> residual_compensation{};
  };
  RelaxationScratch relaxation_{};
  real_t candidate_terminal_voltage_{};
  real_t residual_norm_{};
  real_t relaxation_alpha_{ 2.0 / 3.0 };
  bool configured_{};
  bool has_solution_{};

  friend class PackStepper;
};

} // namespace slide::core
