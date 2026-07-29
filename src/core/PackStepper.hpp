/**
 * @file PackStepper.hpp
 * @brief Transactional staggered electrical/thermal advance for compiled SPM packs.
 * @surface api
 */

#pragma once

#include "EulerLegacy.hpp"
#include "ExponentialModal.hpp"
#include "PackSolver.hpp"

#include <span>
#include <vector>

namespace slide::core {

class PackStepper
{
public:
  PackStepper() = default;
  PackStepper(const PackStepper &) = delete;
  PackStepper &operator=(const PackStepper &) = delete;
  PackStepper(PackStepper &&other) noexcept;
  PackStepper &operator=(PackStepper &&other) noexcept;

  [[nodiscard]] slide::Status configure(
    const CompiledPackTopology &topology,
    std::span<SpmBatch *const> batches,
    unsigned workers = 0);

  //!< Advances the pack by `substeps * dt`, NOT by `dt`.
  //!<
  //!< `substeps` does not subdivide `dt`: the loop runs `substeps` iterations and
  //!< each one advances the full `dt`, at `time + substep * dt` (PackStepper.cpp).
  //!< It is the number of inner electrical steps taken per outer solve — the
  //!< electrical current is re-solved once, then held frozen across them (§3.5
  //!< multirate). The name reads like subdivision, so it is spelled out here:
  //!< `step(I, t, 1.0, ..., substeps = 10)` advances **10 seconds**, not one.
  [[nodiscard]] slide::Status step(
    real_t applied_current,
    real_t time,
    real_t dt,
    std::span<const real_t> boundary_temperature = {},
    PackSolveMode mode = PackSolveMode::sparse_newton,
    real_t current_tolerance = 1e-10,
    int substeps = 1);
  [[nodiscard]] slide::Status stepExponential(
    real_t applied_current,
    real_t time,
    real_t dt,
    std::span<const real_t> boundary_temperature = {},
    PackSolveMode mode = PackSolveMode::sparse_newton,
    real_t current_tolerance = 1e-10);

  std::size_t checkpointSize() const { return checkpoint_.size(); }
  [[nodiscard]] slide::Status checkpoint(std::span<real_t> destination) const;
  [[nodiscard]] slide::Status restore(std::span<const real_t> source);

  /** Refresh the electrical solution without advancing any cell state. */
  [[nodiscard]] slide::Status solveElectrical(
    real_t applied_current,
    PackSolveMode mode = PackSolveMode::sparse_newton,
    real_t current_tolerance = 1e-10)
  {
    return solver_.solve(applied_current, mode, current_tolerance);
  }

  const PackSolution &solution() const { return solver_.solution(); }
  const PackSolveDiagnostics &diagnostics() const { return solver_.diagnostics(); }
  std::span<const real_t> cellExternalHeat() const { return cell_external_heat_; }
  std::span<const real_t> boundaryHeat() const { return boundary_heat_; }
  const PackSolver &solver() const { return solver_; }
  unsigned batchWorkerCount() const { return solver_.batchWorkerCount(); }

private:
  [[nodiscard]] slide::Status stepImpl(real_t applied_current,
                                       real_t time,
                                       real_t dt,
                                       std::span<const real_t> boundary_temperature,
                                       PackSolveMode mode,
                                       real_t current_tolerance,
                                       int substeps,
                                       bool exponential);
  void gatherStates(std::span<real_t> destination) const;
  void scatterStates(std::span<const real_t> source);
  void saveCheckpoint();
  void restoreCheckpoint();

  CompiledPackTopology topology_{};
  std::vector<SpmBatch *> batches_{};
  std::vector<EulerLegacy> steppers_{};
  std::vector<ExponentialModal> exponential_steppers_{};
  PackSolver solver_{};
  std::vector<std::vector<real_t>> current_density_{};
  std::vector<std::size_t> checkpoint_offsets_{};
  std::vector<real_t> checkpoint_{};
  std::vector<real_t> cell_temperature_{};
  std::vector<real_t> cell_external_heat_{};
  std::vector<real_t> boundary_heat_{};
  PackSolution solver_checkpoint_solution_{};
  PackSolveDiagnostics solver_checkpoint_diagnostics_{};
  std::vector<real_t> cell_external_heat_checkpoint_{};
  std::vector<real_t> boundary_heat_checkpoint_{};
  bool solver_checkpoint_has_solution_{};
  bool configured_{};
};

} // namespace slide::core
