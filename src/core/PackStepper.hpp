/**
 * @file PackStepper.hpp
 * @brief Transactional staggered electrical/thermal advance for compiled SPM packs.
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
  [[nodiscard]] slide::Status configure(
    const CompiledPackTopology &topology,
    std::span<SpmBatch *const> batches);

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

  const PackSolution &solution() const { return solver_.solution(); }
  const PackSolveDiagnostics &diagnostics() const { return solver_.diagnostics(); }
  std::span<const real_t> cellExternalHeat() const { return cell_external_heat_; }
  std::span<const real_t> boundaryHeat() const { return boundary_heat_; }
  PackSolver &solver() { return solver_; }
  const PackSolver &solver() const { return solver_; }

private:
  [[nodiscard]] slide::Status stepImpl(real_t applied_current,
                                       real_t time,
                                       real_t dt,
                                       std::span<const real_t> boundary_temperature,
                                       PackSolveMode mode,
                                       real_t current_tolerance,
                                       int substeps,
                                       bool exponential);
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
  bool configured_{};
};

} // namespace slide::core
