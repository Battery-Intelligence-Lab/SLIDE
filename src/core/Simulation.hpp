/**
 * @file Simulation.hpp
 * @brief Minimal Phase-1 single-batch Simulation façade and solution.
 * @surface api
 */

#pragma once

#include "EulerLegacy.hpp"

#include <cassert>
#include <span>
#include <vector>

namespace slide::core {

struct ConstantCurrentExperiment
{
  real_t current_A{}; //!< positive discharges the cell
  real_t duration{};
  real_t step{};
};

struct SimulationSolution
{
  int n_lanes{};
  std::vector<real_t> time{};
  std::vector<real_t> terminal_voltage{}; //!< sample-major [sample*n_lanes + lane]
  slide::Status termination{ slide::Status::Success };

  std::span<const real_t> voltageAt(std::size_t sample) const
  {
    assert(sample < time.size());
    return std::span<const real_t>{ terminal_voltage }
      .subspan(sample * static_cast<std::size_t>(n_lanes),
               static_cast<std::size_t>(n_lanes));
  }
};

class Simulation
{
public:
  Simulation() = default;
  Simulation(const Simulation &) = delete;
  Simulation &operator=(const Simulation &) = delete;
  Simulation(Simulation &&) noexcept = default;
  Simulation &operator=(Simulation &&) noexcept = default;

  [[nodiscard]] slide::Status build(const SpmFactoryInput &input,
                                    const SpmModelOptions &options,
                                    int n_lanes = 1);
  [[nodiscard]] slide::Status solve(const ConstantCurrentExperiment &experiment,
                                    SimulationSolution &output);

  bool valid() const { return batch_.valid(); }
  SpmBatch &batch() { return batch_; }
  const SpmBatch &batch() const { return batch_; }

private:
  SpmBatch batch_{};
  EulerLegacy stepper_{};
};

} // namespace slide::core
