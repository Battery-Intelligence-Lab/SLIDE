/**
 * @file Simulation.cpp
 * @brief Phase-1 Simulation façade implementation.
 */

#include "Simulation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace slide::core {

slide::Status Simulation::build(const SpmFactoryInput &input,
                                const SpmModelOptions &options,
                                int n_lanes)
{
  Simulation candidate;
  auto status = buildSpmBatch(input, options, n_lanes, candidate.batch_);
  if (status != slide::Status::Success)
    return status;
  status = candidate.stepper_.configure(candidate.batch_);
  if (status != slide::Status::Success)
    return status;
  *this = std::move(candidate);
  return slide::Status::Success;
}

slide::Status Simulation::solve(const ConstantCurrentExperiment &experiment,
                                SimulationSolution &output)
{
  if (!valid() || !is_finite(experiment.current_A)
      || !is_finite(experiment.duration) || experiment.duration < 0.0
      || !is_finite(experiment.step) || !(experiment.step > 0.0))
    return slide::Status::Invalid_parameters;

  double step_ratio = experiment.duration / experiment.step;
  const double nearest_ratio = std::round(step_ratio);
  const double snap_tolerance = 16.0 * std::numeric_limits<double>::epsilon()
                                * std::max(1.0, std::abs(step_ratio));
  if (std::abs(step_ratio - nearest_ratio) <= snap_tolerance)
    step_ratio = nearest_ratio;
  const double raw_steps = experiment.duration == 0.0 ? 0.0
                                                      : std::ceil(step_ratio);
  if (!is_finite(raw_steps)
      || raw_steps > static_cast<double>(std::numeric_limits<std::size_t>::max() - 1))
    return slide::Status::Invalid_parameters;
  const auto steps = static_cast<std::size_t>(raw_steps);
  const auto lanes = static_cast<std::size_t>(batch_.n_lanes());
  if (lanes != 0 && steps + 1 > std::numeric_limits<std::size_t>::max() / lanes)
    return slide::Status::Invalid_parameters;

  SimulationSolution solution;
  solution.n_lanes = batch_.n_lanes();
  solution.time.resize(steps + 1);
  solution.terminal_voltage.resize((steps + 1) * lanes);
  std::vector<real_t> current_density(lanes,
                                      experiment.current_A / batch_.electrode_area());
  const auto &layout = batch_.layout();
  real_t time = batch_.state().at(layout.elapsed_time, 0, 0);
  solution.time[0] = time;
  StepCtx observation_ctx{ .time = time, .dt = 0.0, .i_app = current_density };
  auto status = batch_.terminalVoltage(observation_ctx,
                                       std::span<real_t>{ solution.terminal_voltage }.first(lanes));
  if (status != slide::Status::Success)
    return status;

  std::size_t completed{};
  while (completed < steps) {
    const real_t remaining = experiment.duration
                             - static_cast<real_t>(completed) * experiment.step;
    const real_t dt = std::min(experiment.step, remaining);
    status = stepper_.step(batch_, current_density, time, dt);
    if (status != slide::Status::Success) {
      solution.termination = status;
      solution.time.resize(completed + 1);
      solution.terminal_voltage.resize((completed + 1) * lanes);
      output = std::move(solution);
      return status;
    }
    ++completed;
    time += dt;
    solution.time[completed] = time;
    std::copy(stepper_.terminalVoltage().begin(),
              stepper_.terminalVoltage().end(),
              solution.terminal_voltage.begin()
                + static_cast<std::ptrdiff_t>(completed * lanes));
  }

  solution.termination = slide::Status::Success;
  output = std::move(solution);
  return slide::Status::Success;
}

} // namespace slide::core
