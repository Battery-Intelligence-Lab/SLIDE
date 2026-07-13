/**
 * @file ExponentialModal.hpp
 * @brief Allocation-free exact modal diffusion stepper with Strang slow-physics splitting.
 * @surface api
 */

#pragma once

#include "SpmFactory.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <span>
#include <vector>

namespace slide::core {

class ExponentialModal
{
public:
  ExponentialModal() = default;
  explicit ExponentialModal(const SpmBatch &batch)
  {
    const auto status = configure(batch);
    assert(status == slide::Status::Success);
    (void)status;
  }

  [[nodiscard]] slide::Status configure(const SpmBatch &batch)
  {
    if (!batch.valid())
      return slide::Status::Invalid_parameters;
    lanes_ = batch.n_lanes();
    rows_ = batch.state().n_rows();
    backup_.resize(batch.state().size());
    initial_.resize(batch.state().size());
    full_step_.resize(batch.state().size());
    terminal_voltage_.resize(static_cast<std::size_t>(lanes_));
    return slide::Status::Success;
  }

  [[nodiscard]] slide::Status step(SpmBatch &batch,
                                   std::span<const real_t> current_density,
                                   real_t time,
                                   real_t dt)
  {
    if (!batch.valid() || batch.n_lanes() != lanes_
        || batch.state().n_rows() != rows_
        || static_cast<int>(current_density.size()) != lanes_
        || !is_finite(time) || !is_finite(dt) || !(dt > 0.0))
      return slide::Status::Invalid_parameters;
    for (const real_t current : current_density)
      if (!is_finite(current))
        return slide::Status::Invalid_parameters;

    std::memcpy(backup_.data(), batch.state().raw().data(), batch.state().raw().size_bytes());
    const StepCtx ctx{ .time = time, .dt = dt, .i_app = current_density };
    auto status = batch.exponentialStep(ctx, dt, terminal_voltage_);
    if (status != slide::Status::Success)
      return restore(batch, status);

    const auto &layout = batch.layout();
    for (int lane = 0; lane < lanes_; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const real_t dAh = current_density[i] * batch.electrode_area() * dt / 3600.0;
      batch.state().at(layout.elapsed_time, 0, lane) += dt;
      batch.state().at(layout.charge_throughput, 0, lane) += std::abs(dAh);
      batch.state().at(layout.energy_throughput, 0, lane) +=
        std::abs(dAh * terminal_voltage_[i]);
    }
    status = batch.storeStressHistory(dt);
    if (status != slide::Status::Success)
      return restore(batch, status);
    return slide::Status::Success;
  }

  std::span<const real_t> terminalVoltage() const { return terminal_voltage_; }

  /** Step-doubling error control; accepted state is the two-half-step solution. */
  [[nodiscard]] slide::Status stepAdaptive(
    SpmBatch &batch,
    std::span<const real_t> current_density,
    real_t time,
    real_t proposed_dt,
    real_t absolute_tolerance,
    real_t relative_tolerance,
    real_t &accepted_dt,
    real_t &next_dt)
  {
    if (!batch.valid() || batch.state().size() != initial_.size()
        || !(is_finite(proposed_dt) && proposed_dt > 0.0)
        || !(is_finite(absolute_tolerance) && absolute_tolerance > 0.0)
        || !(is_finite(relative_tolerance) && relative_tolerance > 0.0))
      return slide::Status::Invalid_parameters;
    std::memcpy(initial_.data(), batch.state().raw().data(), batch.state().raw().size_bytes());
    real_t trial_dt = proposed_dt;
    for (int attempt = 0; attempt < 10; ++attempt) {
      std::memcpy(batch.state().raw().data(), initial_.data(), batch.state().raw().size_bytes());
      auto status = step(batch, current_density, time, trial_dt);
      if (status == slide::Status::Success)
        std::memcpy(full_step_.data(), batch.state().raw().data(), batch.state().raw().size_bytes());

      std::memcpy(batch.state().raw().data(), initial_.data(), batch.state().raw().size_bytes());
      if (status == slide::Status::Success)
        status = step(batch, current_density, time, 0.5 * trial_dt);
      if (status == slide::Status::Success)
        status = step(batch, current_density, time + 0.5 * trial_dt, 0.5 * trial_dt);
      if (status != slide::Status::Success) {
        trial_dt *= 0.5;
        continue;
      }

      real_t error{};
      const auto current = batch.state().raw();
      for (std::size_t i = 0; i < current.size(); ++i) {
        const real_t scale = absolute_tolerance
                             + relative_tolerance
                                 * std::max(std::abs(current[i]), std::abs(full_step_[i]));
        error = std::max(error, std::abs(current[i] - full_step_[i]) / scale);
      }
      const real_t factor = error > 0.0
                              ? std::clamp(0.9 / std::sqrt(error), 0.2, 2.0)
                              : 2.0;
      if (error <= 1.0) {
        accepted_dt = trial_dt;
        next_dt = trial_dt * factor;
        return slide::Status::Success;
      }
      trial_dt *= factor;
    }
    std::memcpy(batch.state().raw().data(), initial_.data(), batch.state().raw().size_bytes());
    return slide::Status::Numerical_failure;
  }

  /** Clamp a proposed step to the next known discontinuity without crossing it. */
  static real_t alignToEvent(real_t time, real_t proposed_dt, real_t event_time)
  {
    if (!(is_finite(time) && is_finite(proposed_dt) && proposed_dt > 0.0
          && is_finite(event_time) && event_time > time))
      return proposed_dt;
    return std::min(proposed_dt, event_time - time);
  }

private:
  slide::Status restore(SpmBatch &batch, slide::Status status)
  {
    std::memcpy(batch.state().raw().data(), backup_.data(), batch.state().raw().size_bytes());
    return status;
  }

  std::vector<real_t> backup_{};
  std::vector<real_t> initial_{};
  std::vector<real_t> full_step_{};
  std::vector<real_t> terminal_voltage_{};
  int lanes_{};
  int rows_{};
};

} // namespace slide::core
