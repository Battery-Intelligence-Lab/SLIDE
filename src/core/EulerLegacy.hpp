/**
 * @file EulerLegacy.hpp
 * @brief Allocation-free forward-Euler stepper for Phase-1 parity mode.
 */

#pragma once

#include "SpmFactory.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <span>
#include <vector>

namespace slide::core {

class EulerLegacy
{
public:
  EulerLegacy() = default;
  explicit EulerLegacy(const SpmBatch &batch)
  {
    const auto status = configure(batch);
    assert(status == slide::Status::Success);
    (void)status;
  }

  /** Allocate rollback and voltage scratch once, outside the stepping loop. */
  [[nodiscard]] slide::Status configure(const SpmBatch &batch)
  {
    if (!batch.valid())
      return slide::Status::Invalid_parameters;
    backup_.resize(batch.state().size());
    terminal_voltage_.resize(static_cast<std::size_t>(batch.n_lanes()));
    return slide::Status::Success;
  }

  [[nodiscard]] slide::Status step(SpmBatch &batch,
                                   std::span<const real_t> current_density,
                                   real_t time,
                                   real_t dt)
  {
    if (!batch.valid() || batch.state().size() != backup_.size()
        || static_cast<int>(current_density.size()) != batch.n_lanes()
        || terminal_voltage_.size() != current_density.size()
        || !is_finite(time) || !is_finite(dt) || !(dt > 0.0))
      return slide::Status::Invalid_parameters;
    for (const real_t current : current_density)
      if (!is_finite(current))
        return slide::Status::Invalid_parameters;

    std::copy(batch.state().raw().begin(), batch.state().raw().end(), backup_.begin());
    const StepCtx ctx{ .time = time, .dt = dt, .i_app = current_density };
    auto status = batch.evaluate(ctx);
    if (status != slide::Status::Success)
      return status;

    const auto roles = batch.roles();
    for (int row = 0; row < batch.state().n_rows(); ++row) {
      if (roles[static_cast<std::size_t>(row)] != StateRole::ode)
        continue;
      auto state_row = batch.state().row(row);
      const auto derivative_row = batch.derivative().row(row);
      for (int lane = 0; lane < batch.n_lanes(); ++lane)
        state_row[static_cast<std::size_t>(lane)] += dt * derivative_row[static_cast<std::size_t>(lane)];
    }

    const StepCtx accepted_ctx{ .time = time + dt,
                                .dt = 0.0,
                                .i_app = current_density };
    status = batch.terminalVoltage(accepted_ctx, terminal_voltage_);
    if (status != slide::Status::Success) {
      std::copy(backup_.begin(), backup_.end(), batch.state().raw().begin());
      return status;
    }

    const auto &layout = batch.layout();
    for (int lane = 0; lane < batch.n_lanes(); ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const real_t dAh = current_density[i] * batch.electrode_area() * dt / real_t{ 3600 };
      batch.state().at(layout.elapsed_time, 0, lane) += dt;
      batch.state().at(layout.charge_throughput, 0, lane) += std::abs(dAh);
      batch.state().at(layout.energy_throughput, 0, lane) += std::abs(dAh * terminal_voltage_[i]);
    }

    status = batch.storeStressHistory(dt);
    if (status != slide::Status::Success) {
      std::copy(backup_.begin(), backup_.end(), batch.state().raw().begin());
      return status;
    }
    return slide::Status::Success;
  }

  std::span<const real_t> terminalVoltage() const { return terminal_voltage_; }

private:
  std::vector<real_t> backup_{};
  std::vector<real_t> terminal_voltage_{};
};

} // namespace slide::core
