/**
 * @file EulerLegacy.hpp
 * @brief Allocation-free forward-Euler stepper for Phase-1 parity mode.
 * @surface api
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
    ode_rows_.clear();
    rollback_rows_.clear();
    const bool ageing = batch.composition() == SpmComposition::isothermal_ageing
                        || batch.composition() == SpmComposition::thermal_ageing;
    const auto roles = batch.roles();
    for (int row = 0; row < batch.state().n_rows(); ++row) {
      const auto role = roles[static_cast<std::size_t>(row)];
      if (role == StateRole::ode)
        ode_rows_.push_back(row);
      if (ageing || role == StateRole::ode || role == StateRole::cumulative)
        rollback_rows_.push_back(row);
    }
    lanes_ = batch.n_lanes();
    active_lanes_ = batch.trustedLanePeriod();
    state_rows_ = batch.state().n_rows();
    backup_.resize(rollback_rows_.size() * static_cast<std::size_t>(active_lanes_));
    terminal_voltage_.resize(static_cast<std::size_t>(batch.n_lanes()));
    return slide::Status::Success;
  }

  [[nodiscard]] slide::Status step(SpmBatch &batch,
                                   std::span<const real_t> current_density,
                                   real_t time,
                                   real_t dt)
  {
    if (!batch.valid() || batch.n_lanes() != lanes_
        || batch.state().n_rows() != state_rows_
        || static_cast<int>(current_density.size()) != batch.n_lanes()
        || terminal_voltage_.size() != current_density.size()
        || !is_finite(time) || !is_finite(dt) || !(dt > 0.0))
      return slide::Status::Invalid_parameters;
    for (const real_t current : current_density)
      if (!is_finite(current))
        return slide::Status::Invalid_parameters;

    backup(batch);
    const StepCtx ctx{ .time = time, .dt = dt, .i_app = current_density };
    slide::Status status;
    bool voltage_ready = false;
    if (batch.hasFusedEuler()) {
      status = batch.fusedEuler(ctx, dt, terminal_voltage_);
      voltage_ready = status == slide::Status::Success;
    } else {
      status = batch.evaluate(ctx);
      if (status == slide::Status::Success) {
        for (const int row : ode_rows_) {
          auto state_row = batch.state().row(row);
          const auto derivative_row = batch.derivative().row(row);
          for (int lane = 0; lane < batch.n_lanes(); ++lane)
            state_row[static_cast<std::size_t>(lane)] += dt * derivative_row[static_cast<std::size_t>(lane)];
        }
      }
    }
    if (status != slide::Status::Success) {
      restore(batch);
      return status;
    }

    if (!voltage_ready) {
      const StepCtx accepted_ctx{ .time = time + dt,
                                  .dt = 0.0,
                                  .i_app = current_density };
      status = batch.terminalVoltage(accepted_ctx, terminal_voltage_);
      if (status != slide::Status::Success) {
        restore(batch);
        return status;
      }
    }

    const auto &layout = batch.layout();
    for (int lane = 0; lane < active_lanes_; ++lane) {
      const auto i = static_cast<std::size_t>(lane);
      const real_t dAh = current_density[i] * batch.electrode_area() * dt / real_t{ 3600 };
      batch.state().at(layout.elapsed_time, 0, lane) += dt;
      batch.state().at(layout.charge_throughput, 0, lane) += std::abs(dAh);
      batch.state().at(layout.energy_throughput, 0, lane) += std::abs(dAh * terminal_voltage_[i]);
    }
    for (const auto slice : { layout.elapsed_time,
                              layout.charge_throughput,
                              layout.energy_throughput }) {
      auto row = batch.state().row(slice.row_begin);
      for (int lane = active_lanes_; lane < lanes_; ++lane)
        row[static_cast<std::size_t>(lane)] = row[static_cast<std::size_t>(lane % active_lanes_)];
    }

    status = batch.storeStressHistory(dt);
    if (status != slide::Status::Success) {
      restore(batch);
      return status;
    }
    return slide::Status::Success;
  }

  std::span<const real_t> terminalVoltage() const { return terminal_voltage_; }

private:
  void backup(const SpmBatch &batch)
  {
    std::size_t offset{};
    for (const int row : rollback_rows_) {
      const auto source = batch.state().row(row);
      std::copy_n(source.begin(), active_lanes_, backup_.begin() + static_cast<std::ptrdiff_t>(offset));
      offset += static_cast<std::size_t>(active_lanes_);
    }
  }

  void restore(SpmBatch &batch)
  {
    std::size_t offset{};
    for (const int row : rollback_rows_) {
      auto destination = batch.state().row(row);
      std::copy(backup_.begin() + static_cast<std::ptrdiff_t>(offset),
                backup_.begin() + static_cast<std::ptrdiff_t>(offset + active_lanes_),
                destination.begin());
      for (int lane = active_lanes_; lane < lanes_; ++lane)
        destination[static_cast<std::size_t>(lane)] = destination[static_cast<std::size_t>(lane % active_lanes_)];
      offset += static_cast<std::size_t>(active_lanes_);
    }
  }

  std::vector<real_t> backup_{};
  std::vector<real_t> terminal_voltage_{};
  std::vector<int> ode_rows_{};
  std::vector<int> rollback_rows_{};
  int lanes_{};
  int active_lanes_{};
  int state_rows_{};
};

} // namespace slide::core
