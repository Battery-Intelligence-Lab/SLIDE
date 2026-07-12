/**
 * @file CyclerV2.cpp
 * @brief Event-aligned experiment runner and control implementations.
 *
 * MC-3 ownership contract: all CyclerV2 method definitions and runner-only
 * transaction helpers live in this translation unit. Parser mechanics remain
 * in Experiment.cpp; shared semantics live in detail/ExperimentSemantics.hpp.
 */

#include "Experiment.hpp"
#include "detail/ExperimentSemantics.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace slide::core {
namespace {

  [[nodiscard]] slide::Status allocationFailureStatus() noexcept
  {
    return slide::Status::Numerical_failure;
  }

  class RunActivity
  {
  public:
    explicit RunActivity(bool &active) noexcept : active_{ active }
    {
      assert(!active_);
      active_ = true;
    }
    ~RunActivity() noexcept { active_ = false; }

    RunActivity(const RunActivity &) = delete;
    RunActivity &operator=(const RunActivity &) = delete;

  private:
    bool &active_;
  };

  void restoreBatchSnapshot(SpmBatch &batch,
                            std::span<const real_t>
                              state_backup,
                            std::span<const real_t>
                              derivative_backup) noexcept
  {
    const auto state = batch.state().raw();
    const auto derivative = batch.derivative().raw();
    assert(state_backup.size() == state.size());
    assert(derivative_backup.size() == derivative.size());
    if (state_backup.size() != state.size()
        || derivative_backup.size() != derivative.size())
      return;
    std::memcpy(state.data(), state_backup.data(), state.size_bytes());
    std::memcpy(derivative.data(),
                derivative_backup.data(),
                derivative.size_bytes());
  }

} // namespace

slide::Status CyclerV2::configure(SpmBatch &batch,
                                  CyclerIntegrator integrator)
{
  if (in_run_ || !batch.valid() || batch.n_lanes() != 1
      || !(integrator == CyclerIntegrator::euler_legacy
           || integrator == CyclerIntegrator::exponential))
    return slide::Status::Invalid_parameters;
  static_assert(std::is_nothrow_move_assignable_v<EulerLegacy>);
  static_assert(std::is_nothrow_move_assignable_v<ExponentialModal>);
  static_assert(
    std::is_nothrow_move_assignable_v<std::vector<real_t>>);
  try {
    EulerLegacy candidate_euler;
    ExponentialModal candidate_exponential;
    auto status = candidate_euler.configure(batch);
    if (status != slide::Status::Success)
      return status;
    status = candidate_exponential.configure(batch);
    if (status != slide::Status::Success)
      return status;
    std::vector<real_t> candidate_density(1, 0.0);
    std::vector<real_t> candidate_event_backup(batch.state().size(), 0.0);
    std::vector<real_t> candidate_event_derivative_backup(
      batch.derivative().size(), 0.0);
    std::vector<real_t> candidate_run_state_backup(batch.state().size(), 0.0);
    std::vector<real_t> candidate_run_derivative_backup(
      batch.derivative().size(), 0.0);

    euler_ = std::move(candidate_euler);
    exponential_ = std::move(candidate_exponential);
    density_ = std::move(candidate_density);
    event_backup_ = std::move(candidate_event_backup);
    event_derivative_backup_ = std::move(candidate_event_derivative_backup);
    run_state_backup_ = std::move(candidate_run_state_backup);
    run_derivative_backup_ = std::move(candidate_run_derivative_backup);
    integrator_ = integrator;
    batch_ = &batch;
    return slide::Status::Success;
  } catch (const std::bad_alloc &) {
    return allocationFailureStatus();
  } catch (const std::length_error &) {
    return allocationFailureStatus();
  }
}

slide::Status CyclerV2::registerDriveCycle(const DriveCycle &cycle)
try {
  if (in_run_)
    return slide::Status::Invalid_parameters;
  if (cycle.name.empty() || cycle.time.size() != cycle.current.size()
      || cycle.time.size() < 2 || cycle.time.front() != 0.0)
    return slide::Status::Invalid_parameters;
  for (std::size_t i = 0; i < cycle.time.size(); ++i)
    if (!is_finite(cycle.time[i]) || !is_finite(cycle.current[i])
        || (i > 0 && cycle.time[i] <= cycle.time[i - 1]))
      return slide::Status::Invalid_parameters;
  if (findDriveCycle(cycle.name) != nullptr)
    return slide::Status::Invalid_parameters;
  static_assert(std::is_nothrow_move_constructible_v<DriveCycle>);
  drive_cycles_.push_back(cycle);
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  return allocationFailureStatus();
} catch (const std::length_error &) {
  return allocationFailureStatus();
}

slide::Status CyclerV2::voltageAt(real_t current, real_t &voltage)
{
  density_[0] = current / batch_->electrode_area();
  std::array<real_t, 1> output{};
  const auto status = batch_->terminalVoltage({ .i_app = density_ }, output);
  voltage = output[0];
  return status;
}

slide::Status CyclerV2::currentForVoltage(real_t target, real_t &current)
{
  std::array<real_t, 1> input{}, intercept{}, resistance{};
  for (int iteration = 0; iteration < 12; ++iteration) {
    input[0] = current;
    const auto status = batch_->linearizeThevenin(input, intercept, resistance);
    if (status != slide::Status::Success)
      return status;
    const real_t next = (intercept[0] - target) / resistance[0];
    if (!is_finite(next))
      return slide::Status::Invalid_states;
    if (std::abs(next - current) <= 1e-10) {
      current = next;
      return slide::Status::Success;
    }
    current = next;
  }
  return slide::Status::Numerical_failure;
}

slide::Status CyclerV2::currentForPower(real_t target_power,
                                        Direction direction,
                                        real_t &current)
{
  const real_t signed_power = static_cast<real_t>(direction) * target_power;
  if (!(target_power > 0.0) || direction == Direction::none)
    return slide::Status::Invalid_parameters;
  if (current * static_cast<real_t>(direction) <= 0.0)
    current = signed_power / 3.7;
  for (int iteration = 0; iteration < 12; ++iteration) {
    std::array<real_t, 1> input{ current }, intercept{}, resistance{};
    auto status = batch_->linearizeThevenin(input, intercept, resistance);
    if (status != slide::Status::Success)
      return status;
    const real_t voltage = intercept[0] - resistance[0] * current;
    const real_t residual = voltage * current - signed_power;
    const real_t derivative = voltage - resistance[0] * current;
    if (!(is_finite(derivative) && std::abs(derivative) > 1e-12))
      return slide::Status::Numerical_failure;
    const real_t next = current - residual / derivative;
    if (std::abs(next - current) <= 1e-10) {
      current = next;
      return slide::Status::Success;
    }
    current = next;
  }
  return slide::Status::Numerical_failure;
}

slide::Status CyclerV2::evaluateFunction(const ExperimentFunction &function,
                                         real_t time,
                                         real_t local_time,
                                         real_t voltage,
                                         real_t current,
                                         real_t &value) const
{
  if (!function)
    return slide::Status::Invalid_parameters;
  try {
    value = function({ .time = time,
                       .local_time = local_time,
                       .voltage = voltage,
                       .current = current,
                       .power = voltage * current });
  } catch (const std::bad_alloc &) {
    throw;
  } catch (const std::length_error &) {
    throw;
  } catch (...) {
    return slide::Status::Invalid_parameters;
  }
  return is_finite(value) ? slide::Status::Success
                          : slide::Status::Invalid_states;
}

slide::Status CyclerV2::currentForCustom(const ExperimentSegment &segment,
                                         real_t time,
                                         real_t local_time,
                                         real_t &current)
{
  if (!segment.custom_control)
    return slide::Status::Invalid_parameters;
  const real_t maximum_current = 100.0 * batch_->capacity_Ah();
  current = std::clamp(current, -maximum_current, maximum_current);
  for (int iteration = 0; iteration < 24; ++iteration) {
    real_t voltage{};
    auto status = voltageAt(current, voltage);
    if (status != slide::Status::Success)
      return status;
    real_t residual{};
    status = evaluateFunction(segment.custom_control,
                              time,
                              local_time,
                              voltage,
                              current,
                              residual);
    if (status != slide::Status::Success)
      return status;
    if (std::abs(residual) <= 1e-10)
      return slide::Status::Success;

    const real_t h = std::sqrt(std::numeric_limits<real_t>::epsilon())
                     * std::max(real_t{ 1.0 }, std::abs(current));
    const real_t plus_current = std::min(maximum_current, current + h);
    const real_t minus_current = std::max(-maximum_current, current - h);
    assert(plus_current > minus_current);
    real_t plus_voltage{}, minus_voltage{};
    status = voltageAt(plus_current, plus_voltage);
    if (status != slide::Status::Success)
      return status;
    status = voltageAt(minus_current, minus_voltage);
    if (status != slide::Status::Success)
      return status;
    real_t plus_residual{}, minus_residual{};
    status = evaluateFunction(segment.custom_control,
                              time,
                              local_time,
                              plus_voltage,
                              plus_current,
                              plus_residual);
    if (status != slide::Status::Success)
      return status;
    status = evaluateFunction(segment.custom_control,
                              time,
                              local_time,
                              minus_voltage,
                              minus_current,
                              minus_residual);
    if (status != slide::Status::Success)
      return status;
    const real_t derivative = (plus_residual - minus_residual)
                              / (plus_current - minus_current);
    if (!(is_finite(derivative) && std::abs(derivative) > 1e-14))
      return slide::Status::Numerical_failure;
    const real_t next = std::clamp(current - residual / derivative,
                                   -maximum_current,
                                   maximum_current);
    if (std::abs(next - current) <= 1e-10) {
      current = next;
      return slide::Status::Success;
    }
    current = next;
  }
  return slide::Status::Numerical_failure;
}

slide::Status CyclerV2::advance(real_t current, real_t time, real_t dt)
{
  density_[0] = current / batch_->electrode_area();
  return integrator_ == CyclerIntegrator::exponential
           ? exponential_.step(*batch_, density_, time, dt)
           : euler_.step(*batch_, density_, time, dt);
}

const DriveCycle *CyclerV2::findDriveCycle(const std::string &name) const
{
  const auto found = std::find_if(
    drive_cycles_.begin(), drive_cycles_.end(), [&](const auto &cycle) {
      return detail::normalizedEqual(cycle.name, name);
    });
  return found == drive_cycles_.end() ? nullptr : &*found;
}

real_t CyclerV2::driveCurrent(const DriveCycle &cycle, real_t local_time) const
{
  if (local_time <= 0.0)
    return cycle.current.front();
  if (local_time >= cycle.time.back())
    return cycle.current.back();
  const auto upper = std::upper_bound(cycle.time.begin(), cycle.time.end(), local_time);
  const auto i = static_cast<std::size_t>(upper - cycle.time.begin() - 1);
  const real_t fraction = (local_time - cycle.time[i])
                          / (cycle.time[i + 1] - cycle.time[i]);
  return cycle.current[i] + fraction * (cycle.current[i + 1] - cycle.current[i]);
}

slide::Status CyclerV2::run(const Experiment &experiment,
                            real_t sample_step,
                            ExperimentSolution &output)
try {
  if (in_run_ || batch_ == nullptr || experiment.segments.empty()
      || !is_finite(sample_step) || !(sample_step > 0.0))
    return slide::Status::Invalid_parameters;
  RunActivity active_run{ in_run_ };
  run_snapshot_ready_ = false;
  const bool has_schedule = std::any_of(
    experiment.segments.begin(), experiment.segments.end(), [](const auto &segment) {
      return is_finite(segment.scheduled_start) && segment.scheduled_start >= 0.0;
    });
  if (has_schedule
      && (!(is_finite(experiment.segments.front().scheduled_start)
            && experiment.segments.front().scheduled_start >= 0.0)
          || experiment.segments.front().scheduled_start != 0.0))
    return slide::Status::Invalid_parameters;
  real_t previous_scheduled_start{};
  for (const auto &segment : experiment.segments) {
    if (!detail::validSegment(segment))
      return slide::Status::Invalid_parameters;
    if (segment.mode == ControlMode::drive_cycle
        && findDriveCycle(segment.drive_cycle) == nullptr)
      return slide::Status::Invalid_parameters;
    if (is_finite(segment.scheduled_start) && segment.scheduled_start >= 0.0) {
      if (segment.scheduled_start < previous_scheduled_start)
        return slide::Status::Invalid_parameters;
      previous_scheduled_start = segment.scheduled_start;
    }
  }

  const auto state = batch_->state().raw();
  const auto derivative = batch_->derivative().raw();
  if (run_state_backup_.size() != state.size()
      || run_derivative_backup_.size() != derivative.size())
    return slide::Status::Numerical_failure;
  std::memcpy(run_state_backup_.data(), state.data(), state.size_bytes());
  std::memcpy(run_derivative_backup_.data(),
              derivative.data(),
              derivative.size_bytes());
  run_snapshot_ready_ = true;
  static_assert(std::is_nothrow_move_assignable_v<ExperimentSolution>);
  ExperimentSolution solution;
  real_t time = batch_->state().at(batch_->layout().elapsed_time, 0, 0);
  const real_t schedule_origin = time;
  real_t current{};
  real_t voltage{};
  auto status = voltageAt(0.0, voltage);
  if (status != slide::Status::Success)
    return status;
  solution.time.push_back(time);
  solution.voltage.push_back(voltage);
  solution.current.push_back(0.0);
  solution.sample_segment.push_back(0);

  auto appendRestUntil = [&](real_t target, std::size_t segment_index) {
    while (time < target) {
      const real_t dt = std::min(sample_step, target - time);
      const auto rest_status = advance(0.0, time, dt);
      if (rest_status != slide::Status::Success)
        return rest_status;
      time += dt;
      real_t rest_voltage{};
      const auto voltage_status = voltageAt(0.0, rest_voltage);
      if (voltage_status != slide::Status::Success)
        return voltage_status;
      solution.time.push_back(time);
      solution.voltage.push_back(rest_voltage);
      solution.current.push_back(0.0);
      solution.sample_segment.push_back(segment_index);
    }
    current = 0.0;
    return slide::Status::Success;
  };

  for (std::size_t segment_index = 0;
       segment_index < experiment.segments.size();
       ++segment_index) {
    const auto &segment = experiment.segments[segment_index];
    const real_t segment_step = is_finite(segment.sample_period)
                                    && segment.sample_period > 0.0
                                  ? segment.sample_period
                                  : sample_step;
    if (!(segment_step > 0.0)) {
      status = slide::Status::Invalid_parameters;
      break;
    }
    if (is_finite(segment.scheduled_start) && segment.scheduled_start >= 0.0) {
      status = appendRestUntil(schedule_origin + segment.scheduled_start,
                               segment_index);
      if (status != slide::Status::Success)
        break;
    }
    const auto *cycle = segment.mode == ControlMode::drive_cycle
                          ? findDriveCycle(segment.drive_cycle)
                          : nullptr;
    if (segment.mode == ControlMode::drive_cycle && cycle == nullptr) {
      status = slide::Status::Invalid_parameters;
      break;
    }
    const real_t duration = segment.mode == ControlMode::drive_cycle
                              ? cycle->time.back()
                              : segment.duration;
    real_t horizon = duration > 0.0 ? duration : 7.0 * 24.0 * 3600.0;
    bool cut_by_schedule{};
    for (std::size_t next = segment_index + 1;
         next < experiment.segments.size();
         ++next) {
      const real_t next_start = experiment.segments[next].scheduled_start;
      if (!(is_finite(next_start) && next_start >= 0.0))
        continue;
      const real_t until_next = std::max(real_t{},
                                         schedule_origin + next_start - time);
      if (until_next < horizon) {
        horizon = until_next;
        cut_by_schedule = true;
      }
      break;
    }
    real_t local_time{};
    bool event_reached{};
    std::string event_name;
    while (local_time < horizon) {
      const real_t dt = std::min(segment_step, horizon - local_time);
      if (segment.mode == ControlMode::current) {
        const real_t magnitude = segment.value_is_c_rate
                                   ? segment.value * batch_->capacity_Ah()
                                   : segment.value;
        current = static_cast<real_t>(segment.direction) * magnitude;
      } else if (segment.mode == ControlMode::voltage) {
        status = currentForVoltage(segment.value, current);
      } else if (segment.mode == ControlMode::power) {
        status = currentForPower(segment.value, segment.direction, current);
      } else if (segment.mode == ControlMode::rest) {
        current = 0.0;
      } else if (segment.mode == ControlMode::drive_cycle) {
        current = driveCurrent(*cycle, local_time);
      } else if (segment.mode == ControlMode::custom_explicit) {
        status = voltageAt(current, voltage);
        if (status == slide::Status::Success)
          status = evaluateFunction(segment.custom_control,
                                    time,
                                    local_time,
                                    voltage,
                                    current,
                                    current);
      } else if (segment.mode == ControlMode::custom_implicit) {
        status = currentForCustom(segment, time, local_time, current);
      } else {
        status = voltageAt(current, voltage);
        real_t derivative{};
        if (status == slide::Status::Success)
          status = evaluateFunction(segment.custom_control,
                                    time,
                                    local_time,
                                    voltage,
                                    current,
                                    derivative);
        if (status == slide::Status::Success) {
          current += derivative * dt;
          if (!is_finite(current))
            status = slide::Status::Invalid_states;
        }
      }
      if (status != slide::Status::Success)
        break;
      status = voltageAt(current, voltage);
      if (status != slide::Status::Success)
        break;
      // Algebraic controls apply instantaneously. Match PyBaMM's solution
      // convention by reporting the first sample under the first control,
      // rather than a synthetic open-circuit point at the same state.
      if (segment_index == 0 && local_time == 0.0) {
        solution.voltage.front() = voltage;
        solution.current.front() = current;
      }
      auto indicator = [&](real_t event_time,
                           real_t event_local_time,
                           real_t event_voltage,
                           real_t event_current,
                           real_t &value,
                           std::string &name) {
        value = std::numeric_limits<real_t>::max();
        if (segment.voltage_limit > 0.0) {
          if (segment.direction == Direction::none)
            return slide::Status::Invalid_parameters;
          value = segment.direction == Direction::charge
                    ? segment.voltage_limit - event_voltage
                    : event_voltage - segment.voltage_limit;
          name = "voltage cut-off";
        }
        if (segment.current_cutoff > 0.0) {
          const real_t cutoff = segment.cutoff_is_c_rate
                                  ? segment.current_cutoff * batch_->capacity_Ah()
                                  : segment.current_cutoff;
          const real_t current_value = std::abs(event_current) - cutoff;
          if (current_value < value) {
            value = current_value;
            name = "current cut-off";
          }
        }
        for (const auto &termination : segment.custom_terminations) {
          real_t custom_value{};
          const auto custom_status = evaluateFunction(termination.indicator,
                                                      event_time,
                                                      event_local_time,
                                                      event_voltage,
                                                      event_current,
                                                      custom_value);
          if (custom_status != slide::Status::Success)
            return custom_status;
          if (custom_value < value) {
            value = custom_value;
            name = termination.name;
          }
        }
        return slide::Status::Success;
      };
      real_t before_indicator{};
      std::string before_event_name;
      status = indicator(time,
                         local_time,
                         voltage,
                         current,
                         before_indicator,
                         before_event_name);
      if (status != slide::Status::Success)
        break;
      if (before_indicator <= 0.0) {
        event_reached = true;
        event_name = std::move(before_event_name);
        solution.voltage.back() = voltage;
        solution.current.back() = current;
        solution.segment = segment_index;
        break;
      }
      const auto event_state = batch_->state().raw();
      const auto event_derivative = batch_->derivative().raw();
      assert(event_backup_.size() == event_state.size());
      assert(event_derivative_backup_.size() == event_derivative.size());
      std::memcpy(event_backup_.data(),
                  event_state.data(),
                  event_state.size_bytes());
      std::memcpy(event_derivative_backup_.data(),
                  event_derivative.data(),
                  event_derivative.size_bytes());
      const auto restore_event_state = [&]() noexcept {
        std::memcpy(event_state.data(),
                    event_backup_.data(),
                    event_state.size_bytes());
        std::memcpy(event_derivative.data(),
                    event_derivative_backup_.data(),
                    event_derivative.size_bytes());
      };
      status = advance(current, time, dt);
      if (status != slide::Status::Success) {
        restore_event_state();
        break;
      }
      real_t after_voltage{};
      status = voltageAt(current, after_voltage);
      real_t after_current = current;
      if (status == slide::Status::Success && segment.mode == ControlMode::voltage)
        status = currentForVoltage(segment.value, after_current);
      else if (status == slide::Status::Success
               && segment.mode == ControlMode::power)
        status = currentForPower(segment.value, segment.direction, after_current);
      else if (status == slide::Status::Success
               && segment.mode == ControlMode::custom_explicit)
        status = evaluateFunction(segment.custom_control,
                                  time + dt,
                                  local_time + dt,
                                  after_voltage,
                                  after_current,
                                  after_current);
      else if (status == slide::Status::Success
               && segment.mode == ControlMode::custom_implicit)
        status = currentForCustom(segment,
                                  time + dt,
                                  local_time + dt,
                                  after_current);
      if (status == slide::Status::Success && after_current != current)
        status = voltageAt(after_current, after_voltage);
      if (status != slide::Status::Success) {
        restore_event_state();
        break;
      }
      real_t after_indicator{};
      std::string after_event_name;
      status = indicator(time + dt,
                         local_time + dt,
                         after_voltage,
                         after_current,
                         after_indicator,
                         after_event_name);
      if (status != slide::Status::Success) {
        restore_event_state();
        break;
      }
      real_t accepted_dt = dt;
      if (before_indicator > 0.0 && after_indicator <= 0.0) {
        real_t low{}, high = dt;
        for (int iteration = 0; iteration < 45; ++iteration) {
          const real_t middle = 0.5 * (low + high);
          restore_event_state();
          status = advance(current, time, middle);
          if (status != slide::Status::Success)
            break;
          real_t middle_voltage{};
          status = voltageAt(current, middle_voltage);
          real_t middle_current = current;
          if (status == slide::Status::Success
              && segment.mode == ControlMode::voltage)
            status = currentForVoltage(segment.value, middle_current);
          else if (status == slide::Status::Success
                   && segment.mode == ControlMode::power)
            status = currentForPower(segment.value, segment.direction, middle_current);
          else if (status == slide::Status::Success
                   && segment.mode == ControlMode::custom_explicit)
            status = evaluateFunction(segment.custom_control,
                                      time + middle,
                                      local_time + middle,
                                      middle_voltage,
                                      middle_current,
                                      middle_current);
          else if (status == slide::Status::Success
                   && segment.mode == ControlMode::custom_implicit)
            status = currentForCustom(segment,
                                      time + middle,
                                      local_time + middle,
                                      middle_current);
          if (status == slide::Status::Success && middle_current != current)
            status = voltageAt(middle_current, middle_voltage);
          if (status != slide::Status::Success)
            break;
          real_t middle_indicator{};
          std::string middle_event_name;
          status = indicator(time + middle,
                             local_time + middle,
                             middle_voltage,
                             middle_current,
                             middle_indicator,
                             middle_event_name);
          if (status != slide::Status::Success)
            break;
          if (middle_indicator > 0.0)
            low = middle;
          else
            high = middle;
        }
        if (status != slide::Status::Success) {
          restore_event_state();
          break;
        }
        restore_event_state();
        status = advance(current, time, high);
        if (status != slide::Status::Success) {
          restore_event_state();
          break;
        }
        accepted_dt = high;
        status = voltageAt(current, after_voltage);
        if (segment.mode == ControlMode::voltage)
          status = currentForVoltage(segment.value, after_current);
        else if (segment.mode == ControlMode::power)
          status = currentForPower(segment.value, segment.direction, after_current);
        else if (segment.mode == ControlMode::custom_explicit)
          status = evaluateFunction(segment.custom_control,
                                    time + high,
                                    local_time + high,
                                    after_voltage,
                                    after_current,
                                    after_current);
        else if (segment.mode == ControlMode::custom_implicit)
          status = currentForCustom(segment,
                                    time + high,
                                    local_time + high,
                                    after_current);
        if (status == slide::Status::Success && after_current != current)
          status = voltageAt(after_current, after_voltage);
        if (status == slide::Status::Success) {
          real_t final_indicator{};
          status = indicator(time + high,
                             local_time + high,
                             after_voltage,
                             after_current,
                             final_indicator,
                             event_name);
        }
        event_reached = status == slide::Status::Success;
      }
      if (status != slide::Status::Success) {
        restore_event_state();
        break;
      }
      time += accepted_dt;
      local_time += accepted_dt;
      solution.time.push_back(time);
      solution.voltage.push_back(after_voltage);
      solution.current.push_back(after_current);
      solution.sample_segment.push_back(segment_index);
      solution.segment = segment_index;
      if (event_reached)
        break;
    }
    if (status != slide::Status::Success) {
      solution.reason = TerminationReason::error;
      solution.status = status;
      output = std::move(solution);
      return status;
    }
    if (event_reached) {
      solution.reason = TerminationReason::event;
      solution.termination_name = std::move(event_name);
    } else if (duration <= 0.0 && !cut_by_schedule) {
      solution.reason = TerminationReason::limit;
      solution.status = slide::Status::Numerical_failure;
      solution.termination_name = "maximum step duration";
      output = std::move(solution);
      return slide::Status::Numerical_failure;
    } else {
      solution.reason = TerminationReason::final_time;
      solution.termination_name = cut_by_schedule ? "next scheduled start"
                                                  : "final time";
    }
  }
  if (status != slide::Status::Success) {
    solution.reason = TerminationReason::error;
    solution.status = status;
    output = std::move(solution);
    return status;
  }
  solution.status = slide::Status::Success;
  output = std::move(solution);
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  if (run_snapshot_ready_ && batch_ != nullptr)
    restoreBatchSnapshot(
      *batch_, run_state_backup_, run_derivative_backup_);
  return allocationFailureStatus();
} catch (const std::length_error &) {
  if (run_snapshot_ready_ && batch_ != nullptr)
    restoreBatchSnapshot(
      *batch_, run_state_backup_, run_derivative_backup_);
  return allocationFailureStatus();
}

} // namespace slide::core
