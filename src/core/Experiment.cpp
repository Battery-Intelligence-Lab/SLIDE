/**
 * @file Experiment.cpp
 * @brief PyBaMM-style experiment grammar and event-aligned core cycler.
 */

#include "Experiment.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <regex>

namespace slide::core {
namespace {

  std::string normalized(std::string text)
  {
    std::string result;
    result.reserve(text.size());
    bool space = true;
    for (const unsigned char c : text) {
      if (std::isspace(c)) {
        if (!space)
          result.push_back(' ');
        space = true;
      } else {
        result.push_back(static_cast<char>(std::tolower(c)));
        space = false;
      }
    }
    if (!result.empty() && result.back() == ' ')
      result.pop_back();
    return result;
  }

  bool parseNumber(const std::string &text, real_t &value)
  {
    char *end{};
    value = std::strtod(text.c_str(), &end);
    return end == text.c_str() + static_cast<std::ptrdiff_t>(text.size())
           && is_finite(value);
  }

  bool parseDuration(const std::string &text, real_t &seconds)
  {
    static const std::regex pattern{
      R"(^([0-9]+(?:\.[0-9]+)?) ?(s|sec|secs|second|seconds|min|mins|minute|minutes|h|hr|hrs|hour|hours)$)"
    };
    std::smatch match;
    if (!std::regex_match(text, match, pattern) || !parseNumber(match[1].str(), seconds))
      return false;
    const auto unit = match[2].str();
    if (unit == "min" || unit == "mins" || unit.starts_with("minute"))
      seconds *= 60.0;
    else if (unit == "h" || unit == "hr" || unit == "hrs"
             || unit.starts_with("hour"))
      seconds *= 3600.0;
    return seconds > 0.0;
  }

  struct Quantity
  {
    real_t value{};
    enum class Unit : unsigned char { ampere,
                                      c_rate,
                                      volt,
                                      watt } unit{};
  };

  bool parseQuantity(std::string text, Quantity &quantity)
  {
    text.erase(std::remove(text.begin(), text.end(), ' '), text.end());
    if (text.starts_with("c/")) {
      real_t denominator{};
      if (!parseNumber(text.substr(2), denominator) || !(denominator > 0.0))
        return false;
      quantity = { 1.0 / denominator, Quantity::Unit::c_rate };
      return true;
    }
    real_t scale = 1.0;
    Quantity::Unit unit;
    std::size_t suffix{};
    if (text.ends_with("ma")) {
      unit = Quantity::Unit::ampere;
      scale = 1e-3;
      suffix = 2;
    } else if (text.ends_with('a')) {
      unit = Quantity::Unit::ampere;
      suffix = 1;
    } else if (text.ends_with('c')) {
      unit = Quantity::Unit::c_rate;
      suffix = 1;
    } else if (text.ends_with('v')) {
      unit = Quantity::Unit::volt;
      suffix = 1;
    } else if (text.ends_with('w')) {
      unit = Quantity::Unit::watt;
      suffix = 1;
    } else {
      return false;
    }
    real_t value{};
    if (!parseNumber(text.substr(0, text.size() - suffix), value))
      return false;
    quantity = { scale * value, unit };
    return is_finite(quantity.value);
  }

  bool parseTermination(const std::string &keyword,
                        const std::string &value,
                        ExperimentSegment &segment)
  {
    if (keyword == "for")
      return parseDuration(value, segment.duration);
    Quantity quantity;
    if (!parseQuantity(value, quantity))
      return false;
    if (quantity.unit == Quantity::Unit::volt) {
      segment.voltage_limit = quantity.value;
      return quantity.value > 0.0;
    }
    if (quantity.unit == Quantity::Unit::ampere
        || quantity.unit == Quantity::Unit::c_rate) {
      segment.current_cutoff = std::abs(quantity.value);
      segment.cutoff_is_c_rate = quantity.unit == Quantity::Unit::c_rate;
      return segment.current_cutoff > 0.0;
    }
    return false;
  }

  bool parseConditions(const std::string &text, ExperimentSegment &segment)
  {
    if (text.starts_with("for ")) {
      constexpr std::string_view separator{ " or until " };
      const auto event = text.find(separator);
      const auto duration = event == std::string::npos
                              ? text.substr(4)
                              : text.substr(4, event - 4);
      if (!parseDuration(duration, segment.duration))
        return false;
      return event == std::string::npos
             || parseTermination("until",
                                 text.substr(event + separator.size()),
                                 segment);
    }
    if (text.starts_with("until "))
      return parseTermination("until", text.substr(6), segment);
    return false;
  }

  bool splitControlledStep(const std::string &text,
                           std::string &control,
                           std::string &conditions)
  {
    const auto for_position = text.find(" for ");
    const auto until_position = text.find(" until ");
    const auto split = std::min(for_position, until_position);
    if (split == std::string::npos)
      return false;
    control = text.substr(0, split);
    conditions = text.substr(split + 1);
    return !control.empty() && !conditions.empty();
  }

  bool positiveControl(const Quantity &quantity)
  {
    return is_finite(quantity.value) && quantity.value > 0.0;
  }

} // namespace

slide::Status Experiment::parse(std::span<const std::string> steps,
                                Experiment &output,
                                ParseDiagnostic &diagnostic)
{
  diagnostic = {};
  if (steps.empty()) {
    diagnostic.message = "experiment must contain at least one step";
    return slide::Status::Invalid_parameters;
  }
  std::vector<ExperimentSegment> parsed;
  parsed.reserve(steps.size());
  static const std::regex rest_pattern{ R"(^rest for (.+)$)" };
  static const std::regex run_pattern{ R"(^run (.+) \(a\)$)" };

  for (std::size_t index = 0; index < steps.size(); ++index) {
    auto text = normalized(steps[index]);
    std::size_t repetitions = 1;
    constexpr std::string_view repeat_separator{ " * " };
    const auto repeat_position = text.rfind(repeat_separator);
    if (repeat_position != std::string::npos) {
      const auto count = text.substr(repeat_position + repeat_separator.size());
      text.resize(repeat_position);
      std::size_t parsed_count{};
      bool count_valid = !count.empty();
      for (const char digit : count) {
        count_valid = count_valid && digit >= '0' && digit <= '9';
        if (!count_valid || parsed_count > 100'000U) {
          count_valid = false;
          break;
        }
        parsed_count = parsed_count * 10U
                       + static_cast<std::size_t>(digit - '0');
      }
      if (!count_valid || parsed_count < 1U || parsed_count > 1'000'000U) {
        diagnostic.step = index;
        diagnostic.message = "invalid experiment repetition count: " + steps[index];
        return slide::Status::Invalid_parameters;
      }
      repetitions = parsed_count;
    }
    ExperimentSegment segment{ .source = steps[index] };
    std::smatch match;
    bool valid{};
    if (text.starts_with("charge at ") || text.starts_with("discharge at ")) {
      const bool charge = text.starts_with("charge at ");
      const auto prefix = charge ? std::string_view{ "charge at " }
                                 : std::string_view{ "discharge at " };
      segment.direction = charge
                            ? Direction::charge
                            : Direction::discharge;
      std::string control, conditions;
      Quantity quantity;
      valid = splitControlledStep(text.substr(prefix.size()), control, conditions)
              && parseQuantity(control, quantity) && positiveControl(quantity)
              && parseConditions(conditions, segment);
      if (valid && (quantity.unit == Quantity::Unit::ampere || quantity.unit == Quantity::Unit::c_rate)) {
        segment.mode = ControlMode::current;
        segment.value = std::abs(quantity.value);
        segment.value_is_c_rate = quantity.unit == Quantity::Unit::c_rate;
      } else if (valid && quantity.unit == Quantity::Unit::watt) {
        segment.mode = ControlMode::power;
        segment.value = std::abs(quantity.value);
      } else {
        valid = false;
      }
    } else if (text.starts_with("hold at ")) {
      std::string control, conditions;
      Quantity quantity;
      valid = splitControlledStep(text.substr(8), control, conditions)
              && parseQuantity(control, quantity) && positiveControl(quantity)
              && quantity.unit == Quantity::Unit::volt
              && parseConditions(conditions, segment);
      if (valid) {
        segment.mode = ControlMode::voltage;
        segment.value = quantity.value;
      }
    } else if (std::regex_match(text, match, rest_pattern)) {
      segment.mode = ControlMode::rest;
      segment.direction = Direction::none;
      valid = parseDuration(match[1].str(), segment.duration);
    } else if (std::regex_match(text, match, run_pattern)) {
      segment.mode = ControlMode::drive_cycle;
      segment.drive_cycle = match[1].str();
      valid = !segment.drive_cycle.empty();
    }
    if (!valid) {
      diagnostic.step = index;
      diagnostic.offset = 0;
      diagnostic.message = "malformed or unsupported experiment step: " + steps[index];
      return slide::Status::Invalid_parameters;
    }
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
      parsed.push_back(segment);
  }
  output.segments = std::move(parsed);
  return slide::Status::Success;
}

slide::Status CyclerV2::configure(SpmBatch &batch,
                                  CyclerIntegrator integrator)
{
  if (!batch.valid() || batch.n_lanes() != 1)
    return slide::Status::Invalid_parameters;
  auto status = euler_.configure(batch);
  if (status != slide::Status::Success)
    return status;
  status = exponential_.configure(batch);
  if (status != slide::Status::Success)
    return status;
  batch_ = &batch;
  integrator_ = integrator;
  density_.assign(1, 0.0);
  event_backup_.assign(batch.state().size(), 0.0);
  return slide::Status::Success;
}

slide::Status CyclerV2::registerDriveCycle(DriveCycle cycle)
{
  if (cycle.name.empty() || cycle.time.size() != cycle.current.size()
      || cycle.time.size() < 2 || cycle.time.front() != 0.0)
    return slide::Status::Invalid_parameters;
  for (std::size_t i = 0; i < cycle.time.size(); ++i)
    if (!is_finite(cycle.time[i]) || !is_finite(cycle.current[i])
        || (i > 0 && cycle.time[i] <= cycle.time[i - 1]))
      return slide::Status::Invalid_parameters;
  if (findDriveCycle(cycle.name) != nullptr)
    return slide::Status::Invalid_parameters;
  drive_cycles_.push_back(std::move(cycle));
  return slide::Status::Success;
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

slide::Status CyclerV2::advance(real_t current, real_t time, real_t dt)
{
  density_[0] = current / batch_->electrode_area();
  return integrator_ == CyclerIntegrator::exponential
           ? exponential_.step(*batch_, density_, time, dt)
           : euler_.step(*batch_, density_, time, dt);
}

const DriveCycle *CyclerV2::findDriveCycle(const std::string &name) const
{
  const auto found = std::find_if(drive_cycles_.begin(), drive_cycles_.end(), [&](const auto &cycle) { return normalized(cycle.name) == normalized(name); });
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
{
  if (batch_ == nullptr || experiment.segments.empty()
      || !is_finite(sample_step) || !(sample_step > 0.0))
    return slide::Status::Invalid_parameters;
  ExperimentSolution solution;
  real_t time = batch_->state().at(batch_->layout().elapsed_time, 0, 0);
  real_t current{};
  real_t voltage{};
  auto status = voltageAt(0.0, voltage);
  if (status != slide::Status::Success)
    return status;
  solution.time.push_back(time);
  solution.voltage.push_back(voltage);
  solution.current.push_back(0.0);
  solution.sample_segment.push_back(0);

  for (std::size_t segment_index = 0;
       segment_index < experiment.segments.size();
       ++segment_index) {
    const auto &segment = experiment.segments[segment_index];
    const auto *cycle = segment.mode == ControlMode::drive_cycle
                          ? findDriveCycle(segment.drive_cycle)
                          : nullptr;
    if (segment.mode == ControlMode::drive_cycle && cycle == nullptr)
      return slide::Status::Invalid_parameters;
    const real_t duration = segment.mode == ControlMode::drive_cycle
                              ? cycle->time.back()
                              : segment.duration;
    const real_t horizon = duration > 0.0 ? duration : 7.0 * 24.0 * 3600.0;
    real_t local_time{};
    bool event_reached{};
    while (local_time < horizon) {
      const real_t dt = std::min(sample_step, horizon - local_time);
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
      } else {
        current = driveCurrent(*cycle, local_time);
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
      auto indicator = [&](real_t event_voltage, real_t event_current) {
        if (segment.voltage_limit > 0.0)
          return segment.direction == Direction::charge
                   ? segment.voltage_limit - event_voltage
                   : event_voltage - segment.voltage_limit;
        if (segment.current_cutoff > 0.0) {
          const real_t cutoff = segment.cutoff_is_c_rate
                                  ? segment.current_cutoff * batch_->capacity_Ah()
                                  : segment.current_cutoff;
          return std::abs(event_current) - cutoff;
        }
        // Keep the no-event sentinel finite: Release uses -Ofast, under which
        // infinities are explicitly outside the compiler's floating-point model.
        return std::numeric_limits<real_t>::max();
      };
      const real_t before_indicator = indicator(voltage, current);
      if (before_indicator <= 0.0) {
        event_reached = true;
        solution.voltage.back() = voltage;
        solution.current.back() = current;
        solution.segment = segment_index;
        break;
      }
      std::memcpy(event_backup_.data(), batch_->state().raw().data(), batch_->state().raw().size_bytes());
      status = advance(current, time, dt);
      if (status != slide::Status::Success)
        break;
      real_t after_voltage{};
      status = voltageAt(current, after_voltage);
      real_t after_current = current;
      if (status == slide::Status::Success && segment.mode == ControlMode::voltage)
        status = currentForVoltage(segment.value, after_current);
      else if (status == slide::Status::Success
               && segment.mode == ControlMode::power)
        status = currentForPower(segment.value, segment.direction, after_current);
      if (status == slide::Status::Success && after_current != current)
        status = voltageAt(after_current, after_voltage);
      if (status != slide::Status::Success)
        break;
      const real_t after_indicator = indicator(after_voltage, after_current);
      real_t accepted_dt = dt;
      if (before_indicator > 0.0 && after_indicator <= 0.0) {
        real_t low{}, high = dt;
        for (int iteration = 0; iteration < 45; ++iteration) {
          const real_t middle = 0.5 * (low + high);
          std::memcpy(batch_->state().raw().data(), event_backup_.data(), batch_->state().raw().size_bytes());
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
          if (status == slide::Status::Success && middle_current != current)
            status = voltageAt(middle_current, middle_voltage);
          if (status != slide::Status::Success)
            break;
          if (indicator(middle_voltage, middle_current) > 0.0)
            low = middle;
          else
            high = middle;
        }
        if (status != slide::Status::Success)
          break;
        std::memcpy(batch_->state().raw().data(), event_backup_.data(), batch_->state().raw().size_bytes());
        status = advance(current, time, high);
        if (status != slide::Status::Success)
          break;
        accepted_dt = high;
        status = voltageAt(current, after_voltage);
        if (segment.mode == ControlMode::voltage)
          status = currentForVoltage(segment.value, after_current);
        else if (segment.mode == ControlMode::power)
          status = currentForPower(segment.value, segment.direction, after_current);
        if (status == slide::Status::Success && after_current != current)
          status = voltageAt(after_current, after_voltage);
        event_reached = status == slide::Status::Success;
      }
      if (status != slide::Status::Success)
        break;
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
    if (event_reached)
      solution.reason = TerminationReason::event;
    else if (duration <= 0.0) {
      solution.reason = TerminationReason::limit;
      solution.status = slide::Status::Numerical_failure;
      output = std::move(solution);
      return slide::Status::Numerical_failure;
    } else
      solution.reason = TerminationReason::final_time;
  }
  solution.status = slide::Status::Success;
  output = std::move(solution);
  return slide::Status::Success;
}

} // namespace slide::core
