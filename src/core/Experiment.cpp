/**
 * @file Experiment.cpp
 * @brief Allocation-bounded parser for the PyBaMM-style experiment grammar.
 *
 * MC-3 ownership contract: this translation unit owns parser mechanics only.
 * Segment semantics and text normalisation are single-sourced in
 * detail/ExperimentSemantics.hpp; execution belongs to CyclerV2.cpp.
 */

#include "Experiment.hpp"
#include "detail/ExperimentSemantics.hpp"

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cmath>
#include <new>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace slide::core {
namespace {

  void assignDiagnosticNoThrow(std::string &target,
                               std::string_view message) noexcept
  {
    try {
      target.assign(message);
    } catch (...) {
      target.clear();
    }
  }

  slide::Status parserAllocationFailure(ParseDiagnostic &diagnostic,
                                        std::string_view message) noexcept
  {
    assignDiagnosticNoThrow(diagnostic.message, message);
    return slide::Status::Numerical_failure;
  }

  bool parseNumber(std::string_view text, real_t &value)
  {
    if (text.empty())
      return false;
    const auto parsed = std::from_chars(
      text.data(), text.data() + text.size(), value, std::chars_format::general);
    return parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size()
           && is_finite(value);
  }

  bool parseDuration(std::string_view text, real_t &seconds)
  {
    std::size_t number_end{};
    while (number_end < text.size() && text[number_end] >= '0'
           && text[number_end] <= '9')
      ++number_end;
    if (number_end == 0)
      return false;
    if (number_end < text.size() && text[number_end] == '.') {
      const std::size_t fraction = ++number_end;
      while (number_end < text.size() && text[number_end] >= '0'
             && text[number_end] <= '9')
        ++number_end;
      if (number_end == fraction)
        return false;
    }
    const auto number = text.substr(0, number_end);
    if (number_end < text.size() && text[number_end] == ' ')
      ++number_end;
    const auto unit = text.substr(number_end);
    if (!parseNumber(number, seconds)
        || !(unit == "s" || unit == "sec" || unit == "secs"
             || unit == "second" || unit == "seconds" || unit == "min"
             || unit == "mins" || unit == "minute" || unit == "minutes"
             || unit == "h" || unit == "hr" || unit == "hrs"
             || unit == "hour" || unit == "hours"))
      return false;
    if (unit == "min" || unit == "mins" || unit.starts_with("minute"))
      seconds *= 60.0;
    else if (unit == "h" || unit == "hr" || unit == "hrs"
             || unit.starts_with("hour"))
      seconds *= 3600.0;
    return is_finite(seconds) && seconds > 0.0;
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
try {
  constexpr std::size_t max_expanded_segments = 10'000;
  constexpr std::size_t max_expanded_text_bytes = 4U * 1024U * 1024U;
  constexpr std::size_t max_drive_cycle_name_bytes = 1024;
  constexpr std::size_t max_step_bytes = 65'536;
  diagnostic = {};
  if (steps.empty()) {
    diagnostic.message = "experiment must contain at least one step";
    return slide::Status::Invalid_parameters;
  }
  if (steps.size() > max_expanded_segments) {
    diagnostic.message = "experiment expanded segment limit (10000) exceeded";
    return slide::Status::Invalid_parameters;
  }
  std::vector<ExperimentSegment> parsed;
  parsed.reserve(steps.size());
  std::size_t expanded_text_bytes{};
  for (std::size_t index = 0; index < steps.size(); ++index) {
    if (steps[index].size() > max_step_bytes) {
      diagnostic.step = index;
      diagnostic.message = "experiment step exceeds 65536 bytes";
      return slide::Status::Invalid_parameters;
    }
    auto text = detail::normalized(steps[index]);
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
    } else if (text.starts_with("rest for ")) {
      segment.mode = ControlMode::rest;
      segment.direction = Direction::none;
      valid = parseDuration(
        std::string_view{ text }.substr(std::string_view{ "rest for " }.size()),
        segment.duration);
    } else if (text.starts_with("run ") && text.ends_with(" (a)")) {
      segment.mode = ControlMode::drive_cycle;
      constexpr std::size_t prefix_size = std::string_view{ "run " }.size();
      constexpr std::size_t suffix_size = std::string_view{ " (a)" }.size();
      if (text.size() > prefix_size + suffix_size) {
        segment.drive_cycle = text.substr(
          prefix_size, text.size() - prefix_size - suffix_size);
        valid = true;
      }
    }
    if (!valid) {
      diagnostic.step = index;
      diagnostic.offset = 0;
      diagnostic.message = "malformed or unsupported experiment step: " + steps[index];
      return slide::Status::Invalid_parameters;
    }
    if (segment.drive_cycle.size() > max_drive_cycle_name_bytes) {
      diagnostic.step = index;
      diagnostic.message = "drive-cycle name exceeds 1024 bytes";
      return slide::Status::Invalid_parameters;
    }
    if (!detail::validSegment(segment)) {
      diagnostic.step = index;
      diagnostic.offset = 0;
      diagnostic.message = "semantically invalid experiment step: " + steps[index];
      return slide::Status::Invalid_parameters;
    }
    if (repetitions > max_expanded_segments - parsed.size()) {
      diagnostic.step = index;
      diagnostic.message = "experiment expanded segment limit (10000) exceeded";
      return slide::Status::Invalid_parameters;
    }
    // The independently enforced 65536-byte step and 1024-byte drive-cycle
    // limits make this addition representable and well below the 4 MiB cap.
    assert(segment.source.size()
           <= max_expanded_text_bytes - segment.drive_cycle.size());
    const std::size_t retained_text = segment.source.size()
                                      + segment.drive_cycle.size();
    if (retained_text > 0
        && repetitions
             > (max_expanded_text_bytes - expanded_text_bytes) / retained_text) {
      diagnostic.step = index;
      diagnostic.message = "experiment expanded text limit (4194304 bytes) exceeded";
      return slide::Status::Invalid_parameters;
    }
    expanded_text_bytes += repetitions * retained_text;
    parsed.reserve(parsed.size() + repetitions);
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
      parsed.push_back(segment);
  }
  output.segments = std::move(parsed);
  return slide::Status::Success;
} catch (const std::bad_alloc &) {
  return parserAllocationFailure(
    diagnostic, "experiment parser allocation failed");
} catch (const std::length_error &) {
  return parserAllocationFailure(
    diagnostic, "experiment parser size is not representable");
}

} // namespace slide::core
