/**
 * @file ExperimentSemantics.hpp
 * @brief Single-source normalisation and validation shared by experiment TUs.
 *
 * MC-3 contract: these inline definitions are the only parser/runner
 * implementations of normalisation, normalised name equality, direction
 * validation, and segment validation. Keep this header dependency-light and
 * self-contained; do not add parser or runner state here.
 */

#pragma once

#include "../Experiment.hpp"

#include <cctype>
#include <cstddef>
#include <string>
#include <string_view>

namespace slide::core::detail {

inline std::string normalized(std::string_view text)
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

inline bool normalizedEqual(std::string_view left,
                            std::string_view right) noexcept
{
  struct Cursor
  {
    std::string_view text;
    std::size_t position{};
    bool emitted{};
    bool pending_space{};

    bool next(char &output) noexcept
    {
      while (position < text.size()) {
        const auto value = static_cast<unsigned char>(text[position++]);
        if (std::isspace(value)) {
          pending_space = emitted;
          continue;
        }
        if (pending_space) {
          pending_space = false;
          --position;
          output = ' ';
          return true;
        }
        output = static_cast<char>(std::tolower(value));
        emitted = true;
        return true;
      }
      return false;
    }
  } left_cursor{ left }, right_cursor{ right };

  while (true) {
    char left_value{}, right_value{};
    const bool has_left = left_cursor.next(left_value);
    const bool has_right = right_cursor.next(right_value);
    if (has_left != has_right)
      return false;
    if (!has_left)
      return true;
    if (left_value != right_value)
      return false;
  }
}

inline bool validDirection(Direction direction)
{
  return direction == Direction::charge || direction == Direction::none
         || direction == Direction::discharge;
}

inline bool validSegment(const ExperimentSegment &segment)
{
  if (!validDirection(segment.direction) || !is_finite(segment.value)
      || !is_finite(segment.duration) || segment.duration < 0.0
      || !is_finite(segment.voltage_limit) || segment.voltage_limit < 0.0
      || !is_finite(segment.current_cutoff) || segment.current_cutoff < 0.0
      || !is_finite(segment.scheduled_start)
      || !is_finite(segment.sample_period) || segment.sample_period == 0.0)
    return false;

  const bool custom_mode = segment.mode == ControlMode::custom_explicit
                           || segment.mode == ControlMode::custom_implicit
                           || segment.mode == ControlMode::custom_differential;
  if (custom_mode != static_cast<bool>(segment.custom_control))
    return false;
  for (const auto &termination : segment.custom_terminations)
    if (termination.name.empty() || !termination.indicator)
      return false;

  switch (segment.mode) {
  case ControlMode::current:
  case ControlMode::power:
    if (!(segment.value > 0.0)
        || !(segment.direction == Direction::charge
             || segment.direction == Direction::discharge))
      return false;
    break;
  case ControlMode::voltage:
    if (!(segment.value > 0.0))
      return false;
    break;
  case ControlMode::rest:
    if (segment.direction != Direction::none)
      return false;
    break;
  case ControlMode::drive_cycle:
    if (segment.drive_cycle.empty())
      return false;
    break;
  case ControlMode::custom_explicit:
  case ControlMode::custom_implicit:
  case ControlMode::custom_differential:
    break;
  default:
    return false;
  }

  if (segment.voltage_limit > 0.0 && segment.direction == Direction::none)
    return false;
  const bool has_event = segment.voltage_limit > 0.0
                         || segment.current_cutoff > 0.0
                         || !segment.custom_terminations.empty();
  return segment.mode == ControlMode::drive_cycle
         || segment.duration > 0.0 || has_event;
}

} // namespace slide::core::detail
