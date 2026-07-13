/**
 * @file CoreSpmTestHarness.hpp
 * @brief M0.8 shared mechanics for successful SPM build, observation, stepping,
 *        and transparent voltage-error measurement in core unit tests.
 *
 * The caller owns every batch, stepper, input span, scratch span, and output.
 * Unit tags prevent amperes from being passed where A/m2 is required.
 *
 * Ownership invariant: `requireSpmBatch` yields a prvalue, so name it as a
 * local before configuring anything against it. `Recorder` and `CyclerV2` store
 * a raw `SpmBatch *` at configure time, and a configured batch must therefore
 * never be moved, reseated, or reallocated inside a container.
 */

#pragma once

#include "../../src/core/ExponentialModal.hpp"
#include "../../src/core/Numeric.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>

namespace slide::test_support {

namespace detail {

  constexpr std::uint64_t sign_mask = UINT64_C(0x8000000000000000);
  constexpr std::uint64_t exponent_mask = UINT64_C(0x7ff0000000000000);
  constexpr std::uint64_t magnitude_mask = UINT64_C(0x7fffffffffffffff);

  [[nodiscard]] inline core::real_t magnitude(
    const core::real_t &value) noexcept
  {
    const auto bits = std::bit_cast<std::uint64_t>(value) & magnitude_mask;
    return std::bit_cast<core::real_t>(bits);
  }

  /** Conservative preflight: an accepted subtraction cannot overflow. */
  [[nodiscard]] inline bool canEvaluateAbsoluteDifference(
    const core::real_t &left,
    const core::real_t &right) noexcept
  {
    if (!(core::is_finite(left) && core::is_finite(right)))
      return false;
    const auto left_bits = std::bit_cast<std::uint64_t>(left);
    const auto right_bits = std::bit_cast<std::uint64_t>(right);
    if (((left_bits ^ right_bits) & sign_mask) == 0)
      return true;

    const core::real_t left_magnitude = magnitude(left);
    const core::real_t right_magnitude = magnitude(right);
    if (left_magnitude == 0.0 || right_magnitude == 0.0)
      return true;
    const core::real_t smaller = std::min(left_magnitude, right_magnitude);
    const core::real_t larger = std::max(left_magnitude, right_magnitude);
    return larger < std::numeric_limits<core::real_t>::max() - smaller;
  }

  /** Preflight a nonnegative square without evaluating the product. */
  [[nodiscard]] inline bool canSquare(const core::real_t &value) noexcept
  {
    if (!(core::is_finite(value) && value >= 0.0))
      return false;
    if (value <= 1.0)
      return true;
    return value < std::numeric_limits<core::real_t>::max() / value;
  }

  /** Preflight `sum + value * value` without evaluating either risky operation. */
  [[nodiscard]] inline bool canAddSquare(const core::real_t &sum,
                                         const core::real_t &value) noexcept
  {
    if (!(core::is_finite(sum) && sum >= 0.0
          && core::is_finite(value) && value >= 0.0))
      return false;
    const auto value_bits = std::bit_cast<std::uint64_t>(value);
    if ((value_bits & magnitude_mask) == 0)
      return true;

    const core::real_t room =
      std::numeric_limits<core::real_t>::max() - sum;
    if (value <= 1.0)
      return room > value;
    return value < room / value;
  }

  /** Preflight finite `numerator / positive normal denominator` without division. */
  [[nodiscard]] inline bool canDivideByPositiveNormal(
    const core::real_t &numerator,
    const core::real_t &positive_denominator) noexcept
  {
    if (!(core::is_finite(numerator)
          && core::is_finite(positive_denominator)))
      return false;
    const auto denominator_bits =
      std::bit_cast<std::uint64_t>(positive_denominator);
    if ((denominator_bits & sign_mask) != 0
        || (denominator_bits & exponent_mask) == 0)
      return false;

    const core::real_t numerator_magnitude = magnitude(numerator);
    if (numerator_magnitude == 0.0 || positive_denominator >= 1.0)
      return true;
    return numerator_magnitude
           < std::numeric_limits<core::real_t>::max()
               * positive_denominator;
  }

  /**
   * Mutable harness buffers must not overlap any span read after publication.
   * `std::less<T*>` supplies the portable strict total order for unrelated
   * arrays which built-in relational pointer operators do not provide.
   */
  [[nodiscard]] inline bool spansOverlap(
    std::span<const core::real_t> left,
    std::span<const core::real_t> right) noexcept
  {
    if (left.empty() || right.empty())
      return false;
    const auto less = std::less<const core::real_t *>{};
    const auto *const left_end = left.data() + left.size();
    const auto *const right_end = right.data() + right.size();
    return less(left.data(), right_end) && less(right.data(), left_end);
  }

} // namespace detail

struct CurrentA
{
  explicit CurrentA(std::span<const core::real_t> lane_values) noexcept
    : lane_values{ lane_values }
  {}

  std::span<const core::real_t> lane_values;
};

struct CurrentDensityApm2
{
  explicit CurrentDensityApm2(
    std::span<const core::real_t> lane_values) noexcept
    : lane_values{ lane_values }
  {}

  std::span<const core::real_t> lane_values;
};

struct VoltageError
{
  core::real_t maximum_absolute_V;
  core::real_t rms_V;
};

[[nodiscard]] inline core::SpmBatch requireSpmBatch(
  const core::SpmFactoryInput &input,
  const core::SpmModelOptions &options,
  int n_lanes)
{
  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(input, options, n_lanes, batch)
          == Status::Success);
  return batch;
}

inline void requireTerminalVoltage(core::SpmBatch &batch,
                                   CurrentDensityApm2 current_density,
                                   const core::real_t &time_s,
                                   std::span<core::real_t>
                                     terminal_voltage_V)
{
  REQUIRE(batch.valid());
  REQUIRE(batch.n_lanes() > 0);
  const std::size_t lane_count = batch.n_lanes();
  REQUIRE(current_density.lane_values.size() == lane_count);
  REQUIRE(terminal_voltage_V.size() == lane_count);
  REQUIRE_FALSE(
    detail::spansOverlap(current_density.lane_values, terminal_voltage_V));
  REQUIRE(core::is_finite(time_s));
  for (const core::real_t &density_Apm2 : current_density.lane_values)
    REQUIRE(core::is_finite(density_Apm2));

  const core::StepCtx observation{
    .time = time_s,
    .dt = 0.0,
    .i_app = current_density.lane_values,
  };
  REQUIRE(batch.terminalVoltage(observation, terminal_voltage_V)
          == Status::Success);
  for (const core::real_t &voltage_V : terminal_voltage_V)
    REQUIRE(core::is_finite(voltage_V));
}

inline void requireTerminalVoltage(
  core::SpmBatch &batch,
  CurrentA current,
  const core::real_t &time_s,
  std::span<core::real_t>
    current_density_scratch_Apm2,
  std::span<core::real_t>
    terminal_voltage_V)
{
  REQUIRE(batch.valid());
  REQUIRE(batch.n_lanes() > 0);
  const std::size_t lane_count = batch.n_lanes();
  REQUIRE(current.lane_values.size() == lane_count);
  REQUIRE(current_density_scratch_Apm2.size() == lane_count);
  REQUIRE(terminal_voltage_V.size() == lane_count);
  REQUIRE_FALSE(detail::spansOverlap(current.lane_values,
                                     current_density_scratch_Apm2));
  REQUIRE_FALSE(
    detail::spansOverlap(current.lane_values, terminal_voltage_V));
  REQUIRE_FALSE(detail::spansOverlap(current_density_scratch_Apm2,
                                     terminal_voltage_V));
  REQUIRE(core::is_finite(time_s));
  const core::real_t area = batch.electrode_area();
  REQUIRE(core::is_finite(area));
  REQUIRE(area > 0.0);

  for (const core::real_t &current_A : current.lane_values) {
    REQUIRE(core::is_finite(current_A));
    REQUIRE(detail::canDivideByPositiveNormal(current_A, area));
  }
  for (std::size_t lane = 0; lane < lane_count; ++lane) {
    const core::real_t density_Apm2 =
      current.lane_values[lane] / area;
    REQUIRE(core::is_finite(density_Apm2));
    current_density_scratch_Apm2[lane] = density_Apm2;
  }

  requireTerminalVoltage(batch,
                         CurrentDensityApm2{ current_density_scratch_Apm2 },
                         time_s,
                         terminal_voltage_V);
}

inline void requireConstantCurrentTrace(
  core::SpmBatch &batch,
  core::ExponentialModal &stepper,
  CurrentDensityApm2 current_density,
  std::span<const core::real_t>
    sample_time_s,
  std::span<core::real_t>
    sample_major_terminal_voltage_V)
{
  REQUIRE(batch.valid());
  REQUIRE(batch.n_lanes() > 0);
  const std::size_t lane_count = batch.n_lanes();
  REQUIRE(current_density.lane_values.size() == lane_count);
  REQUIRE_FALSE(sample_time_s.empty());
  REQUIRE(sample_time_s.size()
          <= std::numeric_limits<std::size_t>::max() / lane_count);
  REQUIRE(sample_major_terminal_voltage_V.size()
          == sample_time_s.size() * lane_count);
  REQUIRE_FALSE(detail::spansOverlap(current_density.lane_values,
                                     sample_major_terminal_voltage_V));
  REQUIRE_FALSE(detail::spansOverlap(sample_time_s,
                                     sample_major_terminal_voltage_V));

  for (const core::real_t &density_Apm2 : current_density.lane_values)
    REQUIRE(core::is_finite(density_Apm2));
  REQUIRE(core::is_finite(sample_time_s.front()));
  for (std::size_t sample = 1; sample < sample_time_s.size(); ++sample)
    REQUIRE((core::is_finite(sample_time_s[sample])
             && sample_time_s[sample] > sample_time_s[sample - 1]));

  requireTerminalVoltage(
    batch,
    current_density,
    sample_time_s.front(),
    sample_major_terminal_voltage_V.first(lane_count));
  REQUIRE(stepper.configure(batch) == Status::Success);

  for (std::size_t sample = 1; sample < sample_time_s.size(); ++sample) {
    const core::real_t start_time_s = sample_time_s[sample - 1];
    const core::real_t dt_s = sample_time_s[sample] - start_time_s;
    REQUIRE(stepper.step(batch,
                         current_density.lane_values,
                         start_time_s,
                         dt_s)
            == Status::Success);
    const auto voltage_V = stepper.terminalVoltage();
    REQUIRE(voltage_V.size() == lane_count);
    std::ranges::copy(
      voltage_V,
      sample_major_terminal_voltage_V.subspan(sample * lane_count,
                                              lane_count)
        .begin());
  }
  REQUIRE(std::ranges::all_of(sample_major_terminal_voltage_V,
                              [](const core::real_t &voltage_V) {
                                return core::is_finite(voltage_V);
                              }));
}

[[nodiscard]] inline bool computeVoltageError(
  std::span<const core::real_t> actual_voltage_V,
  std::span<const core::real_t> expected_voltage_V,
  VoltageError &output) noexcept
{
  if (actual_voltage_V.empty()
      || actual_voltage_V.size() != expected_voltage_V.size())
    return false;

  if (!detail::canEvaluateAbsoluteDifference(actual_voltage_V[0],
                                             expected_voltage_V[0]))
    return false;
  core::real_t maximum_absolute_V =
    std::abs(actual_voltage_V[0] - expected_voltage_V[0]);
  if (!(core::is_finite(maximum_absolute_V)
        && detail::canSquare(maximum_absolute_V)))
    return false;
  core::real_t square_sum_V2 =
    maximum_absolute_V * maximum_absolute_V;
  if (!core::is_finite(square_sum_V2))
    return false;

  for (std::size_t sample = 1; sample < actual_voltage_V.size(); ++sample) {
    const core::real_t &actual_V = actual_voltage_V[sample];
    const core::real_t &expected_V = expected_voltage_V[sample];
    if (!detail::canEvaluateAbsoluteDifference(actual_V, expected_V))
      return false;

    const core::real_t error_V = std::abs(actual_V - expected_V);
    if (!core::is_finite(error_V))
      return false;
    maximum_absolute_V = std::max(maximum_absolute_V, error_V);
    if (!(detail::canSquare(error_V)
          && detail::canAddSquare(square_sum_V2, error_V)))
      return false;
    square_sum_V2 += error_V * error_V;
    if (!core::is_finite(square_sum_V2))
      return false;
  }

  const core::real_t rms_V =
    std::sqrt(square_sum_V2
              / static_cast<core::real_t>(actual_voltage_V.size()));
  if (!core::is_finite(rms_V))
    return false;
  output = { maximum_absolute_V, rms_V };
  return true;
}

} // namespace slide::test_support
