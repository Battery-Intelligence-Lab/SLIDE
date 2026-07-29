/**
 * @file SnapshotIndexing.hpp
 * @brief Shared SoA-to-SnapshotView indexing for eager recording readers.
 *
 * Owns: `snapshotView`. Implements PLAN.md section 3.7.
 * Cold: allocation-free eager-reader indexing only.
 * @surface internal
 */

#pragma once

#include "../Recorder.hpp"

#include <cstddef>
#include <cstdint>
#include <span>

namespace slide::core::detail {

[[nodiscard]] inline SnapshotView snapshotView(
  std::size_t index,
  std::size_t lanes,
  std::size_t state_values,
  std::span<const std::uint64_t> steps,
  std::span<const real_t> times,
  std::span<const real_t> currents,
  std::span<const real_t> states) noexcept
{
  return { .accepted_step = steps[index],
           .time = times[index],
           .current_density = currents.subspan(index * lanes, lanes),
           .state = states.subspan(
             index * state_values, state_values) };
}

} // namespace slide::core::detail
