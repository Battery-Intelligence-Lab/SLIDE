/**
 * @file RecordedBits.hpp
 * @brief Stable whole-trace fingerprints for pre-refactor floating-point fixtures.
 *
 * Every double contributes all 64 representation bits in a specified byte order.  Two
 * independent 64-bit recurrences and the exact value count make these fixtures compact
 * enough to review while still detecting any changed output bit in a recorded trace.
 */

#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

namespace slide::test_support {

struct RecordedBits
{
  std::uint64_t fnv1a{ UINT64_C(14695981039346656037) };
  std::uint64_t mixed{ UINT64_C(0x6a09e667f3bcc909) };
  std::size_t values{};

  void append(std::span<const double> trace) noexcept
  {
    // Framing makes [a][b,c] distinct from [a,b][c].
    absorb(static_cast<std::uint64_t>(trace.size()));
    for (const double value : trace) {
      absorb(std::bit_cast<std::uint64_t>(value));
      ++values;
    }
  }

private:
  void absorb(std::uint64_t word) noexcept
  {
    constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);
    for (unsigned shift = 0; shift < 64; shift += 8) {
      fnv1a ^= (word >> shift) & UINT64_C(0xff);
      fnv1a *= fnv_prime;
    }

    // A second, structurally different full-width recurrence.  Unsigned overflow is
    // intentional and defined; rotations make word position significant.
    mixed ^= word + UINT64_C(0x9e3779b97f4a7c15) + std::rotl(mixed, 17);
    mixed *= UINT64_C(0xbf58476d1ce4e5b9);
    mixed ^= mixed >> 29;
  }
};

} // namespace slide::test_support
