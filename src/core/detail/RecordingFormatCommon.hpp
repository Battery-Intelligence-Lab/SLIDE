/**
 * @file RecordingFormatCommon.hpp
 * @brief Shared checksum, byte-order, and allocation-failure facts for recording.
 *
 * Owns: `endian_marker`, `crc32`, `headerCrc`, and `allocationFailureStatus`.
 * Implements PLAN.md section 3.7. Cold: file-format checks and allocation translation only.
 * @surface internal
 */

#pragma once

#include "../../types/Status.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>

namespace slide::core::detail {

inline constexpr std::uint32_t endian_marker = 0x01020304U;

inline std::uint32_t crc32(std::span<const std::byte> bytes)
{
  std::uint32_t crc = 0xffffffffU;
  for (const auto byte : bytes) {
    crc ^= std::to_integer<std::uint8_t>(byte);
    for (int bit = 0; bit < 8; ++bit)
      crc = (crc >> 1U) ^ (0xedb88320U & (0U - (crc & 1U)));
  }
  return ~crc;
}

template <class Header>
[[nodiscard]] inline std::uint32_t headerCrc(Header header)
{
  static_assert(std::is_trivially_copyable_v<Header>);
  header.header_crc32 = 0;
  return crc32(std::as_bytes(std::span{ &header, 1 }));
}

inline slide::Status allocationFailureStatus() noexcept
{
  return slide::Status::Numerical_failure;
}

} // namespace slide::core::detail
