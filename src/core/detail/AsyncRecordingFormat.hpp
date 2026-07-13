/**
 * @file AsyncRecordingFormat.hpp
 * @brief The compressed-recording file format: magic, headers, CRC, and the compression bound.
 *
 * Owns: the on-disk constants and header structs shared by the writer (`AsyncRecorder`) and the
 * reader (`CompressedRecording`), plus their CRCs and the codec's worst-case bound.
 * Implements PLAN.md §3.7. Cold: written once per block, read once per open.
 *
 * The writer and the reader must agree byte for byte, so the format has exactly one definition
 * rather than one copy per translation unit (MC-1, PC-10 in spirit).
 * @surface internal
 */

#pragma once

#include "../AsyncRecorder.hpp"

#include <array>
#include <cstdint>
#include <span>

#if defined(SLIDE_WITH_ZSTD)
#include <zstd.h>
#endif

namespace slide::core::detail {

constexpr std::array<char, 8> compressed_magic{ 'S', 'L', 'I', 'D', 'E', 'C', 'M', 'P' };
constexpr std::uint32_t block_magic = 0x314b4c42U; // BLK1
constexpr std::uint32_t endian_marker = 0x01020304U;
constexpr std::uint16_t format_major = 1;
constexpr std::uint16_t format_minor = 0;

#pragma pack(push, 1)
struct CompressedFileHeader
{
  std::array<char, 8> magic{};
  std::uint16_t major{};
  std::uint16_t minor{};
  std::uint32_t endian{};
  std::uint32_t header_size{};
  std::uint32_t header_crc32{};
  std::uint32_t rows{};
  std::uint32_t lanes{};
  std::uint32_t stride{};
  std::uint32_t codec{};
  std::uint64_t snapshots{};
  std::uint64_t file_size{};
  std::uint64_t reserved{};
};

struct CompressedBlockHeader
{
  std::uint32_t magic{};
  std::uint32_t header_size{};
  std::uint32_t codec{};
  std::uint32_t flags{};
  std::uint64_t accepted_step{};
  real_t time{};
  std::uint64_t raw_bytes{};
  std::uint64_t payload_bytes{};
  std::uint32_t raw_crc32{};
  std::uint32_t payload_crc32{};
  std::uint32_t header_crc32{};
  std::uint32_t reserved{};
};
#pragma pack(pop)

static_assert(sizeof(CompressedFileHeader) == 64);
static_assert(sizeof(CompressedBlockHeader) == 64);

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

inline std::uint32_t fileHeaderCrc(CompressedFileHeader header)
{
  header.header_crc32 = 0;
  return crc32(std::as_bytes(std::span{ &header, 1 }));
}

inline std::uint32_t blockHeaderCrc(CompressedBlockHeader header)
{
  header.header_crc32 = 0;
  return crc32(std::as_bytes(std::span{ &header, 1 }));
}



inline slide::Status allocationFailureStatus() noexcept
{
  return slide::Status::Numerical_failure;
}

inline std::size_t compressionBound(CompressionCodec codec, std::size_t raw_bytes)
{
  if (codec == CompressionCodec::none)
    return raw_bytes;
#if defined(SLIDE_WITH_ZSTD)
  if (codec == CompressionCodec::zstd)
    return ZSTD_compressBound(raw_bytes);
#else
  (void)raw_bytes;
#endif
  return 0;
}


} // namespace slide::core::detail
