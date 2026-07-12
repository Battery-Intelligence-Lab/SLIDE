/**
 * @file BoundedFileReader.hpp
 * @brief Exact, allocation-bounded binary file reads shared by cold parsers.
 */

#pragma once

#include "../types/Status.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <istream>
#include <string>
#include <utility>

namespace slide::core::detail {

enum class BoundedFileRead : unsigned char {
  success,
  open_failed,
  too_large,
  io_failed,
};

/** Map the closed internal outcome set to the public parser contract. */
inline slide::Status boundedFileStatus(BoundedFileRead result) noexcept
{
  constexpr std::array statuses{
    slide::Status::Success,
    slide::Status::Invalid_parameters,
    slide::Status::Invalid_parameters,
    slide::Status::Numerical_failure,
  };
  const auto index = static_cast<std::size_t>(result);
  assert(index < statuses.size());
  return statuses[index];
}

/**
 * Read until EOF without trusting a racy size probe. At most `limit` bytes are
 * retained, and `output` is published only after a clean EOF.
 */
inline BoundedFileRead readBoundedStream(std::istream &input,
                                         std::size_t limit,
                                         std::string &output)
{
  constexpr std::size_t chunk_bytes = 8192;
  std::array<char, chunk_bytes> buffer{};
  std::string candidate;
  candidate.reserve(std::min(limit, chunk_bytes));

  for (;;) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const auto count = input.gcount();
    if (input.bad())
      return BoundedFileRead::io_failed;
    if (count < 0)
      return BoundedFileRead::io_failed;
    const auto bytes = static_cast<std::size_t>(count);
    if (bytes > limit - candidate.size())
      return BoundedFileRead::too_large;
    candidate.append(buffer.data(), bytes);
    if (input.eof()) {
      output = std::move(candidate);
      return BoundedFileRead::success;
    }
    if (input.fail() || bytes == 0)
      return BoundedFileRead::io_failed;
  }
}

inline BoundedFileRead readBoundedFile(const std::filesystem::path &path,
                                       std::size_t limit,
                                       std::string &output)
{
  std::ifstream input(path, std::ios::binary);
  if (!input)
    return BoundedFileRead::open_failed;
  return readBoundedStream(input, limit, output);
}

} // namespace slide::core::detail
