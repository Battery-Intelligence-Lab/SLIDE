/**
 * @file CheckedArithmetic.hpp
 * @brief Overflow-checked integer arithmetic for the recording layouts.
 *
 * Owns: `checkedAdd`, `checkedMultiply`, `align64`. Implements PLAN.md §3.7 -- every recording
 * extent is computed from untrusted sizes, so each product and sum is checked before it is used
 * to allocate or to seek. Cold: layout computation only.
 *
 * These lived twice, once in `Recorder.cpp` and once in `AsyncRecorder.cpp`, with identical
 * semantics and different parameter names. One concept, one definition (MC-1).
 * @surface internal
 */

#pragma once

#include <cstdint>
#include <limits>

namespace slide::core::detail {

template <class T>
bool checkedAdd(T left, T right, T &result)
{
  if (left > std::numeric_limits<T>::max() - right)
    return false;
  result = left + right;
  return true;
}

template <class T>
bool checkedMultiply(T left, T right, T &result)
{
  if (left != 0 && right > std::numeric_limits<T>::max() / left)
    return false;
  result = left * right;
  return true;
}

/** Round up to the next 64-byte boundary, refusing to wrap. */
inline bool align64(std::uint64_t value, std::uint64_t &result)
{
  std::uint64_t enlarged{};
  if (!checkedAdd(value, std::uint64_t{ 63 }, enlarged))
    return false;
  result = enlarged & ~std::uint64_t{ 63 };
  return true;
}

} // namespace slide::core::detail
