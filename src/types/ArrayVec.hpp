/**
 * @file ArrayVec.hpp
 * @brief A small class to store variable size arrays. Not working yet.
 * @author Volkan Kumtepeli
 * @date 05 Apr 2022
 */

#pragma once

#include <vector>
#include <array>
#include <cstdlib>
#include <algorithm>
#include <span>

namespace slide {

/**
 * @brief A container class for storing variable-size arrays
 *
 * @tparam Tdata Type of data to store
 *
 * This class provides a way to store multiple variable-size arrays
 * in a single contiguous memory block while maintaining access to
 * individual arrays through spans.
 *
 * @note This class is currently not fully implemented and may not work correctly.
 */
template <typename Tdata>
class ArrayVec
{
  std::vector<std::span<Tdata>> data_span; //!< Vector of spans pointing to data segments

public:
  std::vector<Tdata> data; //!< Contiguous storage for all data

  /**
   * @brief Default constructor
   */
  ArrayVec() = default;

  /**
   * @brief Add a span of data to the container
   * @param spn Span of data to add
   */
  void push_back(std::span<const Tdata> spn)
  {
    const auto new_begin = data.end();
    data.insert(new_begin, spn.begin(), spn.end());
    data_span.emplace_back(new_begin, data.end());
  }

  /**
   * @brief Add a single element to the container
   * @param x Element to add
   */
  void push_back(const Tdata &x)
  {
    const auto new_begin = data.end();
    data.push_back(x);
    data_span.emplace_back(new_begin, new_begin + 1);
  }

  [[nodiscard]] constexpr auto begin() noexcept { return data_span.begin(); }
  [[nodiscard]] constexpr auto end() noexcept { return data_span.end(); }

  [[nodiscard]] constexpr auto cbegin() noexcept { return data_span.cbegin(); }
  [[nodiscard]] constexpr auto cend() noexcept { return data_span.cend(); }
};
} // namespace slide