/*
 * ArrayVec.hpp
 *
 *  A small class to store variable size arrays
 *  Created on: 05 Apr 2022
 *  //!< not working yet.
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <vector>
#include <array>
#include <cstdlib>
#include <algorithm>
#include <span>

namespace slide {
template <typename Tdata>
class ArrayVec
{
  std::vector<std::span<Tdata>> data_span;

public:
  std::vector<Tdata> data;
  ArrayVec() = default;

  void push_back(std::span<const Tdata> spn)
  {
    const auto new_begin = data.end();
    data.insert(new_begin, spn.begin(), spn.end());
    data_span.emplace_back(new_begin, data.end());
  }

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