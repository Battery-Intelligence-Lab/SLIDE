
/*
 * PtrInterval.hpp
 *
 *  Created on: 27 May 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <vector>
#include <array>
#include <span>
#include <iterator>

//#include <type_traits>

namespace slide {
template <typename T>
class PtrInterval //!< std::span is same size. It is now postponed.
{
  T *begin_ptr{ nullptr };
  T *end_ptr{ nullptr };

public:
  [[nodiscard]] constexpr T *begin() noexcept { return begin_ptr; }
  [[nodiscard]] constexpr T *end() noexcept { return end_ptr; }

  size_t size() { return std::distance(begin_ptr, end_ptr); }

  PtrInterval(T &x)
  {
    begin_ptr = &x;
    end_ptr = begin_ptr + 1;
  }

  PtrInterval(std::vector<T> &x)
  {
    begin_ptr = &x[0];
    end_ptr = &x.back() + 1;
  }

  template <size_t N>
  PtrInterval(std::array<T, N> &x)
  {
    begin_ptr = &x[0];
    end_ptr = &x.back() + 1;
  }

  PtrInterval(std::span<T> x)
  {
    begin_ptr = &x[0];
    end_ptr = &x.back() + 1;
  }
};
} // namespace slide