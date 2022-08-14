/*
 * Interval.hpp
 *
 * A small class for returning non-owning view of some storages.
 * It stores a pointer and an interval.
 *
 *  Created on: 05 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <vector>

namespace slide {

template <typename Tcontainer>
class Interval
{
  Tcontainer *root{ nullptr };
  int beg{ 0 }, en{ 0 };

public:
  Interval() = default;

  Interval(Tcontainer &data) : root(&data), en(data.size() - 1) {}

  Interval(Tcontainer &data, int begin_, int end_) : root(&data), beg(begin_), en(end_) {}

  [[nodiscard]] constexpr auto begin() noexcept { return std::begin(*root) + beg; }
  [[nodiscard]] constexpr auto end() noexcept { return std::begin(*root) + en; }

  [[nodiscard]] constexpr auto cbegin() const noexcept { return std::advance(std::cbegin(*root), beg); }
  [[nodiscard]] constexpr auto cend() const noexcept { return std::advance(std::cbegin(*root), en); }
};
} // namespace slide