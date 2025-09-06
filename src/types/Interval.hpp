/**
 * @file Interval.hpp
 * @brief A small class for returning non-owning view of some storages.
 * It stores a pointer and an interval.
 * @author Volkan Kumtepeli
 * @date 05 Apr 2022
 */

#pragma once

#include <vector>

namespace slide {

/**
 * @brief A non-owning view class for accessing intervals of container data
 *
 * @tparam Tcontainer Type of container to provide interval access to
 *
 * This class provides a non-owning view of a specific interval within
 * a container, allowing efficient access to sub-ranges without copying data.
 */
template <typename Tcontainer>
class Interval
{
  Tcontainer *root{ nullptr }; //!< Pointer to the underlying container
  int beg{ 0 }, en{ 0 };       //!< Begin and end indices of the interval

public:
  /**
   * @brief Default constructor
   */
  Interval() = default;

  /**
   * @brief Constructor for full container interval
   * @param data Reference to the container
   */
  Interval(Tcontainer &data) : root(&data), en(data.size() - 1) {}

  /**
   * @brief Constructor for specific interval
   * @param data Reference to the container
   * @param begin_ Begin index of the interval
   * @param end_ End index of the interval
   */
  Interval(Tcontainer &data, int begin_, int end_) : root(&data), beg(begin_), en(end_) {}

  [[nodiscard]] constexpr auto begin() noexcept { return std::begin(*root) + beg; }
  [[nodiscard]] constexpr auto end() noexcept { return std::begin(*root) + en; }

  [[nodiscard]] constexpr auto cbegin() const noexcept { return std::advance(std::cbegin(*root), beg); }
  [[nodiscard]] constexpr auto cend() const noexcept { return std::advance(std::cbegin(*root), en); }
};
} // namespace slide