/**
 * @file VecState.hpp
 * @brief A small class for allocating states.
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 05 Apr 2022
 */

#pragma once

#include <vector>
#include <array>
#include <cstdlib>
#include <algorithm>

namespace slide {
class VecStates
{
  static std::vector<double> st, dst, bst; //!< states, dstates, backup_states
  static std::vector<std::array<int, 2>> locs{ { 0, 0 } };
  static size_t current_id{ 0 };

  size_t id{ 0 };

public:
  VecStates(int n) : id{ current_id++ }
  {
    if (locs.empty()) {
      locs.emplace_back(0, n);
      std::fill_n(std::back_inserter(st), n, 0);
    } else {
    }
  };

  VecStates() : VecStates(0) {};

  //!< auto begin() { return std::begin(*root) + beg; }
  //!< auto end() { return std::begin(*root) + en; }

  //!< auto cbegin() { return std::advance(std::cbegin(*root), beg); }
  //!< auto cend() { return std::advance(std::cbegin(*root), en); }
};
} // namespace slide