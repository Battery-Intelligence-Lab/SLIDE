/*
 * StateStorage.hpp
 *
 *  Created on: 24 May 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <vector>
#include <span>
#include <type_traits>
#include <algorithm>

namespace slide {
class StateStorage
{
  std::vector<std::span<double>> state_span;
  std::vector<double> states;

public:
  size_t index{ 0 };
  constexpr void clear() noexcept
  {
    state_span.clear();
    states.clear();
    index = 0;
  }

  template <typename T>
  void add(T &x)
  {
    add({ &x, 1 });
  }

  template <typename T>
  void add(std::span<T> x)
  {
    state_span.emplace_back(x);
    states.insert(states.end(), x.begin(), x.end()); //!< states.reserve(states.size() + x.size());
  }

  [[nodiscard]] constexpr auto begin() noexcept { return states.begin(); }
  [[nodiscard]] constexpr auto end() noexcept { return states.end(); }

  void restore() //!< Restore states to spans.
  {
    size_t i{ 0 };

    for (auto spn : state_span) {
      std::copy_n(states.begin() + i, spn.size(), spn.begin());
      i += spn.size();
    }
  }
};
} // namespace slide