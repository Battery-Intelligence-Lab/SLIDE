/**
 * @file Pair.hpp
 * @brief Just iterable pair.
 * @author Volkan Kumtepeli
 * @date 11 Aug 2024
 */

#pragma once

#include <array>

namespace slide {

template <typename T = double>
using Pair = std::array<T, 2>;

using DPair = Pair<double>;
using IPair = Pair<int>;

} // namespace slide