/*
 * slide_algorithms.cpp
 *
 * Some algorithms to make easier frequently used operations:
 *
 *  Created on: 24 Oct 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <numeric>
#include <functional>

namespace slide {

auto transform_sum(const auto &SUs, auto &function)
{
  return std::transform_reduce(std::cbegin(SUs), std::cend(SUs), 0.0, std::plus<>(), function);
}

auto transform_max(const auto &SUs, auto &function)
{
  return std::transform_reduce(std::cbegin(SUs), std::cend(SUs), 0.0, std::max<>(), function);
}

auto transform_min(const auto &SUs, auto &function)
{
  return std::transform_reduce(std::cbegin(SUs), std::cend(SUs), 0.0, std::min<>(), function);
}

auto transform_mean(const auto &SUs, auto &function)
{
  if (SUs.empty())
    return 0.0;

  return transform_sum(SUs, function) / SUs.size();
}

} // namespace slide