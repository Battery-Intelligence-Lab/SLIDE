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
  using out_type = decltype(function(SUs[0]));
  return std::transform_reduce(std::cbegin(SUs), std::cend(SUs), out_type(0), std::plus<>(), function);
}

auto transform_max(const auto &SUs, auto &function)
{

  if (SUs.empty())
    return 0.0;

  auto max = function(SUs[0]);

  for (auto iter = std::cbegin(SUs) + 1; iter != std::cend(SUs); ++iter)
    max = std::max(max, function(*iter));

  return max;
}

auto transform_min(const auto &SUs, auto &function)
{
  if (SUs.empty())
    return 0.0;

  auto min = function(SUs[0]);

  for (auto iter = std::cbegin(SUs) + 1; iter != std::cend(SUs); ++iter)
    min = std::min(min, function(*iter));

  return min;
}

auto transform_mean(const auto &SUs, auto &function)
{
  if (SUs.empty())
    return 0.0;

  return transform_sum(SUs, function) / SUs.size();
}

} // namespace slide