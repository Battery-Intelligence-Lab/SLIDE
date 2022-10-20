/*
 * timing.hpp
 *
 * Some utility functions for timing.

 *  Created on: 16 Oct 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */


#pragma once

#include <ctime>
#include <iostream>
#include <cmath>

namespace slide {
struct Clock
{
  std::clock_t tstart{ std::clock() };
  Clock() = default;
  auto now() const { return std::clock(); }
  auto start() const { return tstart; }
  double duration() const { return (now() - start()) / static_cast<double>(CLOCKS_PER_SEC); }
};

inline std::ostream &operator<<(std::ostream &ofs, const Clock clk)
{
  const auto duration = clk.duration();
  ofs << std::floor(duration / 60) << ":"
      << duration - std::floor(duration / 60) * 60
      << " min:sec";

  return ofs;
}
} // namespace slide