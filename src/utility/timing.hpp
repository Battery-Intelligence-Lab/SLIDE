/**
 * @file timing.hpp
 * @brief Some utility functions for timing.
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 16 Oct 2022
 */

#pragma once

#include <ctime>
#include <iostream>
#include <cmath>
#include <chrono>

namespace slide {

/**
 * @brief A simple clock class for measuring elapsed time
 *
 * This class provides functionality to measure elapsed time since construction
 * or a specific start time. It uses std::chrono::steady_clock for reliable
 * timing measurements.
 */
struct Clock
{
  std::chrono::time_point<std::chrono::steady_clock> tstart{ std::chrono::steady_clock::now() };

  /**
   * @brief Default constructor that initializes the start time
   */
  Clock() = default;

  /**
   * @brief Get the current time point
   * @return Current time point
   */
  auto now() const { return std::chrono::steady_clock::now(); }

  /**
   * @brief Get the start time point
   * @return Start time point
   */
  auto start() const { return tstart; }

  /**
   * @brief Calculate the elapsed time since start
   * @return Elapsed time in seconds as a double
   */
  double duration() const
  {
    std::chrono::duration<double> elapsed_seconds = now() - start();
    return elapsed_seconds.count();
  }
};

/**
 * @brief Stream output operator for Clock class
 *
 * @param ofs Output stream
 * @param clk Clock object to output
 * @return Reference to the output stream
 *
 * Outputs the elapsed time in "minutes:seconds" format.
 */
inline std::ostream &operator<<(std::ostream &ofs, const Clock &clk)
{
  const auto duration = clk.duration();
  ofs << std::floor(duration / 60) << ":"
      << duration - std::floor(duration / 60) * 60
      << " min:sec";

  return ofs;
}
} // namespace slide