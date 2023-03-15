/**
 * @file tests_util.hpp
 * @brief Utility functions for tests. (Not really replacing GTest but enough functionality)
 * @author Jorn Reniers, Volkan Kumtepeli
 * @date 12 Nov 2022
 */

#pragma once

#include <cmath>
#include <cassert>
#include <iostream>

namespace slide::tests {

constexpr bool DEBUG_TESTS = false;

/**
 * @brief Test function wrapper.
 * @param fun Function to be tested.
 * @param fun_name Name of the function to be tested.
 * @return true if the test passes, false if it fails or an exception is caught.
 */
bool TEST(auto &&fun, auto &&fun_name)
{

  if (DEBUG_TESTS)
    std::cout << fun_name << " test is started!" << std::endl;

  try {
    return fun();
  } catch (...) {
    return false;
  }
  return true;
}

/**
 * @brief Check if two values are nearly equal within a given absolute error.
 * @param x1 First value to compare.
 * @param x2 Second value to compare.
 * @param abs_error The allowed absolute error between the two values (default: 1e-15).
 * @return true if the values are nearly equal, false otherwise.
 */
bool NEAR(auto x1, auto x2, double abs_error = 1e-15)
{
  const auto abs_diff = std::abs(x1 - x2);
  if (abs_diff < abs_error)
    return true;
  else {
    std::cerr << "Abs diff: " << abs_diff << " x1: " << x1
              << " x2: " << x2 << '\n';
    return false;
  }
}

/**
 * @brief Check if two doubles are equal (wrapper for NEAR function).
 * @param x1 First double value to compare.
 * @param x2 Second double value to compare.
 * @return true if the values are equal, false otherwise.
 */
bool EQ(double x1, double x2)
{
  return NEAR(x1, x2);
}


} // namespace slide::tests
