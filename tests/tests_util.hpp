/*
 * tests_util.hpp
 *
 * Utility functions for tests. (Not really replacing GTest but enough functionality)
 *
 *  Created on: 12 Nov 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <cmath>
#include <cassert>
#include <iostream>

namespace slide::tests {

constexpr bool DEBUG_TESTS = false;

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

bool EQ(double x1, double x2)
{
  return NEAR(x1, x2);
}


} // namespace slide::tests
