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

namespace slide::tests {

bool TEST(auto &&fun, auto &&fun_name)
{
  try {
    return fun();
  } catch (...) {
    return false;
  }
  return true;
}

bool NEAR(auto x1, auto x2, double abs_error = 1e-15)
{
  return std::abs(x1 - x2) < abs_error;
}

bool EQ(double x1, double x2)
{
  return NEAR(x1, x2);
}


} // namespace slide::tests
