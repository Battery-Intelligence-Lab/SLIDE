/*
 * predicate_functions.cpp
 *
 * Predicate functions to help to write everything shorter.
 *
 *  Created on: 07 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

namespace slide::util {
template <typename T>
constexpr auto is_zero(T x) //!< No reference it is probably already very small object.
{
  return x == 0;
}

template <typename T>
constexpr int sign(T x) //!< No reference it is probably already very small object.
{
  return (x > T(0)) - (x < T(0));
}

} // namespace slide::util