/*
 * get description of variables.cpp
 *
 * Free functions to help to write everything shorter.
 *
 * Created on: 05 Apr 2022
 *  Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <string_view>

template <typename T>
consteval int find_description(std::string_view x)
{
  for (size_t i = 0; i < T::description.size(); i++) {
    if (x == T::description[i])
      return i;
  }

  return -1;
}