/**
 * @file descriptions.hpp
 * @brief get description of variables
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 05 Apr 2022
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