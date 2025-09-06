/**
 * @file descriptions.hpp
 * @brief get description of variables
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 05 Apr 2022
 */

#pragma once

#include <string_view>

/**
 * @brief Find the index of a description string in a type's description array
 *
 * @tparam T Type that has a static description array
 * @param x String view to search for
 * @return int Index of the matching description, or -1 if not found
 *
 * This function searches through the static description array of type T
 * to find the index of the provided string view. It's used for finding
 * variable descriptions at compile time.
 */
template <typename T>
consteval int find_description(std::string_view x)
{
  for (size_t i = 0; i < T::description.size(); i++) {
    if (x == T::description[i])
      return i;
  }

  return -1;
}