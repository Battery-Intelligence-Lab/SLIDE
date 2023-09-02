/*
 * DataFrame.hpp
 *
 * Basic data frame class to keep named data. 
 *
 *  Created on: 02 Sep 2023
 *   Author(s): Volkan Kumtepeli
 */

#pragma once

#include <vector> 
#include <string>

namespace slide
{
template <typename T = double>
struct DataFrame : public std::vector<T>
{
    std::vector<std::string> names{};
    DataFrame() = default;

    // Templated constructor for various types
    template <typename NameType>
    explicit DataFrame(std::initializer_list<NameType> nameList)
      : names(std::begin(nameList), std::end(nameList))
    {
    }

    // Function to display the names for demonstration
    void displayNames() const
    {
      std::cout << "Names: ";
      for (const auto &name : names) {
        std::cout << name << ' ';
      }
      std::cout << std::endl;
    }
};

}

