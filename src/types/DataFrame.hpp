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
#include <iostream>
#include <fstream>

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

    // Function to write DataFrame data to CSV
    void to_csv(std::ostream &os, bool header = true) const
    {
      const auto colSize = names.size(); //<! Column size
      if (header) {
        // Write the header
        for (size_t i = 0; i < colSize; ++i) {
          if (i != 0) os << ',';
          os << names[i];
        }
        os << '\n';
      }

      // Write the data
      for(size_t i{}; i < this->size(); ++i)
      {
        if (i != 0 && i % colSize == 0) os << '\n';
        os << (*this)[i];
      }
    }
};

}

