/**
 * @file DataFrame.hpp
 * @brief Basic data frame class to keep named data.
 * @author Volkan Kumtepeli
 * @date 02 Sep 2023
 */

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <string_view>
#include <iomanip> //!< Required for std::setprecision

namespace slide {
template <typename T = double>
struct DataFrame : public std::vector<T>
{
  std::streamsize precision{ 12 };
  std::vector<std::string> header{};

  DataFrame() = default;

  // Templated constructor for various types
  template <typename NameType>
  explicit DataFrame(std::initializer_list<NameType> nameList)
    : header(std::begin(nameList), std::end(nameList))
  {
  }

  // Function to display the header for demonstration
  void displayNames() const
  {
    std::cout << "Names: ";
    for (const auto &name : header) {
      std::cout << name << ' ';
    }
    std::cout << std::endl;
  }

  // Function to write DataFrame data to CSV
  void to_csv(std::ostream &os, bool includeHeader = true) const
  {
    const auto colSize = header.size(); //<! Column size
    if (includeHeader) {
      // Write the header
      for (size_t i = 0; i < colSize; ++i) {
        if (i != 0) os << ',';
        os << header[i];
      }
      os << '\n';
    }

    os << std::fixed << std::setprecision(precision);

    // Write the data
    for (size_t i{}; i < this->size(); ++i) {
      if (i != 0) {
        if (i % colSize == 0)
          os << '\n';
        else
          os << ',';
      }
      os << (*this)[i];
    }
  }

  void to_csv(std::string_view str, bool includeHeader = true)
  {
    std::ofstream out{ str, std::ios::out };
    to_csv(out, header);
  }

  void to_binary(std::filesystem::path pth, bool includeHeader = true)
  {
    io::binary_writer(pth, *this, header, std::ios::out);
  }
};

} // namespace slide
