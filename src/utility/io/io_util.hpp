/**
 * @file io_util.hpp
 * @brief Utility functions for IO.
 * @author Volkan Kumtepeli
 * @date 04 Aug 2024
 */

#pragma once
#include <iostream>
#include <fstream>

namespace slide::io {

inline void ignoreBOM(std::ifstream &in)
{
  char c = '.';

  do {
    in >> c;
  } while (c < 45); //!< Ignore byte order mark (BOM) at the beginning

  in.putback(c);
}


} // namespace slide::io