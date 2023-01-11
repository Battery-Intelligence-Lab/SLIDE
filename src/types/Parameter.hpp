/*
 * Parameter.hpp
 *
 * A small class for named parameters.
 *  Created on: 04 Mar 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <string_view>

struct Parameter
{
  //!< Name should not be copied between all classes. It may be a template argument maybe?
  //!< It should return its value using () argument.
  //!< It we should be able to make it a function.
  //!< Maybe we can use unions between function pointers and doubles.

  double value{};

  constexpr static std::string_view name{}; //!< Name should be shared by all classes.

  auto &operator()() //!< An operator to get its value.
  {
    return value;
  }
};
