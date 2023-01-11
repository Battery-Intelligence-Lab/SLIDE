/*
 * CellDataStoragePolicy.hpp
 *
 * This class is created to generic interface to store data.
 *  Created on: 08 Sep 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */


#pragma once

#include "../../settings/enum_definitions.hpp"

#include <type_traits>

namespace slide::policy {

template <settings::cellDataStorageLevel N>


using CellDataStoragePolicy = std::conditional_t<
  settings::cellDataStorageLevel == 'I',
  First, std::conditional_t<T == 'D', Second, Third>>;


}