/*
 * DataStoragePolicy.hpp
 *
 * This class is created to generic interface to store data.
 *  Created on: 28 Aug 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */


#pragma once

#include <type_traits>

namespace slide::policy {

template <typename T>
inline void no_data(T &&) {} // Dont do anything.


} // namespace slide::policy