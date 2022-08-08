/*
 * util_debug.hpp
 *
 *  Created on: 24 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <iostream>

namespace slide::debug
{
    inline void getLine() // Not working.
    {
        std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    }
}