/*
 * drive_cycles.hpp
 *
 *  Example estimation functions;
 *
 *  Created on: 09 Nov 2022
 *   Author(s): Volkan Kumtepeli
 */

#pragma once

#include "../src/slide.hpp"

#include <string>
#include <memory>

namespace slide::examples {

inline void estimatingOCVparameters()
{
  slide::estimateOCVparameters();
}
} // namespace slide::examples