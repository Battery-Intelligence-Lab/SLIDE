/**
 * @file estimation.hpp
 * @brief Example estimation functions.
 * @author Volkan Kumtepeli
 * @date 09 Nov 2022
 */

#pragma once

#include "../src/slide.hpp"

#include <string>
#include <memory>

namespace slide::examples {

inline void estimatingOCVparameters()
{
  slide::estimateOCVparameters(); // Solution Kokam 2.71 [sp, sn, AMp, AMn] = [0.3853; 0.5645; 3.5997e-06; 6.1982e-06];
}
} // namespace slide::examples