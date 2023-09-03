/*
 * cell_data.hpp
 *
 *  Created on: 07 Feb 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../State.hpp"
#include "../Histogram.hpp"

#include <iostream>

namespace slide {

struct CellCommonHist
{
  Histogram<> I, V, T; //!< histograms for current, voltage, temperature
};

using ThroughputData = State<0, 3>;

} // namespace slide