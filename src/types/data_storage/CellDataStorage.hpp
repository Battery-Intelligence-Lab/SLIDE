/*
 * CellDataStorage.hpp
 *
 *  Created on: 10 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../Histogram.hpp"
#include "cell_data.hpp"
#include "../../settings/enum_definitions.hpp"

#include <utility>
#include <iostream>
#include <array>
#include <vector>

namespace slide {

template <settings::cellDataStorageLevel N>
struct CellDataStorage
{
  template <typename Cell_t>
  inline void initialise(Cell_t &) {} //!< Do nothing.

  template <typename Cell_t>
  inline void storeData(Cell_t &) {} //!< Do nothing.
};

template <>
struct CellDataStorage<settings::cellDataStorageLevel::storeHistogramData> //!< Store as histogram.
{
  std::array<Histogram<>, 3> data;

  template <typename Cell_t>
  inline void initialise(Cell_t &cell) //!< Initialise the histograms.
  {
    data[0] = Histogram<>(-cell.Cap(), cell.Cap()); //!< 1C charge/discharge
    data[1] = Histogram<>(cell.Vmin(), cell.Vmax());
    data[2] = Histogram<>(cell.Tmin(), cell.Tmax());
  }

  template <typename Cell_t>
  inline void storeData(Cell_t &cell)
  {
    data[0].add(cell.I());
    data[1].add(cell.V());
    data[2].add(cell.T());
  }
};

template <>
struct CellDataStorage<settings::cellDataStorageLevel::storeTimeData>
{
  std::vector<double> data; //!< Common data

  template <typename Cell_t>
  inline void initialise(Cell_t &) {} //!< Do nothing.

  template <typename Cell_t>
  inline void storeData(Cell_t &cell)
  {
    const auto throughputs = cell.getThroughputs();
    // #TODO just write all states, throughputs will be included.
    data.push_back(cell.I());
    data.push_back(cell.V());
    data.push_back(cell.SOC());
    data.push_back(cell.T());
    data.push_back(throughputs.time()); // time
    data.push_back(throughputs.Ah());   // Ah
    data.push_back(throughputs.Wh());   // Wh
  }
};
} // namespace slide