/*
 * CellDataStorage.hpp
 *
 *  Created on: 10 Apr 2022
 *   Author(s): Volkan Kumtepeli, Jorn Reniers
 */

#pragma once

#include "../types/Histogram.hpp"
#include "../settings/enum_definitions.hpp"

#include <utility>
#include <iostream>
#include <array>
#include <vector>

namespace slide {

template <settings::CellDataStorageLevel N>
struct CellDataStorage
{
  template <typename Cell_t>
  inline void initialise(Cell_t &) {} //!< Do nothing.

  template <typename Cell_t>
  inline void storeData(Cell_t &) {} //!< Do nothing.
};

template <>
struct CellDataStorage<settings::CellDataStorageLevel::storeHistogramData> //!< Store as histogram.
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
struct CellDataStorage<settings::CellDataStorageLevel::storeTimeData>
{
  std::vector<double> data; //!< Common data

  template <typename Cell_t>
  inline void initialise(Cell_t &) {} //!< Do nothing.

  template <typename Cell_t>
  inline void storeData(Cell_t &cell)
  {
    const auto throughputs = cell.getThroughputs();
    // #TODO just write all states, throughputs will be included.
    data.assign(data.end(),
                { cell.I(), cell.V(), cell.SOC(), cell.T(), throughputs.time(), throughputs.Ah(), throughputs.Wh() });
  }
};

#include "CellDataWriter.hpp"

template <settings::CellDataStorageLevel N>
struct CellData : public CellDataStorage<N>
{
  auto writeData(auto &cell, const std::string &prefix)
  {
    CellDataWriter<N>::writeData(cell, prefix, *this);
  }
};
} // namespace slide