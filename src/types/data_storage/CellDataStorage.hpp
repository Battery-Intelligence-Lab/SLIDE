/*
 * CellDataStorage.hpp
 *
 *  Created on: 10 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <utility>
#include <iostream>

#include "../Histogram.hpp"
#include "cell_data.hpp"
#include "../../settings/enum_definitions.hpp"

namespace slide {

template <settings::cellDataStorageLevel N>
struct CellDataStorage
{
  template <typename Tcell>
  inline void initialise(Tcell &) {} //!< Do nothing.

  template <typename Tcell>
  inline void storeData(Tcell &)
  {
    if constexpr (settings::printBool::printCrit)
      std::cout << "ERROR in Cell::storeData, the settings in constant.h are forbidding from storing data.\n";

  } //!< Do nothing.
};

template <>
struct CellDataStorage<settings::cellDataStorageLevel::storeHistogramData> //!< Store as histogram.
{
  CellCommonHist data;

  template <typename Tcell>
  inline void initialise(Tcell &cell) //!< Initialise the histograms.
  {
    data.I = Histogram<>(-cell.Cap(), cell.Cap()); //!< 1C charge/discharge
    data.V = Histogram<>(cell.Vmin(), cell.Vmax());
    data.T = Histogram<>(cell.Tmin(), cell.Tmax());
  }

  template <typename Tcell>
  inline void storeData(Tcell &cell)
  {
    data.I.add(cell.I());
    data.V.add(cell.V());
    data.T.add(cell.T());
  }
};

template <>
struct CellDataStorage<settings::cellDataStorageLevel::storeTimeData>
{
  std::vector<double> data; //!< Common data

  template <typename Tcell>
  inline void initialise(Tcell &) {} //!< Do nothing.

  template <typename Tcell>
  inline void storeData(Tcell &cell)
  {
    data.push_back(cell.I());
    data.push_back(cell.V());
    data.push_back(cell.SOC());
    data.push_back(cell.T());
    data.push_back(0); // time
    data.push_back(0); // Ah
    data.push_back(0); // Wh
  }
};
} // namespace slide