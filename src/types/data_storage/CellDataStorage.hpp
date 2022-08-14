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

namespace slide {
template <int N>
struct CellDataStorage
{
  template <typename T>
  inline void initialise(T *) {} // Do nothing.

  template <typename... T>
  inline void storeCumulativeData(T...) {} // Do nothing.

  template <typename... T>
  inline void storeInstantenousData(T...)
  {
    if constexpr (settings::printBool::printCrit)
      std::cout << "ERROR in Cell::storeData, the settings in constant.h are forbidding from storing data.\n";

  } // Do nothing.
  inline CellCumulativeData getThroughputData() { return {}; }
  inline void setThroughputData(CellCumulativeData) {}
};

template <>
struct CellDataStorage<1> // Store as histogram.
{
  CellCumulativeData tData; // Throughput Data
  CellCommonHist hist;

  template <typename T>
  inline void initialise(T *ths) // Initialise the histograms.
  {
    hist.I = Histogram<>(-ths->Cap(), ths->Cap()); // 1C charge/discharge
    hist.V = Histogram<>(ths->Vmin(), ths->Vmax());
    hist.T = Histogram<>(ths->Tmin(), ths->Tmax());
  }

  inline void storeCumulativeData(auto dt, auto dAh, auto dWh)
  {
    tData.Time += dt;
    tData.Ah += dAh;
    tData.Wh += dWh;
  }

  inline void storeInstantenousData(auto I, auto V, auto SOC, auto T)
  {
    hist.I.add(I);
    hist.V.add(V);
    hist.T.add(T);
  }

  inline CellCumulativeData getThroughputData() { return tData; }
  inline void setThroughputData(CellCumulativeData tData_)
  {
    tData = std::move(tData_);
  }
};

template <>
struct CellDataStorage<2>
{
  CellCumulativeData tData;      // Throughput Data
  std::vector<CommonData> cData; // Common data

  template <typename T>
  inline void initialise(T *) {} // Do nothing.

  inline void storeCumulativeData(auto dt, auto dAh, auto dWh)
  {
    tData.Time += dt;
    tData.Ah += dAh;
    tData.Wh += dWh;
  }

  inline void storeInstantenousData(auto I, auto V, auto SOC, auto T)
  {
    cData.push_back(CommonData(
      I, V, SOC, T, tData.Ah, tData.Wh, tData.Time));
  }

  inline CellCumulativeData getThroughputData() { return tData; }

  inline void setThroughputData(CellCumulativeData tData_)
  {
    tData = std::move(tData_);
  }
};
} // namespace slide