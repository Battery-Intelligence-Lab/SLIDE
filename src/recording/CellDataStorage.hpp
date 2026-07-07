/**
 * @file CellDataStorage.hpp
 * @brief Template classes for cell data storage with different storage levels
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 10 Apr 2022
 */

#pragma once

#include "../types/Histogram.hpp"
#include "../settings/enum_definitions.hpp"
#include "CellDataWriter.hpp" //!< must be included at global scope (it pulls in <variant> etc.)

#include <utility>
#include <iostream>
#include <array>
#include <vector>

namespace slide {

/**
 * @brief Template class for cell data storage based on storage level
 *
 * @tparam N The cell data storage level that determines what data is stored
 *
 * This is the base template that does nothing by default. Specializations
 * are provided for different storage levels to store various types of data.
 */
template <settings::CellDataStorageLevel N>
struct CellDataStorage
{
  /**
   * @brief Initialize the data storage for a cell
   * @param cell Reference to the cell object
   */
  template <typename Cell_t>
  inline void initialise(Cell_t &) {} //!< Do nothing.

  /**
   * @brief Store data from a cell
   * @param cell Reference to the cell object
   */
  template <typename Cell_t>
  inline void storeData(Cell_t &) {} //!< Do nothing.
};

/**
 * @brief Specialization for storing cell data as histograms
 *
 * This specialization stores current, voltage, and temperature data
 * in histogram format for statistical analysis.
 */
template <>
struct CellDataStorage<settings::CellDataStorageLevel::storeHistogramData>
{
  std::array<Histogram<>, 3> data; //!< Array of histograms for I, V, T data

  /**
   * @brief Initialize the histograms with appropriate ranges
   * @param cell Reference to the cell object to get capacity and limits
   */
  template <typename Cell_t>
  inline void initialise(Cell_t &cell)
  {
    data[0] = Histogram<>(-cell.Cap(), cell.Cap());  //!< Current histogram: ±1C range
    data[1] = Histogram<>(cell.Vmin(), cell.Vmax()); //!< Voltage histogram: cell voltage range
    data[2] = Histogram<>(cell.Tmin(), cell.Tmax()); //!< Temperature histogram: cell temperature range
  }

  /**
   * @brief Store current cell data in histograms
   * @param cell Reference to the cell object
   */
  template <typename Cell_t>
  inline void storeData(Cell_t &cell)
  {
    data[0].add(cell.I()); //!< Add current to histogram
    data[1].add(cell.V()); //!< Add voltage to histogram
    data[2].add(cell.T()); //!< Add temperature to histogram
  }
};

/**
 * @brief Specialization for storing cell data as time series
 *
 * This specialization stores detailed time-series data including
 * current, voltage, SOC, temperature, and throughput information.
 */
template <>
struct CellDataStorage<settings::CellDataStorageLevel::storeTimeData>
{
  std::vector<double> data; //!< Vector to store time series data

  /**
   * @brief Initialize the data storage (no initialization needed for vector)
   * @param cell Reference to the cell object
   */
  template <typename Cell_t>
  inline void initialise(Cell_t &) {} //!< Do nothing.

  /**
   * @brief Store current cell state and throughput data
   * @param cell Reference to the cell object
   *
   * Stores: I, V, SOC, T, time, Ah, Wh in sequence
   */
  template <typename Cell_t>
  inline void storeData(Cell_t &cell)
  {
    const auto throughputs = cell.getThroughputs();
    // #TODO just write all states, throughputs will be included.
    //!< APPEND this timestep's record; assign() would REPLACE the whole history
    //!< (and the iterator + initializer-list form was ill-formed).
    data.insert(data.end(),
                { cell.I(), cell.V(), cell.SOC(), cell.T(), throughputs.time(), throughputs.Ah(), throughputs.Wh() });
  }
};

/**
 * @brief Combined cell data storage and writing class
 *
 * @tparam N The cell data storage level
 *
 * This class combines data storage functionality with data writing
 * capabilities by inheriting from CellDataStorage and providing
 * a writeData method that uses CellDataWriter.
 */
template <settings::CellDataStorageLevel N>
struct CellData : public CellDataStorage<N>
{
  /**
   * @brief Write stored data to file
   * @param cell Reference to the cell object
   * @param prefix Prefix for the output filename
   */
  auto writeData(auto &cell, const std::string &prefix)
  {
    CellDataWriter<N>::writeData(cell, prefix, *this);
  }
};
} // namespace slide