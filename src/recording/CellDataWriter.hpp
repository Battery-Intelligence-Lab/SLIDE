/**
 * @file CellDataWriter.hpp
 * @brief Header file for cell data writing functionality
 *
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 10 Apr 2022
 */

#pragma once

#include "../settings/enum_definitions.hpp"
#include "../utility/free_functions.hpp"
#include "../types/Histogram.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <span>
#include <cstdlib>
#include <array>
#include <span>
#include <variant>

namespace slide {

/**
 * @brief Writes histogram data to a file
 *
 * @param file Output file stream to write to
 * @param histograms Span of histograms to write
 */
inline void writeData(std::ofstream &file, std::span<Histogram<>> histograms)
{
  for (auto &hist : histograms)
    file << hist << "\n\n";
}

/**
 * @brief Implementation of data writing based on storage level
 *
 * @tparam N The cell data storage level
 * @param file Output file stream to write to
 * @param cell Reference to the cell object
 * @param dataStorage Reference to the data storage object
 */
template <settings::CellDataStorageLevel N>
void writeDataImpl(std::ofstream &file, auto &cell, auto &dataStorage)
{
  if constexpr (N >= settings::CellDataStorageLevel::storeHistogramData)
    free::write_data(file, dataStorage.data, 7);
  //!< else write nothing.
}

/**
 * @brief Template class for writing cell data to CSV files
 *
 * @tparam N The cell data storage level that determines what data is written
 *
 * The storage levels determine what data is written:
 * - 0: Nothing is written
 * - 1: General cell info and usage statistics in xxx_cellStats.csv
 * - 2: Cycling data (I, V, T at every time step) in xxx_cellData.csv
 */
template <settings::CellDataStorageLevel N>
struct CellDataWriter
{
  /**
   * @brief Writes cell data to a CSV file
   *
   * The filename is constructed by appending the cell's identification string
   * to the provided prefix, followed by "cellData.csv"
   *
   * @param cell Reference to the cell object
   * @param prefix Prefix for the output filename
   * @param storage Reference to the data storage object
   */
  inline static void writeData(auto &cell, const std::string &prefix, auto &storage)
  {
    constexpr auto suffix = "cellData.csv";
    auto file = free::openFile(cell, PathVar::results, prefix, suffix);
    writeDataImpl<N>(file, cell, storage);
    file.close();
  }
};

} // namespace slide