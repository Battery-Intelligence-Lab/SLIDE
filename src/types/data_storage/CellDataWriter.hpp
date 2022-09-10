/*
 * CellDataWriter.hpp
 *
 *  Created on: 10 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "cell_data.hpp"
#include "../../settings/enum_definitions.hpp"
#include "../../utility/free_functions.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <span>
#include <cstdlib>
#include <array>
#include <span>

namespace slide {

inline void writeData(std::ofstream &file, std::span<Histogram<>> histograms)
{
  for (auto &hist : histograms)
    file << hist << "\n\n";
}

inline void writeVarAndStates(std::ofstream &file, std::span<double> varstate, std::span<double> state)
{
  //!< write throughput data, cell-to-cell variations and the battery state:
  for (auto var : varstate)
    file << var << ',';
  file << '\n';

  for (const auto st_i : state) // Time and Throughput data is written here if available.
    file << st_i << ',';
  file << "\n\n\n";
}

inline void writeData(std::ofstream &file, std::vector<double> &data)
{
  for (size_t i{}; i < data.size(); i++) {
    if (i % 7 == 0) {
      if (i != 0)
        file << '\n';
    } else
      file << ',';

    file << data[i];
  }

  data.clear(); //!< reset the index to 0 since we can overwrite the stored data
}

template <settings::cellDataStorageLevel N>
void writeDataImpl(std::ofstream &file, auto &cell, auto &dataStorage)
{
  if constexpr (N == settings::cellDataStorageLevel::storeCumulativeData)
    writeVarAndStates(file, cell.viewStates(), cell.viewVariations());

  if constexpr (N >= settings::cellDataStorageLevel::storeHistogramData)
    writeData(file, dataStorage.data);
  //!< else write nothing.
}

template <settings::cellDataStorageLevel N>
struct CellDataWriter
{
  /*
   * Writes cell data to a csv file.
   * The name of the csv file starts with the value of prefix,
   * after which the identification string of this cell is appended
   *
   * Depending on the value of DATASTORE_CELL, different things are written
   * 	0 	nothing
   * 	1 	general info about the cell and usage statistics in file xxx_cellStats.csv
   * 	2 	cycling data (I, V, T at every time step) in file xxx_cellData.csv
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