/*
 * CellDataWriter.hpp
 *
 *  Created on: 10 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <span>
#include <cstdlib>

#include "cell_data.hpp"

namespace slide {
inline void writeData(std::ofstream &file, CellCommonHist &hist)
{
  file << hist.I << "\n\n";
  file << hist.V << "\n\n";
  file << hist.T << "\n\n\n";
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

template <int N, typename SU_t, typename storage_t>
void writeDataImpl(std::ofstream &file, SU_t *ths, storage_t &storage)
{

  if constexpr (N == 1) {
    /*
     * first line gives the total utilistion of the cell (time, charge and energy throughput)
     * second line gives the parameters cell-to-cell variation (i.e. what were the values for this cell)
     * third line gives the full battery state (which for an SPM cell indicates which degradation mechanism was active)
     * 4-5 is empty
     * 6-106 give the histogram values (number of data points in each bin), for I, V and T
     * 107-108 is empty
     * 109-208 gives the edges of the bins of the histogram: edge(i-1) < bin(i) < edge(i)
     */
    writeData(file, storage.tData);
    writeVarAndStates(file, ths->viewStates(), ths->viewVariations());
    writeData(file, storage.hist);
  } else if (N == 2) {
    writeData(file, storage.cData);
  } //!< else write nothing.
}

template <int N>
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

  inline std::string getName(auto &cell, const std::string &prefix)
  {
    //!< name of the file, start with the full hierarchy-ID to identify this cell
    return prefix + "_" + cell.getFullID() + "_cellData.csv";
  }

  inline std::ofstream openFile(auto &cell, const std::string &prefix)
  {
    //!< store histograms and degradation state of cell utilisation
    std::string name = getName(cell, prefix); //!< name of the file
    std::ofstream file(name, std::ios_base::app);

    if (!file.is_open()) {
      std::cerr << "ERROR in Cell::writeData, could not open file "
                << name << '\n';
      throw 11;
    }

    return file;
  }

  inline void writeData(auto &cell, const std::string &prefix, auto &storage)
  {
    auto file = openFile(cell, prefix);
    writeDataImpl<N>(file, cell, storage);
    file.close();
  }
};

} // namespace slide