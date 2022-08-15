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
  {
    auto I_bins = hist.I.viewBinValues();
    auto V_bins = hist.V.viewBinValues();
    auto T_bins = hist.T.viewBinValues();

    //!< write bin values
    for (size_t i = 0; i < I_bins.size(); i++)
      file << I_bins[i] << ',' << V_bins[i] << ',' << T_bins[i] << '\n'; //!< #CHECK change how you write the values to remove repetition.
    file << "\n\n";
  }

  {
    auto I_edges = hist.I.getEdgeValues();
    auto V_edges = hist.V.getEdgeValues();
    auto T_edges = hist.T.getEdgeValues();

    //!< write bin edges
    for (size_t i = 0; i < I_edges.size(); i++)
      file << I_edges[i] << ',' << V_edges[i] << ',' << T_edges[i] << '\n';
    //!< data is in bin i if edge[i-1] <= data < edge[i]
    //!< except bin[0] if data < edge[0]
    //!<    and bin[N] if data >= edge[N-1] #CHECK if it is still true.
    file << "\n\n\n";
  }
}

inline void writeVarAndStates(std::ofstream &file, std::span<double> varstate, std::span<double> state)
{
  //!< write throughput data, cell-to-cell variations and the battery state:
  for (auto var : varstate)
    file << var << ',';
  file << '\n';

  for (const auto st_i : state)
    file << st_i << ',';
  file << "\n\n\n";
}

inline void writeData(std::ofstream &file, CellCumulativeData tData)
{
  //!< write throughput data, cell-to-cell variations and the battery state:
  file << tData.Time << ',' << tData.Ah << ',' << tData.Wh << '\n';
}

inline void writeData(std::ofstream &file, std::vector<CommonData> &cData)
{
  for (const auto &c : cData)
    file << c.I << ',' << c.V << ',' << c.SOC
         << ',' << c.T << ',' << c.Ah << ',' << c.Wh
         << ',' << c.Time << '\n';

  //!< reset the index to 0 since we can overwrite the stored data
  cData.clear();
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

  inline std::string getName(auto *ths, const std::string &prefix)
  {
    //!< name of the file, start with the full hierarchy-ID to identify this cell
    return prefix + "_" + ths->getFullID() + "_cellData.csv";
  }

  inline std::ofstream openFile(auto *ths, const std::string &prefix)
  {
    //!< store histograms and degradation state of cell utilisation
    std::string name = getName(ths, prefix); //!< name of the file
    std::ofstream file(name, std::ios_base::app);

    if (!file.is_open()) {
      std::cerr << "ERROR in Cell::writeData, could not open file "
                << name << '\n';
      throw 11;
    }

    return file;
  }

  inline void writeData(auto *ths, const std::string &prefix, auto &storage)
  {
    auto file = openFile(ths, prefix);
    writeDataImpl<N>(file, ths, storage);
    file.close();
  }
};

} // namespace slide