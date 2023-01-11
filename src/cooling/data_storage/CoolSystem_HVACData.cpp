/*
 * CoolSystem_HVACData.hpp
 *
 *  Created on: 02 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "cool_data.hpp"
#include "CoolSystem_HVACData.hpp"
#include "../CoolSystem_HVAC.hpp"

#include <string>
#include <fstream>

namespace slide {
void CoolSystem_HVACData::initialise(CoolSystem_HVAC &cs, double Qac_per_cell)
{
  initialise(cs, tData, Qac_per_cell);
}

void CoolSystem_HVACData::storeData(CoolSystem_HVAC &cs)
{
  if (cs.coolData.cData.time != 0)
    storeData(cs, tData);
} //!< Do nothing.

void CoolSystem_HVACData::writeData(CoolSystem_HVAC &cs, const std::string &prefix)
{
  writeData(cs, prefix, tData);
}

//!< N=1
void CoolSystem_HVACData::initialise(CoolSystem_HVAC &cs, CoolSystem_HVACHist_t &data, double Qac_per_cell)
{
  constexpr double Qacmin = 0;
  const double Qacmax = Qac_per_cell * 1.1;

  data.Qac = Histogram<>(Qacmin, Qacmax);
  data.Eac = Histogram<>(Qacmin * cs.COP, Qacmax * cs.COP);
}

//------------------------
void CoolSystem_HVACData::storeData(CoolSystem_HVAC &cs, CoolSystem_HVACHist_t &data)
{

  const double Eacmean = cData.Eac / cs.coolData.cData.time;

  data.Eac.add(Eacmean / cs.Ncells);
  data.Qac.add(cs.Q_ac / cs.Ncells);

  //!< reset variables to calculate the mean
  cData.Eac = 0;
}

void CoolSystem_HVACData::storeData(CoolSystem_HVAC &cs, CoolSystem_HVACInst_t &data)
{
  const double Eacmean = cData.Eac / cs.coolData.cData.time;

  data.push_back({ cs.Q_ac, Eacmean }); //!< flow rate
}

//!< ----- writeData ------

auto CoolSystem_HVACData::openFile(const std::string &prefix)
{
  std::string name = prefix;
  if constexpr (settings::DATASTORE_COOL == 1)
    name += "_coolSystemStats.csv";
  else if constexpr (settings::DATASTORE_COOL == 2)
    name += "_coolSystemData.csv";
  //!< name of the file, start with the full hierarchy-ID to identify this cell
  /*
   * 1 gives the number of cells
   * 2, 3 are empty
   *
   * 4-104 give the histogram values (number of data points in each bin), for I, V and T
   * 105-106 is empty
   * 107-206 gives the edges of the bins of the histogram: edge(i-1) < bin(i) < edge(i)
   */

  //!< append the new data to the existing file
  std::ofstream file(name, std::ios_base::app);

  if (!file.is_open()) {
    std::cerr << "ERROR in CoolSystem_HVAC::writeData, could not open file " << name << '\n';
    throw 11;
  }

  return file;
}

void CoolSystem_HVACData::writeData(CoolSystem_HVAC &cs, const std::string &prefix, CoolSystem_HVACHist_t &data)
{
  auto file = openFile(prefix);
  //!< write the number of cells and the total operating energy [J]
  file << data.Qac << '\n'
       << data.Eac << '\n';

  file << "\n\n";

  file.close();
}

void CoolSystem_HVACData::writeData(CoolSystem_HVAC &, const std::string &prefix, CoolSystem_HVACInst_t &data)
{
  auto file = openFile(prefix);
  for (auto &bd : data)
    file << bd.Qac << ',' << bd.Eac << '\n';

  file.close();
  data.clear();
}
} // namespace slide