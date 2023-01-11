/*
 * CoolSystemData.hpp
 *
 *  Created on: 02 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "cool_data.hpp"
#include "CoolSystemData.hpp"
#include "../CoolSystem.hpp"
#include "../utility/utility.hpp"

#include <string>
#include <fstream>

namespace slide {
void CoolSystemData::initialise(CoolSystem &cs)
{
  initialise(cs, tData);
} //!< Do nothing.

void CoolSystemData::storeData(CoolSystem &cs)
{
  if (cData.time != 0)
    storeData(cs, tData);
} //!< Do nothing.

void CoolSystemData::writeData(CoolSystem &cs, const std::string &prefix)
{
  writeData(cs, prefix, tData);
}

//!< N=1
void CoolSystemData::initialise(CoolSystem &cs, CoolSystemHist_t &data)
{
  constexpr auto flr_max = 1.1 * settings::cool::flowrate_perCell;
  const auto v_max = flr_max / cs.Across;
  const auto E_max = cs.fluid_rho * cs.Across * v_max * v_max * v_max / cs.eta;

  data.flr = Histogram<>(0, flr_max); //!< 0 to 0.02 m3/s per cell
  data.Q = Histogram<>(0, 5);
  data.T = Histogram<>(0_degC, 100_degC);
  data.E = Histogram<>(0, E_max); //!< #TODO if it is really correct?
}

//------------------------
void CoolSystemData::storeData(CoolSystem &cs, CoolSystemHist_t &data)
{
  const double Qmean_percell = cData.Qevac / cData.time / cs.Ncells;
  const double Emean_percell = cData.E / cData.time / cs.Ncells;
  const double flr_percell = cs.getFlr() / cs.Ncells;

  data.E.add(Emean_percell);
  data.T.add(cs.T());
  data.Q.add(Qmean_percell);
  data.flr.add(flr_percell);

  //!< reset variables to calculate the mean
  cData.E = 0;
  cData.Qevac = 0;
  cData.time = 0;
}

void CoolSystemData::storeData(CoolSystem &cs, CoolSystemInst_t &data)
{
  data.push_back({ cData.time_life,          //!< total time
                   cData.Qevac / cData.time, //!< cooling power
                   cData.E / cData.time,     //!< operating energy
                   cs.getFlr() });           //!< flow rate
}

//!< ----- writeData ------

auto CoolSystemData::openFile(const std::string &prefix)
{
  std::string name = prefix + "_coolSystemStats.csv";
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
    std::cerr << "ERROR in CoolSystem::writeData, could not open file " << name << '\n';
    throw 11;
  }

  return file;
}

void CoolSystemData::writeData(CoolSystem &cs, const std::string &prefix, CoolSystemHist_t &data)
{
  auto file = openFile(prefix);
  //!< write the number of cells and the total operating energy [J]
  file << cs.Ncells << '\n';
  file << cData.E << '\n';
  file << "\n\n";

  file << data.T << '\n'
       << data.Q << '\n'
       << data.flr << '\n'
       << data.E << '\n';

  file << "\n\n";

  file.close();
}

void CoolSystemData::writeData(CoolSystem &, const std::string &prefix, CoolSystemInst_t &data)
{
  auto file = openFile(prefix);
  for (auto &bd : data)
    file << bd.time << ',' << bd.Q << ','
         << bd.E << ',' << bd.flr << '\n';

  file.close();
  data.clear();
}
} // namespace slide