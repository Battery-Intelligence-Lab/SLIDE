/*
 * CoolSystemData.hpp
 *
 *  Created on: 02 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "cool_data.hpp"

#include <string>
#include <fstream>

namespace slide {
class CoolSystem;

struct CoolSystemData
{
  CoolSystemCumulative cData;                   //!< Cumulative variables.
  [[no_unique_address]] CoolSystemData_t tData; //!< Time data (histogram or basic data)

public:
  void initialise(CoolSystem &cs); //!< Do nothing.

  void storeCumulativeData(auto Qchildren, auto Qtotal, auto t, auto E)
  {
    cData.Qevac_life += Qchildren;
    cData.Qevac += Qchildren; //!< only until next data collection interval (reset to 0 in storeData)
    cData.Qabs_life += Qtotal;
    cData.t_life += t;
    cData.time += t;      //!< only until next data collection interval (reset to 0 in storeData)
    cData.time_life += t; //!< time since start of simulations
    cData.E += E;         //!< total energy required to run the fans of the cooling system in the past data-collection period [J]
    cData.Eoperate += E;  //!< total energy required to operate the coolsystem [J]
  }

  void storeData(CoolSystem &cs); //!< Do nothing.

  void writeData(CoolSystem &cs, const std::string &prefix);

private:
  void initialise(CoolSystem &, auto &) {} //!< Do nothing.

  //!< N=1
  void initialise(CoolSystem &cs, CoolSystemHist_t &data);
  //------------------------
  void storeData(CoolSystem &, auto &)
  {
    //!< if constexpr (settings::printBool::printCrit)
    //!<     std::cout << "ERROR in Cell::storeData, the settings in settings.hpp are "
    //!<                  "forbidding from storing data.\n ";

  } //!< Do nothing.

  void storeData(CoolSystem &cs, CoolSystemHist_t &data);
  void storeData(CoolSystem &cs, CoolSystemInst_t &data);

  //!< ----- writeData ------
  auto openFile(const std::string &prefix);

  void writeData(CoolSystem &, auto &, auto &)
  { //!< store nothing
    if constexpr (settings::printBool::printNonCrit)
      std::cerr << "ERROR in CoolSystem::writeData, the settings in constants.hpp are "
                << "forbidding from storing data. DATASTORE_COOL = " << settings::DATASTORE_COOL
                << '\n';
  }

  void writeData(CoolSystem &cs, const std::string &prefix, CoolSystemHist_t &data);
  void writeData(CoolSystem &cs, const std::string &prefix, CoolSystemInst_t &data);
};
} // namespace slide