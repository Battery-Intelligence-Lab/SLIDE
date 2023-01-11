/*
 * CoolSystem_HVACData.hpp
 *
 *  Created on: 02 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "cool_data.hpp"

#include <string>
#include <fstream>

namespace slide {
class CoolSystem_HVAC;

struct CoolSystem_HVACData
{
  CoolSystem_HVACCumulative cData;              //!< Cumulative variables.
  [[no_unique_address]] CoolSystemData_t tData; //!< Time data (histogram or basic data)

public:
  void initialise(CoolSystem_HVAC &cs, double Qac_per_cell); //!< Do nothing.

  //!< void storeCumulativeData(auto Qchildren, auto Qtotal, auto t, auto E)
  //!< {
  //!<     cData.Qevac_life += Qchildren;
  //!<     cData.Qevac += Qchildren; //!< only until next data collection interval (reset to 0 in storeData)
  //!<     cData.Qabs_life += Qtotal;
  //!<     cData.t_life += t;
  //!<     cData.time += t;      //!< only until next data collection interval (reset to 0 in storeData)
  //!<     cData.time_life += t; //!< time since start of simulations
  //!<     cData.E += E;         //!< total energy required to run the fans of the cooling system in the past data-collection period [J]
  //!<     cData.Eoperate += E;  //!< total energy required to operate the coolsystem [J]
  //!< }

  void storeData(CoolSystem_HVAC &cs); //!< Do nothing.

  void writeData(CoolSystem_HVAC &cs, const std::string &prefix);

private:
  void initialise(CoolSystem_HVAC &, auto &, double Qac_per_cell) {} //!< Do nothing.

  //!< N=1
  void initialise(CoolSystem_HVAC &cs, CoolSystem_HVACHist_t &data, double Qac_per_cell);
  //------------------------
  void storeData(CoolSystem_HVAC &, auto &)
  {
    //!< if constexpr (settings::printBool::printCrit)
    //!<     std::cout << "ERROR in Cell::storeData, the settings in settings.hpp are "
    //!<                  "forbidding from storing data.\n ";

  } //!< Do nothing.

  void storeData(CoolSystem_HVAC &cs, CoolSystem_HVACHist_t &data);
  void storeData(CoolSystem_HVAC &cs, CoolSystem_HVACInst_t &data);

  //!< ----- writeData ------
  auto openFile(const std::string &prefix);

  void writeData(CoolSystem_HVAC &, auto &, auto &)
  { //!< store nothing
    if constexpr (settings::printBool::printNonCrit)
      std::cerr << "ERROR in CoolSystem_HVAC::writeData, the settings in constants.hpp are "
                << "forbidding from storing data. DATASTORE_COOL = " << settings::DATASTORE_COOL
                << '\n';
  }

  void writeData(CoolSystem_HVAC &cs, const std::string &prefix, CoolSystem_HVACHist_t &data);
  void writeData(CoolSystem_HVAC &cs, const std::string &prefix, CoolSystem_HVACInst_t &data);
};
} // namespace slide