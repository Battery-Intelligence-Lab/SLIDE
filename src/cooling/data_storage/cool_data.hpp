/*
 * cool_data.hpp
 *
 *  Created on: 02 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../../types/Histogram.hpp"
#include "../../settings/settings.hpp"

#include <type_traits>
#include <vector>

namespace slide {
struct Empty
{
};

struct CoolSystemCumulative
{
  double Qevac{ 0 };      //!< heat evacuated from the child SUs in the past data-collection period [J]
  double Qevac_life{ 0 }; //!< variable for unit testing, total heat evacuated from child SUs over the entire lifetime [J]
  double Qabs_life{ 0 };  //!< variable for unit testing, total heat absorbed by heating up the coolant [J]
  double t_life{ 0 };     //!< total time in the module's life [s]
  double E{ 0 };          //!< total energy required to run the fans of the cooling system in the past data-collection period [J]
  double Eoperate{ 0 };   //!< total energy required to operate the coolsystem [J]
  double time{ 0 };       //!< total time  in the past data-collection period [s]
  double time_life{ 0 };  //!< total time since start of simulations
};

struct CoolSystemInst //!< Instantenous data.
{
  double time{}; //!< total time
  double Q{};    //!< cooling power
  double E{};    //!< operating energy
  double flr{};  //!< flow rate
};

struct CoolSystemHist
{
  //!< histograms for coolant temperature, cooling power to child SUs [W/cell],
  //!< flow rate [m^3 s^-1 per cell], power to operate the fans of the cooling system [W/cell]
  Histogram<> Q, flr, E, T;
};

//!< Data storage
using CoolSystemInst_t = std::vector<CoolSystemInst>;
using CoolSystemHist_t = CoolSystemHist;

using CoolSystemData_t =
  std::conditional_t<settings::DATASTORE_COOL == 1, CoolSystemHist_t, std::conditional_t<settings::DATASTORE_COOL == 2, CoolSystemInst_t, Empty>>;
//-------------------------------------------------------------------------

struct CoolSystem_HVACCumulative
{
  double Eac{ 0 };     //!< total operating power for the AC unit in the present data collection time interval
  double QcoolAC{ 0 }; //!< total heat evacuated by the AC system [J]
};

struct CoolSystem_HVACHist
{
  //!< histograms for coolant temperature, cooling power to child SUs [W/cell],
  //!< flow rate [m^3 s^-1 per cell], power to operate the fans of the cooling system [W/cell]
  Histogram<> Qac, Eac;
};

struct CoolSystem_HVACInst //!< #TODO -> in future try to combine.
{
  double Qac; //!< cooling power of the AC unit
  double Eac; //!< operating energy of the AC unit
};

using CoolSystem_HVACHist_t = CoolSystem_HVACHist;
using CoolSystem_HVACInst_t = std::vector<CoolSystem_HVACInst>;

using CoolSystem_HVACData_t =
  std::conditional_t<settings::DATASTORE_COOL == 1, CoolSystem_HVACHist_t, std::conditional_t<settings::DATASTORE_COOL == 2, CoolSystem_HVACInst_t, Empty>>;

} // namespace slide
