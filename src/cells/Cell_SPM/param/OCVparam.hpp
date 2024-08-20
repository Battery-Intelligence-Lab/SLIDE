/**
 * @file OCVparam.hpp
 * @brief OCV parameter class
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 10 Jul 2022
 */

#pragma once


#include "Pair.hpp"
#include <string>

namespace slide {
//!< Define a struct with the parameters of the OCV curve
//!< These parameters are calculated by the functions in determineOCV.cpp
struct OCVparam
{
  double elec_surf; //!< electrode surface
  DPair e;          //!< volume fraction of active material in the cathode /anode
  DPair thick;      //!< thickness of the cathode/anode

  Pair<std::string> name; //!< name of the CSV file with the cathode/anode OCV curve
  IPair n;                //!< number of points in the cathode/anode OCV curve

  DPair lifracini; //!< lithium fraction in the cathode/anode at 50% soC
  DPair cmax;      //!< maximum lithium concentration cathode/anode [mol m-3]

  double cap;  //!< the capacity of the cell [Ah]
  double Vmax; //!< maximum voltage of the cell [V]
  double Vmin; //!< minimum voltage of the cell [V]
};
} // namespace slide