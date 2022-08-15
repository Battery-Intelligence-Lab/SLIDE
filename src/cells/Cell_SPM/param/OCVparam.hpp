/*
 * OCVparam.hpp
 *  Created on: 10 Jul 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <string>

namespace slide {
//!< Define a struct with the parameters of the OCV curve
//!< These parameters are calculated by the functions in determineOCV.cpp
struct OCVparam
{
  double elec_surf; //!< electrode surface
  double ep;        //!< volume fraction of active material in the cathode
  double en;        //!< volume fraction of active material in the anode
  double thickp;    //!< thickness of the cathode
  double thickn;    //!< thickness of the anode

  std::string namepos; //!< name of the CSV file with the cathode OCV curve
  std::string nameneg; //!< name of the CSV file with the anode OCV curve
  int np;              //!< number of points in the cathode OCV curve
  int nn;              //!< number of points in the anode OCV curve

  double lifracpini; //!< lithium fraction in the cathode at 50% soC
  double lifracnini; //!< lithium fraction in the anode at 50% SOC
  double cmaxp;      //!< maximum lithium concentration in the cathode [mol m-3]
  double cmaxn;      //!< maximum lithium concentration in the anode [mol m-3]
  double cap;        //!< the capacity of the cell [Ah]
  double Vmax;       //!< maximum voltage of the cell [V]
  double Vmin;       //!< minimum voltage of the cell [V]
};
} // namespace slide