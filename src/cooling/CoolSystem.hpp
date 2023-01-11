/*
 * CoolSystem.hpp
 *
 *  Created on: 29 Apr 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../settings/settings.hpp"
#include "../StorageUnit.hpp"
#include "../types/Histogram.hpp"
#include "../types/data_storage/cell_data.hpp"

#include "data_storage/CoolSystemData.hpp"

#include <string>
#include <cstdlib>

namespace slide {
class CoolSystem
{
protected:
  //!< state
  double Tcoolant{ settings::T_ENV }; //!< temperature of the coolant [K]
  double flowrate{};                  //!< flow rate of cooling fluid [m3/s]

  //!< Coolant properties
  double fluid_rho; //!< density of cooling fluid kg / m3
  double fluid_cp;  //!< heat capcity of cooling fluid J / kG / K
  double fluid_V;   //!< total volume of coolant available for this coolsystem m3

  //!< number of cells this system needs to cool
  size_t Ncells{ 0 };

  //!< cooling system properties
  int control_strategy;      //!< integer indicating how the cooling is done
                             //	1 always on
                             //!< 	2 on/off depending on T of hottest child SU [can't be coolant since that would never heat up]
                             //!< 	3 on/off depending on T of hottest cell
                             //!< 	4 proportional to T of hottest child SU, with a minimum of 20% flow rate (else it is off) [can't be coolant since that would never heat up]
                             //!< 	5 proportional to T of hottest cell, with a minimum of 20% flow rate (else it is off)
  double control_onoff_Ton;  //!< T when the system switches on
  double control_onoff_Toff; //!< T when the system switches off
  double control_onoff_flr;  //!< flow rate when on
  double control_prop_T;     //!< T to which the proportional control cools the children
  double control_prop_gain;  //!< gain of the proportional controller
  double Across;             //!< cross section of the fluid [m2], used to convert flow rate to flow speed
  double eta{ 0.5 };         //!< efficiency

  friend struct CoolSystemData;

public:
  CoolSystemData coolData; //!< # Check this should be protected.
  CoolSystem();
  CoolSystem(size_t Ncells, int control);
  virtual ~CoolSystem() = default;

  double T() { return Tcoolant; }
  void setT(double Tnew);
  virtual double dstate(double Etot, double Echildren, double t); //!< calculate the new coolant temperature from a heat exchange of Etot

  virtual double getH();
  auto getFlr() { return flowrate; }
  auto getControl() { return control_strategy; }
  auto getEoperation() { return coolData.cData.Eoperate; }
  void reset_Eoperation() { coolData.cData.Eoperate = 0; }

  void setControl(int control) { control_strategy = control; }
  virtual void control(double Thot_local, double Thot_global);

  virtual void storeData(size_t Ncells);
  virtual void writeData(const std::string &prefix);

  auto getHeatEvac() { return coolData.cData.Qevac_life; }    //!< for unit testing, total heat evacuated from children over entire lifetime
  auto getHeatabsorbed() { return coolData.cData.Qabs_life; } //!< for unit testing, total heat evacuated from children over entire lifetime (heat capacity not constant -> cannot convert Tend-T1 to energy)
  auto getTotalTime() { return coolData.cData.t_life; }       //!< total time this coolsystem has existed for [s]

  virtual CoolSystem *copy() { return new CoolSystem(*this); }
};

} // namespace slide
