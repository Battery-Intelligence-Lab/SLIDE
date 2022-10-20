/*
 * CoolSystemHVAC.hpp
 *
 * This is the coolsystem of the entire battery container.
 * It has the conventional coolsystem-stuff to cool its children (e.g. racks).
 * But additionally, it has an AC unit which exchanged heat with the environment.
 * 	i.e. this coolsystem is not cooled by the coolsystem of its parent module like others
 * 	but it actively cools itself.
 *
 *  Created on: 2 Jun 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "CoolSystem.hpp"
#include "../types/Histogram.hpp"

#include "data_storage/CoolSystem_HVACData.hpp"

namespace slide {
class CoolSystem_HVAC : public CoolSystem
{
protected:
  double Q_ac; //!< cooling power from the AC system [W]. This value will be controlled similarly to the flowrate inside the container
  double COP;  //!< coefficient of performance of the AC system [-]

  double controlAC_onoff_Ton;  //!< T when the AC unit switches on
  double controlAC_onoff_Toff; //!< T when the AC unit switches off, ie minimum temperature
  double controlAC_onoff_Q;    //!< cooling power when on
  double controlAC_prop_T;     //!< T to which the proportional control cools the container
  double controlAC_prop_gain;  //!< gain of the proportional controller for the AC unit

  double getACoperatingPower(double Qac, double t); //!< calculate the operating power of the AC unit

  //!< Data storage
  friend struct CoolSystem_HVACData;

public:
  CoolSystem_HVACData HVACdata; //!< #TODO -> needs to be protected.

  CoolSystem_HVAC();
  CoolSystem_HVAC(size_t Ncells, int control, double Q0);

  double getQcoolAC_tot() { return HVACdata.cData.QcoolAC; };

  double dstate(double Etot, double Echildren, double t) override; //!< calculate the new coolant temperature from a heat exchange of Etot
  void control(double Thot_local, double Thot_global) override;

  virtual void storeData(size_t Ncells) override;
  virtual void writeData(const std::string &prefix) override;

  CoolSystem_HVAC *copy() override { return new CoolSystem_HVAC(*this); }
};

} // namespace slide