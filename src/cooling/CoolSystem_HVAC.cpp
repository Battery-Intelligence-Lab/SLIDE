/*
 * CoolSystemHVAC.cpp
 *
 *  Created on: 2 Jun 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "CoolSystem_HVAC.hpp"
#include "../settings/settings.hpp"
#include "../utility/utility.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

namespace slide {
CoolSystem_HVAC::CoolSystem_HVAC()
{

  //!< Nominal power of the AC system (active cooling of the container)
  //!< For a battery of 2700 cells (continous 1C cycling) we need to evacuate about 4 kW continuously (at the begin of life)
  //!< or about 1.5 W per cell (this is the heat evacuated by the top level module from its children, during equilibrium equal to total heat generation of all cells)
  //!< which is 2.5% of total cell energy (cell power ~ 4V * 15A = 60 W)
  //!< During degradation, R can easily double so we expect to have to cool about 3 W per cell continuous
  //!< So the maximum power should be about double so we can cool down (3 W heating + 7.5 W cooling = cool down at 4.5 W or 50% faster than the cells heated up)
  constexpr double Qac_per_cell = 7.5; //!< cool 7.5 W per cell (5 times the heating load at begin of life)
  Q_ac = Qac_per_cell * Ncells;
  COP = 3; //!< value from Schimpe's paper of thermal model of battery but double check

  controlAC_onoff_Ton = 25.0_degC;
  controlAC_onoff_Toff = 20_degC; //!< if this value is < settings::T_ENV, then coolsystems can cool cells below their initial temperature (which is settings::T_ENV) resulting in negative total heat energy absorbed in cells. See unit tests for module_s
                                  //!< and then some unit tests might fail because they demand positive heating energy (cells must heat up due to cycling)
                                  //!< therefore, it is not recommended to set this T below settings::T_ENV (20 degrees at the moment)
  controlAC_onoff_Q = Q_ac;
  controlAC_prop_T = 20_degC;
  controlAC_prop_gain = 1.0 / (controlAC_onoff_Ton - controlAC_prop_T); //!< 1 at the T where the on/off control would go on, 0 at the target T

  //!< Data storage
  HVACdata.initialise(*this, Qac_per_cell);
}

CoolSystem_HVAC::CoolSystem_HVAC(size_t Ncells, int control, double Q0)
  : CoolSystem(Ncells, control)
{
  /*
   * IN
   * Q0 	ancillary losses [W] independent on the cell power. E.g. constant losses on the converter
   */

  constexpr double Qac_per_cell = 7.5; //!< cool 20 W per cell (cell power ~ 4V * 15A = 60 W so worst case 30% losses)
  Q_ac = Qac_per_cell * Ncells + Q0;
  COP = 3; //!< value from Schimpe's paper of thermal model of battery but double check

  //!< Increase the inertia of the coolsystem if Q0 >>>>  Qac_per_cell * Ncells
  //!< since the inertia in CoolSystem is based on Ncells only, and designed for 7.5W per cell
  fluid_V = fluid_V * Q_ac / (Qac_per_cell * Ncells); //!< increase proportional to the fraction of Q0 of the total heating

  controlAC_onoff_Ton = 25_degC;
  controlAC_onoff_Toff = 20_degC;
  if (control == 3) //!< todo control hotspot -> 20 to 30
    controlAC_onoff_Ton = 30_degC;
  controlAC_onoff_Q = Q_ac;
  controlAC_prop_T = 20_degC;      //!< todo control the local T to 20 degrees
  if (control == 5) {              //!< todo control the hotspot to 25 degrees
    controlAC_prop_T = 25_degC;    //!< todo do something similar to cool3
    controlAC_onoff_Ton = 30_degC; //!< todo
  }
  controlAC_prop_gain = 1.0 / (controlAC_onoff_Ton - controlAC_prop_T); //!< 1 at the T where the on/off control would go on, 0 at the target T

  //!< Data storage
  HVACdata.initialise(*this, Qac_per_cell);
}

double CoolSystem_HVAC::dstate(double Qtotal, double Qchildren, double t)
{
  /*
   * Calculate the new coolant temperature resulting from a heat exchange of Etot [J] over t [s]
   * Note that the new temperature is not set for consistency (see module)
   *
   * IN
   * Qtot 		total heat exchange, with child SUs, neighbour SUs and parent module [J]
   * Qchildren 	heat exchange with child SUs only [J]
   * 				must be identical to Qtot since the HVAC system sits on top so its module can't have neighbours or parents
   * t			time over which this heat has been exchanged [s]
   *
   * THROWS
   * 10 			total heat exchange not the same as the children's heat exchange.
   */

  //!< Ensure Qtot == Qchildren for consistency
  if (std::abs(Qtotal - Qchildren) > 1e-12) {
    std::cerr << "ERROR in coolSystem_HVAC: the total heat exchanged is " << Qtotal
              << " and heat exchange from the children is " << Qchildren
              << ". Both must be the same for an HVAC system, difference is "
              << std::abs(Qtotal - Qchildren) << '\n';
    throw 10;
  }

  //!< calculate total heat exchange, including the cooling from the HVAC
  const double Q = Qtotal - Q_ac * t; //!< - since it is cooling down

  //!< calculate the new temperature
  const double Tnew = T() + Q / (fluid_rho * fluid_cp * fluid_V);

  //!< increase the power required to run the cooling system
  //!< power scales like speed ^ 3 and rho ^1 (kinetic energy/mass * flow rate = rho * A * v * v^2 / 2)
  //!< 		see also https://fluidflowinfo.com/fan-performance-and-fan-laws/, and my undergrad course on wind energy [fan = opposite wind turbine]
  const double v = flowrate / Across;                            //!< speed of the fluid
  const double E = (fluid_rho * Across * cube(v) / 2) / eta * t; //!< this is the energy to speed up air from stationary to the required speed to cool the child SUs [W]
  const double Eac = getACoperatingPower(Q_ac, t);               //!< this is the energy to operate the AC unit

  //!< update the mean cooling power
  coolData.cData.Qevac_life += Qchildren;
  coolData.cData.Qabs_life += Q;
  coolData.cData.t_life += t;
  coolData.cData.Qevac += Qchildren;  //!< heat extracted from children, only until next data collection interval (reset to 0 in storeData)
  coolData.cData.time += t;           //!< only until next data collection interval (reset to 0 in storeData)
  coolData.cData.time_life += t;      //!< time since start of simulations
  coolData.cData.E += E;              //!< increase the total operating energy for the fan
  coolData.cData.Eoperate += E + Eac; //!< increase the total lifetime operating poer

  HVACdata.cData.Eac += Eac;          //!< increase the total operating energy for the AC unit
  HVACdata.cData.QcoolAC += Q_ac * t; //!< heat extracted from the entire coolsystem by the AC system

  //!< cout<<"\t HVAC system: Total heat energy "<<Qtotal<<", HVAC cooling "<<Q_ac*t<<" giving net energy "<<Q<<" and new T "<<Tnew<<" after "<<t<<" seconds"<<endl;
  return Tnew;
}

double CoolSystem_HVAC::getACoperatingPower(double Qac, double t)
{
  /*
   * function to calculate the operating power of the AC system.
   * If the environment is hotter than the target temperature, we have to use active cooling
   * 	Note that active cooling requires a lot of power: if you have 10% losses, you will need 30% of power to cool them away (COP = 3), so overall efficiency = 60%
   * Else, we can simply suck in cold air from the environment with a fan
   *
   * IN
   * Qac 	total heat power which needs to be removed [W]
   * t 	time over which this heat needs to be removed [s]
   */

  using settings::T_ENV;

  if (Qac == 0 || t == 0)
    return 0;

  //!< if the environment is too hot, use active cooling with the AC system
  if (T_ENV > (T() - 5))
    return COP * Qac * t;

  //!< else we can suck in cold air
  const double flr = Qac / (fluid_rho * fluid_cp * (T() - T_ENV)); //!< flow rate of environment air we need (cooling energy = flow rate * rho * cp * dT)
  const double Acr = Across * 2;                                   //!< assume the surface area of the fan to suck in air is twice the surface area of the fan to cool the cells
  const double v = flr / Acr;                                      //!< speed of the air
  const double E = (fluid_rho * Acr * cube(v) / 2) / eta * t;      //!< energy to operate the fan

  return E;
}

void CoolSystem_HVAC::control(double Thot_local, double Thot_global)
{
  /*
   * Control this cool system by adapting the flow rate of the coolant by controlling the fan.
   * Additionally, we also adapt the cooling power of the AC system.
   *
   * depending on the value of control_strategy
   * 		1	fan: 	always on
   * 			AC: 	on if local T > minimum temperature of container
   * 		2 	fan: 	on / off depending on T of the hottest child of the SU (i.e. local hot-spot)
   * 			AC: 	on if local T > Tmax, off if local T < minimum temperature of container
   * 		3	fan: 	on / off depending on T of the hottest cell ultimately connected to this SU (i.e. global hot-spot)
   * 			AC: 	on if T_hotspot > Tmax, off if local T < minimum temperature of container
   * 						note: same 'off' criterium as 2 but on based on hot spot T instead of T of container
   * 		4	fan: 	proportional to T of hottest child SU, with a minimum of 20% flow rate (else it is off)
   * 			AC: 	on proportional to local T, off if local T < minimum temperature of container
   * 		5	fan: 	proportional to T of hottest cell, with a minimum of 20% flow rate (else it is off)
   * 			AC: 	on proportional to hot spot temperature, off if local T < minimum temperature of container
   * 						note: same 'off' criterium as 4 but on based on hot spot T instead of T of container
   *
   * Note that here we can use the local temperature of the container as on-criterium
   * the control of the flowrate is based on the child-SUs. So if they get hot, the fan inside the container will activate (flowrate > 0)
   * but there is no AC yet, so the container will heat up
   * and that will trigger the AC system
   *
   * IN
   * Thot_local 	the local hottest temperature (i.e. T of the hottest child-module) [K]
   * Thot_global 	the global hot spot of the module or battery [K]
   */

  //!< Call the base class function to control the fan inside the container to cool the racks
  CoolSystem::control(Thot_local, Thot_global);

  //!< variables
  double Terr;

  //!< control the cooling power of the AC system
  //!< 	note that the off criterium is always the same, never cool below the minimum temperature
  //!< 	this is because we are controlling an active system, so in theory it would keep cooling to freezing, -100 and even -Inf degrees if you wouldn't stop it
  switch (control_strategy) {
  case 1: //!< always on unless T of the container < minimum temperature
    if (T() < controlAC_onoff_Toff)
      Q_ac = 0;
    else
      Q_ac = controlAC_onoff_Q;
    break;
  case 2: //!< on/off depending on child T -> on and off based on T of container
    if (T() <= controlAC_onoff_Toff)
      Q_ac = 0;
    else if (T() >= controlAC_onoff_Ton)
      Q_ac = controlAC_onoff_Q;
    break;

  case 3: //!< on/off depending on global hot spot temperature -> on based on T of hottest cell, off based on T of container

    if (T() <= controlAC_onoff_Toff)
      Q_ac = 0;
    else if (Thot_global >= controlAC_onoff_Ton)
      Q_ac = controlAC_onoff_Q;
    break;

  case 4: //!< proportional to T of hottest child SU -> proportional to T of container

    //!< error on the local temperature
    Terr = T() - controlAC_prop_T;
    Terr = std::max(Terr, 0.0); //!< ensure we never go negative

    if (T() <= controlAC_onoff_Toff)
      Q_ac = 0;
    else if (control_prop_gain * Terr > 0.2)                                //!< avoid very small cooling powers
      Q_ac = controlAC_onoff_Q * std::min(controlAC_prop_gain * Terr, 1.0); //!< can't do more than the design-power
    else
      Q_ac = 0;
    break;

  case 5: //!< proportional to global hot spot temperature -> on based on T of hottest cell, off based on T of container

    //!< error on the global temperature
    Terr = Thot_global - controlAC_prop_T;
    Terr = std::max(Terr, 0.0); //!< ensure we never go negative

    if (T() <= controlAC_onoff_Toff)
      Q_ac = 0;
    else if (control_prop_gain * Terr > 0.2)                                //!< avoid very small cooling powers
      Q_ac = controlAC_onoff_Q * std::min(controlAC_prop_gain * Terr, 1.0); //!< can't do more than the design-power
    else
      Q_ac = 0;

    break;

  default: //!< always on
    break;
  }
}

void CoolSystem_HVAC::storeData(size_t Ncellsi)
{
  /*
   * Add another data point in the statistics.
   *
   * IN
   * Ncells 	number of cells connected to the module of this CoolSystem
   * 			some usage statistics are stored 'per cell' such that the edges of the bins stay the same for all CoolSystems
   * 			even though cool systems which have to cool 1000s of cells will obviously have much higher powers compared to a system of just 5 cells
   *
   * THROWS
   * 10 		wrong number of cells
   */

  if constexpr (settings::DATASTORE_COOL == 1) {
    if (Ncellsi != Ncells) {
      std::cerr << "Error in CoolSystem_HVAC. This system has " << Ncells
                << " cells connected, but now you are storing data for a system of "
                << Ncellsi << ". throwing an error\n";
      throw 10;
    }

    if (Ncells == 0) {
      std::cerr << "Error in CoolSystem_HVAC. This system has " << Ncells
                << " cells connected, so you cannot store values per cell connected. Throwing an error\n";
      throw 10;
    }
  }

  //!< Store additional data of the AC system
  HVACdata.storeData(*this);

  //!< Call parent function to store conventional data
  CoolSystem::storeData(Ncellsi); //!< sets ttot to 0 so do this after wards
}

void CoolSystem_HVAC::writeData(const std::string &prefix)
{
  /*
   * Writes data to a csv file.
   * The name of the csv file starts with the value of prefix, after which the identification string of this cell is appended
   *
   *
   */

  CoolSystem::writeData(prefix);

  //!< store histograms and degradation state of cell utilisation
  HVACdata.writeData(*this, prefix); //!< #TODO -> since we are doing append we cannot write like this IMPORTANT!!!!!
}
} // namespace slide