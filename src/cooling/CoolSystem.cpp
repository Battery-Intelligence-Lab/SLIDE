/*
 * CoolSystem.cpp
 *
 *  Created on: 29 Apr 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "CoolSystem.hpp"
#include "../modules/Module_p.hpp"
#include "../modules/Module_s.hpp"
#include "../settings/settings.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

namespace slide {
CoolSystem::CoolSystem()
{
  /*
   * Coolant properties:
   *
   * WATER
   * 	rho = 1000 kg m-3
   * 	cp = 4200						see https://en.wikipedia.org/wiki/Table_of_specific_heat_capacities
   * 	k = 0.5918 [W m-1 K-1]
   *
   * AIR
   * 	rho = 1.225 [kg m-3]
   * 	cp = 1.005*1000 [J kg-1 K-1] 	values from Schimpe's paper on energy efficiency evaluation
   * 	k = 0.025 [W m-1 K-1]			https://en.wikipedia.org/wiki/List_of_thermal_conductivities
   * 	Tcoolant heats up by 1 to 2 degrees for a module consisting of 3 cells
   *
   * SIZING
   * 	our orange fan ('the jet') has the following properties
   * 		power: 		500 W
   * 		flow rate: 	65 m3 / min
   * 		diameter: 	0.3 m
   * 		surface: 	0.07 m3
   * 		flow speed:	14 m/s
   * 		rpm:		2800
   * 		pressure:	385 Pa
   * If we use the following values
   * 		flow rate per cell: 	0.0005 m3 / s per cell
   * 		surface area per cell: 	0.000025 m2 per cell
   * then for the big battery (2700 cells) we get properties similar to the jet
   * 		flow rate: 	1.35 m3/s
   * 		surface: 	0.0675 m2
   * 		speed: 		20 m/s
   * 		power: 		661 W (assuming 50% efficiency)
   *
   */

  //!< Coolant properties
  fluid_rho = 1.225;  //!< density  kg / m3
  fluid_cp = 1.005e3; //!< heat capacity  [J / kg K]
  Across = 25e-6;     //!< per cell value (in next constructor this is multiplied by nrCells)
                      //!< 	= 0.0675 m2 for 2700 cells. Our orange jet in the lab has A = 0.07 m2
  eta = 0.5;          //!< assume 50% efficiency

  //!< System controls
  constexpr double V_perCell{ 1.0 }; //!< this might be unrealistically high, but we need a high number for the numerical stability of the system. See dstate

  constexpr auto n_modules = settings::MODULE_NSUs_MAX; //!< Because we do not know number of modules.

  flowrate = settings::cool::flowrate_perCell * n_modules; //!< #TODO Why multiply with Module NSUs_MAX;
  fluid_V = V_perCell * n_modules;
  control_strategy = 1;                      //!< #TODO -> magic number to enum.
  control_onoff_Ton = PhyConst::Kelvin + 35; //!< on/off control: go on at 35 degrees
  const auto t1 = settings::T_ENV + 5;
  control_onoff_Toff = std::max(25.0_degC, t1); //!< on/off control: go off at 25 degrees, or 5 degrees above environmental temperature
  control_onoff_flr = flowrate;
  control_prop_T = PhyConst::Kelvin + 25;
  control_prop_gain = 1.0 / (control_onoff_Ton - control_prop_T); //!< 1 at the T where the on/off control would go on

  //!< Data storage
  coolData.initialise(*this);
}

CoolSystem::CoolSystem(size_t Ncellsi, int control) : CoolSystem()
{
  /*
   * Make a cooling system with the specified number of cells
   *
   * Increase flow rate and cross section.
   * This will keep the flow speed (so h) the same, but reduce the temperature increase of the coolant
   * since there is more mass to cool
   *
   * IN
   * Ncells 	number of cells which this cool system will ultimately have to cool (number of cells, not number of child SUs!)
   * control 	integer indicating the control settings for this cool system
   * 				1 	always on
   * 				2 	on/off depending on T of hottest child SU
   * 						note: cannot be depending on T of coolant, since off -> heat exchange = 0 -> coolant does not heat up -> would never start again
   * 				3 	on/off depending on T of hottest cell
   * 				4	proportional to T of hottest child SU, with a minimum of 20% flow rate (else it is off)
   * 				5	proportional to T of hottest cell, with a minimum of 20% flow rate (else it is off)
   */
  Ncells = Ncellsi;
  Across = 0.000025 * Ncells; //!< values: see previous constructor
  flowrate = settings::cool::flowrate_perCell * Ncells;
  fluid_V = 1.0 * Ncells; //!< this might be unrealistically high, but we need a high number for the numerical stability of the system. See dstate

  control_strategy = control;
  control_onoff_flr = flowrate;
}

void CoolSystem::setT(double Tnew)
{

  //!< Check the new temperature is valid
  if (Tnew < PhyConst::Kelvin || Tnew > PhyConst::Kelvin + 75.0 || std::isnan(Tnew)) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in CoolSystem::setT, the new temperature of "
                << Tnew << " is outside the allowed range from (273+0) K to (273+75) K.\n";
    throw 99;
  }

  //!< set the new temperature
  Tcoolant = Tnew;
}

double CoolSystem::getH()
{
  /*
   * Return the convective heat transfer constant for the cooling system.
   * This is a function of the speed of the fluid (a linear function of the flow speed)
   *
   * Rayleigh number https://en.wikipedia.org/wiki/Rayleigh_number
   * 	Ra = diffusion / convection = (l^2 / a) / (l / u) = u * l / a
   * 		where 	a = k / (rho*cp) [m^2 / s].
   * 					Value for water = 1.4e-7
   * 					value for air = 2.01 e-5
   * 				l = length scale [m],
   * 					assuming one cell it is about 0.2 m
   * 				u = speed [m / s], proportional to the flow rate (constant cross area)
   * 					assume 1
   * Ra for water = 1.4286 e6
   * Ra for air = 9.95 e3
   * and value is linear for flow rate (and coolant speed)
   *
   * Prandtl number https://en.wikipedia.org/wiki/Prandtl_number
   * 	Pr = nu / a
   * 		where nu = mu / rho = viscosity / density
   * 		value for air ~ 0.71
   * 		value for water ~ 7.56
   *
   * h for forced external convection is given by charts
   * https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html
   * 		air 10 - 100,
   * 			h = 12.12 - 1.16*v + 11.6 * sqrt(v); with units of W m-2 K-1
   * 		water 50 - 10000
   */

  if (flowrate == 0) //!< if the coolsystem is off (v == 0), h must be 0 as well
    return 0;        //!< #TODO -> double should not be compared to zero.

  //!< speed of thw fluid

  const double v = flowrate / Across;
  const double h = 12.12 - 1.16 * v + 11.6 * std::sqrt(v);

  return h;
  //!< speed 		h (approximate)
  //!< 0			12
  //!< 2			26
  //!< 5			32
  //!< 10			37
  //!< 20			41
}

double CoolSystem::dstate(double Qtotal, double Qchildren, double t)
{
  /*
   * Calculate the new coolant temperature resulting from a heat exchange of Etot [J] over t [s]
   * Note that the new temperature is not set for consistency (see module)
   *
   * IN
   * Qtot 		total heat exchange, with child SUs, neighbour SUs and parent module [J]
   * Qchildren 	heat exchange with child SUs only [J]
   * t			time over which this heat has been exchanged [s]
   *
   * Notes:
   * Even if flowrate == 0, there is still thermal stuff going on:
   * 		conductive cooling to the first and last cell in a stack
   * 		this module might get cooling from the module above it
   * 		therefore, we cannot say the volume of coolant which should heat up is flowrate * time
   * 		so we must assume there is some fixed volume of coolant for this module which is always at a uniform temperature
   * 		 	this volume is pumped around to cool the cells, but the coolant passing the cells has the same T as the rest of the coolant of this cell
   * 		 	i.e. there is perfect mixing in the reservoir of coolant
   * the volume needs to be quite large for numerical stability
   * 		we only resolve the thermal model ever ~ 20 seconds, at which point it needs to cool a decent amount of energy
   * 		if V is too small, this will result in numerical instability (T increases too much, next step it decreases even more, etc)
   * 		V is the inertia of the system, so the larger the inertia, the more stable the system and the larger time steps we can take
   * 		Note that this is not necessarily unrealistic, in a real battery the air will be mixed in the container, so the heat capacity of modules is ever larger than assumed here
   */

  //!< calculate the new temperature
  const double Tnew = T() + Qtotal / (fluid_rho * fluid_cp * fluid_V);

  //!< increase the power required to run the cooling system
  //!< power scales like speed ^ 3 and rho ^1 (kinetic energy/mass * flow rate = rho * A * v * v^2 / 2)
  //!< 		see also https://fluidflowinfo.com/fan-performance-and-fan-laws/, and my undergrad course on wind energy [fan = opposite wind turbine]
  const double v = flowrate / Across; //!< speed of the fluid
  const double E = (fluid_rho * Across * v * v * v / 2) / eta;
  //!< this is the power to speed up air from stationary to the required speed with effiency eta [W]
  //!< with given parameters, this results in 661 W (our orange jet from the lab consumes 500 W)

  //!< update the mean cooling power
  coolData.storeCumulativeData(Qchildren, Qtotal, t, E * t); //!< increase the total energy
  return Tnew;
}

void CoolSystem::control(double Thot_local, double Thot_global)
{
  /*
   * Control this cool system by adapting the flow rate of the coolant.
   *
   * depending on the value of control_strategy
   * 		1	always on
   * 		2 	on / off depending on T of the hottest child of the SU (i.e. local hot-spot)
   * 		3	on / off depending on T of the hottest cell ultimately connected to this SU (i.e. global hot-spot)
   * 		4	proportional to T of hottest child SU, with a minimum of 20% flow rate (else it is off)
   * 		5	proportional to T of hottest cell, with a minimum of 20% flow rate (else it is off)
   *
   * IN
   * Thot_local 	the local hottest temperature (i.e. T of the hottest child-module) [K]
   * Thot_global 	the global hot spot of the module or battery [K]
   */

  //!< variables
  double Terr;

  switch (control_strategy) {
  case 1: //!< always on
    break;
  case 2: //!< on/off depending on child T

    //!< determine whether the cooling system is on or off
    if (Thot_local >= control_onoff_Ton)
      flowrate = control_onoff_flr;
    if (Thot_local <= control_onoff_Toff)
      flowrate = 0;
    break;

  case 3: //!< on/off depending on global hot spot temperature

    //!< determine whether the cooling system is on or off
    if (Thot_global >= control_onoff_Ton)
      flowrate = control_onoff_flr;
    if (Thot_global <= control_onoff_Toff)
      flowrate = 0;
    break;

  case 4: //!< proportional to T of hottest child SU
    //!< proportionally increase flow rate
    Terr = Thot_local - control_prop_T;
    Terr = std::max(Terr, 0.0);         //!< ensure we never go negative
    if (control_prop_gain * Terr > 0.2) //!< avoid very small flow rates (which would give extremely large dT)
      flowrate = control_onoff_flr * control_prop_gain * Terr;
    else
      flowrate = 0;
    break;

  case 5: //!< proportional to global hot spot temperature
    //!< proportionally increase flow rate
    Terr = Thot_global - control_prop_T;
    Terr = std::max(Terr, 0.0);         //!< ensure we never go negative
    if (control_prop_gain * Terr > 0.2) //!< avoid very small flow rates (which would give extremely large dT)
      flowrate = control_onoff_flr * control_prop_gain * Terr;
    else
      flowrate = 0;
    break;

  default: //!< always on
    break;
  }
}

void CoolSystem::storeData(size_t Ncellsi)
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
      std::cerr << "Error in CoolSystem. This system has " << Ncells
                << " cells connected, but now you are storing data for a system of "
                << Ncellsi << ". Throwing an error.\n";
      throw 10;
    }

    if (Ncells == 0) {
      std::cerr << "Error in CoolSystem. This system has " << Ncells
                << " cells connected, so you cannot store values per cell connected."
                << " Throwing an error.\n";
      throw 10;
    }
  }

  coolData.storeData(*this);
}

void CoolSystem::writeData(const std::string &prefix)
{
  /*
   * Writes data to a csv file.
   * The name of the csv file starts with the value of prefix, after which the identification string of this cell is appended
   *
   * Depending on the value of settings::DATASTORE_CELL, different things are written
   * 	0 	nothing
   * 	1 	histogram
   * 	2 	time data
   *
   */
  coolData.writeData(*this, prefix);
}

} // namespace slide
