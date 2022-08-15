/*
 * free_functions.cpp
 *
 * Free functions to help to write everything shorter.
 *
 *  Created on: 05 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../settings/settings.hpp"
#include "../types/Status.hpp"

namespace slide::free {
template <typename T>
size_t getNcells(T const &SU)
{
  /*  return the number of cells connected to this module
   * 	e.g. if this module has 3 child-modules, each with 2 cells.
   * 	then getNSUs = 3 but getNcells = 6
   */
  return SU->getNcells();
}

template <typename T>
auto getV(T const &SU)
{
  return SU.V();
}

template <typename T>
auto getVmin(T const &SU)
{
  //!< return SU.V();
}

template <bool Print = settings::printBool::printCrit>
auto inline check_SOC(double SOCnew, double SOC_min = 0, double SOC_max = 1)
{
  if (SOCnew < SOC_min || SOCnew > SOC_max) //!< check that the input argument is valid
  {
    if constexpr (Print)
      std::cerr << "ERROR in some Cell, illegal input value of SOC: "
                << SOCnew << " while the minimum is " << SOC_min
                << " and maximum is " << SOC_max << ".\n";
    return false;
  }
  return true;
}

template <bool Print = settings::printBool::printCrit>
auto inline check_Cell_states(auto &su, bool checkV)
{
  if (!su.validStates(Print)) //!< check if valid state
  {
    if constexpr (Print)
      std::cerr << "ERROR in " << su.getID() << "::setStates, illegal State.\n"; //!< #CHECK here add some type id.

    return Status::Invalid_states;
  }

  //!< check the voltage if desired
  if (checkV) {
    double v;
    return su.checkVoltage(v, Print); //!< get the voltage Does not throw anymore!
  }

  return Status::Success;
}

template <bool Print = true>
auto inline check_voltage(double &v, auto &su) //!< Check voltage.
{
  /*
   * -2 	V < VMIN				discharging and outside safety range, cycling should be halted asap to avoid numerical errors
   * -1	VMIN <= V < Vmin		discharging and outside range of cell, but not yet safety limit
   * 0 	Vmin <= V <= Vmax 		valid range
   * 1 	Vmax < V <= VMAX		charging and ...
   * 2 	VMAX < V 				charging and ...
   *
   * note: allow a small margin on Vmin and Vmax for when we are doing a CV phase
   * then occasionally, V can be slightly above or below Vlimit
   */

  constexpr bool printCrit = (settings::printBool::printCrit) && Print;       //!< print if the (global) verbose-setting is above the threshold
  constexpr bool printNonCrit = (settings::printBool::printNonCrit) && Print; //!< print if the (global) verbose-setting is above the threshold

  try {
    v = su.V(Print);
  } catch (int err) {
    std::cout << "We could not calculate voltage!!!\n";
    return Status::V_not_calculated;
  }

  if (v <= 0) {
    if (printCrit)
      std::cout << "Error in Cell::setStates when getting the voltage, which is "
                << v << " restoring the old states.\n";

    return Status::V_not_calculated;
  }

  //!< lower limits
  if (v > (su.Vmin() - settings::MODULE_P_V_ABSTOL) && v < (su.Vmax() + settings::MODULE_P_V_ABSTOL)) {
    return Status::Success; //!< Just to put this here, probably first if will be executed most of the time.
  } else if (v < (su.Vmin() - settings::MODULE_P_V_ABSTOL) && su.isDischarging()) {
    if (printNonCrit)
      std::cout << "The voltage of cell " << su.getFullID() << " is " << v
                << "V which is below its minimum voltage of "
                << su.Vmin() << ", cell temperature = " << K_to_Celsius(su.T())
                << " centigrade and I = " << su.I() << '\n';

    return Status::Vmin_violation;
  } else if (v > (su.Vmax() + settings::MODULE_P_V_ABSTOL) && su.isCharging()) //!< #CHECK
  {
    if (printNonCrit)
      std::cout << "The voltage of cell " << su.getFullID() << " is " << v
                << "V which is above its maximum voltage of "
                << su.Vmax() << ", cell temperature = " << K_to_Celsius(su.T())
                << " centigrade and I = " << su.I() << '\n';

    return Status::Vmax_violation;
  } else if (v < su.VMIN() && su.isDischarging()) //#Check we dont actually look at the is discharging
  {
    if (printCrit)
      std::cout << "The voltage of cell " << su.getFullID() << " is " << v
                << "V which is below its safety limit of "
                << su.VMIN() << ", cell temperature = " << K_to_Celsius(su.T())
                << " centigrade and I = " << su.I() << '\n';

    return Status::VMIN_violation;
  } else if (v > su.VMAX() && su.isCharging()) //!< upper limits
  {
    if (printCrit)
      std::cout << "The voltage of cell " << su.getFullID() << " is " << v
                << "V which is above its safety limit of "
                << su.VMAX() << ", cell temperature = " << K_to_Celsius(su.T())
                << " centigrade and I = " << su.I() << '\n';

    return Status::VMAX_violation;
  }

  return Status::Unknown_problem; //!< We don't know what happened...
}

template <bool Print = settings::printBool::printCrit>
auto inline check_safety(double vi, auto &cyc)
{
  const auto SafetyVmin = cyc.getSafetyVmin(); //!< #CHECK this requires some calculations!
  const auto SafetyVmax = cyc.getSafetyVmax();

  if (vi < SafetyVmin) //!< #CHECK this requires some calculations!
  {
    if constexpr (Print)
      std::cout << "Error in Cycler::??, the voltage of " << vi
                << " V is smaller than the minimum safety voltage of the cycler of "
                << SafetyVmin << " V." << '\n';

    return Status::VMINsafety_violation;
  } else if (vi > SafetyVmax) //!< #CHECK this requires some calculations!
  {
    if constexpr (Print)
      std::cout << "Error in Cycler::??, the voltage of " << vi
                << " V is larger than the maximum safety voltage of the cycler of "
                << SafetyVmax << " V." << '\n';
    return Status::VMAXsafety_violation;
  }

  return Status::SafeVoltage;
}

template <bool Print = true>
auto inline check_current(bool checkV, auto &su) //!< Check voltage.
{

  //!< TBC
}

} // namespace slide::free