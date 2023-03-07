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

#include <fstream>
#include <string>
#include <numeric>
#include <vector>
#include <fstream>
#include <iostream>

namespace slide::free {

template <typename Tsu>
inline Status setVoltage_iterative(Tsu *su, double Vset)
{
  // New redistributeCurrent without PI control:
  //!< get cell voltages
  constexpr int maxIteration = 50;
  double Ia{ su->I() };

  for (int iter{ 0 }; iter < maxIteration; iter++) {
    const auto Va = su->V();
    if (std::abs(Vset - Va) < 1e-6)
      return Status::Success;

    const auto Ib = (Vset > Va) ? Ia - 0.01 * su->Cap() : Ia + 0.01 * su->Cap();
    su->setCurrent(Ib);
    const auto Vb = su->V();

    const double slope = (Ia - Ib) / (Va - Vb);
    Ia = Ib - (Vb - Vset) * slope; //!< False-Position method.
    su->setCurrent(Ia);
  }

  return Status::RedistributeCurrent_failed; // #TODO change to setVoltage failed.
}

inline void write_data(std::ofstream &file, std::vector<double> &data, size_t N = 1)
{
  for (size_t i{}; i < data.size(); i++) {
    if (i % N == 0) {
      if (i != 0)
        file << '\n';
    } else
      file << ',';

    file << data[i];
  }

  data.clear(); //!< reset the index to 0 since we can overwrite the stored data
}

template <typename T>
size_t get_Ncells(T const &SU)
{
  /*  return the number of cells connected to this module
   * 	e.g. if this module has 3 child-modules, each with 2 cells.
   * 	then getNSUs = 3 but getNcells = 6
   */
  return SU->getNcells();
}

template <typename T>
auto get_V(T const &SU) { return SU.V(); }

template <typename T>
auto get_T(T const &SU) { return SU->T(); }


template <typename T>
auto get_Vmin(const T &SU) { return SU->Vmin(); }

template <typename T>
auto get_VMIN(const T &SU) { return SU->VMIN(); }

template <typename T>
auto get_Vmax(const T &SU) { return SU->Vmax(); }

template <typename T>
auto get_VMAX(const T &SU) { return SU->VMAX(); }

template <typename T>
auto get_Cap(const T &SU) { return SU->Cap(); }

template <typename T>
auto get_OCV(const T &SU) { return SU->getOCV(); }

template <typename T>
auto get_I(const T &SU) { return SU->I(); }


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
      std::cerr << "ERROR in " << su.getID() << "::setStates, illegal State.\n"; //!< #TODO here add some type id.

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
    v = su.V();
  } catch (int) {
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
    return Status::Success;                       //!< Just to put this here, probably first if will be executed most of the time.
  } else if (v < su.VMIN() && su.isDischarging()) // #Check we dont actually look at the is discharging
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
  } else if (v < (su.Vmin() - settings::MODULE_P_V_ABSTOL) && su.isDischarging()) {
    if (printNonCrit)
      std::cout << "The voltage of cell " << su.getFullID() << " is " << v
                << "V which is below its minimum voltage of "
                << su.Vmin() << ", cell temperature = " << K_to_Celsius(su.T())
                << " centigrade and I = " << su.I() << '\n';

    return Status::Vmin_violation;
  } else if (v > (su.Vmax() + settings::MODULE_P_V_ABSTOL) && su.isCharging()) //!< #TODO
  {
    if (printNonCrit)
      std::cout << "The voltage of cell " << su.getFullID() << " is " << v
                << "V which is above its maximum voltage of "
                << su.Vmax() << ", cell temperature = " << K_to_Celsius(su.T())
                << " centigrade and I = " << su.I() << '\n';

    return Status::Vmax_violation;
  }

  return Status::Success; //!< #TODO should we also send Vmin/Vmax violations even we are not charging/discharging?
}

template <bool Print = settings::printBool::printCrit>
auto inline check_safety(double vi, auto &cyc)
{
  const auto SafetyVmin = cyc.getSafetyVmin(); //!< #TODO this requires some calculations!
  const auto SafetyVmax = cyc.getSafetyVmax();

  if (vi < SafetyVmin) //!< #TODO this requires some calculations!
  {
    if constexpr (Print)
      std::cout << "Error in Cycler::??, the voltage of " << vi
                << " V is smaller than the minimum safety voltage of the cycler of "
                << SafetyVmin << " V." << '\n';

    return Status::VMINsafety_violation;
  } else if (vi > SafetyVmax) //!< #TODO this requires some calculations!
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

inline std::ofstream openFile(auto &SU, const auto &folder, const std::string &prefix, const std::string &suffix)
{
  const auto name = PathVar::results / (prefix + "_" + SU.getFullID() + "_" + suffix);

  //  std::string name = getName(cell, prefix); //!< name of the file
  std::ofstream file(name, std::ios_base::app); // #TODO app-> initially open then append.

  if (!file.is_open()) {
    std::cerr << "ERROR in Cell::writeData, could not open file "
              << name << '\n';
    throw 11;
  }

  return file;
}

} // namespace slide::free