/*
 * Cell_ECM.hpp
 *
 *  Created on: 17 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../Cell.hpp"
#include "State_ECM.hpp"
#include "../../utility/utility.hpp"
#include "../../settings/settings.hpp"

#include <cstring>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <array>

namespace slide {
class Cell_ECM : public Cell
{

protected:
  State_ECM st{ 0, 0, 0.5, 288 }; //!< States I, Ir, SOC, T;
  //!< parameters:
  double Rp{ 15.8e-3 }, Cp{ 38e3 }; //!< parallel resistance and capacitance
  XYdata_ff OCV;                    //!< SOC vs voltage curve.
  double Rdc{ 2e-3 };               //!< DC resistance [Ohm]

public:
  Cell_ECM();
  Cell_ECM(double capin, double SOCin);

  inline double I() override { return st.I(); }
  inline double getIr() { return st.Ir(); } //!< current through the parallel resistance
  inline double SOC() override { return st.SOC(); }
  inline double T() override { return st.T(); }

  //!< overwrite from Cell
  std::span<double> viewStates() override { return std::span<double>(st.begin(), st.end()); }
  void getStates(getStates_t s) override { s.insert(s.end(), st.begin(), st.end()); }

  auto &getStateObj() { return st; }

  double V(bool print = true) override; //!< crit is an optional argument
  Status setStates(setStates_t s, bool checkStates = true, bool print = true) override;

  double getRtot() override { return Rdc; }          //!< Return the total resistance, V = OCV - I*Rtot
  double getThermalSurface() override { return 0; }; //!< Not implemented?
  double getOCV(bool print = true) override;         //

  Status setSOC(double SOCnew, bool checkV = true, bool print = true) override;
  Status setCurrent(double Inew, bool checkV = true, bool print = true) override;
  inline void setT(double Tnew) override { st.T() = Tnew; }

  virtual bool validStates(bool print = true) override;
  void timeStep_CC(double dt, int steps = 1) override;

  CellThroughputData getThroughputs() { return { st.time(), st.Ah(), st.Wh() }; }


  Cell_ECM *copy() override { return new Cell_ECM(*this); }
};

//!< Implementation:

inline Cell_ECM::Cell_ECM()
{
  ID = "Cell_ECM";
  capNom = 16;
  //!< OCV curve, dummy linear curve with 11 points from 2.0V to 4.4V
  OCV.x = slide::linspace_fix(0.0, 1.0, 11);
  OCV.y = slide::linspace_fix(VMIN(), VMAX(), 11);

  OCV.check_is_fixed();
  cellData.initialise(*this);
}

inline Cell_ECM::Cell_ECM(double capin, double SOCin) : Cell_ECM()
{

  //!< check that the input argument is valid
  if (!free::check_SOC(SOCin))
    throw 10;

  st.SOC() = SOCin;
  setCapacity(capin);
}

inline double Cell_ECM::getOCV(bool print)
{
  /*
   * print 	controls the printing of error messages, (default = true)
   * 			if the SOC is out of the range, an error is always thrown (1)
   * 				if crit is true an error message is printed (if verbose is above the critical level)
   * 				else no error messages are printed
   * 			this is an optional argument. If no value is given, it is assumed to be true
   * 				this can be overwritten by giving an argument of false
   *
   * THROWS
   * 1 	if SOC is outside the allowed range
   * 			passed on from linear interpolation
   */
  return OCV.interp(st.SOC(), (settings::printBool::printCrit && print));
}

inline Status Cell_ECM::setCurrent(double Inew, bool checkV, bool print)
{
  /*
   * sets the current
   *
   * checkV	true, the voltage is checked after setting the current
   * 				if it is outside the safety limits of the cell, error 3 is thrown and the old current is restored
   * 				if it is outside the valid limits of the cell, error 2 is thrown but the new current is kept
   * 				if inside allowed Vrange, it returns the voltage
   * 			false, the voltage is not checked (function returns 0, no errors are thrown)
   *  		if no value of checkV is given, it is set to true
   * print 	controls the printing of error messages
   * 			if true, error messages are printed (if the global printing variable is high enough)
   * 			if false, no messages are printed, but the errors are still thrown
   * 			if no value, the default is true
   *
   * returns the voltage if checkV = true, else it returns 0
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
   * 			and current is in the correct direction. I.e. if charging and V > Vmax or discharging and V < Vmin
   * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
   * 			and current is in the correct direction. I.e. if charging and V > VMAX or discharging and V < VMAX
   * 		if currents are in the 'wrong' direction (e.g. charging but V < Vmin or V < VMIN) then don't throw errors
   * 			since this current is helping to rectify the situation
   */

  const double Iold = I();
  st.I() = Inew;

  const auto status = checkCurrent(checkV, print);

  if (isStatusBad(status))
    st.I() = Iold;

  return status;
}

inline Status Cell_ECM::setSOC(double SOCnew, bool checkV, bool print) //!< Also not used except test functions.
{
  /*
   * checkV	true, the voltage is checked after setting the current
   * 				if it is outside the safety limits of the cell, error 3 is thrown and the old SOC is restored
   * 				if it is outside the valid limits of the cell, error 2 is thrown but the new SOC is kept
   * 				if inside allowed Vrange, it returns the voltage
   * 			false, the voltage is not checked (function returns 0, no errors are thrown)
   *  		if no value of checkV is given, it is set to true
   * print 	controls the printing of error messages
   * 			if true, error messages are printed (if the global printing variable is high enough)
   * 			if false, no messages are printed, but the errors are still thrown
   * 			if no value, the default is true
   *
   * THROWS
   * 10 	illegal SOC, values must be between 0 and 1
   */

  if (!free::check_SOC(SOCnew))
    return Status::SOC_limits_violation;

  const double SOCold = st.SOC();

  st.SOC() = SOCnew;

  if (checkV) {
    double v;
    const auto status = checkVoltage(v, print); //!< get the voltage Does not throw anymore!

    if (isStatusBad(status))
      st.SOC() = SOCold; //!< Restore states here.

    return status;
  }

  return Status::Success;
}

inline double Cell_ECM::V(bool print)
{
  /*
   * print 	controls the printing of error messages, (default = true)
   * 			if the SOC is out of the range, an error is always thrown (1)
   * 				if crit is true an error message is printed (if verbose is above the critical level)
   * 				else no error messages are printed
   * 			this is an optional argument. If no value is given, it is assumed to be true
   * 				this can be overwritten by giving an argument of false
   *
   * THROWS
   * 1 	if SOC is outside the allowed range
   * 			passed on from linear interpolation
   */

  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold
  double ocv;

  try {
    ocv = getOCV(print);
    return ocv - Rp * st.Ir() - Rdc * st.I();
  } catch (int e) {
    if (verb)
      std::cerr << "ERROR in Cell_ECM::getV when getting the OCV.\n";
    return 0;
  }
}

inline Status Cell_ECM::setStates(setStates_t s, bool checkV, bool print)
{
  /*
   */
  auto st_old = st; //!< Back-up values.

  std::copy(s.begin(), s.begin() + st.size(), st.begin()); //!< Copy states.
  s = s.last(s.size() - st.size());                        //!< Remove first Nstates elements from span.

  const Status status = free::check_Cell_states(*this, checkV);

  if (isStatusBad(status))
    st = st_old; //!< Restore states here.

  return status;
}

inline bool Cell_ECM::validStates(bool print)
{
  /*
   * note: does NOT check the voltage, only whether all fields are in the allowed range
   * throws
   * 10 invalid array (wrong length)
   */

  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  //!< check if each value is in the allowed range

  bool range = free::check_SOC(SOC());

  if (T() < Tmin() || T() > Tmax()) {
    if (verb)
      std::cerr << "ERROR in Cell_ECM::validState, T is outside of the range, "
                << Tmin() << " <= T <= " << Tmax()
                << ", value is " << T() << '\n';
    range = false;
  }
  //!< there is no range on the current (Ir or I)

  return range;
}

inline void Cell_ECM::timeStep_CC(double dt, int nstep)
{
  /*
   *	take a time step of dt seconds while keeping the current constant
   */
  if (dt < 0) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Cell_ECM::timeStep_CC, the time step dt must be "
                << "0 or positive, but has value " << dt << '\n';
    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw 10;
  }

  const auto dth = dt / 3600.0;
  //!< take the specified number of time steps
  for (int t = 0; t < nstep; t++) {
    //!< Using forward Euler time integration.
    const auto dAh = st.I() * dth;
    st.SOC() -= dAh / Cap();
    st.Ir() -= dt * (st.Ir() + st.I()) / (Rp * Cp);
    //	dIr/dt =-1/RC Ir - 1/RC I

    //!< increase the cumulative variables of this cell
    if constexpr (settings::data::storeCumulativeData) {
      st.time() += dt;
      st.Ah() += std::abs(dAh);
      st.Wh() += std::abs(dAh * V());
    }
  }
}
} // namespace slide