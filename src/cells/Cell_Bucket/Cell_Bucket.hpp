/*
 * Cell_Bucket.hpp
 *
 *  Created on: 22 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "State_Bucket.hpp"
#include "../Cell.hpp"
#include "../../settings/settings.hpp"
#include "../../utility/utility.hpp"

#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <span>

namespace slide {
class Cell_Bucket : public Cell
{
protected:
  State_Bucket st{ 0, 0.5, settings::T_ENV }; //!< I, T, SOC
  XYdata_ff OCV;                              //!< SOC vs voltage curve.
  double Rdc{ 2e-3 };                         //!< DC resistance [Ohm]

public:
  Cell_Bucket();
  Cell_Bucket(std::string IDi, double capin, double SOCin);

  inline double SOC() override { return st.SOC(); }
  inline double I() const override { return st.I(); }
  double V(bool print = true) override;

  void getStates(getStates_t s) override { s.insert(s.end(), st.begin(), st.end()); }         //!< returns the states of the cell collectively.
  std::span<double> viewStates() override { return std::span<double>(st.begin(), st.end()); } //!< returns the individual states.

  auto &getStateObj() { return st; }

  double getOCV(bool print = true) override; //!< crit is an optional argument
  //!< virtual int getNstates() { return S.size() + 1; } //!< +1 for current

  double getRtot() override { return Rdc; } //!< Return the total resistance, V = OCV - I*Rtot

  Status setCurrent(double Inew, bool checkV = true, bool print = true) override;
  Status setSOC(double SOCnew, bool checkV = true, bool print = true) override;
  Status setStates(setStates_t s, bool checkV = true, bool print = true) override;

  //!< thermal model
  inline double T() override { return st.T(); }
  inline double getThotSpot() override { return T(); }
  double getThermalSurface() override { return 0; }; //!< Not implemented?
  inline void setT(double Tnew) override { st.T() = Tnew; }

  bool validStates(bool print = true) override;
  void timeStep_CC(double dt, int steps = 1) override;

  ThroughputData getThroughputs() { return { st.time(), st.Ah(), st.Wh() }; }


  Cell_Bucket *copy() override { return new Cell_Bucket(*this); }

  //!< dataStorage
  //!< virtual void storeData();
  //!< virtual void writeData(std::string prefix){}; //!< #TODO implement.
};

inline Cell_Bucket::Cell_Bucket()
{
  ID = "Cell_Bucket";
  //!< OCV curve, dummy linear curve with 5 points from 2.0V to 4.4V
  OCV.x = slide::linspace_fix(0.0, 1.0, 3);
  OCV.y = slide::linspace_fix(VMIN(), VMAX(), 3);

  OCV.check_is_fixed();

  cellData.initialise(*this);
}

inline Cell_Bucket::Cell_Bucket(std::string IDi, double capin, double SOCin) : Cell_Bucket()
{
  if (!free::check_SOC(SOCin))
    throw 10; //!< #TODO we need error codes.

  //!< #TODO also check capacity if negative? Use bool instead of throwing?
  ID = std::move(IDi);
  setCapacity(capin);
  st.SOC() = SOCin;
}

inline double Cell_Bucket::getOCV(bool print)
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

inline Status Cell_Bucket::setCurrent(double Inew, bool checkV, bool print)
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

  const auto status = checkCurrent(checkV, print); //!< #TODO this pattern is repeated in all cells.

  if (isStatusBad(status))
    st.I() = Iold;

  return status;
}

inline double Cell_Bucket::V(bool print)
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
  try {
    const double ocv = getOCV(print);
    return ocv - Rdc * I();
  } catch (int e) {
    if (verb)
      std::cerr << "ERROR in Cell_Bucket::getV when getting the OCV.\n";
    return 0;
  }
}

inline Status Cell_Bucket::setSOC(double SOCnew, bool checkV, bool print) //!< Also not used except test functions.
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

inline Status Cell_Bucket::setStates(setStates_t s, bool checkV, bool print)
{
  /*
   */
  auto st_old = st; //!< Back-up values.
  std::cout << "Breakpoint-setStates-1" << std::endl;


  std::copy(s.begin(), s.begin() + st.size(), st.begin()); //!< Copy states.
  std::cout << "Breakpoint-setStates-2" << std::endl;

  s = s.last(s.size() - st.size()); //!< Remove first Nstates elements from span.

  std::cout << "Breakpoint-setStates-3" << std::endl;


  const Status status = free::check_Cell_states(*this, checkV);

  std::cout << "Breakpoint-setStates-4" << std::endl;

  if (isStatusBad(status))
    st = st_old; //!< Restore states here.

  std::cout << "Breakpoint-setStates-5" << std::endl;

  return status;
}

inline bool Cell_Bucket::validStates(bool print)
{
  /*
   * note: does NOT check the voltage, only whether all fields are in the allowed range
   */

  const bool verb = print && settings::printBool::printCrit; //!< print if the (global) verbose-setting is above the threshold

  //!< Check if all fields are present & extract their values
  //!< are all in the allowed range? #TODO change to some error codes.

  bool range = free::check_SOC(SOC());

  if (T() < Tmin() || T() > Tmax()) {
    if (verb)
      std::cerr << "ERROR in Cell_Bucket::validState, T is outside of the range, " << Tmin()
                << " <= T <= " << Tmax() << ", value is " << T() << '\n';
    range = false;
  }
  //!< there is no range on the current
  return range;
}

inline void Cell_Bucket::timeStep_CC(double dt, int nstep)
{
  /*
   *	take a time step of dt seconds while keeping the current constant
   */
  if (dt < 0) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Cell_Bucket::timeStep_CC, the time step dt must be "
                << "0 or positive, but has value " << dt << '\n';
    throw 10;
  }

  const auto dth = dt / 3600.0;

  //!< take the specified number of time steps
  for (int t = 0; t < nstep; t++) {
    //!< Using forward Euler time integration.
    //!< Currently, only the SOC of the cell will change

    const auto dAh = st.I() * dth;

    st.SOC() -= dAh / Cap();

    if constexpr (settings::data::storeCumulativeData) {
      st.time() += dt;
      st.Ah() += std::abs(dAh);
      st.Wh() += std::abs(dAh * V());
    }
  }
}
} // namespace slide
