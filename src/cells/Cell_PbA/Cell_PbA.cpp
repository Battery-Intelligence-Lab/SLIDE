/*
 * Cell.cpp
 *
 *  Created on: 22 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Cell_PbA.hpp"
#include "../../settings/settings.hpp"
#include "../../utility/utility.hpp"

#include <string>
#include <cstring>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>

namespace slide {
Cell_PbA::Cell_PbA()
{
  ID = "Cell_PbA";
  capNom = 20; //!< [Ah] Nominal capacity (data sheet)

  k_OCVp.setCurve(PathVar::data / "Cell_PbA/Lander_corrosion_speed.csv"); // #TODO loads data in constructor.
}

Cell_PbA::Cell_PbA(std::string IDi, double capin, double SOCin) : Cell_PbA()
{
  if (!free::check_SOC(SOCin))
    throw 10;

  //!< #TODO also check capacity if negative? Use bool instead of throwing?
  ID = std::move(IDi);
  setCapacity(capin);
  st.SOC() = SOCin;
}

double Cell_PbA::getOCV()
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
  return OCV.interp(st.SOC(), settings::printBool::printCrit);
}

double Cell_PbA::OCVp()
{
  return OCV.interp(st.SOC(), settings::printBool::printCrit);
}

Status Cell_PbA::setCurrent(double Inew, bool checkV, bool print)
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

  Inew = -Inew; //!< For sign convention.

  const double Iold = st.I();
  st.I() = Inew;

  const auto status = checkCurrent(checkV, print);

  if (isStatusBad(status))
    st.I() = Iold;

  return status;
}

double Cell_PbA::V()
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
  const bool verb = settings::printBool::printCrit; //!< print if the (global) verbose-setting is above the threshold
  try {
    if (isCharging())
      return getOCV() - g_OCV * DOD() + (rho_c / Cap()) * st.I() * (1 + M_c * SOC() / (C_c - st.SOC()));
    else
      return getOCV() - g_OCV * DOD() + (rho_d / Cap()) * st.I() * (1 + M_d * DOD() / (C_d - DOD())); //!< #TODO if are making any problems.
  } catch (int e) {
    if (verb)
      std::cerr << "ERROR in Cell_PbA::getV when getting the OCV.\n";
    return 0;
  }
}

Status Cell_PbA::setSOC(double SOCnew, bool checkV, bool print) //!< Also not used except test functions.
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

Status Cell_PbA::setStates(setStates_t s, bool checkV, bool print)
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

bool Cell_PbA::validStates(bool print)
{
  /*
   * note: does NOT check the voltage, only whether all fields are in the allowed range
   */

  const bool verb = print && settings::printBool::printCrit; //!< print if the (global) verbose-setting is above the threshold

  //!< Check if all fields are present & extract their values

  bool range = free::check_SOC(SOC()); //!< Are we in allowed range?

  if (T() < Tmin() || T() > Tmax()) {
    if (verb)
      std::cerr << "ERROR in Cell_PbA::validState, T is outside of the range, " << Tmin()
                << " <= T <= " << Tmax() << ", value is " << T() << '\n';
    range = false;
  }
  //!< there is no range on the current
  return range;
}

double Cell_PbA::C_deg()
{
  return C_EOL * std::exp(-c_Z * (1 - st.Zw() / (1.6 * Z_IEC)));
}

void Cell_PbA::timeStep_CC(double dt, int nstep)
{
  /*
   *	take a time step of dt seconds while keeping the current constant
   */
  if (dt < 0) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Cell_PbA::timeStep_CC, the time step dt must be "
                << "0 or positive, but has value " << dt << '\n';
    throw 10;
  }

  //!< take the specified number of time steps
  for (int t = 0; t < nstep; t++) {
    //!< Using forward Euler time integration.

    st.SOC() += (st.I() - I_gas()) * dt / (3600.0 * Cap());

    if (Ucorr() < 1.74) {
      const double ks = k_s();
      const double x = std::pow(st.Delta_W() / ks, 1.0 / 0.6) + dt;
      st.Delta_W() = ks * std::pow(x, 0.6);
    } else
      st.Delta_W() += k_s() * dt;

    if (isDischarging()) {
      st.Znom() += std::abs(st.I()) * dt / (3600.0 * Cap());
      st.Zw() += std::abs(st.I()) * f_SOC() * f_acid() * dt / (3600.0 * Cap());
    }

    if (isFullyCharged()) //!< fully charged:
    {
      st.SOC_min() = 1;    //!< Reset minimum SOC.
      st.Delta_tSOC() = 0; //!< Reset elapsed time.
      st.n_bad() = 0;      //!< 0.9999>SOC reset.
    } else {
      st.Delta_tSOC() += dt;
      st.SOC_min() = std::min(st.SOC_min(), st.SOC());
    }

    st.f_stratification() += (f_plus() - f_minus()) * dt;
    st.f_stratification() = std::max(st.f_stratification(), 0.0); //!< Cannot be less than 0.

    if constexpr (settings::data::storeCumulativeData) {
      st.time() += dt;
      st.Ah() += std::abs(dAh);
      st.Wh() += std::abs(dAh * V());
    }
  }
}

double Cell_PbA::I_gas()
{
  return (Cap() / 100) * Igas_0 * std::exp(c_u * (V() - Ugas_0) + c_T * (T() - Tgas_0));
}

double Cell_PbA::k_s()
{

  return k() * std::exp(ks_T * (st.T() - Tcorr_0));
}

double Cell_PbA::k()
{ // Depends on OCVp

  return k_OCVp.interp(OCVp());
}

double Cell_PbA::f_SOC()
{
  return 1 + (cSOC_0 + cSOC_min * (1 - st.SOC_min())) * f_I() * st.Delta_tSOC() / 3600.0;
}

double Cell_PbA::Delta_n()
{
  if (st.SOC() < 0.9)
    return 0;
  else
    return 1 - std::pow((SOC_ref - SOC_max), 2) / 0.0025;

  return;
}
double Cell_PbA::f_acid()
{
  return 1 + f_stratification() * std::sqrt(Iref / (eps + std::abs(st.I())));
}

double Cell_PbA::f_minus_gassing()
{
  return c_minus * std::sqrt(100 / Cap()) *
}

double Cell_PbA::f_minus_diffusion()
{
  return (8 * D_H2SO4 / h_batt * h_batt) * f_stratification() * std::pow(2.0, (st.T() - 20.0_degC) / 10.0);
}

double Cell_PbA::rho_empty()
{
                const auto temp = (rhp.d + rhp.e*rho_nom)*m - rho.f*Cap());
                return rhp.a + std::sqrt(rhp.b + rho.c * rho_nom * (()));
}
double Cell_PbA::Ucorr()
{
                if (isCharging())
                  return Ucorr_0 - SOC_infl * g_OCV * DOD() + 0.5 * (rho_c / Cap()) * st.I() * (1 + M_c * st.SOC() / (C_c - st.SOC()));
                else
                  return Ucorr_0 - SOC_infl * g_OCV * DOD() + 0.5 * (rho_d / Cap()) * st.I() * (1 + M_d * DOD() / (C_d - DOD()));
}

//!< void Cell_PbA::backupStates()
//!< {
//!< 	st_backup.push_back(st);
//!< }

//!< void Cell_PbA::restoreStates()
//!< {
//!< 	st = st_backup.back();
//!< 	st_backup.pop_back();
//!< }

} // namespace slide