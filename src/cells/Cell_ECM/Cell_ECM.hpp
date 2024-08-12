/**
 * @file Cell_ECM.hpp
 * @brief ECM cell.
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 17 Dec 2019
 */

#pragma once

#include "../../settings/settings.hpp"
#include "../Cell.hpp"
#include "State_ECM.hpp"
#include "../../utility/utility.hpp"

#include <cstring>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <array>
#include <span>

namespace slide {

double ocv_eval(const std::vector<double> &ocv_coefs, double z);

template <size_t N_RC = 1>
class Cell_ECM : public Cell
{
public:
  struct R_C_pair
  {
    double R{ 1 }, C{ 1 };
  };

  struct R_Tau_pair
  {
    double R{ 1 }, Tau{ 1 };
  };


private:
  void createDFheaders();
  std::vector<double> ocv_coefs = { -5147.78793933142, 28772.2793299200, -68881.8548465757 };


protected:
  State_ECM<N_RC> st{ settings::T_ENV, 0.5 }; //!< States T, SOC, , I, Ir, ... ;
  //!< parameters:

  std::array<double, N_RC> Rp{}, Tau{}; // Tau = 1/(RC). All initialised zero.
  XYdata_vv OCV;                        //!< SOC vs voltage curve.
  double Rdc{ 2e-3 };                   //!< DC resistance [Ohm]

public:
  Cell_ECM();
  Cell_ECM(double capin, double SOCin);
  Cell_ECM(double capin, double SOCin, double Rdc_, std::span<R_C_pair> R_C_pair_spn);
  Cell_ECM(std::string IDi, double capin, double SOCin, double Rdc_, std::span<R_C_pair> R_C_pair_spn);
  Cell_ECM(double capin, double SOCin, double Rdc_, std::span<R_Tau_pair> R_Tau_pair_spn);
  Cell_ECM(std::string IDi, double capin, double SOCin, double Rdc_, std::span<R_Tau_pair> R_Tau_pair_spn);
  Cell_ECM(std::string IDi, double capin, double SOCin);
  Cell_ECM(std::string IDi);

  double I() const override { return st.I(); }
  double getIr() const { return st.Ir(); } //!< current through the parallel resistance
  double SOC() override { return st.SOC(); }
  double T() override { return st.T(); }

  double getVr() const { return getIr() * getRp(0); }

  //!< overwrite from Cell
  std::span<double> viewStates() override { return std::span<double>(st.begin(), st.end()); }
  void getStates(getStates_t s) override { s.insert(s.end(), st.begin(), st.end()); }
  auto &getStateObj() { return st; }
  double V() override; //!< crit is an optional argument
  Status setStates(setStates_t s, bool checkV = true, bool print = true) override;

  double getRtot() override { return Rdc; } //!< Return the total resistance, V = OCV - I*Rtot
  double getThotSpot() override { return T(); }
  double getThermalSurface() override { return 0; }               //!< #TODO Not implemented?
  double getOCV() override { return ocv_eval(ocv_coefs, SOC()); } //  { return OCV.interp(st.SOC(), settings::printBool::printCrit); }

  // Linear interpolation #TODO add a OCV model. OCV.interp(st.SOC(), settings::printBool::printCrit);
  double getRp(size_t i) const { return Rp[i]; }
  double getTau(size_t i) const { return Tau[i]; }
  double getC(size_t i) const { return getTau(i) / getRp(i); }

  Status setSOC(double SOCnew, bool checkV = true, bool print = true) override;
  Status setCurrent(double Inew, bool checkV = true, bool print = true) override;
  Status setVoltage(double Vnew, bool checkI = true, bool print = true) override;

  inline void setT(double Tnew) override { st.T() = Tnew; }

  virtual bool validStates(bool print = true) override;
  void timeStep_CC(double dt, int steps = 1) override;

  ThroughputData getThroughputs() override { return { st.time(), st.Ah(), st.Wh() }; }

  void storeData() override;
  void writeData(const std::string &prefix) override;

  Cell_ECM<N_RC> *copy() override { return new Cell_ECM<N_RC>(*this); }

  // --- do not use ---

  void set_ocv_coefs(std::vector<double> ocv_coefs_) { ocv_coefs = ocv_coefs_; }
};

// Implementation:

/**
 * Default constructor for Cell_ECM class template.
 */
template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM()
{
  ID = "Cell_ECM_" + std::to_string(N_RC) + "_";
  capNom = 16;
  OCV.x = { -0.1, 0.0, 1.0, 1.1 }; /// OCV curve, dummy linear curve with 3 points from 2.0V to 4.4V
  OCV.y = { VMIN(), Vmin(), Vmax(), VMAX() };

  if constexpr (N_RC >= 1) {
    constexpr double Cp0 = 38e3; // first parallel capacitance
    Rp[0] = 15.8e-3;             // fist parallel (polarisation) resistance default value.
    Tau[0] = (Rp[0] * Cp0);
  }

  if constexpr (N_RC == 2) {
    Rp[1] = 2.5e-3; // second parallel (polarisation) resistance default value.
    Tau[1] = 100.0;
  }

  OCV.check_is_fixed();

  createDFheaders();

  // cellData.initialise(*this);
}

/**
 * Constructor for Cell_ECM class template with given capacity and state of charge.
 * @param capin Capacity input.
 * @param SOCin State of charge input.
 */
template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM(double capin, double SOCin) : Cell_ECM()
{
  /// check that the input argument is valid
  if (!free::check_SOC(SOCin)) throw 10;

  st.SOC() = SOCin;
  setCapacity(capin);
}

/**
 * Constructor for Cell_ECM class template with given capacity, state of charge, and other parameters.
 * @param capin Capacity input.
 * @param SOCin State of charge input.
 * @param Rdc_ DC resistance.
 * @param Rp_ Array of parallel resistance values.
 * @param inv_tau_ Array of inverse time constant values.
 */
template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM(double capin, double SOCin, double Rdc_,
                                std::span<R_C_pair> spn)
  : Cell_ECM(capin, SOCin)
{
  Rdc = Rdc_;

  size_t i{};
  for (auto [Rp_, Cp_] : spn) {
    Rp[i] = Rp_;
    Tau[i] = Rp_ * Cp_;
    i++;
  }
}

template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM(double capin, double SOCin, double Rdc_,
                                std::span<R_Tau_pair> spn)
  : Cell_ECM(capin, SOCin)
{
  Rdc = Rdc_;

  size_t i{};
  for (auto [Rp_, Tau_] : spn) {
    Rp[i] = Rp_;
    Tau[i] = Tau_;
    i++;
  }
}


template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM(std::string IDi, double capin, double SOCin, double Rdc_,
                                std::span<R_C_pair> spn)
  : Cell_ECM(capin, SOCin, Rdc_, spn)
{
  ID = std::move(IDi);
}

template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM(std::string IDi, double capin, double SOCin, double Rdc_,
                                std::span<R_Tau_pair> spn)
  : Cell_ECM(capin, SOCin, Rdc_, spn)
{
  ID = std::move(IDi);
}

template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM(std::string IDi, double capin, double SOCin)
  : Cell_ECM(capin, SOCin)
{
  ID = std::move(IDi);
}

template <size_t N_RC>
inline Cell_ECM<N_RC>::Cell_ECM(std::string IDi)
  : Cell_ECM()
{
  ID = std::move(IDi);
}


/**
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
template <size_t N_RC>
inline Status Cell_ECM<N_RC>::setCurrent(double Inew, bool checkV, bool print)
{
  const double Iold = I();
  st.I() = Inew;

  const auto status = checkCurrent(checkV, print);

  if (isStatusBad(status)) st.I() = Iold;

  return status;
}

/**
 * Sets the voltage of the cell, updates the current accordingly, and checks the current if specified.
 * @param Vnew New voltage value.
 * @param checkI If true, checks the current after setting the voltage (default is true).
 * @param print If true, prints error messages (default is true).
 * @return The status of the operation.
 */
template <size_t N_RC>
inline Status Cell_ECM<N_RC>::setVoltage(double Vnew, bool checkI, bool print)
{
  const double Iold = st.I();
  // #TODO check if V is sensible here.

  const double ocv = getOCV();
  double v_now = ocv - Vnew;

  for (size_t i{}; i < N_RC; i++)
    v_now -= Rp[i] * st.Ir(i);

  const auto Inew = v_now / Rdc;

  st.I() = Inew;

  const auto status = checkCurrent(checkI, print);

  if (isStatusBad(status))
    st.I() = Iold;

  return status;
}

/**
 * Sets the state of charge (SOC) of the cell and checks the voltage if specified.
 * @note This function is mainly used for testing purposes.
 * @param SOCnew New SOC value (must be between 0 and 1).
 * @param checkV If true, the voltage is checked after setting the SOC (default is true).
 *               - If the voltage is outside the safety limits, an error is thrown and the old SOC is restored.
 *               - If the voltage is outside the valid limits, an error is thrown but the new SOC is kept.
 *               - If the voltage is within the allowed range, the function returns the voltage.
 * @param print If true, error messages are printed based on the global printing variable (default is true).
 *              If false, no messages are printed, but errors are still thrown.
 * @return The status of the operation.
 * @throws 10 If SOCnew is illegal (values must be between 0 and 1).
 */
template <size_t N_RC>
inline Status Cell_ECM<N_RC>::setSOC(double SOCnew, bool checkV, bool print) //!< Also not used except test functions.
{
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

/**
 * Calculates the cell voltage based on the current state of the cell.
 * @note Throws an error if the SOC is outside the allowed range.
 * @tparam N_RC The number of RC-elements in the ECM model.
 * @return The calculated cell voltage.
 * @throws 1 If the SOC is outside the allowed range (passed on from linear interpolation).
 */
template <size_t N_RC>
inline double Cell_ECM<N_RC>::V()
{
  const bool verb = settings::printBool::printCrit; //!< print if the (global) verbose-setting is above the threshold
  try {
    const double ocv = getOCV();
    double v_now = ocv - Rdc * st.I();

    for (size_t i{}; i < N_RC; i++)
      v_now -= Rp[i] * st.Ir(i);

    return v_now;
  } catch (...) {
    if (verb)
      std::cerr << "ERROR in Cell_ECM::getV when getting the OCV.\n";
    return 0;
  }
}

template <size_t N_RC>
inline Status Cell_ECM<N_RC>::setStates(setStates_t s, bool checkV, bool print)
{
  /*
   */
  const auto st_old = st; //!< Back-up values.

  std::copy(s.begin(), s.begin() + st.size(), st.begin()); //!< Copy states.
  s = s.last(s.size() - st.size());                        //!< Remove first Nstates elements from span.

  const Status status = free::check_Cell_states(*this, checkV);

  if (isStatusBad(status))
    st = st_old; //!< Restore states here.

  return status;
}

template <size_t N_RC>
inline bool Cell_ECM<N_RC>::validStates(bool print)
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

/**
 * @brief Perform a time step at a constant current.
 */
template <size_t N_RC>
inline void Cell_ECM<N_RC>::timeStep_CC(double dt, int nstep)
{
  assert(dt > 0); // Check if the time step is valid (only in debug mode)

  const auto dth = dt / 3600.0;
  //!< take the specified number of time steps
  for (int t = 0; t < nstep; t++) {
    //!< Using forward Euler time integration.
    const auto dAh = st.I() * dth;
    st.SOC() -= dAh / Cap();

    for (size_t i{}; i < N_RC; i++) // dIr/dt = (I - Ir)/(RC)
      st.Ir(i) += dt * (st.I() - st.Ir(i)) / Tau[i];

    //!< increase the cumulative variables of this cell
    if constexpr (settings::data::storeCumulativeData) {
      st.time() += dt;
      st.Ah() += std::abs(dAh);
      st.Wh() += std::abs(dAh * V());
    }
  }
}


template <size_t N_RC>
void Cell_ECM<N_RC>::createDFheaders()
{
  if (dataStorageLevel == settings::CellDataStorageLevel::storeMostStates) {
    timeData = DataFrame<double>{ "Voltage [V]", "Temperature [K]", "SOC [-]", "Current [A]", "time [s]", "Current throughput [Ah]", "Energy throughput [Wh]" };
  } else {
    throw "dataStorageLevel other than storeMostStates is not implemented yet!\n";
  }
}

template <size_t N_RC>
void Cell_ECM<N_RC>::storeData()
{
  if (dataStorageLevel == settings::CellDataStorageLevel::storeMostStates) {
    timeData.insert(timeData.end(), { V(), T(), SOC(), I(), st.time(), st.Ah(), st.Wh() });
  } else {
    throw "dataStorageLevel other than storeMostStates is not implemented yet!\n";
  }
}

template <size_t N_RC>
void Cell_ECM<N_RC>::writeData(const std::string &prefix)
{

  auto path = PathVar::results / (prefix + "_" + getFullID() + ".slide");
  timeData.to_binary(path);
}


using Cell_Bucket = Cell_ECM<0>;

} // namespace slide