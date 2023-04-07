/*
 * Cycler.cpp
 *
 *  Created on: 19 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Cycler.hpp"
#include "../settings/settings.hpp"
#include "../utility/utility.hpp"
#include "../cells/cells.hpp"

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <span>
#include <algorithm>

namespace slide {
Cycler &Cycler::initialise(StorageUnit *sui, const std::string &IDi)
{
  su = sui;
  ID = IDi;
  return *this;
}

Cycler &Cycler::setDiagnostic(bool newDia)
{
  /*
   * enable or disable diagnostic mode.
   * If diagnostic mode is on, voltage limits are checked at every level of the connected storage unit.
   * For instance, if one cell of a series-string reaches its maximum voltage,
   * this is detected even if the voltage of the entire series string is still within its voltage limits.
   */
  diagnostic = newDia;
  return *this;
}

int Cycler::writeData()
{
  su->writeData(ID);
  index = 0; //!< reset the index.

  return EXIT_SUCCESS;
}

int Cycler::storeData()
{
  if (index >= settings::CELL_NDATA_MAX) // #TODO this one should not be here if we want to specialise every cell.
    writeData();

  su->storeData();
  index++;

  return EXIT_SUCCESS;
}

/**
 * @brief Rest the storage unit for a given amount of time, ignoring voltage limits.
 *
 * This function allows the storage unit to rest for a specified time duration, tlim. Voltage limits
 * are not considered during this period. Data points are stored according to the specified
 * ndt_data parameter.
 *
 * @param[in] tlim Time (in seconds) for which the storage unit should rest.
 * @param[in] dt Time step (in seconds) to be used for time integration.
 * @param[in] ndt_data Integer indicating after how many time steps a data point should be stored.
 *                     If <= 0, no data is stored. Otherwise, a data point is guaranteed at the
 *                     beginning and end of this function, and between those, a data point is added
 *                     every ndt_data*dt seconds.
 * @param[out] th throughput data.
 * @return Status for the reason why the function stopped.
 * @throws 100 Time integration was stopped without the time or voltage limit being reached.
 */
Status Cycler::rest(double tlim, double dt, int ndt_data, ThroughputData &th)
{
  const bool boolStoreData = ndt_data > 0;
  //!< store data point
  if (boolStoreData) storeData();

  //!< Variables
  double dti = dt;   //!< length of time step i
  double ttot = 0;   //!< total time done
  int idat = 0;      //!< consecutive number of time steps done without storing data
  int nOnce = 1;     //!< number of time steps we take at once, will change dynamically
  int nOnceMax = 10; //!< allow maximum this number of steps to be taken at once
                     //!< careful with thermal stability. Thermal model only calculated every nOnce*dt
                     //!< 	so that can be in the unstable region for large batteries with cooling systems
                     //!< so don't have nOnceMax above 10 (even though once everything is in equilibrium, you could take much larger steps)

  bool ninc = true; //!< can nonce increase this iteration?
  if (boolStoreData)
    nOnceMax = std::min(nOnceMax, ndt_data); //!< if we store data, never take more than the interval at which you want to store the voltage

  Status succ = Status::Unknown_problem;
  //!< apply current
  while (ttot < tlim) {
    //!<  set current. It makes the code more stable to do this every iteration.
    //!< 	parallel module::setCurrent will reset all currents to 0 and only allow minor variations to equalise the voltage
    //!< 	if you don't cell setCurrent, then timeStep will call redistributeCurrent, which can incrementally increase cell currents
    //!< 	causing them to diverge. (e.g. I1 = 3A and I2 = -3A)
    //	because we take up to 100 steps at once, we cannot tolerate large currents or the voltage limits will be exceeded
    succ = su->setCurrent(0, false, true); //!< do not check voltagel limits, but print error warnings
                                           //!< std::cout << getStatusMessage(succ) << '\n';
                                           // #TODO actually we should not do this.


    if (isStatusBad(succ)) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "Error in Cycler::rest when setting the current to 0: "
                  << getStatusMessage(succ) << '\n';

      return succ;
    }

    //!< change length of the time step in the last iteration to get exactly tlim seconds
    if (nOnce * dti > (tlim - ttot)) //!< we are close to the total time -> set nOnce to 1
      nOnce = 1;

    if (dti > tlim - ttot) //!< the last time step, ensure we end up exactly at the right time
      dti = tlim - ttot;

    //!< take a number of time steps
    try { // #TODO timeStep_CC to return Status
      su->timeStep_CC(dti, nOnce);
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in Cycler::rest when advancing in time with nOnce = "
                  << nOnce << ". Passing the error on.\n";
      throw e;
    }

    //!< Increase the throughput
    ttot += dti * nOnce;
    idat += nOnce;

    //!< Store a data point if needed
    if (boolStoreData && idat >= ndt_data) {
      storeData();
      idat = 0;
    }

    //!< adapt the number of time steps we take at once
    if (tlim - ttot < 2 * (nOnce * dti)) //!< getting close to the end
      nOnce = static_cast<int>(std::floor(nOnce / 2));
    else {
      //!< increase every other iteration, else thermal model might freak out
      if (ninc) { //!< e.g. rest after cycle -> cell needs to cool down before you can take very long time steps
        nOnce++;  //!< #TODO Increasing nOnce will cause not saving data regularly. We are probably not saving data every ndata_step
        ninc = false;
      } else
        ninc = true;
    }
    nOnce = std::clamp(nOnce, 1, nOnceMax); //!< respect min and maximum

  } //!< end time integration

  //!< check why we stopped the time integration
  if (ttot >= tlim && std::abs(ttot - tlim) < 1) {
    succ = su->setCurrent(0, false, true); // #TODO do it one last time otherwise after timeStep_CC currents change!
    succ = Status::ReachedTimeLimit;

  } else {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "Error in Cycler::rest, stopped time integration for unclear reason after "
                << ttot << "s, we were running with time limit " << tlim << '\n';
    throw 100;
  }

  if (boolStoreData) storeData();

  th.time() += ttot; //!< there is no throughput during resting

  return succ;
}

/**
 * @brief Sets the current to the connected storage unit, not checking individual cell voltage limits unless in diagnostic mode.
 *
 * @param[in] I Current to be set (in Amperes).
 * @param[in] vlim Voltage limit to be respected (in Volts).
 * @return Status for the reason why the function stopped.
 */
Status Cycler::setCurrent(double I, double vlim)
{
  bool checkCellV = diagnostic; //!< check individual voltage limits of cells in the connected SU
                                //!< since violated voltage limits does not cause a fatal error in the code, this is a noncritical error

  const double Vini{ su->V() }, Iini{ su->I() };

  //!< Try setting the current
  const auto status = su->setCurrent(I, checkCellV, settings::printBool::printNonCrit); //!< throws error if checkCellv == true and voltage out of range of a cell

  if (isStatusVoltageLimitsViolated(status)) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Error caught in Cycler::setCurrent of SU " << su->getFullID() << " when trying to set a current of "
                << I << " A. This current cannot be set without violating the voltage limits of one of the cells" << '\n';
    su->setCurrent(Iini, false, false); //!< reset the initial current (without checking voltage limits again)
    return status;
  }

  const double Vnew = su->V();
  //!< check if the voltage limit was exceeded while setting the current
  if (I < 0 && Vini <= vlim && Vnew > vlim) //!< charging -> exceeded if Vini < vlim  & Vnew > vlim
    return Status::ReachedVoltageLimit;
  else {
    if (Vini >= vlim && Vnew < vlim) return Status::ReachedVoltageLimit; //!< #TODO why Vini check?
  }

  return Status::Success; //!< if none of the above, we could set the current without exceeding any voltage limit
}


/**
 * @brief Apply a constant current to the connected storage unit until either a voltage or time limit is reached.
 * @param I[in] Current [A] to be applied (negative for charging, positive for discharging)
 * @param vlim[in] Voltage to which the cell should be (dis)charged (no check is done to ensure compatibility with the current)
 * @param tlim[in] Time [s] for which the current should be applied
 * @param dt[in] Time step [s] to be used for time integration
 * @param ndt_data[in] Integer indicating after how many time steps a data point should be stored (if <= 0, no data is stored)
 * @param th[out] ThroughputData object to store the charge and energy throughput
 * @return Status indicating the reason for stopping the function
 * @throws 100 if time integration was stopped without the time or voltage limit reached
 */
Status Cycler::CC(double I, double vlim, double tlim, double dt, int ndt_data, ThroughputData &th)
{
  /*
   * Apply a constant current to the connected storage unit until either a voltage or time limit.
   * The load is stopped when their first limit is reached (or when the safety limit of the cycler is reached)
   *
   * NOTE: This function has some sort of adaptive time stepping. The diffusion models are always solved every dt
   * but the thermal and degradation models, as well as module-constraints (V equalisation in parallel) are solved every N*dt seconds.
   * N can go up to 10, which significantly speeds up the computation (5 times faster).
   * N is reduced as you approach the target voltage to ensure you don't overshoot.
   * However, it introduces a problem
   * 		- robustness: if cells have degraded a lot, their effective C-rate will go much higher. So for the same N, their voltage changes more
   * 			this gives problems in the large battery simulation where the code will crash at some point (e.g. 50% degradation)
   * 			because you are still far from Vmax or Vmin, but then after N steps you can suddenly overcharge (negative anode li-fraction) or overdischarge (anode li-fraction > 1)
   */
  const bool boolStoreData = ndt_data > 0;

  if (boolStoreData) storeData();

  auto succ = setCurrent(I, vlim);
  if (!isStatusSuccessful(succ)) return succ; //!< stop if we could not successfully set the current

  //!< Variables
  double dti = dt; //!< length of time step i
  double ttot{};   //!< total time done

  int idat = 0;                    //!< consecutive number of time steps done without storing data
  int nOnce = 1;                   //!< number of time steps we take at once, will change dynamically
  constexpr double sfactor = 10.0; //!< increase nOnce if the voltage headroom is bigger than sfactor*dV (dV = change in this iteration)
  int nOnceMax = 2;                //!< allow maximum this number of steps to be taken at once
  if (boolStoreData)               //!< if we store data, never take more than the interval at which you want to store the voltage
    nOnceMax = std::min(nOnceMax, ndt_data);
  bool vtot = false;              //!< boolean indicating if vlim has been reached
  double vo{ su->V() }, vi{ vo }; //!< voltage in the previous/this iteration
  bool allowUp = true;            //!< do we allow nOnce to increase?

  //!< check the voltage limit
  if ((I < 0 && vi > vlim) || (I > 0 && vi < vlim)) //!< charging -> exceeded if Vnew > vlim
    return Status::ReachedVoltageLimit;

  //!< apply current
  while (ttot < tlim) {
    auto succ = setCurrent(I, vlim); // #TODO this was not here I added to get nice results from
    //!< change length of the time step in the last iteration to get exactly tlim seconds
    if (nOnce * dti > tlim - ttot) //!< we are close to the total time -> set nOnce to 1
      nOnce = 1;

    if (dti > tlim - ttot) //!< the last time step, ensure we end up exactly at the right time
      dti = tlim - ttot;

    //!< take a number of time steps
    try {
      su->timeStep_CC(dti, nOnce);
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in Cycler::CC when advancing in time with nOnce = "
                  << nOnce << ", V = " << vi << " and Vlim " << vlim << ". The error is "
                  << e << "on.\n";
      return Status::timeStep_CC_failed;
    }

    //!< get the voltage
    if (diagnostic) {
      const auto status = su->checkVoltage(vi, true);
      if (!isStatusSuccessful(status)) { //!< in diagnostic mode and the voltage of one of the cells was violated
        if constexpr (settings::printBool::printNonCrit)
          std::cout << "in Cycler::CC, the voltage of a cell became invalid in time step when the total voltage was "
                    << vi << '\n';
        return status; //!< previously -4
      }
    } else {
      const auto v_prev = vi;
      vi = su->V();
      if (vi <= 0) {
        if constexpr (settings::printBool::printNonCrit)
          std::cout << "error in Cycler::CC in time step when getting the voltage. The previous voltage was "
                    << v_prev << "V. Terminating the CC phase.\n";
        return Status::V_not_calculated;
      }
    }

    auto myVnow = su->V();
    //!< Increase the throughput
    const auto dt_now = dti * nOnce;
    th.time() += dt_now;
    ttot += dt_now;
    th.Ah() += std::abs(I) * dt_now / 3600.0;
    idat += nOnce;
    th.Wh() += std::abs(I) * dt_now / 3600 * vi;

    //!< Store a data point if needed
    if (boolStoreData && idat >= ndt_data) {
      storeData();
      idat = 0;
    }

    //!< check the voltage limit
    if (((I < 0) && (vi > vlim)) || ((I > 0) && (vi < vlim))) { //!< charging -> exceeded if Vnew > vlim
      vtot = true;                                              //!< indicate the voltage limit was reached
      break;
    }

    //!< #TODO this is pretty much unnecessary.
    // Check if the given limit at the beginning inside voltage limits.
    // Otherwise clamp!

    //!< adapt the number of time steps we take at once
    if (I < 0) {                           //!< if charging, voltage is increasing -> increase if dV < Vlim-V
      if (vi - vo < (vlim - vi) / sfactor) //!< still far of, increase nOnce
        nOnce++;
      else if (vi - vo >= (vlim - vi)) //!< very close, immediately reset to 1 to avoid overshoot
        nOnce = 1;
      else //!< fairly close, reduce nOnce
        nOnce--;


      //!< check the highest cell voltage, if it is close, reduce nOnce aggressively
      /*
       * Note: the code below makes the simulation works with more variation and more degraded cells
       * e.g. v > 4.1 n-2 and v > 4.15 n = 1
       * 		unit test procedure_test stress test with large variation: simulate to a capacity of 10%
       * 		large battery simulation coolsystem 2: simulate 3592 cycles instead of 3295
       *
       * however, it also significantly slows down the simulation time.
       * 	e.g. same settings for large batter increase simulation time from 5s per cycle to 12.2s per cycle
       *
       * v>4.15 n-2
       * 	large distribution unit test goes to 30% capacity
       * 	large simulation takes 25 s per cycle
       */
      /*	if(su->getVhigh() >= 4.1){
                      nOnce-= 2;										//!< reduce nOnce
                      allowUp = false;								//!< block nOnce from increasing

              }
              if(su->getVhigh() >= 4.15){
                      nOnce = 1;										//!< reduce nOnce
                      allowUp = false;								//!< block nOnce from increasing

              }*/
    } else //!< on discharge be careful with steep OCV change at low SoC
    {
      if (vo - vi < (vi - vlim) / sfactor && allowUp)
        nOnce++;
      else if (vo - vi >= (vi - vlim) * 2.0) //!< add the *2 to reduce faster on discharge
        nOnce = 1;
      else
        nOnce -= 2; //!< reduce by 2 to reduce faster


      //!< check the lowest cell voltage, and if it is on the steep part, reduce nOnce aggressively
      if (su->getVlow() < 3.1) { //!< this gets directly to the cell level, so we know for sure the value which matters is 3.1
        nOnce -= 2;              //!< reduce nOnce
        allowUp = false;         //!< block nOnce from increasing
      }
      if (su->getVlow() < 3.0) //!< #TODO this are for cells? What about modules?
      {                        //!< this gets directly to the cell level, so we know for sure the value which matters is 3
        nOnce = 1;
        allowUp = false;
      }
    }

    //!< check minimum and maximum of nOnce
    nOnce = std::clamp(nOnce, 1, nOnceMax); //!< respect min and max
    vo = vi;                                //!< update the voltage from the previous time step

  } //!< end time integration

  //!< check why we stopped the time integration
  if (ttot >= tlim)
    succ = Status::ReachedTimeLimit;
  else if (vtot)
    succ = Status::ReachedVoltageLimit;
  else {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "Error in Cycler::CC, stopped time integration for unclear reason after "
                << ttot << "s and voltage " << vi << "V, we were running with time limit "
                << tlim << " and voltage limit " << vlim << "V and current " << I << '\n';
    throw 100;
  }

  if (boolStoreData) storeData();

  return succ;
}

Status Cycler::CV(double Vset, double Ilim, double tlim, double dt, int ndt_data, ThroughputData &th)
{
  // New CV function using setVoltage.
  using slide::util::sign; //!< #TODO normally sign is very sensitive function
  const bool boolStoreData = ndt_data > 0;

  //!< *************************************************************** INITIALISE *************************************************************************
  //!< store data point
  if (boolStoreData)
    storeData();

  //!< check if we are already at the limit
  double vprev{ su->V() }, Ii{ su->I() }, vi{}; //!< voltage and current in the previous and present time step
  double dV = std::abs(vprev - Vset);
  bool Vlimit = (dV < settings::MODULE_P_V_ABSTOL || dV / Vset < settings::MODULE_P_V_RELTOL); // #TODO seems unnecessary.

  if ((std::abs(Ii) < Ilim) && Vlimit)
    return Status::ReachedCurrentLimit;


  //!< *******************************************************  apply voltage  ****************************************************************************

  const auto nt = static_cast<size_t>(std::ceil(tlim / dt)); //!< number of time steps

  double dti = dt;         //!< length of time step i
  double ttot = 0;         //!< total time done
  bool Itot = false;       //!< boolean indicating if Ilim has been reached
  bool Vtolerance = false; //!< boolean indicating if the voltage tolerance was reached
  double dI, a;            //!< change in current we will apply to keep the voltage constant
  bool reach;

  for (size_t i = 0; i < nt; i++) {
    su->setVoltage(Vset);

    //!< change length of the time step in the last iteration to get exactly tlim seconds
    if (i + 1 == nt)
      dti = tlim - dt * (nt - 1);

    //!< take a time step
    try {
      su->timeStep_CC(dti); // #TODO should return status.
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in Cycler::CV of module " << su->getFullID() << " in time step "
                  << i << " when advancing in time. Passing the error on, " << e << '\n';
      return Status::timeStep_CC_failed;
    }

    //!< increase the throughput
    const auto dAh = std::abs(su->I() * dti / 3600.0);
    vi = su->V();
    ttot += dti;
    th.Ah() += dAh;
    th.Wh() += dAh * vi;


    //!< Store a data point if needed
    if (boolStoreData && ((i + 1) % ndt_data == 0))
      storeData();
    const auto safetyStatus = free::check_safety(vi, *this);

    if (safetyStatus != Status::SafeVoltage)
      return safetyStatus;

    //!< check the current limit
    Ii = su->I();
    dV = std::abs(vi - Vset);
    Vlimit = (dV < settings::MODULE_P_V_ABSTOL || dV / Vset < settings::MODULE_P_V_RELTOL);
    if (std::abs(Ii) < Ilim && Vlimit) {
      Itot = true; //!< indicate the current limit was reached
      Vtolerance = true;
      break;
    }

    //!< fail-safe to avoid an eternal loop if the current is extremely small but the voltage still is not correct
    //!< in this case, stop even though the voltage tolerance has not been achieved
    if (std::abs(Ii) < Ilim / 10.0 || std::abs(Ii) < 1e-6) {
      Itot = true;        //!< indicate the current limit was reached
      Vtolerance = false; //!< but without reaching the voltage tolerance
      break;
    }
  } //!< end time integration

  //!< *********************************************************** TERMINATE ******************************************************************************
  //!< check why we stopped the time integration
  Status succ;
  if (ttot >= tlim)
    succ = Status::ReachedTimeLimit; //!< time limit
  else if (Itot)                     //!< #TODO why a bool named Itot? It looks like double.
  {
    if (Vtolerance)
      succ = Status::ReachedCurrentLimit; //!< current limit
    else
      succ = Status::ReachedSmallCurrent; //!< #TODO is this ever possible? current became way too small, but voltage limit still not reached
  } else {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "Error in Cycler::CV, stopped time integration for unclear reason after "
                << ttot << "s and voltage " << vi << "V, we were running with time limit " << tlim
                << " and current limit " << Ilim << "V and set voltage " << Vset << '\n';

    throw 100;
  }

  //!< store data point
  if (boolStoreData)
    storeData();

  return succ;
}

Status Cycler::CCCV(double I, double Vset, double Ilim, double dt, int ndt_data, ThroughputData &th)
{
  //!< #TODO check all input parameters are sensible
  return Cycler::CCCV_with_tlim(I, Vset, Ilim, TIME_INF, dt, ndt_data, th);
}

Status Cycler::CCCV_with_tlim(double I, double Vset, double Ilim, double tlim, double dt, int ndt_data, ThroughputData &th)
{
  /*
   * Function to do both a CC and CV (dis)charge after each other.
   *
   * IN
   * I 	absolute value of the current
   * Vset voltage to which the cell should be (dis)charged
   *
   */

  //!< check all input parameters are sensible
  I = std::abs(I);
  Ilim = std::abs(Ilim);

  if (Vset < su->Vmin() || Vset > su->Vmax()) //!< #TODO su Vmax and Vmin require lots of computation
  {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Cycler::CCCV, the set voltage of " << Vset
                << " is outside the allowed range for this StorageUnit, which is "
                << su->Vmin() << " to " << su->Vmax() << ".\n";

    return Status::Invalid_Vset;
  }

  //!< Check if we need to charge or discharge
  const auto status = su->setCurrent(0, false, false);

  if (isStatusFailed(status)) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "error in Cycler::CCCV when getting the OCV, error "
                << getStatusMessage(status) << ".\n";

    return status;
  }

  I = (su->V() > Vset) ? I : -I; //!< OCV larger than Vset so we need to discharge

  ThroughputData th1{}, th2{};
  auto succ = CC(I, Vset, tlim, dt, ndt_data, th1); //!< do the CC phase

  if constexpr (settings::printBool::printNonCrit)
    if (succ != Status::ReachedVoltageLimit)
      std::cout << "Cycler::CCCV could not complete the CC phase, terminated with "
                << getStatusMessage(succ) << ". Trying to do a CV phase.\n";

  const auto t_remaining = tlim - th1.time();
  succ = CV(Vset, Ilim, t_remaining, dt, ndt_data, th2);
  if constexpr (settings::printBool::printNonCrit)
    if (succ != Status::ReachedCurrentLimit)
      std::cout << "Cycler::CCCV could not complete the CV phase, terminated with "
                << getStatusMessage(succ) << '\n';

  th.time() = th1.time() + th2.time();
  th.Ah() = th1.Ah() + th2.Ah();
  th.Wh() = th1.Wh() + th2.Wh();

  return succ;
}

double Cycler::testCapacity(double &Ah, double &ttot)
{
  /*
   * Calculate the charge capacity of the connected SU.
   * note, this function does not affect the Cell, it restores the original states at the end
   *
   * in diagnostic mode, we cannot really do a CV since the small voltage errors during CV will cause some cells to exceed their voltage limit
   * therefore, we have to measure it with a slow (dis)charge
   */
  constexpr double dt = 1; //!< use a 2 second time step for accuracy (probably 5 would be fine as well)
  constexpr double crate = 1.0 / 25.0;
  ThroughputData th1{}, th2{};

  std::vector<double> sini; //!< #TODO if it is avoidable.
  sini.clear();

  su->getStates(sini);

  std::span<double> sini_span{ sini };
  //!< #TODO once we introduce a different temperature, set T to Tref

  //!< *********************************************************** 2 full charge / discharge cycle ***********************************************************************

  //!< fully charge and discharge the cell
  const auto cap = su->Cap();
  //!< fully charge battery to its specified maximum voltage )
  auto status = CC(-crate * cap, su->Vmax(), TIME_INF, dt, 0, th1);
  ttot += th1.time();

  if (status != Status::ReachedVoltageLimit) {
    if constexpr (settings::printBool::printCrit)
      std::cout << "Error in a subfunction of Cycler::getCapacity "
                << getStatusMessage(status) << ".\n";

    //!< restore the original states
    su->setStates(sini_span, false, true);

    return 0;
  }

  //!< fully discharge the battery
  status = CC(crate * cap, su->Vmin(), TIME_INF, dt, 0, th2);
  ttot += th2.time();

  if (status != Status::ReachedVoltageLimit) {
    if constexpr (settings::printBool::printCrit)
      std::cout << "Error in a subfunction of Cycler::getCapacity "
                << getStatusMessage(status) << ".\n";

    //!< restore the original states
    su->setStates(sini_span, false, true);
    return 0;
  }

  Ah = th1.Ah() + th2.Ah();
  return th2.Ah();
}

Status Cycler::Profile(std::span<double> I_vec, double vlim, double tlim, double dt, int ndt_data, double &Ah, double &Wh)
{
  return Status::Unknown_problem;
}
} // namespace slide
