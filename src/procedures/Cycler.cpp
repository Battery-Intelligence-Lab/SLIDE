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
void Cycler::initialise(StorageUnit *sui, const std::string &IDi)
{
  su = sui;
  ID = IDi;
}

void Cycler::setDiagnostic(bool newDia)
{
  /*
   * enable or disable diagnostic mode.
   * If diagnostic mode is on, voltage limits are checked at every level of the connected storage unit.
   * For instance, if one cell of a series-string reaches its maximum voltage,
   * this is detected even if the voltage of the entire series string is still within its voltage limits.
   */
  diagnostic = newDia;
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

  //!< then store the data
  su->storeData();
  index++;

  return EXIT_SUCCESS;
}

Status Cycler::rest(double tlim, double dt, int ndt_data, double &Ah, double &Wh)
{
  /*
   * Rest the storage unit for a given amount of time.
   * Voltage limits are completely ignored during this period
   *
   * IN
   * tlim		time [s] for which the current should be applied.
   * dt 		time step [s] to be used for time integration
   * ndt_data integer indicating after how many time steps a data point should be stored
   * 				if <= 0, no data is stored
   * 				else a data point is guaranteed at the beginning and end of this function
   * 					and between those, a data point is added every ndt_data*dt seconds
   *
   * OUT
   * int 		returns why the function stopped
   * 				-5 error when setting the current
   * 				-4 in diagnostic mode, the voltage limit of one of the connected cells has been violated
   * 				-3 We get an error when getting the cell voltage
   * 				-2  lower safety voltage limit encountered
   * 				-1 	upper safety voltage limit encountered
   * 				1 	voltage limit reached
   * 				2 	time limit reached
   * Ah 		charge throughput
   * Wh		energy throughput
   *
   * THROWS
   * 100 		time integration was stopped without the time or voltage limit reached
   */

  const bool boolStoreData = ndt_data > 0;
  //!< store data point
  if (boolStoreData) storeData();
  //!< cout<<"resting and voltage is "<<su->V()<<", current = "<<su->I()<<endl;

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

  constexpr bool prdet = false; //!< print details of what is happening at every time step

  //!< apply current
  while (ttot < tlim) {
    //!<  set current. It makes the code more stable to do this every iteration.
    //!< 	parallel module::setCurrent will reset all currents to 0 and only allow minor variations to equalise the voltage
    //!< 	if you don't cell setCurrent, then timeStep will call redistributeCurrent, which can incrementally increase cell currents
    //!< 	causing them to diverge. (e.g. I1 = 3A and I2 = -3A)
    //	because we take up to 100 steps at once, we cannot tolerate large currents or the voltage limits will be exceeded
    try {
      auto succ = su->setCurrent(0, false, true); //!< do not check voltagel limits, but print error warnings
                                                  //!< std::cout << getStatusMessage(succ) << '\n';
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "Error in Cycler::rest when setting the current to 0: " << e << ", throwing it on" << '\n';
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      throw e;
    }

    //!< change length of the time step in the last iteration to get exactly tlim seconds
    if (nOnce * dti > (tlim - ttot)) //!< we are close to the total time -> set nOnce to 1
      nOnce = 1;

    if (dti > tlim - ttot) //!< the last time step, ensure we end up exactly at the right time
      dti = tlim - ttot;

    //!< take a number of time steps
    try {
      su->timeStep_CC(dti, nOnce);
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in Cycler::rest when advancing in time with nOnce = "
                  << nOnce << ". Passing the error on.\n";
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
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

    //!< check minimum and maximum of nOnce
    nOnce = std::clamp(nOnce, 1, nOnceMax); //!< respect min and maximum

    if constexpr (prdet)
      std::cout << " \t nOnce at end = " << nOnce << '\n';

  } //!< end time integration

  //!< check why we stopped the time integration
  Status succ = Status::Unknown_problem;
  if (ttot >= tlim && std::abs(ttot - tlim) < 1)
    succ = Status::ReachedTimeLimit;
  else {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "Error in Cycler::rest, stopped time integration for unclear reason after "
                << ttot << "s, we were running with time limit " << tlim << '\n';
    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw 100;
  }

  //!< store data point
  if (boolStoreData)
    storeData();

  //!< there is no throughput during resting
  Wh = 0;
  Ah = 0;

  return succ;
}

Status Cycler::setCurrent(double I, double vlim)
{
  /*
   * Function to set the current to the connected storage unit.
   * Individual cell voltage limits are NOT checked against unless we are running in diagnostic mode
   *
   * IN
   * I 	current to be set [A]
   * vlim voltage limit to be respected [V]
   *
   * OUT
   * int
   * 		-5 error when setting the current
   * 		-4 in diagnostic mode, the voltage limit of one of the connected cells has been violated -> current cannot be set
   * 		-3 We get an error when getting the cell voltage after setting the current -> current cannot be set
   * 		-2 lower safety voltage limit is reached when setting the current -> current cannot be set
   * 		-1 upper safety voltage limit is reached when setting the current -> current cannot be set
   * 		0 the current can be set
   * 		1 the current can be set, but in doing so the voltage limit is crossed
   */

  bool checkCellV = diagnostic; //!< check individual voltage limits of cells in the connected SU
                                //!< since violated voltage limits does not cause a fatal error in the code, this is a noncritical error

  const double Vini = su->V(false);
  const double Iini = su->I();

  //!< Try setting the current
  const auto status = su->setCurrent(I, checkCellV, settings::printBool::printNonCrit); //!< throws error if checkCellv == true and voltage out of range of a cell

  if (isStatusVoltageLimitsViolated(status)) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Error caught in Cycler::setCurrent of SU " << su->getFullID() << " when trying to set a current of "
                << I << " A. This current cannot be set without violating the voltage limits of one of the cells" << '\n';
    su->setCurrent(Iini, false, false); //!< reset the initial current (without checking voltage limits again)
    return status;
  } else if (status == Status::ParallelUnit_failed) {
    if constexpr (settings::printBool::printCrit) //!< this is a more critical error than the others since solving the eqn should never fail
      std::cout << "Error 14 caught in Cycler::setCurrent of SU " << su->getFullID() << ". This means a parallel unit "
                << "did not manage to solve the equation to set the current to the desired level of " << I << "." << '\n';
    su->setCurrent(Iini, false, false); //!< reset the initial current (without checking voltage limits again)

    return status;
  }

  //!< check the safety limits of the cell
  const double Vnew = su->V(settings::printBool::printCrit);

  if (Vnew <= 0) {
    if constexpr (settings::printBool::printCrit)
      std::cout << "Error in Cycler::setCurrent when getting the voltage after trying to set a current of "
                << I << " A." << '\n';

    su->setCurrent(Iini, false, false); //!< reset the initial current (without checking voltage limits again)
    return Status::V_not_calculated;
  } else {
    const auto status = free::check_safety(Vnew, *this);

    if (status != Status::SafeVoltage) //!< #TODO Can also put isStatusSafe but this should be faster
    {
      su->setCurrent(Iini, false, false); //!< reset the initial current (without checking voltage limits again)
      return status;

      //!< 	std::cout << "Error in Cycler::setCurrent, after setting a current of " << I << " A, the voltage of "
      //!<   << Vnew << "V is smaller than the minimum safety voltage of the cycler of " << getSafetyVmin()
      //!<   << " V." << '\n';
    }
  }

  //!< check if the voltage limit was exceeded while setting the current
  if (I < 0) { //!< charging -> exceeded if Vini < vlim  & Vnew > vlim
    if (Vini <= vlim && Vnew > vlim)
      return Status::ReachedVoltageLimit;
  } else {
    if (Vini >= vlim && Vnew < vlim) //!< #TODO why Vini check?
      return Status::ReachedVoltageLimit;
  }

  //!< if none of the above, we could set the current without exceeding any voltage limit
  return Status::Success;
}

Status Cycler::CC(double I, double vlim, double tlim, double dt, int ndt_data, double &Ah, double &Wh, double &ttot)
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
   *
   * IN
   * I 		current [A] to be applied
   * 				negative is charge
   * 				positive is discharge
   * vlim 	voltage to which the cell should be (dis)charged.
   * 				no check is done to ensure compatibility with the current.
   * 				i.e. it is possible to set I > 0  (discharge) but Vlim > OCV (voltage limit is larger than the OCV of the cell)
   * tlim		time [s] for which the current should be applied.
   * dt 		time step [s] to be used for time integration
   * ndt_data integer indicating after how many time steps a data point should be stored
   * 				if <= 0, no data is stored
   * 				else a data point is guaranteed at the beginning and end of this function
   * 					and between those, a data point is added every ndt_data*dt seconds
   *
   * OUT
   * int 		returns why the function stopped
   * 				-5 error when setting the current
   * 				-4 in diagnostic mode, the voltage limit of one of the connected cells has been violated
   * 				-3 We get an error when getting the cell voltage
   * 				-2  lower safety voltage limit encountered
   * 				-1 	upper safety voltage limit encountered
   * 				1 	voltage limit reached
   * 				2 	time limit reached
   * Ah 		charge throughput
   * Wh		energy throughput
   *
   * THROWS
   * 100 		time integration was stopped without the time or voltage limit reached
   */

  //!< cout<<"CC starting with voltage limit"<<vlim<<" set current "<<I<<" and initial voltage "<<su->V()<<endl; //
  //!< store data point

  const bool boolStoreData = ndt_data > 0;

  if (boolStoreData)
    storeData();

  //!< set current
  auto succ = setCurrent(I, vlim);
  if (succ != Status::Success)
    return succ; //!< stop if we could not successfully set the current

  //	std::cout<<"CC with Vlim"<<vlim<<" has set the current, voltage is now "<<su->V()<<endl; //

  //!< Variables
  double dti = dt; //!< length of time step i

  ttot = 0;                        //!< total time done
  int idat = 0;                    //!< consecutive number of time steps done without storing data
  int nOnce = 1;                   //!< number of time steps we take at once, will change dynamically
  constexpr double sfactor = 10.0; //!< increase nOnce if the voltage headroom is bigger than sfactor*dV (dV = change in this iteration)
  int nOnceMax = 10;               //!< allow maximum this number of steps to be taken at once
  if (boolStoreData)               //!< if we store data, never take more than the interval at which you want to store the voltage
    nOnceMax = std::min(nOnceMax, ndt_data);
  bool vtot = false;            //!< boolean indicating if vlim has been reached
  double vo = su->V();          //!< voltage in the previous iteration
  double vi{ vo };              //!< voltage in this iteration
  bool allowUp = true;          //!< do we allow nOnce to increase?
  constexpr bool prdet = false; //!< print details of what is happening at every time step
  Wh = 0;
  Ah = 0;
  int i = 0;

  //!< check the voltage limit
  if (((I < 0) && (vi > vlim)) || ((I > 0) && (vi < vlim))) { //!< charging -> exceeded if Vnew > vlim
    return Status::ReachedVoltageLimit;
  }

  //!< apply current
  while (ttot < tlim) {
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
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
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
    ttot += dti * nOnce;
    Ah = I * ttot / 3600.0;
    i += nOnce;
    idat += nOnce;
    Wh += I * dti * nOnce / 3600 * vi;

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

    const auto safetyStatus = free::check_safety(vi, *this);
    //!< #TODO this is pretty much unnecessary.
    // Check if the given limit at the beginning inside voltage limits.
    // Otherwise clamp!

    if (safetyStatus != Status::SafeVoltage)
      return safetyStatus;

    if constexpr (prdet) {
      std::cout << "Cycler with nOnce = " << nOnce
                << ", vlim = " << vlim << ", I = " << I
                << " with vi = " << vi << ", vo = " << vo << ", dV = "
                << std::abs(vi - vo) << " and headroom = " << std::abs(vlim - vi) << "\t";
    }

    //!< adapt the number of time steps we take at once
    if (I < 0) {                           //!< if charging, voltage is increasing -> increase if dV < Vlim-V
      if (vi - vo < (vlim - vi) / sfactor) //!< still far of, increase nOnce
        nOnce++;
      else if (vi - vo >= (vlim - vi)) { //!< very close, immediately reset to 1 to avoid overshoot
        nOnce = 1;
        if constexpr (prdet)
          std::cout << "- charging and getting too close, set to 1 ";
      } else { //!< fairly close, reduce nOnce
        nOnce--;
        if constexpr (prdet)
          std::cout << "- charging and getting close, reduce by 1 ";
      }

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
                      if(prdet)
                              std::cout<<"- highest V is above 4.15V per cell, reducing nOnce by 2 ";

              }
              if(su->getVhigh() >= 4.15){
                      nOnce = 1;										//!< reduce nOnce
                      allowUp = false;								//!< block nOnce from increasing
                      if(prdet)
                              std::cout<<"- highest V is above 4.17V per cell, setting to 1 ";

              }*/
    } else //!< on discharge be careful with steep OCV change at low SoC
    {
      if (vo - vi < (vi - vlim) / sfactor && allowUp) {
        nOnce++;
        if constexpr (prdet)
          std::cout << "- increasing nOnce by 1 ";
      } else if (vo - vi >= (vi - vlim) * 2.0) { //!< add the *2 to reduce faster on discharge
        nOnce = 1;
        if constexpr (prdet)
          std::cout << "- discharging and getting too close, setting nOnce to 1 ";
      } else {
        nOnce -= 2; //!< reduce by 2 to reduce faster
        if constexpr (prdet)
          std::cout << "- decreasing nOnce by 2 cause we are getting close ";
      }

      //!< check the lowest cell voltage, and if it is on the steep part, reduce nOnce aggressively
      if (su->getVlow() < 3.1) { //!< this gets directly to the cell level, so we know for sure the value which matters is 3.1
        nOnce -= 2;              //!< reduce nOnce
        allowUp = false;         //!< block nOnce from increasing
        if constexpr (prdet)
          std::cout << "- below 3.1V per cell, reducing nOnce by 2 ";
      }
      if (su->getVlow() < 3.0) //!< #TODO this are for cells? What about modules?
      {                        //!< this gets directly to the cell level, so we know for sure the value which matters is 3
        nOnce = 1;
        allowUp = false;
        if constexpr (prdet)
          std::cout << "- lowest cell is below 3V, setting nOnce to 1 ";
      }
    }

    //!< check minimum and maximum of nOnce
    nOnce = std::clamp(nOnce, 1, nOnceMax); //!< respect min and max
    vo = vi;                                //!< update the voltage from the previous time step
    if constexpr (prdet)
      std::cout << " \t nOnce at end = " << nOnce << '\n';

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
    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw 100;
  }

  //!< store data point
  if (boolStoreData)
    storeData();

  return succ;
}

Status Cycler::CV(double Vset, double Ilim, double tlim, double dt, int ndt_data, double &Ah, double &Wh, double &ttot)
{
  /*
   * Apply a constant voltage to the connected storage unit until either a current or time limit.
   * The load is stopped when their first limit is reached (or when the safety limit of the cycler is reached).
   *
   * NOTE: this function should only be called when the cell is already very close to the voltage
   * Else setting the current to reach the voltage in the first time step might fail because the current required is too large
   *
   * IN
   * Vset		voltage [V] to be applied
   * Ilim 	absolute value of the current below which the CV phase is stopped
   * tlim		time [s] for which the current should be applied.
   * dt 		time step [s] to be used for time integration
   * ndt_data integer indicating after how many time steps a data point should be stored
   * 				if <= 0, no data is stored
   * 				else a data point is guaranteed at the beginning and end of this function
   * 					and between those, a data point is added every ndt_data*dt seconds
   *
   * OUT
   * int 		returns why the function stopped
   * 				-5 error when setting the current
   * 				-4 in diagnostic mode, the voltage limit of one of the connected cells has been violated
   * 				-3 We get an error when getting the cell voltage in diagnostic mode
   * 				-2  lower safety voltage limit encountered
   * 				-1 	upper safety voltage limit encountered
   * 				0	the current became too small but the voltage tolerance was not satisfied
   * 				1 	current limit reached while the voltage tolerance is achieved
   * 				2 	time limit reached
   *
   * Ah 		charge throughput
   * Wh 		energy throughput
   *
   * THROWS
   * 100 		time integration was stopped without the time or voltage limit reached
   */
  using slide::util::sign; //!< #TODO normally sign is very sensitive function
  const bool boolStoreData = ndt_data > 0;

  //!< Gain of the PI control (similar to Module_p::setCurrent())
  constexpr double f = 0.2; //!< apply 30% of the estimated error based on the OCV on the current

  int sgn = sign(su->I()); //!< sgn of what we are doing, > 0 for discharge, < for discharge
  if (sgn == 0)            //!< if 0 current, use the OCV
    sgn = (su->V() < Vset) ? -1 : 1;

  //!< *************************************************************** INITIALISE *************************************************************************
  //!< store data point
  if (boolStoreData)
    storeData();

  Ah = 0;
  Wh = 0;

  //!< check if we are already at the limit
  double vprev{ su->V() }, Ii{ su->I() }, vi{}; //!< voltage and current in the previous and present time step
  double dV = std::abs(vprev - Vset);
  bool Vlimit = (dV < settings::MODULE_P_V_ABSTOL || dV / Vset < settings::MODULE_P_V_RELTOL);

  if ((std::abs(Ii) < Ilim) && Vlimit)
    return Status::ReachedCurrentLimit;

  //!< set an initial current.
  double sign2 = sign(vprev - Vset);
  if (std::abs(Ii) < 1e-3) { //!< if no current yet, set a small one
    const auto succ = setCurrent(0.1 * sign2, Vset);
    if (succ > Status::Success)
      return succ; //!< stop if we could not successfully set the initial current
  } else {
    auto succ = setCurrent(0.99 * Ii, Vset); //!< if there was a current running, assume this is CCCV and we need to slighly reduce |I| to keep V constant
    while (succ > Status::Success) {         //!< we could not set the current, try starting with a lower current
      Ii = Ii / 2.0;
      succ = setCurrent(0.5 * Ii, Vset);

      if (std::abs(Ii) < 0.1) //!< if the current is too low, give up and return succ
        return succ;
    }
  }
  vprev = su->V();
  Ii = su->I();

  //!< *******************************************************  apply voltage  ****************************************************************************

  const auto nt = static_cast<size_t>(std::ceil(tlim / dt)); //!< number of time steps

  double dti = dt;          //!< length of time step i
  ttot = 0;                 //!< total time done
  bool Itot = false;        //!< boolean indicating if Ilim has been reached
  bool Vtolerance = false;  //!< boolean indicating if the voltage tolerance was reached
  double R = su->getRtot(); //!< approximate DC resistance
  double dI, a;             //!< change in current we will apply to keep the voltage constant
  bool reach;

  for (size_t i = 0; i < nt; i++) {
    //!< change length of the time step in the last iteration to get exactly tlim seconds
    if (i + 1 == nt)
      dti = tlim - dt * (nt - 1);

    //!< take a time step
    try {
      su->timeStep_CC(dti);
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in Cycler::CV of module " << su->getFullID() << " in time step "
                  << i << " when advancing in time. Passing the error on, " << e << '\n';
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      return Status::timeStep_CC_failed;
    }

    //!< Get the voltage
    if (diagnostic) {
      const auto status = su->checkVoltage(vi, true);
      if (!isStatusSuccessful(status)) { //!< in diagnostic mode and the voltage of one of the cells was violated
        if constexpr (settings::printBool::printNonCrit)
          std::cout << "in Cycler::CV, the voltage of a cell became invalid in time step " << i
                    << " when the total voltage was " << vi << '\n';
        return status; //!< previously -3
      }
    } else {
      vi = su->V();
      if (vi <= 0) {
        if constexpr (settings::printBool::printNonCrit)
          std::cout << "error in Cycler::CV in time step " << i
                    << " when getting the voltage. The previous voltage was "
                    << vi << "V. Terminating the CV phase.\n";
        return Status::V_not_calculated;
      }
    }

    //!< increase the throughput
    const auto dAh = su->I() * dti / 3600.0;
    ttot += dti;
    Ah += dAh;
    Wh += dAh * vi;

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

    //!< Adapt the current to keep the voltage constant
    //!< 	1 OCV effect
    //!< 		v_{t+1} - v_t ~ (v_t - v_{t-1}) * abs[ (I + dI)/I  ]
    //		i.e. if dI = 0, then dv(t) == dv(t-1)
    //!< 			 if dI < 0, then the voltage is going to change less |dv(t)| < dv(t-1)|, i.e. increase and decrease less
    //!< 			 if dI > 0, then the voltage is going to change more |dv(t)| > dv(t-1)|, i.e. increase and decrease more
    //!< 		v_{t+1} ~ v_t + (v_t - v_{t-1}) * abs[ (I + dI)/I  ]
    //!< 			assuming abs(dI) < abs(I), you can drop the absolute values
    //!< 		v_{t+1} ~ v_t + (v_t - v_{t-1}) * (I + dI)/I
    //!< 2 R effect
    //!< 		dV = - R*dI
    //!< 3 both together
    //!< 		v_{t+1} ~ v_t + (v_t - v_{t-1}) * (I + dI)/I - R*dI
    //!< 				~ 2*v_t - v_{t-1} + dI [v_t/I - v_{t-1}/I - R]
    //!< 				~ 2*v_t - v_{t-1} + dI * a
    //!< 		v_{t+1} == Vset
    //!< 			so dI = ( Vset - 2*v_t + v_{t-1} ) / a
    //!< 				notice the numerator is like a central difference scheme
    a = (vi - vprev) / Ii - R;
    dI = Vset - 2 * vi + vprev;
    dI = dI / a;
    dI = f * dI; //!< 	we change a fraction f of this number. The smaller f, the longer it takes to converge but the more stable once we reach it
    if (sgn < 0) {
      dI = std::clamp(dI, 0.0, -Ii / 2.0); //!< charging, so dI should be >= 0 to slowly reduce the current at max we half the current (I <0 and dI > 0,
      if (vi < Vset)
        dI = 0; //!< avoid overshoot
    } else {
      dI = std::clamp(dI, -Ii / 3.0, 0.0); //!< discharging, so dI should be <= 0 to slowly reduce the current
      if (i == 0)                          //!< CV discharge at Vmin often does too much in first iteration
        dI = std::max(dI, -Ii / 5.0);      //!< since CC did one step too many, we have exceeded the Vlimit, so dI will be very large
      if (vi > Vset)                       //!< avoid overshoot
        dI = 0.0;
    } //!< ensure the current can never flip sgn

    //!< Attempt to set the current
    reach = false;
    Ii = Ii + dI;
    while (!reach) {

      const auto status = su->setCurrent(Ii, diagnostic, true); //!< set the new current, can throw errors if not valid

      if (isStatusOK(status))
        reach = true;
      else {
        Ii = 0.99 * Ii;            //!< 	so try a smaller current
        if (std::abs(Ii) < 1e-6) { //!< avoid eternal loop with very small currents by setting 0A, which should always work
          Ii = 0;
          su->setCurrent(Ii, false);
          reach = true;
        }
      }
    }

    //!< cout<<"\t CV with Vset = "<<Vset<<" and Ilim = "<<Ilim<<" iteration "<<i<<" has V = "<<vi<<" and I = "<<Ii<<" and dI = "<<dI<<" and sgn "<<sgn <<endl; //

    //!< this is very inefficient for parallel modules. setCurrent takes ages since it equalises voltages, which is done again by CC
    //!< better make a function dI which applies a change in current, without checking cell voltages
    vprev = vi;

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

    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw 100;
  }

  //!< store data point
  if (boolStoreData)
    storeData();

  return succ;
}

Status Cycler::CCCV(double I, double Vset, double Ilim, double dt, int ndt_data, double &Ah, double &Wh, double &dtime)
{
  /*
   * Function to do both a CC and CV (dis)charge after each other.
   *
   * IN
   * I 	absolute value of the current
   * Vset voltage to which the cell should be (dis)charged
   *
   */

  const double tlim = std::numeric_limits<double>::max(); //!< #TODO why tlim is not good?
  //!< check all input parameters are sensible

  return Cycler::CCCV_with_tlim(I, Vset, Ilim, tlim, dt, ndt_data, Ah, Wh, dtime);
}

Status Cycler::CCCV_with_tlim(double I, double Vset, double Ilim, double tlim, double dt, int ndt_data, double &Ah, double &Wh, double &dtime)
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
  if (I < 0) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Cycler::CCCV is called with I = " << I << " while we need a positive value. "
                << "Recovering by taking the absolute value.\n";
    I = std::abs(I);
  }

  if (Ilim < 0) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Cycler::CCCV is called with Ilim = " << Ilim
                << " while we need a positive value. Recovering by taking the absolute value.\n";
    Ilim = std::abs(Ilim);
  }

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

  const auto sgn = (su->V() > Vset) ? 1 : -1; //!< OCV larger than Vset so we need to discharge
  //!< do the CC phase
  double Ah1, Ah2, Wh1, Wh2, dtime1, dtime2;
  auto succ = CC(sgn * I, Vset, tlim, dt, ndt_data, Ah1, Wh1, dtime1);

  if (succ != Status::ReachedVoltageLimit) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Cycler::CCCV could not complete the CC phase, terminated with "
                << getStatusMessage(succ) << ". Trying to do a CV phase.\n";
  }

  //!< do the CV phase
  succ = CV(Vset, Ilim, (tlim - dtime1), dt, ndt_data, Ah2, Wh2, dtime2);
  if (succ != Status::ReachedCurrentLimit) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Cycler::CCCV could not complete the CV phase, terminated with "
                << getStatusMessage(succ) << '\n';
  }

  Ah = Ah1 + Ah2;
  Wh = Wh1 + Wh2;
  dtime = dtime1 + dtime2;
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
  constexpr double tlim = std::numeric_limits<double>::max();
  constexpr double dt = 1; //!< use a 2 second time step for accuracy (probably 5 would be fine as well)
  constexpr double crate = 1.0 / 25.0;
  double Ah1{ 0 }, cap1{ 0 }, Wh{};

  Ah = 0;

  std::vector<double> sini; //!< #TODO if it is avoidable.
  sini.clear();

  su->getStates(sini);

  std::span<const double> sini_span{ sini };
  //!< todo once we introduce a different temperature, set T to Tref

  //!< *********************************************************** 2 full charge / discharge cycle ***********************************************************************

  //!< fully charge and discharge the cell
  double dtime{};
  const auto cap = su->Cap();
  //!< fully charge battery to its specified maximum voltage )
  auto status = CC(-crate * cap, su->Vmax(), tlim, dt, 0, Ah1, Wh, dtime);
  ttot += dtime;

  if (status != Status::ReachedVoltageLimit) {
    if constexpr (settings::printBool::printCrit)
      std::cout << "Error in a subfunction of Cycler::getCapacity "
                << getStatusMessage(status) << ".\n";

    //!< restore the original states
    su->setStates(sini_span, false, true);

    return 0;
  }

  //!< fully discharge the battery
  status = CC(crate * cap, su->Vmin(), tlim, dt, 0, cap1, Wh, dtime);
  ttot += dtime;

  if (status != Status::ReachedVoltageLimit) {
    if constexpr (settings::printBool::printCrit)
      std::cout << "Error in a subfunction of Cycler::getCapacity "
                << getStatusMessage(status) << ".\n";

    //!< restore the original states
    su->setStates(sini_span, false, true);
    return 0;
  }

  Ah = std::abs(Ah1) + std::abs(cap1);
  return std::abs(cap1);
}

Status Cycler::Profile(std::span<double> I_vec, double vlim, double tlim, double dt, int ndt_data, double &Ah, double &Wh)
{
  return Status::Unknown_problem;
}
} // namespace slide
