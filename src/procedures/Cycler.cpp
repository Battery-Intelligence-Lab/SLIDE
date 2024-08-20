/**
 * @file Cycler.cpp
 * @brief Cycler class
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 19 Dec 2019
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

    dti = std::min(dti, tlim - ttot); //!< the last time step, ensure we end up exactly at the right time

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

    nOnce = std::clamp(nOnce, 1, nOnceMax); //!< respect min and maximum

  } //!< end time integration

  succ = su->setCurrent(0, false, true); // #TODO do it one last time otherwise after timeStep_CC currents change!
  if (isStatusBad(succ)) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "Error in Cycler::rest, stopped time integration for unclear reason after "
                << ttot << "s, we were running with time limit " << tlim << '\n'
                << getStatusMessage(succ) << '\n';
    return succ;
  }

  if (boolStoreData) storeData();

  th.time() += ttot; //!< there is no throughput during resting

  return Status::ReachedTimeLimit;
}

/**
 * @brief Sets the current to the connected storage unit
 *
 * @param[in] I Current to be set (in Amperes).
 * @param[in] vlim Voltage limit to be respected (in Volts).
 * @return Status for the reason why the function stopped.
 */
Status Cycler::setCurrent(double I, double vlim, double &v_now)
{
  const bool checkCellV = true;

  const double Vini{ su->V() }, Iini{ su->I() };

  //!< check the voltage limit
  if ((I < 0 && Vini > vlim) || (I > 0 && Vini < vlim)) //!< charging -> exceeded if Vnew > vlim
    return Status::ReachedVoltageLimit;

  //!< Try setting the current
  const auto status = su->setCurrent(I, checkCellV, settings::printBool::printNonCrit); //!< throws error if checkCellv == true and voltage out of range of a cell

  if (isStatusBad(status)) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "Error caught in Cycler::setCurrent of SU " << su->getFullID() << " when trying to set a current of "
                << I << " A. This current cannot be set without violating the voltage limits of one of the cells" << '\n';
    su->setCurrent(Iini, false, false); //!< reset the initial current (without checking voltage limits again)
    return status;
  }

  const double Vnew = su->V();
  //!< check if the voltage limit was exceeded while setting the current
  //!< charging -> exceeded if Vini < vlim  & Vnew > vlim //!< #TODO why Vini check?
  if (Vnew <= 0) {
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "error in Cycler::setCurrent in time step when getting the voltage. "
                << "The previous voltage was " << Vini << "V. Terminating the corresponding phase.\n";
    return Status::V_not_calculated;
  } else if ((I < 0 && Vnew > vlim) || (I > 0 && Vnew < vlim) || isStatusWarning(status))
    return Status::ReachedVoltageLimit;

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
  /* #TODO
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

  double vi{}; //!< voltage now
  auto succNow = setCurrent(I, vlim, vi);
  if (!isStatusSuccessful(succNow)) return succNow; //!< stop if we could not successfully set the current

  auto succ = Status::ReachedTimeLimit;

  //!< Variables
  double dti = dt; //!< length of time step i
  double ttot{};   //!< total time done

  int idat = 0;      //!< consecutive number of time steps done without storing data
  int nOnce = 1;     //!< number of time steps we take at once, will change dynamically
  int nOnceMax = 2;  //!< allow maximum this number of steps to be taken at once
  if (boolStoreData) //!< if we store data, never take more than the interval at which you want to store the voltage
    nOnceMax = std::min(nOnceMax, ndt_data);

  while (ttot < tlim) {

    succNow = setCurrent(I, vlim, vi); // #TODO this was not here I added to get nice results from
    if (!isStatusSuccessful(succNow))
      return succNow; //!< stop if we could not successfully set the current

    dti = std::min(dti, tlim - ttot); //!< the last time step, ensure we end up exactly at the right time

    //!< take a number of time steps
    try {
      su->timeStep_CC(dti, nOnce); // #TODO this should also have some error codes.
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in Cycler::CC when advancing in time with nOnce = "
                  << nOnce << ", V = " << vi << " and Vlim " << vlim << ". The error is "
                  << e << "on.\n";
      return Status::timeStep_CC_failed;
    }

    //!< Increase the throughput
    const auto dt_now = dti * nOnce;
    th.time() += dt_now;
    ttot += dt_now;
    th.Ah() += std::abs(I) * dt_now / 3600.0;
    idat += nOnce;
    th.Wh() += std::abs(I) * dt_now / 3600 * vi; // #TODO This should be (v_before + v_after)/2

    //!< Store a data point if needed
    if (boolStoreData && idat >= ndt_data) {
      storeData();
      idat = 0;
    }

  } //!< end time integration

  if (boolStoreData) storeData();

  return succ;
}

Status Cycler::CV(double Vset, double Ilim, double tlim, double dt, int ndt_data, ThroughputData &th)
{
  // New CV function using setVoltage.
  const bool boolStoreData = ndt_data > 0;

  //!< ************************** INITIALISE **************************
  //!< store data point
  if (boolStoreData) storeData();

  //!< check if we are already at the limit
  double vprev{ su->V() }, Ii{ su->I() }, vi{}; //!< voltage and current in the previous and present time step

  if ((std::abs(Ii) < Ilim)) return Status::ReachedCurrentLimit;

  int idat = 0;    //!< consecutive number of time steps done without storing data
  double dti = dt; //!< length of time step i
  double ttot = 0; //!< total time done

  auto succ{ Status::ReachedTimeLimit }; //!< time limit;

  while (ttot < tlim) {
    auto succNow = su->setVoltage(Vset);
    if (!isStatusOK(succNow)) return succNow; //!< stop if we could not successfully set the current

    dti = std::min(dti, tlim - ttot); //!< change length of the time step in the last iteration to get exactly tlim seconds

    //!< take a time step
    try {
      su->timeStep_CC(dti); // #TODO should return status.
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "error in Cycler::CV of module " << su->getFullID() << " at time "
                  << ttot << " when advancing in time. Passing the error on, " << e << '\n';
      return Status::timeStep_CC_failed;
    }

    //!< increase the throughput
    const auto dAh = std::abs(su->I() * dti / 3600.0);
    vi = su->V();
    ttot += dti;
    idat++;
    th.Ah() += dAh;
    th.Wh() += dAh * vi;

    //!< Store a data point if needed
    if (boolStoreData && idat >= ndt_data) {
      storeData();
      idat = 0;
    }
    //!< check the current limit
    Ii = su->I();
    if (std::abs(Ii) < Ilim) {
      succ = Status::ReachedCurrentLimit; //!< current limit
      break;
    }
  } //!< end time integration

  if (boolStoreData) storeData();

  return succ;
}

Status Cycler::CCCV(double I, double Vset, double Ilim, double dt, int ndt_data, ThroughputData &th)
{
  //!< #TODO check all input parameters are sensible
  return Cycler::CCCV_with_tlim(I, Vset, Ilim, TIME_INF, dt, ndt_data, th);
}

/**
 * @brief Performs both a CC (Constant Current) and CV (Constant Voltage) (dis)charge sequentially.
 *
 * This function first performs a CC (dis)charge and then a CV (dis)charge.
 *
 * @param[in] I Absolute value of the current.
 * @param[in] Vset Voltage to which the cell should be (dis)charged.
 * @param[in] Ilim Absolute value of the current limit for the CV phase.
 * @param[in] tlim Time limit for the CC and CV phases combined.
 * @param[in] dt Time step for time integration.
 * @param[in] ndt_data Integer indicating after how many time steps a data point should be stored.
 * @param[out] th ThroughputData reference to store the time, charge, and energy throughput.
 * @return Status indicating the success or failure of the operation.
 */
Status Cycler::CCCV_with_tlim(double I, double Vset, double Ilim, double tlim, double dt, int ndt_data, ThroughputData &th)
{
  //!< #TODO check all input parameters are sensible
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

  I = (su->V() > Vset) ? std::abs(I) : -std::abs(I); //!< OCV larger than Vset so we need to discharge

  ThroughputData th1{}, th2{};
  auto succ = CC(I, Vset, tlim, dt, ndt_data, th1); //!< do the CC phase

  if constexpr (settings::printBool::printNonCrit)
    if (succ != Status::ReachedVoltageLimit)
      std::cout << "Cycler::CCCV could not complete the CC phase, terminated with "
                << getStatusMessage(succ) << ". Trying to do a CV phase.\n";

  const auto t_remaining = tlim - th1.time();

  // Due to the fixed time step CC does not exactly stop on Vlimit but slightly passes it like 2.68 instead of 2.7 V limit.
  // So if we passed voltage limit set to the voltage limit anyway.
  // But at the same time since in series/parallel module voltages you may also stop due to an individual cell limit instead of the total limit then you probably
  // Stopped before let's say 2.7 like 2.75, this time we set the CV voltage to 2.75 V as individual cell is probably experiencing a lower voltage limit.
  auto V_CV = su->V();
  if ((I > 0 && V_CV < Vset) || (I < 0 && V_CV > Vset))
    V_CV = Vset;


  succ = CV(V_CV, Ilim, t_remaining, dt, ndt_data, th2);
  if constexpr (settings::printBool::printNonCrit)
    if (succ != Status::ReachedCurrentLimit)
      std::cout << "Cycler::CCCV could not complete the CV phase, terminated with "
                << getStatusMessage(succ) << '\n';

  th.time() = th1.time() + th2.time();
  th.Ah() = th1.Ah() + th2.Ah();
  th.Wh() = th1.Wh() + th2.Wh();

  return succ;
}

/**
 * @brief Calculates the charge capacity of the connected Storage Unit (SU).
 *
 * This function performs a full charge and discharge cycle and calculates the charge capacity.
 * It does not affect the Cell; the original states are restored at the end.
 *
 * @param[out] Ah Charge capacity in ampere-hours.
 * @param[out] ttot Total time spent during the full charge and discharge cycle.
 * @return Charge capacity in ampere-hours.
 */
double Cycler::testCapacity(double &Ah, double &ttot)
{
  constexpr double dt = 1; //!< use a 2 second time step for accuracy (probably 5 would be fine as well)
  constexpr double crate = 1.0 / 25.0;
  ThroughputData th1{}, th2{};

  std::vector<double> sini; //!< #TODO if it is avoidable.
  sini.clear();

  su->getStates(sini);
  int n_sini = 0;
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
    su->setStates(sini, n_sini, false, true);

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
    su->setStates(sini, n_sini, false, true); // #TODO we forgot to reset Ah and others maybe.
    return 0;
  }

  Ah = th1.Ah() + th2.Ah();
  return th2.Ah();
}

Status Cycler::Profile(std::span<double> I_vec, double vlim, double tlim, double dt, int ndt_data, double &Ah, double &Wh)
{
  return Status::NotImplementedYet;
}

} // namespace slide
