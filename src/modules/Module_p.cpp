/**
 * @file Module_p.cpp
 * @brief Implementation of parallel module class.
 * @author Jorn Reniers, Volkan Kumtepeli
 * @date 18 Dec 2019
 */

#include "Module_p.hpp"

#include "../settings/settings.hpp"
#include "../utility/utility.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <array>
#include <algorithm>
#include <ctime>

namespace slide {

/**
 * @brief Calculate the total resistance of the module.
 *
 * Computes the total resistance using the formula for resistances in parallel:
 * 1/Rtot = sum(1/R_i)
 *
 * @return The total resistance of the module.
 */
double Module_p::getRtot() // #TODO -> This function seems to be very expensive.
{
  if (SUs.empty()) return 0; // If there are no cells connected, return 0

  // Start from the cell furthest away
  double rtot = Rcontact.back() + SUs.back()->getRtot();

  // Then iteratively come closer, updating resistance  Rcontact[i] + (Rcell[i] \\ Rtot)
  for (int i = static_cast<int>(SUs.size()) - 2; i >= 0; i--) // #TODO Also check if Rcontact.empty() or in future it will be array.
    rtot = Rcontact[i] + (SUs[i]->getRtot() * rtot) / (SUs[i]->getRtot() + rtot);

  return rtot;
}

/**
 * @brief Calculate the voltage of each SU accounting for contact resistance.
 *
 * Computes the terminal voltage Vt for each SU in the module while considering
 * contact resistances and stores the values in the Vall span.
 *
 * @param Vall A span to store the voltage values.
 * @param print Unused parameter (can be removed if not required).
 */
void Module_p::getVall(std::span<double> Vall, bool print) // #TODO span may not be the best container here.
{
  /*
   * Return the voltage of SU[i] as seen from the terminal while accounting for the contact resistance
   *
   * Contact resistances:
   * 		we assume the terminals of the parallel module are on either side (i.e. next to SUs[0] and SUs[N-1])
   * 		the contact resistances are in the 'horizontal' paths, and the current of all subsequent cells goes through it
   * 			i.e. the current through Rc[0] = i[0] + i[1] + ... + i[N-1}
   * 		The series-resistance of every cell is already included in the cells themselves (as Rdc), and in the cell voltage v[i]
   * 		The terminal voltage Vt must be the same for the paths to all cells, i.e.
   * 			Vt = v[0] - R[0] (I[0] + I[1] + ... + I[N-1]) = v[0] - R[0]*I[0:N-1]
   * 			   = V[1] - R[1]*I[1:N-1] - R[0]*I[0:N-1]
   * 			   = v[i] - sum{ R[j]*I[j..N-1], j=0..i }
   * 			   = v[i] - sum{ R[j]*sum(I[k], k=j..N-1), j=0..i }
   */

  if (SUs.empty()) return;

  if (Vall.size() < SUs.size()) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Module::getVall, container is too small!\n";
    throw 10;
  }

  double I_cumulative{ 0 };
  for (size_t i{}; i < SUs.size(); i++) {
    const auto j = SUs.size() - 1 - i; // Inverse indexing.
    Vall[j] = SUs[j]->V();             // Store the voltage of the current SU

    I_cumulative += SUs[j]->I(); // Update cumulative current

    for (auto k{ j }; k < SUs.size(); k++) // Update the voltage values in Vall considering contact resistances
      Vall[k] -= I_cumulative * Rcontact[j];
  }
}

/**
 * @brief Redistribute the current among the SUs to balance the voltage.
 *
 * Iteratively adjusts the current of each SU to balance the voltages across them.
 *
 * @param checkV Unused parameter (can be removed if not required).
 * @param print Unused parameter (can be removed if not required).
 * @return The status of the operation.
 */
Status Module_p::redistributeCurrent(bool checkV, bool print)
{
  // New redistributeCurrent without PI control:
  //!< get cell voltages
  std::array<double, settings::MODULE_NSUs_MAX> Va, Vb, Ia, Ib; //!< #TODO if we should make them vector.

  //!< voltage and initial current of each cell //!< #TODO it is a constant value SU.
  constexpr int maxIteration = 2500;
  const auto nSU = getNSUs();

  auto StatusNow = Status::RedistributeCurrent_failed;

  if (nSU <= 1) return Status::Success;

  double Itot{ 0 };
  getVall(Va, print);
  for (size_t i = 0; i < nSU; i++) {
    Ia[i] = SUs[i]->I();
    Itot += Ia[i]; // We also need to preserve sum of the currents!
  }

  int iter{ 0 };
  for (; iter < maxIteration; iter++) {
    double error{ 0 };

    const auto Vmean = std::accumulate(Va.begin(), Va.begin() + nSU, 0.0) / nSU; // Compute the mean voltage

    for (size_t i = 0; i < nSU; i++) // Compute the error between mean voltage and individual cell voltages
      error += std::abs(Vmean - Va[i]);


    if (error < 1e-10) // Return success if the error is below the threshold
      return Status::Success;

    // Update the currents based on the difference between mean and individual voltages
    for (size_t i = 0; i < nSU; i++) {
      Ia[i] = Ia[i] - (Vmean - Va[i]) * SUs[i]->Cap();
      SUs[i]->setCurrent(Ia[i]);
    }

    getVall(Va, print);
  }

  if constexpr (settings::printNumIterations)
    std::cout << "redistributeCurrent iterations: " << iter << '\n';

  return StatusNow;
}

Status Module_p::setVoltage(double Vnew, bool checkI, bool print)
{
  constexpr int maxIteration = 50;
  const auto nSU = getNSUs();

  auto StatusNow = Status::RedistributeCurrent_failed;

  std::array<double, settings::MODULE_NSUs_MAX> Iolds, Ia, Va;

  //!< get the old currents so we can revert if needed
  for (size_t i{}; i < SUs.size(); i++)
    Ia[i] = Iolds[i] = SUs[i]->I();

  int iter{ 0 };
  for (; iter < maxIteration; iter++) {
    getVall(Va, print);

    double error{ 0 };
    for (size_t i = 0; i < nSU; i++)
      error += std::abs(Vnew - Va[i]);

    if (error < 1e-10) {
      StatusNow = Status::Success;
      break;
    }

    for (size_t i = 0; i < nSU; i++) {
      Ia[i] += -0.2 * (Vnew - Va[i]) / 0.001;
      SUs[i]->setCurrent(Ia[i]);
    }
  }

  if constexpr (settings::printNumIterations)
    std::cout << "setVoltage iterations: " << iter << '\n';

  return StatusNow; // #TODO add some voltage/current etc. check
}

Status Module_p::setCurrent(double Inew, bool checkV, bool print)
{
  /*
   * Set the current of a parallel module
   * This function takes small steps adapting the current of each connected cell until the total current is reached
   * 	 step 1 	change I by a bit  -> measure V_cell -> derive Rtot
   * 	 step 2 	do 50% of Inew-I() by increasing current proportionally to this resistance
   * 	 step 3   	iteratively change I of the cell with the smallest V (charge) or biggest V (discharge)
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
   * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
   * 15 	after setting the current, the voltage of the cells are too far apart
   */

  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  constexpr int maxIteration = 10550;
  const auto nSU = getNSUs();

  auto StatusNow = Status::Success;

  std::array<double, settings::MODULE_NSUs_MAX> Iolds, Ia, Va;

  double Itot{ 0 };

  //!< get the old currents so we can revert if needed
  for (size_t i{}; i < SUs.size(); i++) {
    Ia[i] = Iolds[i] = SUs[i]->I();
    Itot += Iolds[i];
  }

  const auto dI = (Inew - Itot) / nSU; // #TODO must change if charge/discharge.

  for (size_t i{}; i < SUs.size(); i++) {
    Ia[i] += dI;
    StatusNow = SUs[i]->setCurrent(Ia[i]); // #TODO worse status should be here.
    if (!isStatusOK(StatusNow)) return StatusNow;
  }

  getVall(Va, print);

  int iter{ 0 };
  for (; iter < maxIteration; iter++) {
    double Vmean{ 0 }, error{ 0 };

    for (size_t i = 0; i < nSU; i++)
      Vmean += Va[i];

    Vmean /= nSU;

    for (size_t i = 0; i < nSU; i++)
      error += std::abs(Vmean - Va[i]);

    if (error < 1e-9) {
      StatusNow = Status::Success;
      break;
    }

    for (size_t i = 0; i < nSU; i++) {
      Ia[i] = Ia[i] - (Vmean - Va[i]) * SUs[i]->Cap();
      if (!isStatusSuccessful(SUs[i]->setCurrent(Ia[i])))
        std::cout << "This should not happen!!" << std::endl; // #TODO
    }
    getVall(Va, print);
  }

  //!< voltage of cell i is outside the valid range, but within safety limits
  //!< indicate this happened but continue setting states
  if (isStatusWarning(StatusNow)) { // #TODO maybe we should not need to set current equally immediately?
    if (verb)
      std::cout << "warning in Module_p::setCurrent, the voltage of cell with id "
                << SUs[0]->getFullID() << " is outside the allowed range for Inew = " << Inew / getNSUs()
                << ". Continue for now since we are going to redistribute the current to equalise the voltages.\n";

  } else if (isStatusBad(StatusNow)) {
    if (verb)
      std::cout << "ERROR " << getStatusMessage(StatusNow) << " in Module_p::setCurrent when setting the current of cell "
                << " with id " << SUs[0]->getFullID() << " for Inew = " << Inew / nSU
                << ". Try to recover using the iterative version of setCurrent.\n";
    //!< throw error, the catch statement will use the iterative function
    // #TODO give exact id, temporarily set to SUs[0]
  }

  if constexpr (settings::printNumIterations)
    if (iter != 0) std::cout << "setCurrent iterations: " << iter << '\n';
  // #TODO set old currents back here!
  return StatusNow; //!< #TODO problem
}


/**
 * @brief Perform a time step at a constant current.
 *
 * Takes a time step at a constant current by either explicitly solving the system of equations
 * or letting every cell take a CC time step and checking if the voltage equation is satisfied.
 *
 * @param dt The time step.
 * @param nstep Number of steps.
 */
void Module_p::timeStep_CC(double dt, int nstep)
{
  // Check if the time step is valid
  if (dt < 0) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Module_p::timeStep_CC, the time step dt must be 0 or positive, but has value "
                << dt << '\n';
    throw 10;
  }

  //!< we simply take one CC time step on every cell
  auto task_indv = [&](int i) { SUs[i]->timeStep_CC(dt, nstep); };

  try {
    run(task_indv, getNSUs(), (par ? -1 : 1));
  } catch (int e) {
    std::cout << "Error in Module_p::timeStep_CC with module ID " << getFullID()
              << ". error " << e << ", throwing it on.\n";
    throw e;
  }

  //!< **************************************************** Calculate the thermal model once for the nstep * dt time period *****************************************************************

  //!< update the time since the last update of the thermal model
  if (!blockDegAndTherm) {
    therm.time += nstep * dt;

    //!< Increase the heat from the contact resistances
    double Ii = 0; //!< current through resistor I
    for (size_t i = 0; i < SUs.size(); i++) {
      for (size_t j = i; j < SUs.size(); j++) // #TODO very important! Not calculating current for Rcontact properly!!!!!!!!!!
        Ii += SUs[j]->I();                    //!< resistor i sees the currents through the cells 'behind' them

      therm.Qcontact += Rcontact[i] * sqr(Ii) * nstep * dt;
    }

    //!< If this module has a parent module, this parent will call the thermal model with the correct parameters
    //!< which will include heat exchange with the module's neighbours and cooling from the cooling system of the parent module.

    //!< if there is no parent, this module is the top-level.
    //!< It then directly exchanges heat with the environment at a fixed temperature
    //!< If there is no parent, assume we update T every nstep*dt. So update the temperature now
    if (!parent) {
      //!< double check that we have an HVAC coolsystem (see constructor)
      if (typeid(*getCoolSystem()) != typeid(CoolSystem_HVAC)) {
        std::cerr << "ERROR in module_p::timeStep_CC in module " << getFullID() << ". this is a top-level"
                  << "module but does not have an HVAC coolsystem for active cooling with the environment.\n";
        throw 14;
      }

      //!< Call the thermal model without heat exchanges with neighbours or parents (since this module doesn't have any)
      //!< 	Heat exchange between this module and the environment is done by the AC-part of the HVAC coolsystem
      //!< 	heat exchange between this module and its children is done the conventional way by the HVAC coolsystem
      //!< set the new temperature since we have calculated all the temperatures so there is no risk for inconsistency
      double Tneigh[1], Kneigh[1], Aneigh[1];
      //!< Make the arrays even though they will not be used (length should be 0 but I don't think you can make an array of length 0)
      setT(thermalModel(0, Tneigh, Kneigh, Aneigh, therm.time)); //!< the 0 signals there are no neighbours or parents

      /*
      double Tneigh[1] = {settings::T_ENV};							//!< T of environment
      double Kneigh[1] = {cool->getH()}; 					//!< h to environment is same as to children, since the speed of the coolant is the same
      double Aneigh[1] = {therm.A};						//!< A to the environment is the A of this module
      //!< note that this will have more heat exchange than to children, since A = min(A_this, A_parent)
      //!< 	children will have a smaller A themselves, so the resulting A = A_child < A of this module
      setT(thermalModel(1, Tneigh, Kneigh, Aneigh, therm.time));*/
    }

    //!< control the cooling system
    double Tlocal = 0;
    for (auto &SU : SUs)
      Tlocal = std::max(Tlocal, SU->T());

    cool->control(Tlocal, getThotSpot());
  }

  Vmodule_valid = false; //!< we have changed the SOC/concnetration, so the stored voltage is no longer valid

  //!< check if the cell's voltage is valid #TODO I changed this to make redistribute everytime!
  auto status = redistributeCurrent(false, true); //!< don't check the currents
  if (status != Status::Success) {
    auto status = redistributeCurrent(false, true); //!< don't check the currents
    throw 100000;                                   //!< #TODO
  }
}

} // namespace slide