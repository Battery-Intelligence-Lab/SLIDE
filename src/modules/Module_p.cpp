/*
 * Module_p.cpp
 *
 *  Created on: 18 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
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

double Module_p::getRtot() // #TODO -> This function seems to be very expensive.
{
  /*
   * Return the total resistance
   * 		V(I) = OCV - I*Rtot
   * 		with V and OCV the total values for this module
   *
   * parallel: 1/Rtot = sum( 1/R_i )
   */

  //!< If there are no cells connected, return 0
  if (SUs.empty()) return 0;

  //!< check if there are contact resistances only until SUs size. //!< #TODO why are we checking this?
  const bool noRc = std::all_of(Rcontact.begin(), Rcontact.begin() + SUs.size(), util::is_zero<double>);

  if (noRc) { //!< no contact resistance
    double rtot = 0;
    for (auto &SU : SUs)
      rtot += 1.0 / SU->getRtot();

    return 1.0 / rtot;
  }

  //!< with contact resistance
  //!< start from the cell furthest away
  double rtot = Rcontact.back() + SUs.back()->getRtot();

  //!< then iteratively come closer, every time Rcontact[i] + (Rcell[i] \\ Rtot)
  //!< 				= Rc[i] + Rcell[i]*Rtot / (Rcel[i]*Rtot)
  for (int i = SUs.size() - 2; i >= 0; i--) //!< #TODO bug if there are less than 2 SUs.
    rtot = Rcontact[i] + (SUs[i]->getRtot() * rtot) / (SUs[i]->getRtot() + rtot);

  return rtot;
}


void Module_p::getVall(std::span<double> Vall, bool print)
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
    Vall[j] = SUs[j]->V();

    I_cumulative += SUs[j]->I();

    for (auto k{ j }; k < SUs.size(); k++)
      Vall[k] -= I_cumulative * Rcontact[j];
  }
}

Status Module_p::redistributeCurrent_new(bool checkV, bool print)
{
  // New redistributeCurrent without PI control:
  //!< get cell voltages
  std::array<double, settings::MODULE_NSUs_MAX> Va, Vb, Ia, Ib; //!< #TODO if we should make them vector.

  //!< voltage and initial current of each cell //!< #TODO it is a constant value SU.
  constexpr int maxIteration = 250;
  const auto nSU = getNSUs();

  auto StatusNow = Status::RedistributeCurrent_failed;

  if (nSU <= 1) return Status::Success;

  double Itot{ 0 };
  getVall(Va, print);
  for (size_t i = 0; i < nSU; i++) {
    Ia[i] = SUs[i]->I();
    Itot += Ia[i]; // We also need to preserve sum of the currents!
  }

  for (int iter{ 0 }; iter < maxIteration; iter++) {
    double Vmean{ 0 }, error{ 0 };

    for (size_t i = 0; i < nSU; i++)
      Vmean += Va[i];

    Vmean /= nSU;

    for (size_t i = 0; i < nSU; i++)
      error += std::abs(Vmean - Va[i]);

    if (error < 1e-10)
      return Status::Success;

    for (size_t i = 0; i < nSU; i++) {
      Ia[i] = Ia[i] - 0.2 * (Vmean - Va[i]) / 0.001; // SUs[i]->getRtot();
      SUs[i]->setCurrent(Ia[i]);
    }
    getVall(Va, print);
  }

  return StatusNow;
}

Status Module_p::setVoltage(double Vnew, bool checkI, bool print)
{

  // #TODO check if V is sensible here.

  const double Iold = I();

  for (auto &SU : SUs) {
    const auto status = SU->setVoltage(Vnew, checkI, print);

    if (!isStatusOK(status)) {
      setCurrent(Iold, false, false);
      return status;
    }
  }

  return Status::Success; // #TODO status should not be success but worst status given by setVoltages.
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
  double v;

  //!< get the old currents so we can revert if needed
  std::vector<double> Iolds; //!< #TODO we need to remove vector somehow?
  Iolds.clear();

  for (auto &SU : SUs)
    Iolds.push_back(SU->I());

  //!< allocate the current uniformly
  for (size_t i = 0; i < getNSUs(); i++) {
    auto status = SUs[i]->setCurrent(Inew / getNSUs(), checkV, print);

    //!< voltage of cell i is outside the valid range, but within safety limits
    //!< indicate this happened but continue setting states
    if (isStatusWarning(status)) { // #TODO maybe we should not need to set current equally immediately?
      {
        if (verb)
          std::cout << "warning in Module_p::setCurrent, the voltage of cell " << i << " with id "
                    << SUs[i]->getFullID() << " is outside the allowed range for Inew = " << Inew / getNSUs()
                    << ". Continue for now since we are going to redistribute the current to equalise the voltages.\n";
      }
    } else if (isStatusBad(status)) {
      if (verb)
        std::cout << "ERROR " << getStatusMessage(status) << " in Module_p::setCurrent when setting the current of cell "
                  << i << " with id " << SUs[i]->getFullID() << " for Inew = " << Inew / getNSUs()
                  << ". Try to recover using the iterative version of setCurrent.\n";
      //!< throw error, the catch statement will use the iterative function
    }

  } //!< loop

  //!< Redistribute the current to equalise the voltages
  const auto status = redistributeCurrent_new(checkV, print);

  if (!isStatusOK(status)) {
    for (size_t i = 0; i < getNSUs(); i++) //!< revert to the old current
      SUs[i]->setCurrent(Iolds[i], false, true);
    if (verb)
      std::cout << "warning in Module_p::setCurrent, after redistribute, the voltage of one of the cells is "
                << "outside the allowed but inside the safe range for Inew = " << Inew << ". Continue for now.\n"
                << getStatusMessage(status) << '\n';

    return status;
  }

  return Status::Success; //!< #TODO problem
}

bool Module_p::validSUs(SUs_span_t c, bool print)
{
  /*
   * Checks the cells are a valid combination for a parallel-connected module
   * the voltage differences are within the tolerance
   *
   * If the number of cells is the same as in this module, use the contact resistances
   */
  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  //!< Check the voltage of each cell is valid and within the error tolerance #TODO it is better to supply both module and Rcontact.
  std::array<double, settings::MODULE_NSUs_MAX> Vall;
  getVall(Vall, print);

  auto [V_min_it, V_max_it] = std::minmax_element(Vall.begin(), Vall.begin() + SUs.size());

  const double dV = *V_max_it - *V_min_it; //!< Check that this limit is below the absolute or relative threshold
  if (dV > settings::MODULE_P_V_ABSTOL || dV > settings::MODULE_P_V_RELTOL * (*V_max_it)) {
    if (verb)
      std::cout << "error in Module_p::validSUs for SU = " << getFullID() << ", the maximum voltage is in cell"
                << std::distance(Vall.begin(), V_max_it) << " and is " << *V_max_it
                << " while the minimum voltage is in cell" << std::distance(Vall.begin(), V_min_it) << " and is "
                << *V_min_it << " which is an error of " << dV << " and the allowed absolute tolerance is "
                << settings::MODULE_P_V_ABSTOL << ", the allowed relative tolerance gives an error of "
                << settings::MODULE_P_V_RELTOL * (*V_max_it) << '\n';

    return false;
  } //!< else the voltage is valid

  return true;
}

void Module_p::timeStep_CC(double dt, int nstep)
{
  /*
   * Take a time step at a constant current.
   * There are two ways to do this:
   * 		rootFinding::Current_EQN with dti = dt, which explicitly solves the system of equations
   * 		let every cell take a CC time step, and check if the voltage equation is satisfied
   * 			if not, redistribute the current using setCurrent()
   *
   * The second approach is probably quicker since it only solves the system of equations
   * if the voltage difference becomes too large
   */

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

  //!< check if the cell's voltage is valid
  if (!validSUs(SUs, false)) {
    try {
      auto status = redistributeCurrent_new(false, true); //!< don't check the currents

      if (status != Status::Success)
        throw 100000; //!< #TODO
    } catch (int e) {
      std::cerr << "error in Module_p::timeStep_CC when redistributing the current. Throwing the error on " << e << '\n';
      throw e;
    }
  }
}

Module_p *Module_p::copy()
{
  //!< check the type of coolsystem we have #TODO for a better way.
  int cooltype = 0;
  if (typeid(*getCoolSystem()) == typeid(CoolSystem_HVAC))
    cooltype = 1;
  else if (typeid(*getCoolSystem()) == typeid(CoolSystem_open))
    cooltype = 2;

  Module_p *copied_ptr = new Module_p(getID(), cool->T(), true, par, getNcells(), cool->getControl(), cooltype);

  copied_ptr->Rcontact = Rcontact;
  copied_ptr->setT(T());

  copied_ptr->SUs.clear();
  for (size_t i{ 0 }; i < getNSUs(); i++) {
    copied_ptr->SUs.emplace_back(SUs[i]->copy()); // #TODO remove when we have Module<...>
    copied_ptr->SUs.back()->setParent(copied_ptr);
  }

  return copied_ptr;
}
} // namespace slide