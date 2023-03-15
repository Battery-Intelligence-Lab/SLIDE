/*
 * Module_s.cpp
 *
 * series-connected base Module
 *
 *  Created on: 9 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Module_s.hpp"
#include "../utility/utility.hpp"

// #include "settings/settings.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <thread>
#include <memory>
#include <ctime>
#include <numeric>

namespace slide {

double Module_s::getOCV()
{ //!< sum of the open circuit voltage of all cells
  try {
    return transform_sum(SUs, free::get_OCV<SU_t>);
  } catch (int e) {
    if constexpr (settings::printBool::printCrit) //!< print if the (global) verbose-setting is above the threshold
      std::cout << "Error in Module_s::getOCV when getting the OCV of a SU in module "
                << ID << ", throwing it on.\n";
    throw e;
  }
}

double Module_s::getRtot()
{
  /*
   * Return the total resistance
   * 		V(I) = OCV - I*Rtot
   * 		with V and OCV the total values for this module
   */

  double rtot = 0;
  for (size_t i = 0; i < SUs.size(); i++)
    rtot += SUs[i]->getRtot() + Rcontact[i];

  return rtot;
}

double Module_s::V()
{
  //!< sum of the voltage of all cells
  //!< bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  //!< if the stored value is up to date, return that one
  if (Vmodule_valid)
    return Vmodule;

  //!< else calculate the voltage
  Vmodule = 0;
  for (size_t i = 0; i < getNSUs(); i++) {
    const auto v_i = SUs[i]->V();

    if (v_i <= 0) //!< SU has an invalid voltage.
      return 0;   // #TODO check if this is same for getOCV?

    Vmodule += v_i - Rcontact[i] * SUs[i]->I();
  }

  Vmodule_valid = true;
  return Vmodule;
}

Status Module_s::setCurrent(double Inew, bool checkV, bool print)
{
  /*
   * checkV 	check the module voltage is in the allowed range AND that the voltage of each cell is inside the allowed range
   *
   * returns the voltage (if checkV)
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
   * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
   */

  //!< the current is the same in all cells
  bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  //!< Set the current, if checkVi this also gets the cell voltages
  Vmodule_valid = false; //!< we are changing the current, so the stored voltage is no longer valid

  std::array<double, settings::MODULE_NSUs_MAX> Iolds;

  for (int i = 0; i < getNSUs(); i++) {
    Iolds[i] = SUs[i]->I();
    const auto status = SUs[i]->setCurrent(Inew, checkV, print);

    if (isStatusBad(status)) {
      if (verb)
        std::cerr << "ERROR in Module_s::setCurrent when setting the current of cell " << i
                  << " with id " << SUs[i]->getFullID() << " for Inew = "
                  << Inew << ". Restoring the old currents and throwing on error "
                  << getStatusMessage(status) << '\n';


      for (int j{ i }; j > 0; j--)
        SUs[j]->setCurrent(Iolds[j], false, false); // #TODO this can start from j{i-1} since probably SUs[i] failed to assign any currents.

      return status; // #TODO setCurrent here may be costly.
    }
  }
  //!< Check if the voltage is valid  #TODO
  //!< #TODO Here we need module specific voltage.
  return Status::Success;
}

void Module_s::timeStep_CC(double dt, int nstep)
{
  /*
   * a time step at constant current is simply a time step of every individual cell
   *
   * THROWS
   * 10 	negative time step
   * 13 	in multithreading, an exception was thrown in one of the child threads of timeSetp_CC_i;
   * 14 	this module has no parent (i.e. is top level) but does not have an HVAC coolsystem
   * 			the top level needs an HVAC system for heat exchange with the environment
   * 			lower-level modules (with a parent module) will get cooling from the coolsystems of their parent so they have a regular coolSystem
   */
  if (dt < 0) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Module_s::timeStep_CC, the time step dt must be 0 or positive, but has value " << dt << '\n';
    throw 10;
  }

  //!< we simply take one CC time step on every cell
  auto task_indv = [&](int i) { SUs[i]->timeStep_CC(dt, nstep); };

  try {
    run(task_indv, getNSUs(), (par ? -1 : 1)); // #TODO SU based for_each.
  } catch (int e) {
    std::cout << "Error in Module_s::timeStep_CC with module ID " << getFullID()
              << ". error " << e << ", throwing it on.\n";
    throw e;
  }

  //!< **************************************************** Calculate the thermal model once for the nstep * dt time period *****************************************************************

  //!< update the time since the last update of the thermal model
  if (!blockDegAndTherm) {
    therm.time += nstep * dt;

    //!< Increase the heat from the contact resistances
    for (const auto r : Rcontact)                  // #TODO calculating I() many times.
      therm.Qcontact += r * sqr(I()) * nstep * dt; //!< each resistor sees the total module current

    //!< If this module has a parent module, this parent will call the thermal model with the correct parameters
    //!< which will include heat exchange with the module's neighbours and cooling from the cooling system of the parent module.

    //!< if there is no parent, this module is the top-level.
    //!< It then directly exchanges heat with the environment at a fixed temperature through the HVAC system
    //!< If there is no parent, assume we update T every nstep*dt. So update the temperature now
    if (!parent) {

      //!< double check that we have an HVAC coolsystem (see constructor)
      if (typeid(*getCoolSystem()) != typeid(CoolSystem_HVAC)) {
        std::cerr << "ERROR in module_s::timeStep_CC in module " << getFullID() << ". this is a top-level"
                  << " module but does not have an HVAC coolsystem for active cooling with the environment.\n";
        throw 14;
      }

      //!< Call the thermal model without heat exchanges with neighbours or parents (since this module doesn't have any)
      //!< 	Heat exchange between this module and the environment is done by the AC-part of the HVAC coolsystem
      //!< 	heat exchange between this module and its children is done the conventional way by the HVAC coolsystem
      //!< set the new temperature since we have calculated all the temperatures so there is no risk for inconsistency
      double Tneigh[1], Kneigh[1], Aneigh[1];                    //!< Make the arrays even though they will not be used (length should be 0 but I don't think you can make an array of length 0)
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
    const double Tlocal = transform_max(SUs, free::get_T<SU_t>); // #TODO Battery also has this.
    cool->control(Tlocal, getThotSpot());
  }

  Vmodule_valid = false; //!< we have changed the SOC/concnetration, so the stored voltage is no longer valid
}

} // namespace slide