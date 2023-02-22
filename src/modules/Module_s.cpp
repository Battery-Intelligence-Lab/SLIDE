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

double Module_s::getOCV(bool print)
{ //!< sum of the open circuit voltage of all cells
  try {
    return transform_sum(SUs, free::get_OCV<SU_t>);
  } catch (int e) {
    if (print && (settings::printBool::printCrit)) //!< print if the (global) verbose-setting is above the threshold
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

double Module_s::V(bool print)
{
  //!< sum of the voltage of all cells
  //!< bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  //!< if the stored value is up to date, return that one
  if (Vmodule_valid)
    return Vmodule;

  //!< else calculate the voltage
  Vmodule = 0;
  for (size_t i = 0; i < getNSUs(); i++) {
    const auto v_i = SUs[i]->V(print);

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
  bool validVcell = true; //!< #TODO
  Vmodule_valid = false;  //!< we are changing the current, so the stored voltage is no longer valid
  double Iold = I();

  for (size_t i = 0; i < getNSUs(); i++) {
    try {
      SUs[i]->setCurrent(Inew, checkV, print); //!< note: this will throw 2 or 3 if the voltage of a cell is illegal
    } catch (int e) {
      if (e == 2) {
        //!< voltage of cell i is outside the valid range, but within safety limits
        //!< indicate this happened but continue setting states
        validVcell = false;
        if (verb)
          std::cout << "warning in Module_s::setCurrent, the voltage of cell " << i
                    << " with id " << SUs[i]->getFullID()
                    << " is outside the allowed range for Inew = " << Inew << ". Continue for now\n";
      } else {
        //!< error 10 is illegal state
        //!< error 3 means the voltage is outside the safety limit
        if (verb)
          std::cerr << "ERROR in Module_s::setCurrent when setting the current of cell " << i << " for Inew = "
                    << Inew << ". Restoring the old currents and throwing on error " << e << '\n';
        setCurrent(Iold, false, print); //!< restore the original current without checking validity (they should be valid)
        std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
        throw e;
      }
    }
  }

  //!< check and return the voltage of the module
  //!< Check if the voltage is valid

  //!< #TODO Here we need module specific voltage.

  return Status::Success;
}

bool Module_s::validSUs(SUs_span_t c, bool print)
{
  /*
   * Checks the cells are a valid combination for a series-connected module
   * the current is the same in each cell
   *
   * Note that nothing else is checked. This function should not be relied on on its own to check validity
   * For that, you shoud use setStates
   *
   */

  // #TODO -> Unnecessary function. Check validity when settings SUs or states.
  // Algorithms should not leave anything in an invalid state.

  bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold
  bool result{ true };
  //!< check the currents are the same
  const double Imod = c[0]->I();
  for (size_t i = 1; i < c.size(); i++) {
    const double err = std::abs(Imod - c[i]->I());
    bool val = err < settings::MODULE_P_I_ABSTOL || err < Imod * settings::MODULE_P_I_RELTOL; //!< #TODO should not be || here. Ok got it due to 0 current.
    //!< in complex modules, an s can be made out of p modules, and p modules are allowed to have a small error in their current
    if (!val) {
      if (verb)
        std::cout << "error in Module_s::validCells, the current of cell " << i << " is "
                  << c[i]->I() << "A while the current in the 0th cell is " << Imod << "A.\n";

      result = false;
      break;
    }
  }

  return result; //!< else the cells are valid
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
    for (const auto r : Rcontact)
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

Module_s *Module_s::copy()
{
  //!< check the type of coolsystem we have #TODO for a better way. Also same for both modules.

  int cooltype = 0;

  if (typeid(*getCoolSystem()) == typeid(CoolSystem_HVAC))
    cooltype = 1;
  else if (typeid(*getCoolSystem()) == typeid(CoolSystem_open))
    cooltype = 2;

  Module_s *copied_ptr = new Module_s(getID(), cool->T(), true, par, getNcells(), cool->getControl(), cooltype);

  copied_ptr->Rcontact = Rcontact;
  copied_ptr->setT(T());

  for (size_t i{ 0 }; i < getNSUs(); i++) {
    copied_ptr->SUs.emplace_back(SUs[i]->copy());
    copied_ptr->SUs.back()->setParent(copied_ptr);
  }

  return copied_ptr;
}

} // namespace slide