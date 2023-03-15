/*
 * Module_base.cpp
 *
 *  Created on: 29 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Module.hpp"

#include "../settings/settings.hpp"

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <string_view>

//!< common implementation for all base-modules
namespace slide {

Module::Module(std::string_view ID_, double Ti, bool print, bool pari, int Ncells, int coolControl, int cooltype)
  : StorageUnit(ID_), par{ pari }
{
  /*
   * IN
   * Ti 			temperture of the coolant
   * print 		print error warnings
   * pari 		use multithreaded computation to calculate the time steps of the connected SUs or not
   * Ncells 		number of cells which will be connected to this module.
   * 				This is necessary to properly size the cooling system:
   * 				modules with more cells will have more heat generation, so need better cooling.
   * 				This is implemented by increasing the thermally active surface area of the module,
   * 				and the flow rate (and cross section) of the coolant in the CoolSystem.
   * coolControl 	integer deciding how the cooling system is controlled.
   * cooltype 	determines what type of coolsystem this module has
   * 				0: regular coolsystem
   * 				1: this module is the top level module and has an HVAC coolsystem
   * 					which means it has a fan to cool the child SUs and an air conditioning unit to cool the module from the environment
   * 				2: this module is a mid-level module and has an open coolsystem (pass through between parent module and child modules)
   */

  //!< Set the module temperature
  therm.A = 0.0042 * 10 * Ncells; //!< thermally active surface area. The first number is the thermal active surface area of a cell
  double Q0 = 0;                  //!< constant ancillary losses. There are none since a module only has cells

  if (cooltype == 1)
    cool = make<CoolSystem_HVAC>(Ncells, coolControl, Q0);
  else if (cooltype == 2)
    cool = make<CoolSystem_open>(Ncells, coolControl);
  else
    cool = make<CoolSystem>(Ncells, coolControl);

  cool->setT(Ti);
}

Module::Module(std::string_view ID_, double Ti, bool print, bool pari, int Ncells, CoolSystem_t &&coolControlPtr, int cooltype)
  : StorageUnit(ID_), cool(std::move(coolControlPtr)), par{ pari }
{
  therm.A = 0.0042 * 10 * Ncells; //!< thermally active surface area. The first number is the thermal active surface area of a cell
  cool->setT(Ti);
}


void Module::setSUs(SUs_span_t c, bool checkCells, bool print)
{
  /*
   *	Note: this function can change the number of cells connected to a module
   *	It also sets the parent of the storage units in c to this module
   *	All connection resistances are reset to 0.
   *
   * THROWS
   * 10 	the cells are illegal to be combined into a module and checkStates = true
   * 		or the cells already belong to a different module
   */

  const bool verb = print && (settings::printBool::printCrit);

  //!< check the cells don't have parents yet (unless it is this module)
  for (size_t i = 0; i < c.size(); i++) {
    auto p = c[i]->getParent();
    if (p != nullptr && p != this) //!< #TODO probably cannot be this since it is a unique_ptr
    {
      if (verb)
        std::cerr << "ERROR in Module::setCells, SU " << i << ", already has a parent "
                  << "with full ID: " << parent->getFullID() << ". Throwing 10.\n";
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      throw 10;
    }
  }

  //!< check the cells are valid according to this module layout
  // #TODO checkCells & setSUs removed from here

  //!< connect the cells to this module
  Vmodule_valid = false; //!< we are changing the SUs, so the stored voltage is no longer valid

  SUs.clear();
  Rcontact.clear();
  size_t r{ 0 };
  for (auto &SU : c) {
    r += SU->getNcells();
    SU->setParent(this); //!< Set the parent of all cells/ Does not throw.
    SUs.push_back(std::move(SU));
    Rcontact.push_back(0); //!< Make Rcontact zero? #TODO
  }

  Ncells = r;
}


Status Module::checkVoltage(double &v, bool print) noexcept
{
  /*
   * Calculate the module voltage and check whether it is valid
   * 		it checks both the limits of the module and the limits of each cell individually
   * 		the 'worst' number is returned.
   * 			i.e. the one with the largest absolute value,
   * 			i.e. if one cell is outside the safety limits and the rest is fine, than the one outside the safety limits is returned
   */

  //!< const bool printCrit = print && (settings::printBool::printCrit);		//!< print if the (global) verbose-setting is above the threshold
  //!< const bool printNonCrit = print && (settings::printBool::printNonCrit); //!< print if the (global) verbose-setting is above the threshold

  //!< check the voltage of the module
  v = V(); // -> Hey we dont need to calculate this anymore.
  //!< #TODO here was a not useful chuck of code for repeated checking.
  //!< We may have constant limits for module voltage.

  //!< check the voltages of the connected cells
  double vi;
  auto res = Status::Success;
  for (const auto &SU : SUs)
    res = std::max(res, SU->checkVoltage(vi, print)); // get the worst status.

  return res;
}

double Module::getVhigh()
{
  //!< return the voltage of the cell with the highest voltage
  //!< 	note CELL not child SU
  double Vhigh = Vmin(); // #TODO
  for (const auto &SU : SUs)
    Vhigh = std::max(Vhigh, SU->getVhigh()); //!< will be called recursively to the cell levels

  return Vhigh;
}

double Module::getVlow()
{
  //!< return the voltage of the cell with the lowest voltage note CELL not child SU
  double Vlow = Vmax(); // #TODO
  for (const auto &SU : SUs)
    Vlow = std::min(Vlow, SU->getVlow()); //!< will be called recursively to the cell levels

  return Vlow;
}

void Module::getStates(getStates_t s)
{
  /*
   * Returns the module-states, as well as the states of all cells in one long array
   * [s0 s1 s2 ... sn Tmod]
   * where s0 is the array with the states of the first cell of this module
   */
  for (const auto &SU : SUs)
    SU->getStates(s); //!< pass a vector, the next nsi locations will be automatically filled with the states of cell i

  s.push_back(T()); //!< store the module temperature
}

bool Module::validStates(bool print)
{
  /*
   * function to check if the states are compatible for the cells connected to this series-connected module
   *
   * This is checked by trying to set the states of this module to the given value and checking if the resulting state is valid.
   * Note that the voltage limits of cells are NOT checked, that is only done by setStates(checkV=true)
   * then the original state of this module is restored.
   *
   * 	todo temperatures?
   *
   * This is a slow function, so only execute it if you have to
   */

  const bool verb = print && (settings::printBool::printNonCrit); //!< print if the (global) verbose-setting is above the threshold

  return true; // #TODO here we probably need to check if all submodule states valid!
}

Status Module::setStates(setStates_t s, bool checkV, bool print)
{
  /*
   * The states are first set, and then checked for validity (valid states are always checked, valid voltages only if checkV == true).
   * If they are not valid, the original states are restored and an error is thrown.
   *
   * This is because we cannot check s is valid otherwise since we don't know the exact type of the connectedStorageUnits.
   * Therefore, we don't know if s[0] is the SOC of a Cell, or the first state of a Module
   * So the only way to check if s is valid is by actually setting the states, the polymorphism will select the correct states
   *
   * IN
   * checkV 	if false
   * 				cells check whether their state is valid (e.g. SOC in range) but not if their voltage is valid
   * 				module checks that the module constraints are valid (e.g. same I if series module, same V if parallel module
   * 			if true, cells also check whether their voltage is valid and so does the entire module
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range, new states are kept
   * 3 	checkV is true && the voltage is outside the safety range, old sates are restored
   */

  const bool verb = print && (settings::printBool::printCrit);

  std::vector<double> sorig;
  getStates(sorig);

  std::span<double> spn_orig{ sorig };

  Vmodule_valid = false; //!< we are changing the states, so the stored voltage is no longer valid

  //!< set the new cell states
  for (size_t i = 0; i < getNSUs(); i++) {

    const Status status = SUs[i]->setStates(s, checkV, print); //!<  set the states

    if (verb && isStatusWarning(status))
      std::cout << "warning in Module::setStates, the voltage of cell " << i << " with id "
                << SUs[i]->getFullID() << " is outside the allowed range. Continue for now.\n";
    else if (isStatusBad(status)) {
      if (verb)
        std::cerr << "ERROR in Module::setStates when setting the state of cell " << i
                  << ". Restoring the old states, status: " << getStatusMessage(status) << '\n';

      for (size_t j = 0; j <= i; j++)
        SUs[i]->setStates(spn_orig, false, print); //!< restore the original states without checking validity (they should be valid)

      return status;
    }

  } //!< end loop to set the cell states

  //!< set the module temperature
  assert(s[0] >= PhyConst::Kelvin); //!< #TODO here we are checking but should we?
  setT(s[0]);
  s = s.last(s.size() - 1);

  //!< check that the cells are valid for this module configuration (same I if series, same V if parallel)
  /*
   * Note: only check this if checkV is on, else you can get an eternal loop if the initial state is not valid.
   * (in that case, setState(sorig) will also fail, and call its own setState(sorig), etc.)
   * So we must have a way of restoring the state without checking anything.
   * OR alternatively, at the beginning we check that the initial state is valid.
   */
  if (checkV) { // #TODO
    // bool valCells = validSUs(SUs, print);
    // if (!valCells) {
    //   if (verb)
    //     std::cerr << "ERROR in Module:setStates for SU = " << getFullID() << ", the state is illegal "
    //               << "according to validCells(), returning the old states and throwing 10.\n";

    //   return Status::Invalid_states;
    // }
  }

  return Status::Success; //!< return success.
}

double Module::thermalModel_cell()
{
  /*
   * Calculate the thermal model of each cell individually
   */

  //!< array with the new temperatures of the child SUs
  double Tnew[settings::MODULE_NSUs_MAX];

  //!< dummy arrays
  double Tsu[1], Ksu[1], Asu[1];
  Tsu[0] = 0;
  Ksu[0] = 0;
  Asu[0] = 0;
  double tim = 0;

  //!< calculate cells' temperature
  for (size_t i = 0; i < getNSUs(); i++) {
    try {
      Tnew[i] = SUs[i]->thermalModel(3, Tsu, Ksu, Asu, tim);
    } catch (int e) {
      if constexpr (settings::printBool::printCrit)
        std::cout << "Error in module " << getFullID() << " when calculating the thermal balance of child SU "
                  << i << ", error " << e << ". throwing it on.\n";
      throw e;
    }
  }

  //!< Set all the new temperatures to the children
  for (size_t i = 0; i < getNSUs(); i++)
    SUs[i]->setT(Tnew[i]);

  //!< the temperature of a module doesn't change
  return T();
}

double Module::thermalModel_coupled(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim)
{
  /*
   * Calculate the thermal model of this SU
   *
   * It calls the thermal model of all child-SUs as well
   * It cools the child SUs, and transfers heat between them.
   * Cells at the edges get the same amount of cooling as cells in the middle.
   * But at one edge, they exchange heat with the module (with T of coolant; and heat transfer coefficient similar to cell to cell).
   * Cells in the middle exchange heat with cells on both sides
   *
   * throws
   * 98 	invalid time keeping
   * 99 	invalid module temperature
   */

  //!< Check the time keeping
  if constexpr (settings::printBool::printNonCrit)
    if (therm.time > 15 * 60)
      std::cout << "Warning in Module::thermalModel, the time since this function was called last is very large, "
                << therm.time << " which might lead to excessive temperature variations.\n";

  //!< then check whether our internal time keeping matches up with the external one
  if (std::abs(therm.time - tim) > 1) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Module::thermalModel of SU " << getFullID() << ", according to the cell's internal timing, "
                << therm.time << "s have passed since the last thermal model solution."
                << " The external time provided was " << tim << "s, which is more than 1s difference. Throwing an error.\n";
    throw 98;
  }

  //!< Only calculate the thermal model if time has actually passed. Else the cooling temperature will become NaN since the volume is 0
  if (tim != 0) {
    //!< array with the new temperatures of the child SUs
    double Tnew[settings::MODULE_NSUs_MAX];

    //!< ************************************************************* heat exchange with child SUs *******************************************************************
    //!< Make the arrays for heat exchange to each child SU.
    //!< They have max length 3 (cooling from here, and 2 neighbours)
    double Tsu[3], Ksu[3], Asu[3]; //!< parent, left neighbour, right neighbour

    //!< The first element is the cooling from the parent module to the child SUs
    Tsu[0] = T();
    Ksu[0] = cool->getH();
    Asu[0] = therm.A;

    //!< Loop to cool each child SU
    for (size_t i = 0; i < getNSUs(); i++) {

      //!< left cell
      Ksu[1] = therm.k_cell2cell; //!< conductive heat exchange via long sides of cell

      if (i > 0) {
        Tsu[1] = SUs[i - 1]->T(); //!< left is cell i-1
        Asu[1] = SUs[i - 1]->getThermalSurface();
      } else { //!< left is module
        Asu[1] = getThermalSurface();
        Tsu[1] = T();
      }

      //!< right cell
      Ksu[2] = therm.k_cell2cell;
      if (i + 1 < getNSUs()) //!< Last cell. getNSUs() is unsigned therefore getNSUs() -1 is omitted.
      {
        Asu[2] = SUs[i + 1]->getThermalSurface();
        Tsu[2] = SUs[i + 1]->T();
      } else { //!< right is edge of module
        Asu[2] = getThermalSurface();
        Tsu[2] = T();
      }

      //!< calculate thermal balance of the child SU
      try {
        Tnew[i] = SUs[i]->thermalModel(3, Tsu, Ksu, Asu, tim);
      } catch (int e) {
        if constexpr (settings::printBool::printCrit)
          std::cout << "Error in module " << getFullID() << " when calculating the thermal balance of child SU "
                    << i << ", error " << e << ". throwing it on.\n";
        throw e;
      }
    }

    //!< The cooling fluid in the module heats up from cooling all cells
    double Etot = 0;
    for (size_t i = 0; i < getNSUs(); i++) {
      const double Atherm = std::min(Asu[0], SUs[i]->getThermalSurface());
      Etot += Ksu[0] * Atherm * (SUs[i]->T() - T()) * tim;

      //!< additional cooling to the first and last cell of the stack (which both have 1 edge from the coolsystem)
      if (i == 0)
        Etot += Ksu[1] * Atherm * (SUs[i]->T() - T()) * tim;

      if (i == getNSUs() - 1) //!< #TODO problem.
        Etot += Ksu[2] * Atherm * (SUs[i]->T() - T()) * tim;
    }

    //!< Add up the heat generated in all the contact resistances of this Module
    Etot += therm.Qcontact;

    double Echildren = Etot; //!< cooling energy extracted from the children

    //!< *********************************************************** Heat exchange with neighbours and parent **********************************************************

    //!< The module exchanges heat with its parent and neighbours
    for (int i = 0; i < Nneighbours; i++) {
      const double Atherm = std::min(Aneighb[i], therm.A);
      Etot += Kneighbours[i] * Atherm * (Tneighbours[i] - T()) * tim;
    }

    //!< *************************************************************** New coolant temperature *******************************************************************

    //!< Calculate the new temperature of the coolant
    //!< rho * cp * dT/dt = Qtot / V
    //!< 		where 	Qtot is total power in W
    //!< 				V is the total volume of fluid over the time period, the product of the flow rate and the time (therm_v*tim)
    //!< so integrated over time this is
    //!< rho * cp * (Tnew - Told) = Etot / V
    double Tcool_new = cool->dstate(Etot, Echildren, tim);

    //!< Set all the new temperatures to the children
    for (size_t i = 0; i < getNSUs(); i++)
      SUs[i]->setT(Tnew[i]);

    //!< reset the time since the last update of the thermal model
    therm.time = 0;

    //!< Check the new temperature is valid
    if (Tcool_new < PhyConst::Kelvin || Tcool_new > PhyConst::Kelvin + 75.0 || std::isnan(Tcool_new)) {
      if constexpr (settings::printBool::printCrit) {
        std::cerr << "ERROR in Module::thermalModel of SU " << getFullID() << ", the new temperature of " << Tcool_new << " is outside the allowed range from (273+0) K to (273+75) K";
        std::cerr << ". The time since the last time this function was called is " << tim << '\n';
      }
      throw 99;
    }

    //!< return the new cooling temperature
    return Tcool_new;

  }    //!< if(tim != 0)
  else //!< no time has passed -> no heat generated -> just return present temperature
    return T();
}

/*
 * Calculate the thermal model of this SU and its children
 *
 * throws
 * 98 	invalid time keeping
 * 99 	invalid module temperature
 */
double Module::thermalModel(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim)
{
  double Tout;
  try {
    if constexpr (settings::T_MODEL == 0)
      Tout = T();
    else if constexpr (settings::T_MODEL == 1)
      Tout = thermalModel_cell();
    else if constexpr (settings::T_MODEL == 2)
      Tout = thermalModel_coupled(Nneighbours, Tneighbours, Kneighbours, Aneighb, tim);
  } catch (int e) {
    //!< indicate we have ran the thermal model
    therm.time = 0;
    therm.Qcontact = 0;
    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw e;
  }

  //!< Reset the cumulative thermal variables
  therm.time = 0;
  therm.Qcontact = 0;

  //!< control the cooling system is done at the time integration function
  return Tout;
}

double Module::getCoolingLoad()
{
  /*
   * Calculate and return the energy required to operate the cooling system of this module and all of its children
   * since the last time this function was called.
   * unit: [J]
   *
   * Coolsystems keep track of how much energy they consume cumulatively.
   * This variable can be reset to 0 by coolsystem::reset_Eoperation().
   * After calculating how much energy we have used, this function is called such that all coolsystems are reset.
   * This means we start again from 0, so the next time this function is called it will return the amount of energy consumed since now.
   */

  double Etot = getCoolSystem()->getEoperation(); //!< energy to run coolsystem of this module
  getCoolSystem()->reset_Eoperation();            //!< reset to 0

  //!< loop for the children, and if they are modules, add their cooling energy
  for (auto &SU : SUs)
    if (auto m = dynamic_cast<Module *>(SU.get())) //!< cast to a module pointer instead of storage unit
    {
      Etot += m->getCoolingLoad();            //!< add energy
      m->getCoolSystem()->reset_Eoperation(); //!< reset cooling energy for that child
    }

  return Etot;
}

} // namespace slide