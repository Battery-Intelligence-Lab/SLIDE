/*
 * Battery.cpp
 *
 *  Created on: 11 Jun 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Battery.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <memory>

namespace slide {

Battery::Battery() : StorageUnit("Battery") {}
Battery::Battery(std::string IDi) : Battery() { ID = std::move(IDi); }

void Battery::setSeriesandParallel(unsigned int ser, unsigned int par)
{
  nseries = ser;
  nparallel = par;
}

void Battery::setModule(Deep_ptr<Module> &&module)
{
  /*
   * Connect the given module to this battery.
   * The module must have no parent yet, and it must have a conventional or open coolsystem.
   * 	HVAC coolsystems are not allowed, since if you make a battery, the HVAC is part of the battery and not the module with the cells.
   *
   * THROWS
   * 10 	the module cannot be connected to this battery
   */

  //!< check the Module is compatible
  //!< 		has no parent yet
  //!< 		does not have an HVAC cool system
  if (module->getParent() != nullptr) {
    std::cerr << "ERROR in Battery::setModule. The parent of module " << module->getFullID() << " is "
              << module->getParent()->getFullID() << ". We can only connect modules without parent. Throwing 10.\n";
    throw 10;
  }
  if (module->getCoolSystem() == nullptr) {
    std::cerr << "ERROR in Battery::setModule. Module " << module->getFullID() << " no coolsystem. Throwing 10.\n";
    throw 10;
  } else if (typeid(*module->getCoolSystem()) == typeid(CoolSystem_HVAC)) {
    std::cerr << "ERROR in Battery::setModule. Module " << module->getFullID()
              << " has an HVAC coolsystem so we cannot connect it to the battery. Throwing 10.\n";
    throw 10;
  }

  //!< check that we don't have a module or coolsystem yet
  if (cells != nullptr || cool != nullptr) {
    std::cerr << "ERROR in Battery::setModule. Battery " << getFullID()
              << " already has a module or a coolsystem. Throwing 10.\n";
    throw 10;
  }

  //!< connect here
  cells = std::move(module);

  //!< set the parent of the module to this
  cells->setParent(this);

  //!< make HVAC coolsystem
  /* Average losses in the converter (which are fairly constant independent on power, so a constant ancillary load
   * 		single cell: 1.5 kW 		1500 W  per cell
   * 		10 parallel: 2.2 kW		  220 W   per cell
   * 		10 series: 1.5 kW 		  150 W   per cell
   * 		10s 10p = 2.2 kW 		    22 W    per cell
   * 		100s 10p = 2.2 kW 		  2.2 W   per cell
   * simply estimate based on losses at max voltage and 1C current
   * and then do * 2
   */
  double Q0 = conv.getLosses(Vmax(), Cap()) * 2; //!< idle losses of the converter (due to switching and others)
  cool = make<CoolSystem_HVAC>(getNcells(), cells->getCoolSystem()->getControl(), Q0);

  //!< Scale the power electronic converter correspondingly
  conv.setPower(Cap() * Vmax());
}

Status Battery::checkVoltage(double &v, bool print) noexcept
{
  Status r = cells->checkVoltage(v, print);
  v *= nseries;
  return r;
}

void Battery::getStates(getStates_t s)
{
  //!< States of the cells and the temperature of the battery
  cells->getStates(s);
  s.push_back(T());
}

void Battery::setBlockDegAndTherm(bool block)
{
  blockDegAndTherm = block;
  cells->setBlockDegAndTherm(block);
}

Status Battery::setStates(setStates_t s, bool checkStates, bool print)
{
  auto status = cells->setStates(s, checkStates, print);
  setT(s.back()); //!< #TODO probably here we need to check?
  return status;
}

double Battery::getAndResetConvLosses()
{
  double loss = convlosses;
  convlosses = 0;
  return loss;
}

double Battery::thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim)
{
  /*
   * Calculate the thermal model of the battery.
   * This is the top level, so there are no neighbours or parents. Therefore, we ignore the inputs
   */

  double Tbatt;

  if constexpr (settings::T_MODEL == 0)
    Tbatt = T();
  else {
    //!< Calculate the model of the cells, the battery's HVAC coolsystem will cool the cells
    double Tsu[1], Ksu[1], Asu[1];
    Tsu[0] = T();
    Ksu[0] = cool->getH();
    Asu[0] = getThermalSurface();
    double Tcells = cells->thermalModel(1, Tsu, Ksu, Asu, tim);

    //!< The battery heats up from cooling all cells
    double Echildren = Ksu[0] * Asu[0] * (cells->T() - T()) * tim; //!< cooling energy extracted from the children
    //!< double Qcells = Echildren;

    //!< and from the losses in the converter
    Echildren += conv.getLosses(V(), I()) * tim;

    //!< there are no neighbours or parents, so ignore inputs
    if (Nneighb != 0) {
      std::cerr << "ERROR in Battery::thermalModel, the battery must be a stand-alone unit so it cannot have"
                   "neighbours or parents to exchange heat with.\n";
      throw 10;
    }

    double Etot = Echildren;

    //!< update the battery temperature
    Tbatt = cool->dstate(Etot, Echildren, tim);
    //!< cout<<"Battery thermal balance: Qcells = "<<Qcells<<", converter "<<conv->getLosses(cells->V(), cells->I())*tim<<" resutling in new T "<<Tbatt<<" for battery power "<<cells->V()*cells->I()<<endl;

    //!< Check the new temperature is valid
    if (Tbatt < PhyConst::Kelvin || Tbatt > PhyConst::Kelvin + 75.0 || std::isnan(Tbatt)) //!< #TODO -> 75.0 is magical number.
    {
      if constexpr (settings::printBool::printCrit)
        std::cerr << "ERROR in Battery::thermalModel for battery " << getFullID() << ", the new temperature of "
                  << Tbatt << " is outside the allowed range from (273+0) K to (273+75) K\n";
      throw 99;
    }

    //!< Set all the new temperatures
    cells->setT(Tcells);
    setT(Tbatt);
  }

  //!< return the new cooling temperature
  return Tbatt;
}

double Battery::getCoolingLoad()
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

  //!< heat of the coolsystem of this battery
  double Etot = getCoolSystem()->getEoperation(); //!< energy to run coolsystem of this module
  getCoolSystem()->reset_Eoperation();            //!< reset to 0

  //!< coolsystem of child cells
  Etot += cells->getCoolingLoad();

  return Etot;
}

void Battery::timeStep_CC(double dt, int steps)
{

  //!< integrate in time for the cells
  cells->timeStep_CC(dt, steps);

  //!< increase the losses from the converter
  double l = conv.getLosses(V(), I()) * dt * steps; //!< losses [J] during this period
  convlosses += l;
  convlosses_tot += l;

  //!< Calculate the thermal model
  if (!blockDegAndTherm) {
    //!< Call the thermal model without heat exchanges with neighbours or parents (since this module doesn't have any)
    double Tneigh[1], Kneigh[1], Aneigh[1];              //!< Make the arrays even though they will not be used (length should be 0 but I don't think you can make an array of length 0)
    thermalModel(0, Tneigh, Kneigh, Aneigh, steps * dt); //!< the 0 signals there are no neighbours or parents
  }

  //!< control the cooling system
  double Tlocal = transform_max(cells->getSUs(), free::get_T<Module::SU_t>);
  cool->control(Tlocal, getThotSpot());

//!< data storage
#if DATASTORE_BATT > 1
  timetot += dt * steps;
#endif
}

void Battery::storeData()
{
  //!< Store data of the cells and the cooling system
  cells->storeData();
  cool->storeData(getNcells());

//!< store local data
#if DATASTORE_BATT == 2
  //!< add a new point in all the arrays
  batData.push_back({ I(), V(), T(), conv->getLosses(V(), I()), timetot });
#endif
}
void Battery::writeData(const std::string &prefix)
{

//!< Write data for the cooling system
#if 0 //!< DATASTORE_COOL > 0
		cool->writeData(prefix + '_' + getFullID()); //!< append the ID to the prefix since the cooling system does not have the ID of this module
#endif

  //!< write data for the cells
  cells->writeData(prefix);

//!< write local data
#if DATASTORE_BATT == 2
  std::string name = prefix + "_" + getFullID() + "_battData.csv"; //!< name of the file, start with the full hierarchy-ID to identify this cell

  //!< append the new data to the existing file
  std::ofstream file;
  file.open(name, std::ios_base::app);

  if (!file.is_open()) {
    std::cerr << "ERROR in Cell::writeData, could not open file " << name << '\n';
    throw 11;
  }

  for (const auto &d : batData)
    file << d.Timetot << ',' << d.Icells << ',' << d.Vcells
         << ',' << d.Tbatt << ',' << d.convloss << '\n';

  file.close();

  batData.clear(); //!< reset
#endif
}

} // namespace slide