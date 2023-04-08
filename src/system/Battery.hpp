/*
 * Battery.hpp
 *
 *  Created on: 11 Jun 2020
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../StorageUnit.hpp"
#include "../modules/Module.hpp"
#include "../cooling/CoolSystem_HVAC.hpp"
#include "../power_conversion/Converter.hpp"
#include "../types/data_storage/cell_data.hpp"

#include <memory>
#include <vector>
#include <span>

namespace slide {

class Battery : public StorageUnit
{
protected:
  Deep_ptr<Module> cells{};                  //!< Module with the cells of this battery
  Deep_ptr<CoolSystem_HVAC> cool{};          //!< HVAC system of the battery
  Converter conv{};                          //!< power electronic converter. Dual step DC/DC and DC/AC
  unsigned int nseries{ 1 }, nparallel{ 1 }; //!< number of series/parallel 'copies' of this module

  double convlosses{};     //!< losses in the converter during a given period (set to 0 by reset_convlosses)
  double convlosses_tot{}; //!< total cumulative losses in the converter during the entire lifetime

#if DATASTORE_BATT > 1
  double timetot{ 0 }; //!< time so far [s]
  BatteryData batData;
#endif

public:
  Battery();
  Battery(std::string IDi);
  void setModule(Deep_ptr<Module> &&module);
  void setSeriesandParallel(unsigned int ser, unsigned int par);

  //!< basic getters and setters
  double Cap() const override { return cells->Cap() * nparallel; }
  double Vmin() const override { return cells->Vmin() * nseries; }
  double VMIN() const override { return cells->VMIN() * nseries; }
  double Vmax() const override { return cells->Vmax() * nseries; }
  double VMAX() const override { return cells->VMAX() * nseries; }
  double I() const override { return cells->I() * nparallel; }
  double getRtot() override { return cells->getRtot() * nseries / static_cast<double>(nparallel); }
  size_t getNcells() override { return cells->getNcells() * nseries * nparallel; }

  Module *getCells() { return cells.get(); }
  //!< int getNstates() { return cells->getNstates() + 1; } //!< +1 for the temperature of this battery
  void getStates(getStates_t &s) override; //!< returns one long array with the states
  void setBlockDegAndTherm(bool block);

  Status setCurrent(double Inew, bool checkV = true, bool print = true) override
  {
    return cells->setCurrent(Inew / nparallel, checkV, print);
  }

  Status setStates(setStates_t s, bool checkStates = true, bool print = true) override; //!< opposite of getStates, check the states are valid?
  double getAndResetConvLosses();
  double getConvLosses_total() { return convlosses_tot; }
  void resetConvLosses() { convlosses = 0; }

  //!< voltage
  double getOCV() override { return cells->getOCV() * nseries; }
  double V() override { return cells->V() * nseries; }
  Status checkVoltage(double &v, bool print) noexcept override;      //!< get the voltage and check if it is valid
  double getVhigh() override { return cells->getVhigh() * nseries; } //!< return the voltage of the cell with the highest voltage
  double getVlow() override { return cells->getVlow() * nseries; }   //!< return the voltage of the cell with the lowest voltage

  //!< thermal model
  double T() override { return cool->T(); }
  double getThotSpot() override { return std::max(cells->getThotSpot(), T()); } //!< the T of the hottest element in the SU
  double getThermalSurface() override { return cells->getThermalSurface(); }    //!< return the 'A' for the thermal model of this SU (Q = hA*dT)

  double thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim) override; //!< calculate the thermal model of this SU
  void setT(double Tnew) override { cool->setT(Tnew); }
  auto *getCoolSystem() { return cool.get(); }
  double getCoolingLoad(); //!< return the energy required to run the entire coolingsystem of this module and all its children

  //!< functionality
  bool validStates(bool print = true) override { return cells->validStates(print); }
  //!< checks if a state array is valid
  void timeStep_CC(double dt, int steps = 1) override; //!< take a number of time steps

  //!< Data collection of cycling data (I, V, T, etc. for every cell)
  void storeData() override;
  void writeData(const std::string &prefix) override;

  Battery *copy() override { return new Battery(*this); }
};

} // namespace slide
