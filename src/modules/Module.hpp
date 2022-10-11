/*
 * Module.hpp
 *
 *  Created on: 29 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../settings/settings.hpp"
#include "../StorageUnit.hpp"
#include "../cooling/CoolSystem_HVAC.hpp"
#include "../cooling/CoolSystem_open.hpp"
#include "../types/data_storage/cell_data.hpp"
#include "../types/data_storage/module_data.hpp"
#include "../utility/free_functions.hpp"

#include <vector>
#include <cstdlib>
#include <memory>
#include <iostream>
#include <span>


namespace slide {
struct ModuleThermalParam
{
  double k_cell2cell{ 5 };                             //!< conductive heat transfer coefficient for heat transfer betweenthe child SUs
  double A{ 0.0042 * 10 * settings::MODULE_NSUs_MAX }; //!< thermally active surface area of this module. The first number is the thermal active surface area of a cell
  double Qcontact{ 0 };                                //!< heat energy generated in the contact resistances since the last time the thermal model was solved
  double time{ 0 };                                    //!< time since the last update of the thermal model [s]
};

//!< template <int DATASTORE_MODULE = DATASTORE_MODULE>
class Module : public StorageUnit //, public ModuleDataStorage<DATASTORE_MODULE>
{
protected:
  //!< connected child SUs
  using moduleSUs_t = std::vector<std::unique_ptr<StorageUnit>>;
  using moduleSUs_span_t = std::span<std::unique_ptr<StorageUnit>>;

  moduleSUs_t SUs;
  std::vector<double> Rcontact; //!< array with the contact resistance for cell i

  //!< thermal model
  std::unique_ptr<CoolSystem> cool{ nullptr }; //!< cooling system of this module //!< std::make_unique<CoolSystem>(settings::MODULE_NSUs_MAX, 1);
  ModuleThermalParam therm;

  //!< voltage
  size_t Ncells;               //!< Number of cells this module contains.
  double Vmodule{ 0 };         //!< voltage of the module
  bool Vmodule_valid{ false }; //!< boolean indicating if stored the voltage of the module is valid
  bool par;                    //!< if true, some functions will be calculated parallel using multithreaded computing
                               //!< data storage
#if DATASTORE_MODULE > 0
  CellCumulativeData tData;
#endif

#if DATASTORE_MODULE > 1
  std::vector<CommonData> cData; //!< Common data
#endif
  size_t calculateNcells()
  {
    /*  return the number of cells connected to this module
     * 	e.g. if this module has 3 child-modules, each with 2 cells.
     * 	then getNSUs = 3 but getNcells = 6
     */
    size_t r{ 0 };
    for (const auto &SU : SUs)
      r += SU->getNcells();

    Ncells = r; //!< #TODO this function needs to call parent to update their number of cells if they are module.

    return r;
  }

  double thermalModel_cell();
  double thermalModel_coupled(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim);

public:
  //!< common implementation for all base-modules
  inline size_t getNSUs() { return SUs.size(); } //!< note that these child-SUs can be modules themselves (or they can be cells)
  moduleSUs_t &getSUs() { return SUs; }

  virtual Status checkVoltage(double &v, bool print) noexcept override; //!< get the voltage and check if it is valid
  double getVhigh() override;                                           //!< return the voltage of the cell with the highest voltage
  double getVlow() override;                                            //!< return the voltage of the cell with the lowest voltage
  void getStates(getStates_t &s) override;                              //!< returns one long array with the states: getStates-array of each cell followed by the states of the module (Tmod)
                                                                        //	[s0 s1 s2 ... sn T], where s1 = cells[0].getStates()

  //	int Vstatus

  virtual bool validStates(bool print = true) override; //!< check if a state-array is valid for this module (uses setStates)
  virtual Status setStates(setStates_t s, bool checkV = true, bool print = true) override;
  //!< Set the states of this module to the given one
  //!< note: setStates is the master function to check if states and cells are valid
  //!< if checkV=true, then also the cell and module voltages are checked
  //!< the other functions just call setStates to check validity

  //!< void backupStates();  //!< Back-up states.
  //!< void restoreStates(); //!< restore backed-up states.

  //!< thermal model
  double T() override { return cool->T(); } //!< get the temperature of this module

  double getThermalSurface() override { return therm.A; };
  double getCoolingLoad(); //!< return the energy required to run the entire coolingsystem of this module and all its children
  double thermalModel(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim) override;

  //!< different implementation for series vs parallel modules
  virtual bool validSUs(moduleSUs_span_t c, bool print = true) = 0;         //!< check if the cells in this array are valid for this module
  virtual bool validSUs(bool print = true) { return validSUs(SUs, print); } //!< check if the cells in valid for this module

  virtual void setSUs(moduleSUs_span_t c, bool checkCells = true, bool print = true);
  //!< Sets the cells of this module. Checks module-constraints
  //!< does not check if the states of cells are valid, nor the voltages of the cells and module
  //!< it only checks whether the cells are ok for this module (same current if series, same voltage if parallel)
  //!< note on virtual functions, see StorageUnit.hpp

  CoolSystem *getCoolSystem() { return cool.get(); }

  void setRcontact(std::span<double> Rc) //!< #TODO if ok.
  {
    /*
     * Set the contact resistance of each cell.
     * In parallel modules, this is the 'horizontal' resistance through which all currents of the cells 'behind it' goes.
     */
    assert(Rc.size() == getNSUs());
    Rcontact.resize(Rc.size());
    std::copy(Rc.begin(), Rc.end(), Rcontact.begin());
    Vmodule_valid = false; //!< we are changing the resistance, so the stored voltage is no longer valid
  }

  virtual Module *copy() override = 0;

  void storeData() override
  {
    for (const auto &SU : SUs) //!< Tell all connected cells to store their data
      SU->storeData();

    //!< Store data for the coolsystem
    cool->storeData(getNcells());

#if (DATASTORE_MODULE > 1) //!< Store data of this module
    moduleData.push_back(ModuleData(ahtot, whtot, timetot, I(), V(), T()));
#endif
  }

  void writeData(const std::string &prefix) override
  {
    //!< Tell all connected cells to write their data

    for (const auto &SU : SUs)
      SU->writeData(prefix);

#if (DATASTORE_MODULE == 0)
    if constexpr (settings::printBool::printNonCrit)
      std::cout << "ERROR in Module::writeData, the settings in constants.hpp are forbidding from storing data.\n";
#endif

#if (DATASTORE_MODULE > 1)                                             //!< Write data for this module
    std::string name = prefix + "_" + getFullID() + "_ModuleData.csv"; //!< name of the file, start with the full hierarchy-ID to identify this cell

    //!< append the new data to the existing file
    std::ofstream file(name, std::ios_base::app);

    if (!file.is_open()) {
      std::cerr << "ERROR in Module::writeData, could not open file " << name << '\n';
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      throw 11;
    }

    for (const auto &data : this->moduleData)
      file << data.Itot << ',' << data.Vtot << ',' << data.Ttot << ','
           << data.Ahtot << ',' << data.Whtot << ',' << data.Timetot << '\n';

    file.close();

    this->moduleData.clear();
#endif

//!< Write data for the cooling system
#if 0 //!< (DATASTORE_COOL > 0)
			cool->writeData(prefix + '_' + getFullID()); //!< append the ID to the prefix since the cooling system does not have the ID of this module
#endif
  }

  void setBlockDegAndTherm(bool block)
  {
    for (auto &SU : SUs)
      SU->setBlockDegAndTherm(block);

    blockDegAndTherm = block;
  }

  void setT(double Tnew) override //!< set a module temperature
  {
    /*
     * Set the temperature of the coolant in this module.
     * Note that it does NOT immediately change the temperatures of the child-SUs.
     * That is done by the thermal model over time
     */
    cool->setT(Tnew);
  }

  auto getSUTemperature(size_t i) { return SUs[i]->T(); } //!< return an array with the temperatures of the children of this module. //!< Note that these can be both modules and cells

  size_t getNcells() override
  {
    /*  return the number of cells connected to this module
     * 	e.g. if this module has 3 child-modules, each with 2 cells.
     * 	then getNSUs = 3 but getNcells = 6
     */
    return Ncells;
  }

  double getThotSpot() override //!< get the maximum temperature of the cells or the module
  {
    //!< Return the temperature of the hottest element of the module.
    //!< Note that this will be the T of a cell, since child-modules will pass on this function to their cells
    double Thot = cool->T();
    for (const auto &SU : SUs)
      Thot = std::max(Thot, SU->getThotSpot());

    return Thot;
  }
};
} // namespace slide