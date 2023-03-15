/*
 * Module.hpp
 *
 *  Created on: 29 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../StorageUnit.hpp"
#include "../cooling/cooling.hpp"
#include "../types/State.hpp"
#include "../settings/settings.hpp"
#include "../utility/utility.hpp"


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

class Module : public StorageUnit
{

public:
  //!< connected child SUs
  using SU_t = Deep_ptr<StorageUnit>; // #TODO in future it should store directly storage unit itself.
  using SUs_t = std::vector<SU_t>;
  using SUs_span_t = std::span<SU_t>;
  using CoolSystem_t = Deep_ptr<CoolSystem>;

protected:
  SUs_t SUs;
  std::vector<double> Rcontact; //!< array with the contact resistance for cell i

  //!< thermal model
  CoolSystem_t cool{ nullptr }; //!< cooling system of this module //!< make<CoolSystem>(settings::MODULE_NSUs_MAX, 1);
  ModuleThermalParam therm;

  //!< voltage
  size_t Ncells;               //!< Number of cells this module contains.
  double Vmodule{ 0 };         //!< voltage of the module
  bool Vmodule_valid{ false }; //!< boolean indicating if stored the voltage of the module is valid
  bool par{ true };            //!< if true, some functions will be calculated parallel using multithreaded computing
                               //!< data storage

  State<0, settings::data::N_CumulativeModule> st_module;
  std::vector<double> data; //!< Time data


  size_t calculateNcells() override
  {
    /*  return the number of cells connected to this module
     * 	e.g. if this module has 3 child-modules, each with 2 cells.
     * 	then getNSUs = 3 but getNcells = 6
     */
    //!< #TODO this function needs to call parent to update their number of cells if they are module.
    return Ncells = transform_sum(SUs, free::get_Ncells<SU_t>);
  }

  double thermalModel_cell();
  double thermalModel_coupled(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim);

public:
  Module() : StorageUnit("Module") {}
  Module(std::string_view ID_) : StorageUnit(ID_) {}
  Module(std::string_view ID_, double Ti, bool print, bool pari, int Ncells, int coolControl, int cooltype);
  Module(std::string_view ID_, double Ti, bool print, bool pari, int Ncells, CoolSystem_t &&coolControlPtr, int cooltype);

  //!< common implementation for all base-modules
  size_t getNSUs() { return SUs.size(); } //!< note that these child-SUs can be modules themselves (or they can be cells)
  SUs_t &getSUs() { return SUs; }

  const auto &operator[](size_t i) const { return SUs[i]; }
  auto &operator[](size_t i) { return SUs[i]; }

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

  //!< thermal model
  double T() override { return cool->T(); } //!< get the temperature of this module

  double getThermalSurface() override { return therm.A; };
  double getCoolingLoad(); //!< return the energy required to run the entire coolingsystem of this module and all its children
  double thermalModel(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim) override;

  virtual void setSUs(SUs_span_t c, bool checkCells = true, bool print = true);
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

    if constexpr (settings::DATASTORE_MODULE >= settings::moduleDataStorageLevel::storeTimeData) //!< Store data of this module
      data.insert(data.end(), { st_module.Ah(), st_module.Wh(), st_module.time(), I(), V(), T() });
  }

  void writeData(const std::string &prefix) override
  {
    for (const auto &SU : SUs) //!< Tell all connected cells to write their data
      SU->writeData(prefix);


    if constexpr (settings::DATASTORE_MODULE >= settings::moduleDataStorageLevel::storeTimeData) //!< Write data for this module
    {
      std::string name = prefix + "_" + getFullID() + "_ModuleData.csv"; //!< name of the file, start with the full hierarchy-ID to identify this cell

      //!< append the new data to the existing file
      std::ofstream file(name, std::ios_base::app); // #TODO if first time open + header, if not append.

      if (!file.is_open()) {
        std::cerr << "ERROR in Module::writeData, could not open file " << name << '\n';
        throw 11;
      }

      free::write_data(file, data, 6);

      file.close();
    }

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

  /*****************************************************************
   * Set the temperature of the coolant in this module.
   * Note that it does NOT immediately change the temperatures of the child-SUs.
   * That is done by the thermal model over time
   *****************************************************************/
  void setT(double Tnew) override { cool->setT(Tnew); } //!< set a module temperature

  /*****************************************************************
   * return the number of cells connected to this module
   * e.g. if this module has 3 child-modules, each with 2 cells.
   * then getNSUs = 3 but getNcells = 6
   *****************************************************************/
  size_t getNcells() override { return Ncells; }

  //!< Return the temperature of the hottest element of the module.
  //!< Note that this will be the T of a cell, since child-modules will pass on
  //!< this function to their cells
  double getThotSpot() override //!< get the maximum temperature of the cells or the module
  {
    double Thot = cool->T();
    for (const auto &SU : SUs)
      Thot = std::max(Thot, SU->getThotSpot());

    return Thot;
  }
};
} // namespace slide