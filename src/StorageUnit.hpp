/*
 * StorageUnit.hpp
 *
 * parent class for all cells and modules, interface (ie only defines very general functions)
 *
 *  Created on: 9 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

#include "settings/settings.hpp"
#include "types/Status.hpp"
#include "types/Deep_ptr.hpp"
#include "utility/free_functions.hpp"

#include <string>
#include <cstdlib>
#include <vector>
#include <span>
#include <string_view>
#include <utility>

namespace slide {
class StorageUnit
{
protected:
  std::string ID{ "StorageUnit" };           //!< identification string
  StorageUnit *parent{ nullptr };            //!< pointer to the SU 'above' this one [e.g. the module to which a cell is connected]
  bool blockDegAndTherm{ false };            //!< if true, degradation and the thermal ODE are ignored
  using setStates_t = std::span<double> &;   //!< To pass states to read, non-expandable container.
  using getStates_t = std::vector<double> &; //!< To pass states to save, expandable container.
  using viewStates_t = std::span<double>;
  virtual size_t calculateNcells() { return 0; }

public:
  //!< basic getters and setters
  StorageUnit() = default;
  StorageUnit(std::string_view ID_) : ID(ID_) {}
  StorageUnit(std::string_view ID_, StorageUnit *parent_, bool blockDegAndTherm_)
    : ID(ID_), parent(parent_), blockDegAndTherm(blockDegAndTherm_) {}

  virtual ~StorageUnit() = default;
  const std::string &getID() { return ID; }
  void setID(std::string IDi) { ID = std::move(IDi); }

  //!< Return the full ID string, including the ID of the parent module
  virtual std::string getFullID() { return (parent != nullptr) ? parent->getFullID() + "_" + getID() : getID(); }

  virtual double Cap() const = 0; //

  auto *getParent() { return parent; };

  virtual double I() const = 0; //
  virtual double getRtot() = 0;
  virtual size_t getNcells() = 0; //!< return the number of single cells connected to this SU

  inline bool isCharging() { return I() < 0; }    //!< negative means charge.
  inline bool isDischarging() { return I() > 0; } //!< positive means discharge.

  virtual void getStates(getStates_t s) = 0;       //!< returns one long array with the states
  virtual viewStates_t viewStates() { return {}; } //!< Only for cells to see individual states.
  void setBlockDegAndTherm(bool block) { blockDegAndTherm = block; }
  virtual void setParent(StorageUnit *p) { parent = p; }                             //!< set the parent
  virtual Status setCurrent(double Inew, bool checkV = true, bool print = true) = 0; //
  virtual Status setVoltage(double Vnew, bool checkI = true, bool print = true)
  {
    return free::setVoltage_iterative(this, Vnew);
  }


  virtual Status setStates(setStates_t s, bool checkStates = true, bool print = true) = 0; //!< opposite of getStates, check the states are valid?

  virtual void backupStates() {}  //!< Back-up states.
  virtual void restoreStates() {} //!< restore backed-up states.

  //!< virtual int getNstates() = 0;
  //!< virtual double SOC() = 0;
  //!<  voltage
  virtual double getOCV() = 0;
  virtual double V() = 0;                                          //!< print is an optional argument
  virtual Status checkVoltage(double &v, bool print) noexcept = 0; //!< get the voltage and check if it is valid
  virtual double getVhigh() = 0;                                   //!< return the voltage of the cell with the highest voltage
  virtual double getVlow() = 0;                                    //!< return the voltage of the cell with the lowest voltage

  virtual double Vmin() const = 0; //!< lower voltage limit of the cell
  virtual double VMIN() const = 0; //!< safety cut off
  virtual double Vmax() const = 0; //!< upper voltage limit of the cell
  virtual double VMAX() const = 0; //!< safety cut off

  //!< thermal model
  virtual double T() = 0;
  virtual double getThotSpot() = 0;                                                                               //!< the T of the hottest element in the SU
  virtual double getThermalSurface() = 0;                                                                         //!< return the 'A' for the thermal model of this SU (Q = hA*dT)
  virtual double thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim) = 0; //!< calculate the thermal model of this SU
  virtual void setT(double Tnew) = 0;

  //!< functionality
  virtual bool validStates(bool print = true) = 0;        //!< checks if a state array is valid
  virtual StorageUnit *copy() = 0;                        //!< copy this SU to a new object
  virtual void timeStep_CC(double dt, int steps = 1) = 0; //!< take a number of time steps

  //!< Data collection of cycling data (I, V, T, etc. for every cell)
  virtual void storeData() = 0;
  virtual void writeData(const std::string &prefix) = 0;
};

// Free functions:
// template <typename T, typename... Args>
// auto make_SU(Args &&...args) { return Deep_ptr<StorageUnit>(new T(std::forward<Args>(args)...)); }

template <typename T, typename... Args>
auto make_SU(Args &&...args)
{
  return Deep_ptr<StorageUnit>(std::unique_ptr<T>(new T(std::forward<Args>(args)...)));
}

} // namespace slide
