/*
 * Cell.hpp
 *
 *  Created on: 22 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "cell_limits.hpp"
#include "../settings/settings.hpp"
#include "../StorageUnit.hpp"
#include "../types/Histogram.hpp"
#include "../types/data_storage/CellData.hpp"
#include "../types/Status.hpp"
#include "../utility/utility.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <span>

namespace slide {
class Cell : public StorageUnit
{
protected:
  double capNom{ 16 }; //!< capacity [Ah]

  CellData<settings::DATASTORE_CELL> cellData;

public:
  constexpr static CellLimits limits{ defaultCellLimits }; // #TODO make it changable.

  Cell() : StorageUnit("cell") {}

  Cell(const std::string &ID_) : StorageUnit(ID_) {}
  virtual ~Cell() = default;

  double Cap() const final override { return capNom; }
  void setCapacity(double capacity) { capNom = capacity; }

  constexpr double Vmin() const override { return limits.Vmin; }
  constexpr double VMIN() const override { return limits.VMIN; }
  constexpr double VMAX() const override { return limits.VMAX; }
  constexpr double Vmax() const override { return limits.Vmax; }
  constexpr double Tmax() { return limits.Tmax; }
  constexpr double Tmin() { return limits.Tmin; }

  //!< virtual double V(bool print = true) = 0;   //!< crit is an optional argument
  //!< virtual double getOCV(bool print = true) = 0; //!< crit is an optional argument

  double getVhigh() final { return V(); } //!< return the voltage of the cell with the highest voltage
  double getVlow() final { return V(); }  //!< return the voltage of the cell with the lowest voltage

  virtual double getThotSpot() override { return T(); }

  virtual Status checkCurrent(bool checkV, bool print) noexcept
  {
    double v;
    Status Vstatus = checkV ? checkVoltage(v, print) : Status::Success;
    //!< #TODO Current checking part is missing!
    return Vstatus;
  }

  virtual Status checkVoltage(double &v, bool print) noexcept override
  {
    /*
     * -2 	V < VMIN				discharging and outside safety range, cycling should be halted asap to avoid numerical errors
     * -1	VMIN <= V < Vmin		discharging and outside range of cell, but not yet safety limit
     * 0 	Vmin <= V <= Vmax 		valid range
     * 1 	Vmax < V <= VMAX		charging and ...
     * 2 	VMAX < V 				charging and ...
     *
     * note: allow a small margin on Vmin and Vmax for when we are doing a CV phase
     * then occasionally, V can be slightly above or below Vlimit
     */

    return free::check_voltage(v, *this);
  }

  size_t getNcells() override final { return 1; } //!< this is a single cell

  virtual Status setSOC(double SOCnew, bool checkV = true, bool print = true) = 0;
  virtual double SOC() = 0;

  //!< virtual int getNstates() = 0;

  //!< thermal model
  //!< virtual double getThermalSurface() = 0; //!< todo not implemented
  virtual double thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim) override
  {
    /*
     * Calculate the thermal model of this cell
     */
    //!< todo not implemented #TODO
    //!< need something similar as SPM cell (keep track of heat generation and time)
    //!< and then here you can solve the ODE
    return T();
  }

  //!< dataStorage
  virtual void storeData() override { cellData.storeData(*this); }                                  //!< Add another data point in the array.
  virtual void writeData(const std::string &prefix) override { cellData.writeData(*this, prefix); } // #TODO *this may be Cell not actual type.

  virtual ThroughputData getThroughputs() { return {}; }

  virtual std::array<double, 4> getVariations() const noexcept { return {}; } // #TODO will be deleted in future.
};

} // namespace slide