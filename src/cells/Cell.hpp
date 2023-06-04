/**
 * @file Cell.hpp
 * @brief Cell class definition
 * @author Jorn Reniers, Volkan Kumtepeli
 * @date 22 Nov 2019
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

/**
 * @brief Abstract Class representing a single battery cell.
 */
class Cell : public StorageUnit
{
protected:
  double capNom{ 16 }; //!< capacity [Ah].

  CellData<settings::DATASTORE_CELL> cellData; //!< Cell data storage.

public:
  constexpr static CellLimits limits{ defaultCellLimits }; // Default cell limits. #TODO make it changable.

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

  double getVhigh() final { return V(); } //!< return the voltage of the cell with the highest voltage
  double getVlow() final { return V(); }  //!< return the voltage of the cell with the lowest voltage

  virtual Status setSOC(double SOCnew, bool checkV = true, bool print = true) = 0;
  virtual double SOC() = 0;
  virtual double getThotSpot() override { return T(); }
  size_t getNcells() override final { return 1; } //!< this is a single cell

  virtual Status checkCurrent(bool checkV, bool print) noexcept
  {
    double v;
    Status Vstatus = checkV ? checkVoltage(v, print) : Status::Success;
    //!< #TODO Current checking part is missing!
    return Vstatus;
  }

  /**
   * @brief Check the voltage status of the cell.
   * @param v Reference to the voltage value.
   * @param print Boolean flag to indicate whether to print the result or not.
   * @return Status of the voltage check.
   */
  virtual Status checkVoltage(double &v, bool print) noexcept override { return free::check_voltage(v, *this); }


  //!< thermal model
  //!< virtual double getThermalSurface() = 0; //!< todo not implemented

  /**
   * @brief Calculate the thermal model of the cell.
   * @param Nneighb Number of neighboring cells.
   * @param Tneighb Pointer to an array of neighboring cells' temperatures.
   * @param Kneighb Pointer to an array of neighboring cells' thermal conductivities.
   * @param Aneighb Pointer to an array of neighboring cells' contact areas.
   * @param tim Time duration for the thermal model calculation.
   * @return Cell temperature after thermal model calculation.
   * @note The thermal model is not implemented yet. Requires something similar to an SPM cell
   * (keep track of heat generation and time) and then solve the ODE here.
   */
  virtual double thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim) override
  {
    return T();
  }

  //!< dataStorage
  virtual void storeData() override { cellData.storeData(*this); }                                  //!< Add another data point in the array.
  virtual void writeData(const std::string &prefix) override { cellData.writeData(*this, prefix); } // #TODO *this may be Cell not actual type.

  virtual ThroughputData getThroughputs() { return {}; }

  virtual std::array<double, 4> getVariations() const noexcept { return {}; } // #TODO will be deleted in future.
};

} // namespace slide