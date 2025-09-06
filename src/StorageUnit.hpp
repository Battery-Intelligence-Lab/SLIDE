/**
 * @file StorageUnit.hpp
 * @brief Abstract base class for all cells and modules
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 9 Dec 2019
 *
 * This is the parent class for all cells and modules, providing a common interface
 * that defines very general functions for battery storage units.
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

/**
 * @brief Abstract base class for all storage units (cells and modules)
 *
 * This class provides the common interface for all battery storage units,
 * including cells and modules. It defines virtual functions that must be
 * implemented by derived classes for electrical, thermal, and data management
 * operations.
 */
class StorageUnit
{
protected:
  std::string ID{ "StorageUnit" };             //!< identification string
  StorageUnit *parent{ nullptr };              //!< pointer to the SU 'above' this one [e.g. the module to which a cell is connected]
  bool blockDegAndTherm{ false };              //!< if true, degradation and the thermal ODE are ignored
  using setStates_t = std::span<const double>; //!< To pass states to read, non-expandable container.
  using getStates_t = std::vector<double> &;   //!< To pass states to save, expandable container.
  using viewStates_t = std::span<double>;
  virtual size_t calculateNcells() { return 0; }

public:
  /**
   * @brief Default constructor
   */
  StorageUnit() = default;

  /**
   * @brief Constructor with ID
   * @param ID_ Identification string for this storage unit
   */
  StorageUnit(std::string_view ID_) : ID(ID_) {}

  /**
   * @brief Constructor with ID, parent, and degradation/thermal blocking
   * @param ID_ Identification string for this storage unit
   * @param parent_ Pointer to parent storage unit
   * @param blockDegAndTherm_ Whether to block degradation and thermal calculations
   */
  StorageUnit(std::string_view ID_, StorageUnit *parent_, bool blockDegAndTherm_)
    : ID(ID_), parent(parent_), blockDegAndTherm(blockDegAndTherm_) {}

  /**
   * @brief Virtual destructor
   */
  virtual ~StorageUnit() = default;

  /**
   * @brief Get the ID of this storage unit
   * @return Reference to the ID string
   */
  const std::string &getID() { return ID; }

  /**
   * @brief Set the ID of this storage unit
   * @param IDi New ID string
   */
  void setID(std::string IDi) { ID = std::move(IDi); }

  /**
   * @brief Get the full hierarchical ID string
   * @return Full ID including parent module IDs (e.g., "module1_cell2")
   */
  virtual std::string getFullID() { return (parent != nullptr) ? parent->getFullID() + "_" + getID() : getID(); }

  /**
   * @brief Get the capacity of this storage unit
   * @return Capacity in Ah
   */
  virtual double Cap() const = 0;

  /**
   * @brief Get the parent storage unit
   * @return Pointer to parent storage unit
   */
  auto *getParent() { return parent; };

  /**
   * @brief Get the current through this storage unit
   * @return Current in A (negative for charging, positive for discharging)
   */
  virtual double I() const = 0;

  /**
   * @brief Get the total resistance of this storage unit
   * @return Total resistance in Ohms
   */
  virtual double getRtot() = 0;

  /**
   * @brief Get the number of single cells connected to this storage unit
   * @return Number of cells
   */
  virtual size_t getNcells() = 0;

  /**
   * @brief Check if the storage unit is charging
   * @return True if charging (negative current)
   */
  inline bool isCharging() { return I() < 0; }

  /**
   * @brief Check if the storage unit is discharging
   * @return True if discharging (positive current)
   */
  inline bool isDischarging() { return I() > 0; }

  /**
   * @brief Get all states of this storage unit
   * @param s Vector to store the states
   */
  virtual void getStates(getStates_t s) = 0;

  /**
   * @brief Get a view of the states (only for cells)
   * @return Span view of the states
   */
  virtual viewStates_t viewStates() { return {}; }

  /**
   * @brief Set whether to block degradation and thermal calculations
   * @param block True to block degradation and thermal calculations
   */
  void setBlockDegAndTherm(bool block) { blockDegAndTherm = block; }

  /**
   * @brief Set the parent storage unit
   * @param p Pointer to parent storage unit
   */
  virtual void setParent(StorageUnit *p) { parent = p; }

  /**
   * @brief Set the current through this storage unit
   * @param Inew New current value in A
   * @param checkV Whether to check voltage limits
   * @param print Whether to print status messages
   * @return Status of the operation
   */
  virtual Status setCurrent(double Inew, bool checkV = true, bool print = true) = 0;

  /**
   * @brief Set the voltage of this storage unit
   * @param Vnew New voltage value in V
   * @param checkI Whether to check current limits
   * @param print Whether to print status messages
   * @return Status of the operation
   */
  virtual Status setVoltage(double Vnew, bool checkI = true, bool print = true) { return free::setVoltage_iterative(this, Vnew); }

  /**
   * @brief Set the states of this storage unit
   * @param s Span containing the new states
   * @param n Number of states processed
   * @param checkStates Whether to validate the states
   * @param print Whether to print status messages
   * @return Status of the operation
   */
  virtual Status setStates(setStates_t s, int &n, bool checkStates = true, bool print = true) = 0;

  /**
   * @brief Get the time derivatives of the states
   * @param s Vector to store the derivatives
   */
  virtual void get_dxdt(getStates_t s) {}

  /**
   * @brief Backup the current states
   */
  virtual void backupStates() {}

  /**
   * @brief Restore the backed-up states
   */
  virtual void restoreStates() {}

  /**
   * @brief Get the open circuit voltage
   * @return Open circuit voltage in V
   */
  virtual double getOCV() = 0;

  /**
   * @brief Get the terminal voltage
   * @return Terminal voltage in V
   */
  virtual double V() = 0;

  /**
   * @brief Check if the voltage is within valid limits
   * @param v Reference to store the voltage value
   * @param print Whether to print status messages
   * @return Status of the voltage check
   */
  virtual Status checkVoltage(double &v, bool print) noexcept = 0;

  /**
   * @brief Get the highest voltage among all cells
   * @return Highest voltage in V
   */
  virtual double getVhigh() = 0;

  /**
   * @brief Get the lowest voltage among all cells
   * @return Lowest voltage in V
   */
  virtual double getVlow() = 0;

  /**
   * @brief Get the minimum operating voltage
   * @return Minimum voltage in V
   */
  virtual double Vmin() const = 0;

  /**
   * @brief Get the minimum safety voltage
   * @return Minimum safety voltage in V
   */
  virtual double VMIN() const = 0;

  /**
   * @brief Get the maximum operating voltage
   * @return Maximum voltage in V
   */
  virtual double Vmax() const = 0;

  /**
   * @brief Get the maximum safety voltage
   * @return Maximum safety voltage in V
   */
  virtual double VMAX() const = 0;

  /**
   * @brief Get the temperature of this storage unit
   * @return Temperature in K
   */
  virtual double T() = 0;

  /**
   * @brief Get the temperature of the hottest element
   * @return Hot spot temperature in K
   */
  virtual double getThotSpot() = 0;

  /**
   * @brief Get the thermal surface area for heat transfer
   * @return Surface area in m²
   */
  virtual double getThermalSurface() = 0;

  /**
   * @brief Calculate the thermal model
   * @param Nneighb Number of neighboring elements
   * @param Tneighb Array of neighboring temperatures
   * @param Kneighb Array of thermal conductivities
   * @param Aneighb Array of contact areas
   * @param tim Time step
   * @return Temperature after thermal calculation
   */
  virtual double thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim) = 0;

  /**
   * @brief Set the temperature of this storage unit
   * @param Tnew New temperature in K
   */
  virtual void setT(double Tnew) = 0;

  /**
   * @brief Check if the current states are valid
   * @param print Whether to print error messages
   * @return True if states are valid
   */
  virtual bool validStates(bool print = true) = 0;

  /**
   * @brief Create a copy of this storage unit
   * @return Pointer to the copied storage unit
   */
  virtual StorageUnit *copy() = 0;

  /**
   * @brief Take a time step with constant current
   * @param dt Time step in s
   * @param steps Number of steps to take
   */
  virtual void timeStep_CC(double dt, int steps = 1) = 0;

  /**
   * @brief Store current data for this storage unit
   */
  virtual void storeData() = 0;

  /**
   * @brief Write stored data to file
   * @param prefix Prefix for the output filename
   */
  virtual void writeData(const std::string &prefix) = 0;
};

/**
 * @brief Factory function to create a storage unit with deep pointer management
 *
 * @tparam T Type of storage unit to create (must inherit from StorageUnit)
 * @tparam Args Types of constructor arguments
 * @param args Constructor arguments for the storage unit
 * @return Deep_ptr<StorageUnit> containing the created storage unit
 *
 * This function creates a new storage unit of type T and wraps it in a Deep_ptr
 * for automatic memory management. The created object is owned by the Deep_ptr.
 */
template <typename T, typename... Args>
auto make_SU(Args &&...args)
{
  return Deep_ptr<StorageUnit>(std::unique_ptr<T>(new T(std::forward<Args>(args)...)));
}

} // namespace slide
