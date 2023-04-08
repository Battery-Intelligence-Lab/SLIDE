/*
 * Cycler.hpp
 *
 *  Created on: 19 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../StorageUnit.hpp"
#include "../types/data_storage/cell_data.hpp"

#include <string>
#include <memory>
#include <span>

namespace slide {
class Cycler
{
private:
  std::string ID{ "Cycler" }; //!< identification string of this cycler
  StorageUnit *su{ nullptr }; //!< (pointer to) something of type storage unit which is connected to this cycler

  size_t index{ 0 };        //!< Cycler should keep its on index for data writing.
  bool diagnostic{ false }; //!< are we running in diagnostic mode or not?

  //!< secondary functions
  Status setCurrent(double I, double vlim); //!< sets the current to the connected SU

public:
  Cycler() = default;
  Cycler(StorageUnit *sui) : su(sui) {}
  Cycler(StorageUnit *sui, const std::string &IDi) : ID(IDi), su(sui) {}
  Cycler(Deep_ptr<StorageUnit> &sui, const std::string &IDi) : Cycler(sui.get(), IDi) {}

  Cycler &initialise(StorageUnit *sui, const std::string &IDi);
  Cycler &initialise(Deep_ptr<StorageUnit> &sui, const std::string &IDi) { return initialise(sui.get(), IDi); }

  Cycler &setDiagnostic(bool newDia);
  double getSafetyVmin() { return su->VMIN() * 0.99; } //!< #TODO probably causing many calculations.
  double getSafetyVmax() { return su->VMAX() * 1.01; }

  Status rest(double tlim, double dt, int ndt_data, ThroughputData &th);
  Status CC(double I, double vlim, double tlim, double dt, int ndt_data, ThroughputData &th);
  Status CV(double Vset, double Ilim, double tlim, double dt, int ndt_data, ThroughputData &th);
  Status CCCV(double I, double Vset, double Ilim, double dt, int ndt_data, ThroughputData &th);
  Status CCCV_with_tlim(double I, double Vset, double Ilim, double tlim, double dt, int ndt_data, ThroughputData &th);

  Status Profile(std::span<double> I_vec, double vlim, double tlim, double dt, int ndt_data, double &Ah, double &Wh);

  int storeData();
  int writeData();

  double testCapacity(double &Ah, double &ttot);
  double testCapacity()
  {
    double Ah{}, ttot{};
    return testCapacity(Ah, ttot);
  }


  StorageUnit *getSU() { return su; }
};
} // namespace slide