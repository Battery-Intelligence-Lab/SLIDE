/*
 * Cycler.hpp
 *
 *  Created on: 19 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "../StorageUnit.hpp"

#include <string>
#include <memory>
#include <span>

namespace slide {
class Cycler
{
private:
  std::string ID{ "ArbitraryCycler" }; //!< identification string of this cycler
  StorageUnit *su{ nullptr };          //!< (pointer to) something of type storage unit which is connected to this cycler

  size_t index{ 0 };        //!< Cycler should keep its on index for data writing.
  bool diagnostic{ false }; //!< are we running in diagnostic mode or not?

  //!< secondary functions
  Status setCurrent(double I, double vlim); //!< sets the current to the connected SU

public:
  Cycler() = default;
  Cycler(StorageUnit *sui, const std::string &IDi) : ID(IDi), su(sui) {}
  Cycler(std::unique_ptr<StorageUnit> &sui, const std::string &IDi) : Cycler(sui.get(), IDi) {}

  void initialise(StorageUnit *sui, const std::string &IDi);
  void initialise(std::unique_ptr<StorageUnit> &sui, const std::string &IDi) { initialise(sui.get(), IDi); }

  void setDiagnostic(bool newDia);
  double getSafetyVmin() { return su->VMIN() * 0.99; } //!< #CHECK probably causing many calculations.
  double getSafetyVmax() { return su->VMAX() * 1.01; }

  Status rest(double tlim, double dt, int ndt_data, double &Ah, double &Wh);
  Status CC(double I, double vlim, double tlim, double dt, int ndt_data, double &Ah, double &Wh, double &ttot);
  Status CV(double Vset, double Ilim, double tlim, double dt, int ndt_data, double &Ah, double &Wh, double &ttot);
  Status CCCV(double I, double Vset, double Ilim, double dt, int ndt_data, double &Ah, double &Wh, double &ttot);
  Status CCCV_with_tlim(double I, double Vset, double Ilim, double tlim, double dt, int ndt_data, double &Ah, double &Wh, double &ttot);

  Status Profile(std::span<double> I_vec, double vlim, double tlim, double dt, int ndt_data, double &Ah, double &Wh);

  int storeData();
  int writeData();

  double testCapacity(double &Ah, double &ttot);

  StorageUnit *getSU() { return su; }
};
} // namespace slide