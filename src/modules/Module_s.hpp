/*
 * Module_base_s.hpp
 *
 * series-connected base Module
 *
 *  Created on: 9 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <string_view>
#include <memory>

#include "Module.hpp"
#include "../settings/settings.hpp"
#include "../cells/Cell.hpp"
#include "../types/data_storage/cell_data.hpp"

namespace slide {
class Module_s : public Module //<DATASTORE_MODULE>
{
protected:
#if TIMING //!< std::clock_t tstart;
  TimingData_Module_s timeData{};
#endif
public:
  Module_s();
  Module_s(std::string_view IDi, double Ti, bool print, bool pari, int Ncells, int coolControl, int cooltype);
  //!< Module_s(std::string_view IDi, bool pari, int Ncells, std::unique_ptr<CoolSystem> &&cool_); #TODO

  //!< functions from Module_base
  double Cap() override;  //!< module capacity (sum of cells)
  double Vmin() override; //!< module capacity (sum of cells)
  double VMIN() override;
  double Vmax() override; //!< module capacity (sum of cells)
  double VMAX() override;
  double I() override;                       //!< module capacity (sum of cells)
  double getOCV(bool print = true) override; //!< module voltage (sum of cells), print is an optional argument
  double getRtot() override;
  double V(bool print = true) override; //!< module voltage (sum of cells), print is an optional argument

  //!< bool validSUs(bool print = true);
  bool validSUs(moduleSUs_span_t c, bool print = true) override; //!< check if the cells in this array are valid for this module

  Status setCurrent(double Inew, bool checkV = true, bool print = true) override; //!< set a module current
  void timeStep_CC(double dt, bool addData = false, int steps = 1) override;

  Module_s *copy() override;
  TimingData_Module_s getTimings();
  void setTimings(TimingData_Module_s td);
};
} // namespace slide
