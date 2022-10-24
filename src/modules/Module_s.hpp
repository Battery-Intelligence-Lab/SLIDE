/*
 * Module_base_s.hpp
 *
 * series-connected base Module
 *
 *  Created on: 9 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "Module.hpp"
#include "../settings/settings.hpp"
#include "../cells/Cell.hpp"
#include "../types/data_storage/cell_data.hpp"
#include "../utility/utility.hpp"

#include <string_view>
#include <memory>

namespace slide {
class Module_s : public Module //<DATASTORE_MODULE>
{
protected:
#if TIMING //!< std::clock_t tstart;
  TimingData_Module_s timeData{};
#endif
public:
  Module_s() : Module("moduleS") {} //!< note this constructor should never be used. It can't determine which coolsystem to use
  Module_s(std::string_view ID_, double Ti, bool print, bool pari, int Ncells, int coolControl, int cooltype)
    : Module(ID_, Ti, print, pari, Ncells, coolControl, cooltype) {}
  //!< Module_s(std::string_view IDi, bool pari, int Ncells, std::unique_ptr<CoolSystem> &&cool_); #TODO

  //!< functions from Module_base
  // Since cells are in series following functions are just sum.
  double Vmin() override { return transform_sum(SUs, free::get_Vmin<SU_t>); } //!< SUM all SUs' Vmin
  double VMIN() override { return transform_sum(SUs, free::get_VMIN<SU_t>); } //!< SUM all SUs' VMIN
  double Vmax() override { return transform_sum(SUs, free::get_Vmax<SU_t>); } //!< SUM all SUs' Vmax
  double VMAX() override { return transform_sum(SUs, free::get_VMAX<SU_t>); } //!< SUM all SUs' VMAX

  double I() override { return (SUs.empty() ? 0 : SUs[0]->I()); } //!< the current is the same in all cells, so we can use first value. Return zero if SUs empty.

  double getOCV(bool print = true) override; //!< module voltage (sum of cells), print is an optional argument
  double getRtot() override;
  double V(bool print = true) override; //!< module voltage (sum of cells), print is an optional argument

  //!< bool validSUs(bool print = true);
  bool validSUs(SUs_span_t c, bool print = true) override; //!< check if the cells in this array are valid for this module

  Status setCurrent(double Inew, bool checkV = true, bool print = true) override; //!< set a module current
  void timeStep_CC(double dt, int steps = 1) override;

  Module_s *copy() override;
  TimingData_Module_s getTimings();
  void setTimings(TimingData_Module_s td);
};
} // namespace slide
