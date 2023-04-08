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
class Module_s : public Module
{
public:
  Module_s() : Module("moduleS") {} //!< note this constructor should never be used. It can't determine which coolsystem to use
  Module_s(std::string_view ID_, double Ti, bool print, bool pari, int Ncells_, int coolControl, int cooltype)
    : Module(ID_, Ti, print, pari, Ncells_, coolControl, cooltype) {}
  //!< Module_s(std::string_view IDi, bool pari, int Ncells, std::unique_ptr<CoolSystem> &&cool_); #TODO

  //!< functions from Module_base
  // Since cells are in series following functions are just sum.
  double Vmin() const override { return transform_sum(SUs, free::get_Vmin<SU_t>); } //!< SUM all SUs' Vmin
  double VMIN() const override { return transform_sum(SUs, free::get_VMIN<SU_t>); } //!< SUM all SUs' VMIN
  double Vmax() const override { return transform_sum(SUs, free::get_Vmax<SU_t>); } //!< SUM all SUs' Vmax
  double VMAX() const override { return transform_sum(SUs, free::get_VMAX<SU_t>); } //!< SUM all SUs' VMAX

  double I() const override { return (SUs.empty() ? 0 : SUs[0]->I()); }           //!< the current is the same in all cells, so we can use first value. Return zero if SUs empty.
  double Cap() const override { return transform_min(SUs, free::get_Cap<SU_t>); } //!< module capacity is the capacity of the smallest cell


  double getOCV() override; //!< module voltage (sum of cells), print is an optional argument
  double getRtot() override;
  double V() override; //!< module voltage (sum of cells), print is an optional argument

  Status setCurrent(double Inew, bool checkV = true, bool print = true) override; //!< set a module current
  void timeStep_CC(double dt, int steps = 1) override;

  Module_s *copy() override { return new Module_s(*this); }
};
} // namespace slide
