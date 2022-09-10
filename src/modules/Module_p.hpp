/*
 * Module_base_p.hpp
 *
 *  Created on: 18 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "Module.hpp"
#include "../types/data_storage/cell_data.hpp"

#include <string>

namespace slide {
class Module_p : public Module
{
protected:
#if TIMING //!< std::clock_t tstart;
  TimingData_Module_p timeData{};
#endif
public:
  Module_p();
  Module_p(std::string IDi, double Ti, bool print, bool pari, int Ncells, int coolControl, int cooltype);

  //!< functions from Module_base
  double Cap() override;  //!< module capacity (sum of cells)
  double Vmin() override; //!< module capacity (sum of cells)
  double VMIN() override;
  double VMAX() override;
  double Vmax() override;                    //!< module capacity (sum of cells)
  double I() override;                       //!< module capacity (sum of cells)
  double getOCV(bool print = true) override; //!< module voltage (sum of cells), print is an optional argument
  double getRtot() override;
  double V(bool print = true) override;      //!< module voltage (sum of cells), print is an optional argument
  double getVi(size_t i, bool print = true); //!< get the voltage of SU[i] while accounting for the contact resistance

  Status setI_iterative(double Inew, bool checkV = true, bool print = true);      //!< set a module current iteratively
  Status setCurrent(double Inew, bool checkV = true, bool print = true) override; //!< set a module current
  Status redistributeCurrent(bool checkV = true, bool print = true);              //!< redistribute the total module current to the different cells

  bool validSUs(moduleSUs_span_t c, bool print = true) override; //!< check if the cells in this array are valid for this module

  void timeStep_CC(double dt, int steps = 1) override;

  Module_p *copy() override;
  TimingData_Module_p getTimings();
  void setTimings(TimingData_Module_p td);
};
} // namespace slide
