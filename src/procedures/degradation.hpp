/*
 * Degradation.hpp
 *
 * Header file for the degradation simulations.
 * The exact degradation procedures are defined in the Cycler.
 * In the functions defined here, the functions from the Cycler are called various times with slightly different parameters (e.g. different temperatures, C rates, etc).
 * As such, we can simulate the effect the different parameters have on battery degradation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "Cell.hpp"
#include "settings/settings.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <memory>
//#include <thread>

//!< Definitions of custom datatypes:
struct CycleAgeingConfig
{
  double Vma{ 4.2 };
  double Vmi{ 2.7 };
  double Tc{ 45 };
  double Ccha{ 1 };
  double Cdis{ 1 };
  double SOCma{ 100 };
  double SOCmi{ 0 };

  CycleAgeingConfig(double Vma, double Vmi, double Tc, double Ccha, double Cdis, double SOCma, double SOCmi)
    : Vma(Vma), Vmi(Vmi), Tc(Tc), Ccha(Ccha), Cdis(Cdis), SOCma(SOCma), SOCmi(SOCmi) {}

  double Ti() const { return Tc + PhyConst::Kelvin; }
  std::string get_name(const std::string &pref) const
  {
    //!< Example output: pref + "T45_1C1D_SOC0-100";
    return pref + "T" + std::to_string((int)Tc) + "_" + std::to_string((int)Ccha) + "C" + std::to_string((int)Cdis) + "D" + "_" + "SOC" + std::to_string((int)SOCmi) + "-" + std::to_string((int)SOCma);
  }
};

struct CalendarAgeingConfig
{
  double V{ 4.2 };
  double Tc{ 5 };
  double SOC{ 100 };

  CalendarAgeingConfig(double V, double Tc, double SOC) : V(V), Tc(Tc), SOC(SOC) {}

  double Ti() const { return Tc + PhyConst::Kelvin; }
  std::string get_name(const std::string &pref) const
  {
    //!< Example output: pref + "Cal-T45_SOC100";
    return pref + "Cal-T" + std::to_string((int)Tc) + "_SOC" + std::to_string((int)SOC);
  }
};

struct ProfileAgeingConfig
{
  double Vma{ 4.2 };
  double Vmi{ 2.7 };
  double Tc{ 45 };
  double SOCma{ 100 };
  double SOCmi{ 0 };
  std::string csvName{ "Current Profile drive cycle HWFET.csv" };
  std::string namePrefix{ "prof-HWFET" };

  ProfileAgeingConfig(double Vma, double Vmi, double Tc, double SOCma, double SOCmi, std::string &csvName, std::string &namePrefix)
    : Vma(Vma), Vmi(Vmi), Tc(Tc), SOCma(SOCma), SOCmi(SOCmi), csvName(csvName), namePrefix(namePrefix) {}

  double Ti() const { return Tc + PhyConst::Kelvin; }
  std::string get_name(const std::string &pref) const
  {
    //!< Example output: pref + "T45_1C1D_SOC0-100";
    return pref + namePrefix + "-T" + std::to_string((int)Tc) + "_" + "SOC" + std::to_string((int)SOCmi) + "-" + std::to_string((int)SOCma);
  }
};

//!< Auxiliary functions for multi-threaded simulations
void Calendar_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose, //!< simulate one calendar ageing experiment
                  double V, double Ti, int Time, int mode, int timeCycleData, int timeCheck, struct checkUpProcedure &proc, std::string name);
void Cycle_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose, double Vma, double Vmi, //!< simulate one cycle ageing experiment
               double Ccha, bool CVcha, double Icutcha, double Cdis, bool CVdis, double Icutdis, double Ti, int timeCycleData, int nrCycles, int nrCap, struct checkUpProcedure &proc, std::string name);
void Profile_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose, std::string profName, int n, int limit, //!< simulate one drive cycle ageing experiment
                 double Vma, double Vmi, double Ti, int timeCycleData, int nrProfiles, int nrCap, struct checkUpProcedure &proc, std::string name);

//!< Degradation experiments
void CycleAgeing(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, int verbose);    //!< simulate a range of cycle ageing experiments (different temperatures, SOC windows, currents)
void CalendarAgeing(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, int verbose); //!< simulate a range of calendar ageing experiments (different temperatures, SOC levels)
void ProfileAgeing(const struct slide::Model_SPM &M, std::string pref, const struct DEG_ID &degid, int cellType, int verbose);  //!< simulate a range of drive cycle experiments (different cycles, different temperatures, etc.)

//!< Configuration struct for above-given functions.

void Cycle_one(const struct slide::Model_SPM &M, const struct DEG_ID &degid, int cellType, int verbose, //!< simulate one cycle ageing experiment
               const struct CycleAgeingConfig &cycAgConfig, bool CVcha, double Icutcha, bool CVdis, double Icutdis, int timeCycleData, int nrCycles, int nrCap, struct checkUpProcedure &proc, const std::string &pref);
