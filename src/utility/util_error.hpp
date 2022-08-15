/*
 * Util_error.hpp
 *
 * Utility functions for error. So the code is less verbose.
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "../cells/Cell.hpp"

class Cell;

namespace slide::util::error {
void checkInputParam_CalAge(Cell &c, double V, double Ti, int Time, int timeCheck, int mode);
void checkInputParam_CC_V_CV_I(Cell &c, double Crate, double Vset, double Ccut);
void checkInputParam_CycAge(Cell &c, double Vma, double Vmi, double Ccha, double Ccutcha,
                            double Cdis, double Ccutdis, double Ti, int nrCycles, int nrCap);

inline void checkInputParam_CalAge(Cell &c, double V, double Ti, int Time, int timeCheck, int mode)
{

  //!< Check the input parameters
  bool vmax = V > c.Vmax(); //!< check if the maximum voltage is below the cell maximum voltage
  if (vmax)
    std::cerr << "Error in Cycler::CalendarAgeing. The voltage " << V << " is too high. The maximum value is " << c.Vmax() << ".\n";

  bool vmin = V < c.Vmin(); //!< check if the minimum voltage is above the cell minimum voltage
  if (vmin)
    std::cerr << "Error in Cycler::CalendarAgeing. The voltage " << V << " is too low. The minimum value is " << c.Vmin() << ".\n";

  bool Temin = Ti < settings::Tmin_K; //!< check the temperature is above 0 degrees, TMIN is defined in State.hpp
  if (Temin)
    std::cerr << "Error in Cycler::CalendarAgeing. The temperature " << Ti << "K is too low. The minimum value is 273.\n";

  bool Temax = Ti > settings::Tmax_K; //!< check the temperature is below 60 degrees, TMAX is defined in State.hpp
  if (Temax)
    std::cerr << "Error in Cycler::CalendarAgeing. The temperature " << Ti << " is too high. The maximum value is (273+60).\n";

  bool mod = (mode < 0 || mode > 2); //!< check the setting in 'mode' is allowed
  if (mod)
    std::cerr << "Error in Cycler::CalendarAgeing. The resting mode " << mode << " is illegal. Only values of 0, 1 or 2 are allowed.\n";

  bool cycles = (std::fmod(Time, timeCheck) > 1); //!< check the total resting time is a multiple of the the time between two check-ups
  if (cycles)                                     //!< allow a margin of 1 day (for numerical errors)
    std::cerr << "Error in Cycler::CalendarAgeing. The total resting time " << Time
              << " is not a multiple of the time between two check ups " << timeCheck << ".\n";

  if (vmax || vmin || Temin || Temax || mod || cycles)
    throw 1014;

  //!< print a warning if mode is 2 because it will take very long to simulate
  if (mode == 2)
    std::cout << "Warning in Cycler::CalendarAgeing. You have chosen to float the cell at a constant voltage for a long period. "
                 "It will take long to simulate this.\n";
}

inline void checkInputParam_CC_V_CV_I(Cell &c, double Crate, double Vset, double Ccut)
{

  //!< check the voltage which should be reached
  bool vmax = Vset > c.Vmax(); //!< check if the maximum voltage is below the cell maximum voltage
  if (vmax)
    std::cerr << "Error in BasicCycler::CC_V_CV_I. The voltage " << Vset << " is too high. The maximum value is " << c.Vmax() << ".\n";

  bool vmin = Vset < c.Vmin(); //!< check if the minimum voltage is above the cell minimum voltage
  if (vmin)
    std::cerr << "Error in BasicCycler::CC_V_CV_I. The voltage " << Vset << " is too low. The minimum value is " << c.Vmin() << ".\n";

  if (vmax || vmin)
    throw 1005;

  //!< check the current threshold for the CV phase
  bool currCV = Ccut < 0;
  if (currCV) {
    std::cerr << "Error in BasicCycler::CC_V_CV_I. The cutoff C rate " << Ccut << " is negative. It must be positive.\n";

    throw 1008;
  }

  //!< check the C rate for the current during the CC phase
  bool currCC = Crate < 0;
  if (currCC) {
    std::cerr << "Error in BasicCycler::CC_V_CV_I. The Crate " << Crate << " is negative. It must be positive.\n";
    throw 1010;
  }
}

inline void checkInputParam_CycAge(Cell &c, double Vma, double Vmi, double Ccha, double Ccutcha,
                                   double Cdis, double Ccutdis, double Ti, int nrCycles, int nrCap)
{
  //!< Check the input parameters
  bool vmax = Vma > c.Vmax(); //!< check if the maximum voltage is below the cell maximum voltage
  if (vmax)
    std::cerr << "Error in Cycler::cycleAgeing. The maximum voltage " << Vma << " is too high. The maximum value is " << c.Vmax() << ".\n";

  bool vmin = Vmi < c.Vmin(); //!< check if the minimum voltage is above the cell minimum voltage
  if (vmin)
    std::cerr << "Error in Cycler::cycleAgeing. The minimum voltage " << Vmi << " is too low. The minimum value is " << c.Vmin() << ".\n";

  bool Temin = Ti < settings::Tmin_Cell_K; //!< check the temperature is above 0 degrees, TMIN is defined in State.hpp ->settings::Tmin_K
  if (Temin)
    std::cerr << "Error in Cycler::cycleAgeing. The temperature " << Ti << "K is too low. The minimum value is 273.\n";

  bool Temax = Ti > settings::Tmax_Cell_K; //!< check the temperature is below 60 degrees, TMAX is defined in State.hpp ->settings::Tmax_K
  if (Temax)
    std::cerr << "Error in Cycler::cycleAgeing. The temperature " << Ti << " is too high. The maximum value is (273+60).\n";

  bool cchar = Ccha <= 0; //!< check the charging Crate is positive
  if (cchar)
    std::cerr << "Error in Cycler::cycleAgeing. The charging Crate " << Ccha << " is negative. It must be positive.\n";

  bool cdis = Cdis <= 0; //!< check the discharging Crate is positive
  if (cdis)
    std::cerr << "Error in Cycler::cycleAgeing. The discharging Crate " << Cdis << " is negative. It must be positive.\n";

  bool ccvchar = Ccutcha <= 0; //!< check the charging cutoff current is positive
  if (ccvchar)
    std::cerr << "Error in Cycler::cycleAgeing. The Crate of the cutoff current for the CV charge " << Ccutcha << " is negative. It must be positive.\n";

  bool ccvdis = Ccutdis <= 0; //!< check the discharging cutoff current is positive
  if (ccvdis)
    std::cerr << "Error in Cycler::cycleAgeing. The Crate of the cutoff current for the CV discharge " << Ccutdis << " is negative. It must be positive.\n";

  bool cycles = nrCycles <= nrCap; //!< check the number of cycles between consecutive check-ups is lower than the total number of cycles
  if (cycles)
    std::cerr << "Error in Cycler::cycleAgeing. The number of cycles between two check ups " << nrCap << " is higher than the total number of cycles " << nrCycles << ".\n";

  if (vmax || vmin || Temin || Temax || cchar || cdis || ccvchar || ccvdis || cycles)
    throw 1014;
}

} // namespace slide::util::error