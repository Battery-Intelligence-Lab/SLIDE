/*
 * Module_p.cpp
 *
 *  Created on: 18 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Module_p.hpp"

#include "../settings/settings.hpp"
#include "../utility/utility.hpp"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <array>
#include <algorithm>
#include <ctime>

namespace slide {

double Module_p::Vmin()
{
  //!< the voltage limits are the most constraining limits of all cells
  //!< ie the highest Vmin
  double vm{ 0 };

  for (auto &SU : SUs)
    vm = std::max(vm, SU->Vmin());

  return vm;
}

double Module_p::VMIN()
{
  //!< the voltage limits are the most constraining limits of all cells
  //!< ie the highest Vmin
  double vm{ 0 };
  for (auto &SU : SUs)
    vm = std::max(vm, SU->VMIN());

  return vm;
}

double Module_p::Vmax()
{
  //!< the voltage limits are the most constraining limits of all cells
  //!< ie the lowest Vmax
  double vm = std::numeric_limits<double>::max();
  for (auto &SU : SUs)
    vm = std::min(vm, SU->Vmax());
  return vm;
}

double Module_p::VMAX()
{
  //!< the voltage limits are the most constraining limits of all cells
  //!< ie the lowest Vmax
  double vm = std::numeric_limits<double>::max();
  for (auto &SU : SUs)
    vm = std::min(vm, SU->VMAX());
  return vm;
}

double Module_p::I()
{
  //!< the current is the sum  of the current of each cell. Returns 0 if empty.
  double i{ 0 };
  for (auto &SU : SUs)
    i += SU->I();

  return i;
}

double Module_p::getOCV(bool print)
{
  //!< the voltage is the same across all cells, so just return the V of the first cell
  double v{ 0 };
  for (auto &SU : SUs)
    v += SU->getOCV(print);

  return (v / SUs.size());
}

double Module_p::getRtot()
{
  /*
   * Return the total resistance
   * 		V(I) = OCV - I*Rtot
   * 		with V and OCV the total values for this module
   *
   * parallel: 1/Rtot = sum( 1/R_i )
   */

  //!< If there are no cells connected, return 0
  if (SUs.empty())
    return 0;

  //!< check if there are contact resistances only until SUs size. //!< #TODO why are we checking this?
  const bool noRc = std::all_of(Rcontact.begin(), Rcontact.begin() + SUs.size(), util::is_zero<double>);

  if (noRc) //!< no contact resistance
  {
    double rtot = 0;
    for (auto &SU : SUs)
      rtot += 1.0 / SU->getRtot();

    return 1.0 / rtot;
  }

  //!< with contact resistance
  //!< start from the cell furthest away
  double rtot = Rcontact.back() + SUs.back()->getRtot();

  //!< then iteratively come closer, every time Rcontact[i] + (Rcell[i] \\ Rtot)
  //!< 				= Rc[i] + Rcell[i]*Rtot / (Rcel[i]*Rtot)
  for (int i = getNSUs() - 2; i >= 0; i--) //!< #TODO bug if there are less than 2 SUs.
    rtot = Rcontact[i] + (SUs[i]->getRtot() * rtot) / (SUs[i]->getRtot() + rtot);

  return rtot;
}

double Module_p::getVi(size_t i, bool print)
{
  /*
   * Return the voltage of SU[i] as seen from the terminal while accounting for the contact resistance
   *
   * Contact resistances:
   * 		we assume the terminals of the parallel module are on either side (i.e. next to SUs[0] and SUs[N-1])
   * 		the contact resistances are in the 'horizontal' paths, and the current of all subsequent cells goes through it
   * 			i.e. the current through Rc[0] = i[0] + i[1] + ... + i[N-1}
   * 		The series-resistance of every cell is already included in the cells themselves (as Rdc), and in the cell voltage v[i]
   * 		The terminal voltage Vt must be the same for the paths to all cells, i.e.
   * 			Vt = v[0] - R[0] (I[0] + I[1] + ... + I[N-1]) = v[0] - R[0]*I[0:N-1]
   * 			   = V[1] - R[1]*I[1:N-1] - R[0]*I[0:N-1]
   * 			   = v[i] - sum{ R[j]*I[j..N-1], j=0..i }
   * 			   = v[i] - sum{ R[j]*sum(I[k], k=j..N-1), j=0..i }
   */

  if (i >= getNSUs()) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Module::getVi, you ask the voltage of cell " << i
                << " but the size of the cell array is " << getNSUs() << '\n';
    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw 10;
  }

  double v = SUs[i]->V(print); //!< the voltage of cell i
  if (v == 0)
    return 0;

  for (size_t j = 0; j <= i; j++) {
    double Ij{ 0 };                        //!< the current through parallel resistance j
    for (size_t k = j; k < getNSUs(); k++) //!< the sum of all currents 'behind' this resistance, i.e. from j to the last one
      Ij += SUs[k]->I();
    v -= Rcontact[j] * Ij;
  }

  return v;
}

double Module_p::V(bool print)
{
  //!< Return the mean voltage of all cells, since there can be a small difference between the voltages of cells

  //!< if the stored value is up to date, return that one
  if (Vmodule_valid)
    return Vmodule;

  //!< Else calculate the new voltage
  Vmodule = 0;
  for (size_t i = 0; i < SUs.size(); i++) {
    const auto v = getVi(i, print);
    if (v == 0)
      return 0;

    Vmodule += getVi(i, print); //!< #TODO there is definitely something fishy. getVi already does some calculations.
  }

  Vmodule /= SUs.size();
  Vmodule_valid = true;

  return Vmodule;
}

Status Module_p::redistributeCurrent(bool checkV, bool print)
{
  /*
   * Function which redistributes the current between different cells in the module to achieve the same voltage between all.
   * It does this by taking a bit of current from the cell with the largest voltage, and giving it to the cell with the smallest voltage.
   *
   * Note: Ideally, it works like a PI controller (the change in current is a function of the voltage error, i.e. the difference between
   * the voltage of this cell and the mean voltage of all cells)
   * Unfortunately, that is difficult to get to work correctly for a couple of reasons
   * 		absolute errors depend on how many series-cells there are. For a series of 1000 cells, even small different in
   * 		every cell can accumulate to large errors, which would give large corrections
   * 			this could be fixed by using relative voltages though
   * 		the OCV curve is not linear. Near 0% SOC the OCV is very steep, such that minor difference in SOC or current
   * 		result in huge voltage differences (both absolute and relative)
   * 		Due to the cell-to-cell variations, some cells might already be on this steep bit while others aren't, making
   * 		it difficult to assess what the steepness is
   * 			especially if the parallel module is made up of series or parallel modules, which will smooten out the behaviour of the weakest cell
   * 		The resistance of cells is different (e.g. more degraded cells can have double the resistance). So similar
   * 		corrections in current give different voltage responses
   * 		etc.
   * All these issues could be fixed if enough time and attention were spent on them. I did no such thing.
   * Instead, the correction is not directly proportional to the voltage error.
   * I start with an initial magnitude which is a function of the voltage error.
   * Then over time I reduce the correction (which is the factor f), such that you converge to the good outcome.
   * The numbers used below have been tweaked to ensure it is stable and works.
   * But a proper PI controller which addresses the issues above would be better
   *
   * IN
   * checkV 	check the voltage or not
   * print 	print messages
   *
   * OUT
   * number of iterations
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
   * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
   * 14 	failed to redistribute the current. Maximum number of iterations, or minimum change in current reached
   */

#if TIMING
  Clock clk;
#endif

  //!< variables
  const size_t nIterationsMax = 1000 * getNSUs();                      //!< allow a maximum number of iterations //!< Very big iteration....
  const double dImin = settings::MODULE_P_I_ABSTOL / 10.0 / getNSUs(); //!< allow a minimum change in current
  double dI;                                                           //!< change in current in this iteration
  const bool verb = print && (settings::printBool::printCrit);         //!< print if the (global) verbose-setting is above the threshold
  Vmodule_valid = false;                                               //!< we are changing the current to individual cells, so the stored voltage is no longer valid

  //!< get cell voltages
  std::array<double, settings::MODULE_NSUs_MAX> Vo, Iold; //!< #TODO if we should make them vector.

  //!< voltage and initial current of each cell //!< #TODO it is a constant value SU.
  for (size_t i = 0; i < getNSUs(); i++) {
    Vo[i] = getVi(i, print);
    Iold[i] = SUs[i]->I();
  }

  //!< check if there are contact resistances #TODO if we can do better.
  const bool noRc = std::all_of(Rcontact.begin(), Rcontact.begin() + SUs.size(), util::is_zero<double>);

  //!< variables
  double f;             //!< fraction of the cell current to change per iteration
  size_t imax_prev = 0; //!< index of cell with std::max voltage in previous iteration
  size_t imin_prev = 0;
  double dV;          //!< voltage difference between min and std::max
  double dV_prev = 0; //!< voltage difference with the previous value of f
  size_t k = 0;       //!< iteration with this value of f
  size_t ktot = 0;    //!< total number of iterations

  bool reached = false;

  //!< Initial value of the change in current per iteration
  double Cmean = 0; //!< mean crate of the current in the module
  for (auto &SU : SUs)
    Cmean += std::abs(SU->I()) / SU->Cap();

  Cmean = Cmean / SUs.size();

  if (Cmean < 1e-3) //!< if the current of each cell is 0
    Cmean = 0.1;    //!< then use dI of 0.1C to equalise the currents

  //!< iterate until the voltages are satisfactory close
  while (!reached) //!< #TODO if this algorithm is efficient.
  {

    //!< find the cells with the smallest and largest V
    const auto [minIter, maxIter] = std::minmax_element(std::begin(Vo), std::begin(Vo) + getNSUs());

    size_t imin = std::distance(std::begin(Vo), minIter); //!< indices with extreme voltages.
    size_t imax = std::distance(std::begin(Vo), maxIter);

    //!< stop iterating if this difference is below the threshold
    dV = *maxIter - *minIter;
    if (dV < settings::MODULE_P_V_ABSTOL || dV / Vo[imax] < settings::MODULE_P_V_RELTOL)
      reached = true;

    /* in the very first iteration, initialise f base on dV relative to range of the voltage
     * 	for one cell (Vmax-Vmin = 1.5):
     * 		dV > 1 		f = 0.3
     * 		dV > 0.1 	f = 0.2
     * 		else 		f = 0.05
     * so for N cells in series (Vmax - Vmin = N * 1.5)
     * 		dV > N 		f = 0.3
     * 		dV > 0.1*N	f = 0.2
     * 		else 		f = 0.05
     */
    if (ktot == 0) {
      const int Nseries = static_cast<int>((Vmax() - Vmin()) / 1.5); //!< #TODO this probably takes some time.
      if (dV > 1.0 * Nseries)
        f = 0.3; //!< for 1V, change 30% of the cell current
      else if (dV > 0.1 * Nseries)
        f = 0.2; //!< for 0.1V, change 10% of the cell current
      else
        f = 0.05; //!< else, change 5% of the cell current

      if ((Vo[imin] - SUs[imin]->Vmin() < 0.2) || (SUs[imax]->Vmax() - Vo[imax] < 0.1))
        f = 0.02; //!< if one of the children is close to its voltage limits, change 2% of the cell current
    }

    //!< initialise the voltage difference with this value of f;
    if (k == 0)
      dV_prev = dV;

    //!< If the difference is too large, swap currents
    if (!reached) {

      //!< swap currents
      dI = f * Cmean * SUs[imax]->Cap(); //!< current we swap in this iteration
      //!< dI = f * std::max(std::abs(SUs[imax]->I()), std::abs(SUs[imin]->I())); 	//!< current we swap this iteration. f*Icell, where Icell is smallest of the currents of both cells
      SUs[imax]->setCurrent(SUs[imax]->I() + dI, false, print); //!< don't check the voltages, since at the end of a full (dis)charge, all cells will be just over their allowed limit
      SUs[imin]->setCurrent(SUs[imin]->I() - dI, false, print); //!< and then this would throw 2
      k++;
      ktot++;

      //!< update voltages
      if (noRc) { //!< no contact resistance, so only for cells which change current, does the V change
        try {
          Vo[imax] = getVi(imax, print); //!< #TODO if getVi is throwing.
          Vo[imin] = getVi(imin, print);
        } catch (int e) {
          std::cout << "Error " << e << " in redistributeCurrent iteration " << ktot << " for module "
                    << getFullID() << " when getting the voltage after setting currents of "
                    << SUs[imax]->I() + dI << " and " << SUs[imin]->I() - dI
                    << " with dI = " << dI << ", dV = " << dV << ", Vmax = " << Vo[imax]
                    << ", Vmin = " << Vo[imin] << ". Giving up.\n";

          //!< try to recover by almost reverting the changes
          //!< (If you fully rever, you will get an eternal loop since the next iteration will repeat the same mistake)
          try {
            SUs[imax]->setCurrent(SUs[imax]->I() - 0.95 * dI, false, print);
            SUs[imin]->setCurrent(SUs[imin]->I() + 0.95 * dI, false, print);
            Vo[imax] = getVi(imax, print);
            Vo[imin] = getVi(imin, print);
            std::cout << "The code recovered by resetting 95% of the changes, continue as normal.\n";
          } catch (int e) {
            std::cout << "Even after resetting the changes, we still get error " << e
                      << " so we are restoring the original currents and giving up.\n";
            SUs[imax]->setCurrent(SUs[imax]->I() - 0.05 * dI, false, print); //!< don't check the voltage, it should be fine but we are not sure
            SUs[imin]->setCurrent(SUs[imin]->I() + 0.05 * dI, false, print);
            reached = true;
          }
        }
      } else {
        try {
          for (size_t i = 0; i < getNSUs(); i++)
            Vo[i] = getVi(i, print);
          /*
           * with contactR, the terminal voltage from each cell changes, since a change in I[i] affects the current through the parallel resistances
           * and therefore the voltage to each cell
           */
        } catch (int e) {
          std::cout << "Error " << e << " in redistributeCurrent iteration " << ktot
                    << " for module " << getFullID() << " when getting the voltage after setting currents of "
                    << SUs[imax]->I() + dI << " and " << SUs[imin]->I() - dI << " with dI = " << dI
                    << ", dV = " << dV << ", Vmax = " << Vo[imax] << ", Vmin = " << Vo[imin] << ". Giving up\n";
          //!< try to recover by almost reverting the changes
          try {
            SUs[imax]->setCurrent(SUs[imax]->I() - 0.95 * dI, false, print);
            SUs[imin]->setCurrent(SUs[imin]->I() + 0.95 * dI, false, print);
            for (size_t i = 0; i < getNSUs(); i++)
              Vo[i] = getVi(i, print);
            std::cout << "The code recovered by resetting 95% of the changes, continue as normal.\n";
          } catch (int e) {
            std::cout << "Even after resetting the changes, we still get error " << e
                      << " so we are restoring the original currents and giving up.\n";
            SUs[imax]->setCurrent(SUs[imax]->I() - 0.05 * dI, false, print); //!< don't check the voltage, it should be fine but we are not sure
            SUs[imin]->setCurrent(SUs[imin]->I() + 0.05 * dI, false, print);
            reached = true;
          }
        }
      } //!< end update the voltages

      //!< if dI is too big, we are simply swapping current from one cell to the other and then back
      //!< avoid an eternal loop by reducing f in this case
      if (imax_prev == imin && imin_prev == imax) {
        f /= 2.0;
        k = 0;
      }
      imax_prev = imax;
      imin_prev = imin;

      //!< If dV is decreasing and we have done enough iterations, decrease f slowly
      if (k > getNSUs() && dV < dV_prev / 2.0) {
        f = f / 2.0;
        k = 0;
      }

      //!< avoid eternal loop, if dI is too small, or we have done too many iterations, give up
      if ((ktot > nIterationsMax) || dI < dImin) {
        //!< find the cells with the smallest and largest V
        const auto [min, max] = std::minmax_element(std::begin(Vo), std::begin(Vo) + getNSUs());

        imin = std::distance(std::begin(Vo), min);
        imax = std::distance(std::begin(Vo), max);

        dV = Vo[imax] - Vo[imin];
        if (dV > settings::MODULE_P_V_ABSTOL && dV / Vo[imax] > settings::MODULE_P_V_RELTOL) {
          if (verb) {
            std::cerr << "error in Module_p::redistributeCurrent, the total number if iterations is " << ktot
                      << " and the change in current is " << dI << ", limits for both are " << nIterationsMax << " and "
                      << dImin << ". Stop iterating with a voltage error of " << std::abs(Vo[imax] - Vo[imin])
                      << ". The allowed absolute tolerance is " << settings::MODULE_P_V_ABSTOL
                      << " and relative tolerance " << settings::MODULE_P_V_RELTOL * Vo[imax] << '\n';
          }
          std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
          return Status::RedistributeCurrent_failed;
        }
      }
    }
  }

  std::cout << "Number of iterations in redistributeCurrent() is: " << ktot << '\n';

  if (checkV) {
    double v;
    auto status = checkVoltage(v, print);

    if (isStatusWarning(status)) {
      //!< valid range -> keep old states, but throw 2 to notify something is not as it should be
      if (verb)
        std::cerr << "Error in Module_p::redistributeCurrent(). The voltage of one of the cells is "
                     "outside the allowed limits. allowing it for now.\n";
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      return status;
    } else if (isStatusBad(status)) { //!< outside safety limits (< VMIN and discharge || > VMAX and charge)
      //!< safety limits -> restore old states and throw 3
      if (verb)
        std::cerr << "Error in Module_p::redistributeCurrent(). The voltage of one of the cells is "
                     "outside the safety limits. Restoring the old currents and throwing 3.\n";

      for (size_t i = 0; i < getNSUs(); i++)
        SUs[i]->setCurrent(Iold[i], false, print); //!< don't check the voltage since we restore the original currents, which should be valid

      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      return status;
    }
  }

#if TIMING
  timeData.redistributeCurrent += clk.duration(); //!< time in seconds
#endif
  return Status::Success;
}

Status Module_p::setI_iterative(double Inew, bool checkV, bool print)
{
  /*
   * Set the current of a parallel module
   * This function takes small steps adapting the current of each connected cell until the total current is reached
   * 	 step 1 	change I by a bit  -> measure V_cell -> derive Rtot
   * 	 step 2 	do 50% of Inew-I() by increasing current proportionally to this resistance
   * 	 step 3   	iteratively change I of the cell with the smallest V (charge) or biggest V (discharge)
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
   * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
   * 15 	after setting the current, the voltage of the cells are too far apart
   */

#if TIMING
  Clock clk;
#endif

  //!< Gains in the PI controllers
  //!< fraction of current changes set in various steps
  const double f2 = 0.5; //!< allocated 50% of the change in current proportional to the resistance
  const double f3 = 0.2; //!< allocate 20% of the change in current per cell in each iteration in step 3
                         //!< note, this is always relative to the error, i.e f* (Inew-I())/N
                         //!< so it will automatically become smaller, it is more the gain of a P controller

  //!< check if there are contact resistances #TODO why
  const bool noRc = std::all_of(Rcontact.begin(), Rcontact.begin() + SUs.size(), util::is_zero<double>);

  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold
  double dI = Inew - I();                                      //!< total change in current needed
  double v = 0;                                                //!< end voltage
  double vn;                                                   //!< voltage of cell n

  //!< Change the current if needed
  if (dI != 0) {
    double dIi = dI / getNSUs(); //!< initial guess for change in current per cell (uniform if all cells are equal)
    dIi = dIi / 2.0;             //!< in first step, change half the current
    Vmodule_valid = false;       //!< we are changing the current, so the stored voltage is no longer valid

    //!< STEP 1: CHANGE ICELL BY dIi/2 AND SEE HOW CELL VOLTAGE CHANGES to estimate resistance
    std::array<double, settings::MODULE_NSUs_MAX> Vo, Ri, Iold; //#TODO these probably need to be vectors since NSUs will increase.
    double Rtot = 0;                                            //!< total resistance of the parallel module, 1/Rtot = sum ( 1/R(i) )
    for (size_t i = 0; i < getNSUs(); i++) {
      Vo[i] = getVi(i, print); //!< terminal voltage reached from cell i
      Iold[i] = SUs[i]->I();
      SUs[i]->setCurrent(SUs[i]->I() + dIi, false, print); //!< don't check the voltage since this is an intermediate step
    }

    for (size_t i = 0; i < getNSUs(); i++) {
      vn = getVi(i, print);
      Ri[i] = (Vo[i] - vn) / dIi;
      Ri[i] = std::max(Ri[i], 1e-5); //!< avoid R=0, which would give a NaN. Instead give very small value
      Rtot += 1.0 / Ri[i];
    }

    Rtot = 1.0 / Rtot;

    //!< STEP 2 FROM INITIAL CELL CURRENTS, INCREASE THEIR VALUE PROPORTIONALLY TO THEIR RESISTANCE
    //!< at Iold(i), V(i) == V(j)
    //!< 	if linear, dI(i) will change V(i) by -Ri*dI(i)
    //!< 	current divider rule: dI(i) = (Rtot / R(i)) dItot
    for (size_t i = 0; i < SUs.size(); i++) {
      dIi = f2 * dI * Rtot / Ri[i];
      SUs[i]->setCurrent(Iold[i] + dIi, false, print); //!< don't check the voltage since this is an intermediate step
    }

    for (size_t i = 0; i < SUs.size(); i++) //!< only get voltage once all currents are set, since through contactR I[j] affects V[k]
      Vo[i] = getVi(i, print);

    //!< STEP 3 ITERATE TO INCREASE THE VOLTAGE OF THE MOST-DIFFERENT CELL
    //!< if dI > 0, we are discharging more (or charging less), so the cell with the highest voltage should get a larger discharge (smaller charge) current
    //!< if dI < 0, we are charging more (or discharging less), so the cell with the lowest voltage should get a larger charge (smaller discharge) current
    //!< every iteration, allocate 10% of the remaining current per cell
    bool reached = false; //!< boolean indicating we are exactly at the correct current
    bool almost = false;  //!< boolean indicating wether we are almost there (i.e. one more iteration)
    size_t ind;           //!< index of the cell which should be changed
                          //		int k=0;//
    while (!reached) {
      dIi = f3 * (Inew - I()) / SUs.size(); //!< if uniform, each cell should get (Inew-I())/Ncell. Allocate a fraction f of this amount now

      //!< if the remaining current is smaller than dI, then allocate all of the remaining current
      if (almost || std::abs(I() - Inew) < std::abs(dIi)) {
        almost = true;
        dIi = Inew - I();
      }

      if (dI < 0)
        ind = std::distance(std::begin(Vo), std::min_element(std::begin(Vo), std::begin(Vo) + SUs.size()));
      else
        ind = std::distance(std::begin(Vo), std::max_element(std::begin(Vo), std::begin(Vo) + SUs.size()));

      SUs[ind]->setCurrent(SUs[ind]->I() + dIi, false, print); //!< don't check the voltage since this is an intermediate step

      if (noRc) //!< no Rcontact so only voltage of cell ind changes
        Vo[ind] = getVi(ind, print);
      else
        for (size_t i = 0; i < getNSUs(); i++) //!< voltage of all cells changes
          Vo[i] = getVi(i, print);

      //!< check tolerances
      if (almost) //!< we have done the final step where we allocate all the remaining current
        reached = true;
      if (std::abs(Inew - I()) < settings::MODULE_P_I_ABSTOL) //!< we are very close (absolutely) so next time allocate all remaining current
        almost = true;
      if (std::abs(Inew - I()) < settings::MODULE_P_I_RELTOL * Inew) //!< relative
        almost = true;
    }

    //!< check the voltages are valid
    if (checkV) {
      const auto status = checkVoltage(v, print);
      if (isStatusWarning(status)) {
        //!< valid range -> keep old states, but throw 2 to notify something is not as it should be
        if (verb)
          std::cout << "Error in Module_p::setCurrent(). The voltage of one of the cells is outside the allowed limits. allowing it for now.\n";
        std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
        return status;
      } else if (isStatusBad(status)) { //!< outside safety limits (< VMIN and discharge || > VMAX and charge)
        //!< safety limits -> restore old states and throw 3
        if (verb)
          std::cout << "Error in Module_p::setCurrent(). The voltage of one of the cells is outside the safety limits. Restoring the old currents and throwing 3.\n";

        for (size_t i = 0; i < getNSUs(); i++)
          SUs[i]->setCurrent(Iold[i], false, print); //!< don't check the voltage since we restore the original currents, which should be valid

        std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
        return status;
      }
    }

  } //!< if dI != 0
  else
    v = V(print); //!< the current does not change, so just return the voltage

  //!< check the voltage difference is between the tolerance
  if (checkV) {
    const bool val = validSUs(SUs, print);
    if (!val) {
      if (verb)
        std::cout << "error in Module_p::setCurrent for SU " << getFullID() << ", the cell voltage is not valid "
                  << "after setting the module voltage to reach a total current of " << Inew
                  << ". Trying to recover by redistributing the current.\n";

      //!< redistribute the current, if that is successful, we have valid cell voltages
      try {
        auto status = redistributeCurrent(checkV, print); //!< #TODO this is utterly wrong.
        if (status != Status::Success)
          throw 100000;
      } catch (int) {
        if (verb)
          std::cout << "error in Module_p::setCurrent. We tried redistributing the current to equalise "
                       "the cell voltages after setting the new current, "
                       " but redistributeCurrent() failed too. Give up and throwing 15";
        std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
        throw 15; //!< Why 15?
      }
    }
  }

  //	std::cout<<"setCurrent() for SU "<<getFullID()<<" finishing.\n"; //

#if TIMING
  timeData.setCurrent += clk.duration(); //!< time in seconds
#endif
  return Status::Success;
}

Status Module_p::setCurrent(double Inew, bool checkV, bool print)
{
  /*
   * Set the current of a parallel module
   * This function takes small steps adapting the current of each connected cell until the total current is reached
   * 	 step 1 	change I by a bit  -> measure V_cell -> derive Rtot
   * 	 step 2 	do 50% of Inew-I() by increasing current proportionally to this resistance
   * 	 step 3   	iteratively change I of the cell with the smallest V (charge) or biggest V (discharge)
   *
   * THROWS
   * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
   * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
   * 15 	after setting the current, the voltage of the cells are too far apart
   */

#if TIMING
  std::clock_t tstart = std::clock();
#endif

  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold
  double v;

  //!< get the old currents so we can revert if needed
  std::vector<double> Iolds; //!< #TODO we need to remove vector somehow?
  Iolds.clear();

  for (auto &SU : SUs)
    Iolds.push_back(SU->I());

  //!< allocate the current uniformly
  try {
    for (size_t i = 0; i < getNSUs(); i++) {

      auto status = SUs[i]->setCurrent(Inew / getNSUs(), checkV, print);

      //!< voltage of cell i is outside the valid range, but within safety limits
      //!< indicate this happened but continue setting states
      if (isStatusWarning(status)) {
        if (verb)
          std::cout << "warning in Module_p::setCurrent, the voltage of cell " << i << " with id "
                    << SUs[i]->getFullID() << " is outside the allowed range for Inew = " << Inew / getNSUs()
                    << ". Continue for now since we are going to redistribute the current to equalise the voltages.\n";
      } else {
        if (verb)
          std::cout << "ERROR " << (int)status << " in Module_p::setCurrent when setting the current of cell "
                    << i << " with id " << SUs[i]->getFullID() << " for Inew = " << Inew / getNSUs()
                    << ". Try to recover using the iterative version of setCurrent.\n";
        //!< throw error, the catch statement will use the iterative function
      }

    } //!< loop

    //!< Redistribute the current to equalise the voltages
    try {
      redistributeCurrent(checkV, print); //!< this will check the voltage limits if needed
      v = V();
    } catch (int e) {
      if (e == 2) {
        if (verb)
          std::cout << "warning in Module_p::setCurrent, after redistribute, the voltage of one of the cells is "
                    << "outside the allowed but inside the safe range for Inew = " << Inew << ". Continue for now.\n";
      } else {
        if (verb)
          std::cout << "Error in Mpodule_p::setCurrent when redistributing the current: " << e
                    << ", Try to recover using the iterative version of setCurrent.\n";
        std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
        throw e;
      }
    }
  } //!< try allocating uniform current
  catch (int) {
    //!< revert to the old current
    for (size_t i = 0; i < getNSUs(); i++)
      SUs[i]->setCurrent(Iolds[i], false, true);

    try {
      setI_iterative(Inew, checkV, print); //!< this will check the voltage limits if needed
      v = V();
      if (verb)
        std::cout << "Module_p::setCurrent, We managed to recover using the iterative version. Continue as normal.\n";
    } catch (int e) {
      if (verb)
        std::cerr << "Module_p::setCurrent, Even the iterative version failed with error " << e
                  << ". Giving up, restore the original current and throwing the error on.\n";

      for (size_t i = 0; i < getNSUs(); i++)
        SUs[i]->setCurrent(Iolds[i], false, true);

      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      throw e;
    }
  } //!< catch statement of uniform allocation

#if TIMING
  timeData.setCurrent += (std::clock() - tstart) / static_cast<double>(CLOCKS_PER_SEC); //!< time in seconds
#endif
  return Status::Success; //!< #TODO problem
}

bool Module_p::validSUs(moduleSUs_span_t c, bool print)
{
  /*
   * Checks the cells are a valid combination for a parallel-connected module
   * the voltage differences are within the tolerance
   *
   * If the number of cells is the same as in this module, use the contact resistances
   */
#if TIMING
  std::clock_t tstart = std::clock();
#endif
  const bool verb = print && (settings::printBool::printCrit); //!< print if the (global) verbose-setting is above the threshold

  //!< Check the voltage of each cell is valid and within the error tolerance #TODO it is better to supply both module and Rcontact.
  auto Vi = [N = getNSUs(), c, print, this](size_t i) {
    double v = c[i]->V(print);

    if (c.size() != N)
      return v;

    for (size_t j = 0; j <= i; j++) { //!< account for the contact resistances
      double Ij = 0;
      for (size_t k = j; k < N; k++) //!< the sum of all currents 'behind' this resistance, i.e. from j to the last one
        Ij += c[k]->I();
      v -= Rcontact[j] * Ij;
    }

    return v;
  };

  size_t i_min{}, i_max{};
  double V_min{ 1e16 }, V_max{ 0 };

  for (size_t i = 0; i < c.size(); i++) //!< find the cells with the smallest and largest V
  {
    const auto V_i = Vi(i);
    if (V_i < V_min) {
      V_min = V_i;
      i_min = i; //!< cell with the smallest voltage
    }

    if (V_i > V_max) {
      V_max = V_i;
      i_max = i; //!< cell with the largest voltage
    }
  }

  //!< Check that this limit is below the absolute or relative threshold
  const double dV = V_max - V_min;
  bool result{ true };
  if (dV > settings::MODULE_P_V_ABSTOL && dV > settings::MODULE_P_V_RELTOL * V_max) //!< #TODO if this should be || not &&
  {
    if (verb) {
      std::cout << "error in Module_p::validSUs for SU = " << getFullID() << ", the maximum voltage is in cell"
                << i_max << " and is " << V_max << " while the minimum voltage is in cell" << i_min << " and is "
                << V_min << " which is an error of " << dV << " and the allowed absolute tolerance is "
                << settings::MODULE_P_V_ABSTOL << ", the allowed relative tolerance gives an error of "
                << settings::MODULE_P_V_RELTOL * V_max << '\n';
    }
    result = false;
  } //!< else the voltage is valid

#if TIMING
  timeData.validSUs += (std::clock() - tstart) / static_cast<double>(CLOCKS_PER_SEC); //!< time in seconds
#endif
  return result;
}

void Module_p::timeStep_CC(double dt, int nstep)
{
  /*
   * Take a time step at a constant current.
   * There are two ways to do this:
   * 		rootFinding::Current_EQN with dti = dt, which explicitly solves the system of equations
   * 		let every cell take a CC time step, and check if the voltage equation is satisfied
   * 			if not, redistribute the current using setCurrent()
   *
   * The second approach is probably quicker since it only solves the system of equations
   * if the voltage difference becomes too large
   */

#if TIMING
  std::clock_t tstart = std::clock();
#endif

  if (dt < 0) {
    if constexpr (settings::printBool::printCrit)
      std::cerr << "ERROR in Module_p::timeStep_CC, the time step dt must be 0 or positive, but has value "
                << dt << '\n';
    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw 10;
  }

  //!< we simply take one CC time step on every cell
  auto task_indv = [&](int i) { SUs[i]->timeStep_CC(dt, nstep); };

  try {
    run(task_indv, getNSUs(), (par ? -1 : 1));
  } catch (int e) {
    std::cout << "Error in Module_p::timeStep_CC with module ID " << getFullID()
              << ". error " << e << ", throwing it on.\n";
    std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
    throw e;
  }

  //!< **************************************************** Calculate the thermal model once for the nstep * dt time period *****************************************************************

  //!< update the time since the last update of the thermal model
  if (!blockDegAndTherm) {
    therm.time += nstep * dt;

    //!< Increase the heat from the contact resistances
    double Ii = 0; //!< current through resistor I
    for (size_t i = 0; i < SUs.size(); i++) {
      for (size_t j = i; j < SUs.size(); j++)
        Ii += SUs[j]->I(); //!< resistor i sees the currents through the cells 'behind' them

      therm.Qcontact += Rcontact[i] * std::pow(Ii, 2) * nstep * dt;
    }

    //!< If this module has a parent module, this parent will call the thermal model with the correct parameters
    //!< which will include heat exchange with the module's neighbours and cooling from the cooling system of the parent module.

    //!< if there is no parent, this module is the top-level.
    //!< It then directly exchanges heat with the environment at a fixed temperature
    //!< If there is no parent, assume we update T every nstep*dt. So update the temperature now
    if (!parent) {
      //!< double check that we have an HVAC coolsystem (see constructor)
      if (typeid(*getCoolSystem()) != typeid(CoolSystem_HVAC)) {
        std::cerr << "ERROR in module_p::timeStep_CC in module " << getFullID() << ". this is a top-level"
                  << "module but does not have an HVAC coolsystem for active cooling with the environment.\n";
        std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
        throw 14;
      }

      //!< Call the thermal model without heat exchanges with neighbours or parents (since this module doesn't have any)
      //!< 	Heat exchange between this module and the environment is done by the AC-part of the HVAC coolsystem
      //!< 	heat exchange between this module and its children is done the conventional way by the HVAC coolsystem
      //!< set the new temperature since we have calculated all the temperatures so there is no risk for inconsistency
      double Tneigh[1], Kneigh[1], Aneigh[1];
      //!< Make the arrays even though they will not be used (length should be 0 but I don't think you can make an array of length 0)
      setT(thermalModel(0, Tneigh, Kneigh, Aneigh, therm.time)); //!< the 0 signals there are no neighbours or parents

      /*
      double Tneigh[1] = {settings::T_ENV};							//!< T of environment
      double Kneigh[1] = {cool->getH()}; 					//!< h to environment is same as to children, since the speed of the coolant is the same
      double Aneigh[1] = {therm.A};						//!< A to the environment is the A of this module
      //!< note that this will have more heat exchange than to children, since A = min(A_this, A_parent)
      //!< 	children will have a smaller A themselves, so the resulting A = A_child < A of this module
      setT(thermalModel(1, Tneigh, Kneigh, Aneigh, therm.time));*/
    }

    //!< control the cooling system
    double Tlocal = 0;
    for (auto &SU : SUs)
      Tlocal = std::max(Tlocal, SU->T());

    cool->control(Tlocal, getThotSpot());
  }

  Vmodule_valid = false; //!< we have changed the SOC/concnetration, so the stored voltage is no longer valid

  //!< check if the cell's voltage is valid
  if (!validSUs(SUs, false)) {
    try {
      auto status = redistributeCurrent(false, true); //!< don't check the currents

      if (status != Status::Success)
        throw 100000; //!< #TODO
    } catch (int e) {
      std::cout << "error in Module_p::timeStep_CC when redistributing the current. Throwing the error on " << e << '\n';
      std::cout << "Throwed in File: " << __FILE__ << ", line: " << __LINE__ << '\n';
      throw e;
    }
  }

#if TIMING
  timeData.timeStep += (std::clock() - tstart) / static_cast<double>(CLOCKS_PER_SEC); //!< time in seconds
#endif
}

TimingData_Module_p Module_p::getTimings()
{
#if TIMING
  return timeData;
#else
  return {};
#endif
}

void Module_p::setTimings(TimingData_Module_p td)
{
#if TIMING
  timeData = std::move(td);
#endif
}

Module_p *Module_p::copy()
{
  //!< check the type of coolsystem we have #TODO for a better way.
  int cooltype = 0;
  if (typeid(*getCoolSystem()) == typeid(CoolSystem_HVAC))
    cooltype = 1;
  else if (typeid(*getCoolSystem()) == typeid(CoolSystem_open))
    cooltype = 2;

  Module_p *copied_ptr = new Module_p(getID(), cool->T(), true, par, getNcells(), cool->getControl(), cooltype);

  copied_ptr->Rcontact = Rcontact;
  copied_ptr->setT(T());

  copied_ptr->SUs.clear();
  for (size_t i{ 0 }; i < getNSUs(); i++) {
    copied_ptr->SUs.emplace_back(SUs[i]->copy());
    copied_ptr->SUs.back()->setParent(copied_ptr);
  }

#if TIMING
  copied_ptr->setTimings(T_redistributeCurrent, T_setI, T_validSUs, T_timeStep, T_timeStepi);
#endif

  return copied_ptr;
}
} // namespace slide
