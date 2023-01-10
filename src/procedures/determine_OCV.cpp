/*
 * determineOCV.cpp
 *
 * The functions below can be used to find parameters relating to the OCV of a cell if a measured OCV curve is supplied.
 *
 * As input, the user has to supply the electrode OCV curves and the cell OCV curve.
 * See readOCVinput()
 *
 * Then the functions below calculate
 * 		the 'windows' of the half-cell OCV curve you want to use, i.e. the li-fractions of each electrode at the fully charged and fully discharged state
 * 		the amount of active material of each electrode
 *
 * These values can be used in the Cell-constructors:
 * 		the initial concentration (for a cell starting at 50% SOC)
 * 		the electrode surface, thickness, volume fraction, effective surface
 * 			changing the radius of the particles is not recommended because it means the spatial discretisation changes, and hence the model matrices have to be recalculated in MATLAB
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */
#include "determine_OCV.hpp"
#include "../utility/utility.hpp"
#include "../settings/settings.hpp"

#include <thread>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <ctime>
#include <utility>

namespace slide {
bool validOCV(bool checkRange, slide::XYdata_vv &data)
{
  /*
   * Function to check if the OCV curves supplied by the users have the correct format.
   * the first column must have strictly increasing values starting at 0.
   * 		the electrode OCV curves must end at 1 (they have the li-fractions, so the range must be from 0 to 1)
   * 		the electrode OCV curve can end at any point
   * the second column must have realistic voltages (0 <= V <= 10)
   *
   * If one of the conditions is not satisfied (i.e. the OCV curve is illegal)
   * an error message is printed and the function returns false.
   *
   * IN
   * checkRange 	boolean to indicate if we need to check that the range of the first column goes to 1
   * 				set to true for the electrode OCV curves (where the first column gives the lithium fraction, so from 0 to 1)
   * 				set to false for the cell OCV curve (where the first column gives the discharged charge, so from 0 to some positive value)
   * x 			matrix with the OCV curve.
   * 				the first column gives the 'x-values' (lithium fraction or discharged charge)
   * 				the second column gives the 'y-values' (open circuit voltage)
   *
   * OUT
   * bool 		true if the OCV curve has the correct format
   */

  //!< Check if the range of the first column is valid
  bool range;                                                                 //!< boolean to indicate of the first column has the correct range
  if (checkRange) {                                                           //!< check the first column has the correct range, from 0 to 1
    range = std::abs(data.x[0]) < 1e-5 && std::abs(data.x.back() - 1) < 1e-5; //!< allow a small error of e-5
    if (!range)
      std::cerr << "illegal value in deterimineOCV::validOCV: the range of the first column is wrong. It goes from "
                << data.x[0] << " to " << data.x.back() << " instead of from 0 to 1.\n";
  } else { //!< check the first column has the correct starting point, 0
    range = std::abs(data.x[0]) < 1e-5;
    if (!range)
      std::cerr << "illegal value in deterimineOCV::validOCV: the range of the first column is wrong. It starts at "
                << data.x[0] << " instead of starting at 0.\n";
  }

  //!< Loop through the curve to check that the values of the first column are increasing and the voltages in the second column are realistic.
  bool decreasei = false;   //!< boolean to indicate if this row has a decreasing value in the 1st column
  bool decrease = false;    //!< boolean to indicate if there are decreasing values in the 1st column
  bool unrealistic = false; //!< boolean to indicate if there is an unrealistic value in the 2nd column
  for (size_t i = 0; i < data.size(); i++) {
    //!< check for increasing values
    if (i > 0)
      decreasei = (data.x[i] <= data.x[i - 1]);
    if (decreasei) { //!< print an error message
      std::cerr << "illegal value in deterimineOCV::validOCV: the first column has decreasing or double values at rows "
                << i - 1 << " and " << i << ".\n"
                << "values are " << data.x[i - 1] << " and " << data.x[i] << ".\n";
      decrease = true; //!< we have found a decreasing value
    }

    //!< check for realistic voltages
    if ((data.y[i] < 0) || (data.y[i] > 10)) { //!< print an error message
      std::cerr << "illegal value in deterimineOCV::validOCV: the second column has an illegal value of "
                << data.y[i] << " at row " << i << ". The values have to be between 0 and 10.\n";
      unrealistic = true; //!< we have found an unrealistic value
    }
  }

  //!< return true if the format is correct
  return range && !decrease && !unrealistic; //!< the range must be correct and we cannot have any decreasing or unrealistic values
}


auto discharge_noexcept(const slide::XYdata_vv &OCVp, const slide::XYdata_vv &OCVn, double cap, const double AMp, const double AMn, const double cmaxp,
                        const double cmaxn, double sp, double sn, double Vend, slide::XYdata_vv &OCV, slide::XYdata_vv &OCVanode,
                        slide::XYdata_vv &OCVcathode, double fp[], double fn[])
{
  /*
   * Function to simulate a CC discharge at a very low current, so we can ignore resistances and diffusion.
   *
   * It assumes the concentrations are uniform at all times, so we don't have to account for diffusion.
   * Therefore, we don't need to use the full cell model.
   * Rather, we can simulate two 'buckets of lithium' and we only need to exchange lithium between both.
   *
   * IN
   * OCVp 	cathode OCV curve
   * OCVn 	anode OCV curve
   * cap		capacity of the cell [Ah]
   * AMp 		amount of active material on the cathode [m3], elec_surf * thicp * ep
   * AMn 		amount of active material on the anode [m3], elec_surf * thickn * en
   * cmaxp 	maximum li-concentration in the cathode [mol m-3]
   * cmaxn	maximum li-concentration in the anode [mol m-3]
   * sp 		start li-fraction at the cathode (between 0 and 1)
   * sn 		start li-fraciton at the anode (between 0 and 1)
   * Vend 	voltage when the OCV discharge should end [V]
   *
   * OUT
   * OCV 		OCV curve of the cell simulated with the given parameters
   * 			the first column [i][0] gives the discharged charge in [Ah]
   * 			the second column [i][1] gives the cell OCV in [V]
   * OCVanode OCV of the anode
   * 			the first column [i][0] gives the lithium fraction in the anode
   * 			the second column [i][1] gives the voltage of the anode
   * OCVcathode OCV of the cathode
   * 			the first column [i][0] gives the lithium fraction in the cathode
   * 			the second column [i][1] gives the voltage of the cathode
   * fp		array with the lithium fraction in the cathode at 100% SOC, 50% SOC and 0% SOC
   * fn	 	array with the lithium fraction in the anode at 100% SOC, 50% SOC and 0% SOC
   */

  //!< *********************************************************** 1 variables ***********************************************************************
  using namespace PhyConst;
  //!< Variables
  const double n = 1.0;       //!< number of electrons involved in the reaction
  double Ah = 0;              //!< discharged charge up to now
  const double I = 0.1 * cap; //!< magnitude of the discharge current [A]
  const double dt = 5;        //!< time step to use in the simulation [s]
  const bool bound = true;    //!< check the surface concentration is within the allowed limits when doing linear interpolation.
  bool end = false;           //!< boolean to check if we have reached the end voltage

  //!< If we store a data point for every time step, the output arrays will be very long.
  //!< Therefore, store only one value every 'nstore' time steps
  const size_t nstore = 100; //!< store one data point every 100 time steps

  bool flag{ true };
  double ocvpi, ocvni, v;

  auto interpolate = [&]() {
    bool flag_ocvpi, flag_ocvni;

    std::tie(ocvpi, flag_ocvpi) = slide::linInt_noexcept(bound, OCVp.x, OCVp.y, OCVp.size(), sp);
    std::tie(ocvni, flag_ocvni) = slide::linInt_noexcept(bound, OCVn.x, OCVn.y, OCVn.size(), sn);

    return ocvpi - ocvni; //!< OCV = OCV_cathode - OCV_anode, cell voltage in this time step
  };

  auto store = [&]() {
    OCV.x.push_back(Ah);           //!< discharge charge
    OCV.y.push_back(v);            //!< open circuit voltage
    OCVanode.x.push_back(sn);      //!< anode lithium fraction
    OCVanode.y.push_back(ocvni);   //!< anode potential
    OCVcathode.x.push_back(sp);    //!< cathode lithium fraction
    OCVcathode.y.push_back(ocvpi); //!< cathode potential
  };

  auto step = [&]() {
    for (size_t i = 0; (!end && i < nstore); i++) {
      v = interpolate();

      Ah += I * dt / 3600;
      sp += I * dt / (n * F) / AMp / cmaxp;
      sn -= I * dt / (n * F) / AMn / cmaxn;

      end = v <= Vend || sp <= 0 || sp >= 1 || sn <= 0 || sn >= 1;
      //!< if (sp_n != 1 && sn_n != 1)
      //!< 	std::cout << "V = " << V << " ocvpi = " << ocvpi << " ocvni = " << ocvni << " sp = " << sp << " sn = " << sn << '\n';
    }
  };

  //!< Store the lithium fractions at the start (=100% SOC)
  fp[0] = sp;
  fn[0] = sn;

  //!< *********************************************************** 2 discharge & measure voltage ***********************************************************************
  //!< std::cout << "===========================\n";
  v = interpolate();
  store();
  while (!end) {
    step();
    store();
  }

  //!< //!< *********************************************************** 3 output parameters ***********************************************************************

  //!< Store the lithium fractions at the end (= 0% SOC)
  fp[2] = sp;
  fn[2] = sn;

  //!< Then the lithium fractions at 50% SOC are the average between the fractions at 100% and 0%
  //!< This is true because SOC is defined based on charge throughput, which is directly (linearly) linked with the lithium concentration
  fp[1] = (fp[0] + fp[2]) / 2.0;
  fn[1] = (fn[0] + fn[2]) / 2.0;

  return flag; //!< Flag is never false.
}

auto discharge_noexcept(const slide::XYdata_vv &OCVp, const slide::XYdata_vv &OCVn, double cap, const double AMp, const double AMn, const double cmaxp,
                        const double cmaxn, double sp, double sn, double Vend, slide::XYdata_vv &OCV, double fp[], double fn[])
{
  /*
   * Function to simulate a CC discharge at a very low current, so we can ignore resistances and diffusion.
   *
   * It assumes the concentrations are uniform at all times, so we don't have to account for diffusion.
   * Therefore, we don't need to use the full cell model.
   * Rather, we can simulate two 'buckets of lithium' and we only need to exchange lithium between both.
   *
   * IN
   * OCVp 	cathode OCV curve
   * OCVn 	anode OCV curve
   * cap		capacity of the cell [Ah]
   * AMp 		amount of active material on the cathode [m3], elec_surf * thicp * ep
   * AMn 		amount of active material on the anode [m3], elec_surf * thickn * en
   * cmaxp 	maximum li-concentration in the cathode [mol m-3]
   * cmaxn	maximum li-concentration in the anode [mol m-3]
   * sp 		start li-fraction at the cathode (between 0 and 1)
   * sn 		start li-fraciton at the anode (between 0 and 1)
   * Vend 	voltage when the OCV discharge should end [V]
   *
   * OUT
   * OCV 		OCV curve of the cell simulated with the given parameters
   * 			the first column [i][0] gives the discharged charge in [Ah]
   * 			the second column [i][1] gives the cell OCV in [V]
   * OCVanode OCV of the anode
   * 			the first column [i][0] gives the lithium fraction in the anode
   * 			the second column [i][1] gives the voltage of the anode
   * OCVcathode OCV of the cathode
   * 			the first column [i][0] gives the lithium fraction in the cathode
   * 			the second column [i][1] gives the voltage of the cathode
   * fp		array with the lithium fraction in the cathode at 100% SOC, 50% SOC and 0% SOC
   * fn	 	array with the lithium fraction in the anode at 100% SOC, 50% SOC and 0% SOC
   */

  //!< *********************************************************** 1 variables ***********************************************************************
  using namespace PhyConst;
  //!< Variables
  constexpr double n = 1.0;    //!< number of electrons involved in the reaction
  double Ah = 0;               //!< discharged charge up to now
  const double I = 0.1 * cap;  //!< magnitude of the discharge current [A]
  const double dt = 5;         //!< time step to use in the simulation [s]
  constexpr bool bound = true; //!< check the surface concentration is within the allowed limits when doing linear interpolation.
  bool end = false;            //!< boolean to check if we have reached the end voltage

  double v, ocvpi, ocvni;
  bool flag{ true };

  //!< If we store a data point for every time step, the output arrays will be very long.
  //!< Therefore, store only one value every 'nstore' time steps
  constexpr size_t nstore = 100; //!< store one data point every 100 time steps

  //!< Store the lithium fractions at the start (=100% SOC)
  fp[0] = sp;
  fn[0] = sn;

  //!< *********************************************************** 2 discharge & measure voltage ***********************************************************************

  auto interpolate = [&]() {
    bool status_ocvpi, status_ocvni;
    std::tie(ocvpi, status_ocvpi) = slide::linInt_noexcept(bound, OCVp.x, OCVp.y, OCVp.size(), sp);
    std::tie(ocvni, status_ocvni) = slide::linInt_noexcept(bound, OCVn.x, OCVn.y, OCVn.size(), sn);

    if (status_ocvpi || status_ocvni) {
      std::cout << "Error in deterimineOCV::discharge from linInt when getting the voltage. positive li-fraction is "
                << sp << " negative li-fraction is " << sn << ". Both should be between 0 and 1.\n";
      flag = false;
    }

    v = ocvpi - ocvni; //!< OCV = OCV_cathode - OCV_anode, cell voltage in this time step
  };

  auto store = [&]() {
    OCV.x.push_back(Ah); //!< discharge charge
    OCV.y.push_back(v);  //!< open circuit voltage
  };

  auto step = [&]() {
    for (size_t i = 0; (flag && !end && i < nstore); i++) {
      interpolate();

      Ah += I * dt / 3600;
      sp += I * dt / (n * F) / AMp / cmaxp;
      sn -= I * dt / (n * F) / AMn / cmaxn;

      end = v <= Vend || sp <= 0 || sp >= 1 || sn <= 0 || sn >= 1;
    }
  };

  interpolate();
  while (flag && !end) {
    store();
    step();
  }

  if (!flag)
    return flag;

  store();

  //!< //!< *********************************************************** 3 output parameters ***********************************************************************

  //!< Store the lithium fractions at the end (= 0% SOC)
  fp[2] = sp;
  fn[2] = sn;

  //!< Then the lithium fractions at 50% SOC are the average between the fractions at 100% and 0%
  //!< This is true because SOC is defined based on charge throughput, which is directly (linearly) linked with the lithium concentration
  fp[1] = (fp[0] + fp[2]) / 2.0;
  fn[1] = (fn[0] + fn[2]) / 2.0;

  return flag;
}

void discharge(const slide::XYdata_vv &OCVp, const slide::XYdata_vv &OCVn, double cap, const double AMp, const double AMn, const double cmaxp,
               const double cmaxn, double sp, double sn, double Vend, slide::XYdata_vv &OCV, slide::XYdata_vv &OCVanode,
               slide::XYdata_vv &OCVcathode, double fp[], double fn[])
{
  /*
   * Function to simulate a CC discharge at a very low current, so we can ignore resistances and diffusion.
   * Throwing version of discharge_noexcept, see "discharge_noexcept" for further information.
   */

  //!< *********************************************************** 1 variables ***********************************************************************

  const auto flag = discharge_noexcept(OCVp, OCVn, cap, AMp, AMn, cmaxp, cmaxn, sp, sn, Vend, OCV, OCVanode, OCVcathode, fp, fn);

  if (!flag) {
    throw 1; //!< It can only throw due to IntLin therefore, its error code it 1.
  }
}

void readOCVinput(const std::string &namepos, const std::string &nameneg, const std::string &namecell,
                  slide::XYdata_vv &OCVp, slide::XYdata_vv &OCVn, slide::XYdata_vv &OCVcell)
{
  /*
   * Function to read the csv files witht the OCV curves which the user has to supply as input.
   *
   * As input, the user has to supply the electrode OCV curves.
   * These must come in the form of two csv files, each with two columns.
   * The first column of each file must give the lithium fraction in increasing order (i.e. going from 0 to 1)
   * The second column of each file must give the voltage vs li/li+ of the electrode at the corresponding li fraction
   *
   * Additionally, the user has to supply the cell OCV curve.
   * This must be one csv file with 2 columns.
   * The first column gives the discharged charge [Ah], starting from 0 (i.e. cell is fully charged) to c (i.e. the cell is fully discharged) [c is the capacity which can be discharged]
   * The second column gives the OCV of the cell [V].
   *
   * IN
   * namepos 	name of the CSV file with the OCV curve of the cathode
   * nameneg 	name of the CSV file with the OCV curve of the anode
   * namecell name of the CSV file with the OCV curve of the cell
   * np 		number of data points in the file with the cathode OCV
   * nn 		number of data points in the file with the anode OCV
   * ncell 	number of data points in the file with the cell OCV
   *
   * OUT
   * OCVp 	matrix with the cathode OCV
   * OCVn 	matrix with the anode OCV
   * OCVcell 	matrix with the cell OCV curve
   *
   * THROWS
   * 10000 	one of the files has an illegal format
   */

  try {
    OCVp.setCurve(PathVar::data / namepos);
    OCVn.setCurve(PathVar::data / nameneg);
    OCVcell.setCurve(PathVar::data / namecell);
  } catch (int e) {
    //!< std::cout << "Throw test: " << 80 << '\n';
    std::cerr << "ERROR in determineOCV::readOCVinput, an error " << e << " happened while reading the files. Throwing the error on.\n";
    throw e;
  }

  //!< Check the OCV curves are valid
  const bool valp = validOCV(true, OCVp); //!< check the cathode OCV curve. The first column must go from 0-1
  if (!valp) {
    std::cerr << "ERROR in determineOCV::readOCVinput, the file with the cathode OCV curve " << namepos << " has an illegal format.\n";
    throw 10000;
  }

  const bool valn = validOCV(true, OCVn); //!< check the anode OCV curve. The first column must go from 0-1
  if (!valn) {
    std::cerr << "ERROR in determineOCV::readOCVinput, the file with the anode OCV curve " << nameneg << " has an illegal format.\n";
    throw 10000;
  }

  const bool valcell = validOCV(false, OCVcell); //!< check the cell OCV curve. The first column must start at 0 but can end anywhere
  if (!valcell) {
    std::cerr << "ERROR in determineOCV::readOCVinput, the file with the cell OCV curve " << namecell << " has an illegal format.\n";
    throw 10000;
  }
}

double calculateError(bool bound, slide::XYdata_vv &OCVcell, slide::XYdata_vv &OCVsim)
{
  /*
   * Function to calculate the root mean square error between the OCV curve of the cell supplied by the user and the simulated OCV curve
   *
   * IN
   * bound	boolean deciding what to do if the value of x is out of range of xdat for linear interpolation
   * 			i.e. how to 'extend' OCVsim to the same capacity as OCVcell if OCVsim has a lower capacity
   * 				if true, the value will be set to 0
   * 				if false, the value will be set to the last point of OCVsim (e.g. 2.7)
   * OCVcell 	OCV curve of the cell, 2 columns
   * OCVsim 	simulated OCV curve, 2 columns
   *
   * OUT
   * rmse		RMSE between both curves
   */

  double rmse = 0; //!< root mean square error

  //!< loop through all data points
  for (size_t i = 0; i < OCVcell.size(); i++) {
    //!< Get the simulated OCV at the discharged charge of this point on the measured OCV curve of the cell
    auto [Vsimi, status] = slide::linInt_noexcept(bound, OCVsim.x, OCVsim.y, OCVsim.size(), OCVcell.x[i]);
    //!< if status is not 0 then the simulated voltage is already set to 0;

    const double err = OCVcell.y[i] - Vsimi; //!< Calculate the error
    rmse += err * err;                       //!< sum ( (Vcell[i] - Vsim[i])^2, i=0..ncell )
  }

  //!< Calculate the RMSE
  rmse = std::sqrt(rmse / OCVcell.size());
  return rmse;
}

auto cost_OCV(const slide::XYdata_vv &OCVp, const slide::XYdata_vv &OCVn, const double AMp, const double AMn, double sp, double sn, const double cmaxp,
              const double cmaxn, const slide::XYdata_vv &OCVcell)
{

  using namespace PhyConst;
  //!< Variables
  constexpr double n = 1.0;    //!< number of electrons involved in the reaction
  double Ah = 0;               //!< discharged charge up to now
  constexpr bool bound = true; //!< check the surface concentration is within the allowed limits when doing linear interpolation.
  bool end = false;            //!< boolean to check if we have reached the end voltage

  const auto Vend = OCVcell.y.back();

  double sqr_err = 0; //!< square error
  double v = 0;       //!< Simulation OCV.

  for (size_t i = 0; i < OCVcell.size(); i++) {
    //!< Original equations:
    //!< Ah += I * dt / 3600;  -> I*dt = dAs
    //!< sp += I * dt / (n * F) / AMp / cmaxp;
    //!< sn -= I * dt / (n * F) / AMn / cmaxn;

    if (!end) {

      const auto Ah_cell = OCVcell.x[i];

      const auto dAh = Ah_cell - Ah; //!< Ah need to be added to reach cell's Ah.
      const auto dAs = dAh * 3600.0; //!< Ah += I * dt / 3600;  I*dt = dAs, Amper seconds

      Ah = Ah_cell;
      sp += dAs / (n * F) / AMp / cmaxp;
      sn -= dAs / (n * F) / AMn / cmaxn;

      const auto [ocvpi, status_ocvpi] = slide::linInt_noexcept(bound, OCVp.x, OCVp.y, OCVp.size(), sp); //!< #TODO: convert to interp
      const auto [ocvni, status_ocvni] = slide::linInt_noexcept(bound, OCVn.x, OCVn.y, OCVn.size(), sn);

      if (status_ocvpi || status_ocvni)
        v = 0;
      else
        v = ocvpi - ocvni; //!< OCV = OCV_cathode - OCV_anode, cell voltage in this time step
    } else {
      v = 0;
    }

    const double err = OCVcell.y[i] - v; //!< Calculate the error
    sqr_err += err * err;                //!< sum ( (Vcell[i] - Vsim[i])^2, i=0..ncell )

    end = v <= Vend || sp < 0 || sp > 1 || sn < 0 || sn > 1;
  }

  return sqr_err;
}

void fitAMnAndStartingPoints(int hierarchy, int ap, double AMp, slide::FixedData<double> AMn_space, slide::FixedData<double> sp_space,
                             slide::FixedData<double> sn_space, double cmaxp, double cmaxn, double *err, std::array<double, 4> &par,
                             slide::XYdata_vv &OCVp, slide::XYdata_vv &OCVn, slide::XYdata_vv &OCVcell)
{ /*
   * Function to scan a given search space for 3 parameters
   * (amount of anode active material, starting point for the cathode and starting point for the anode).
   * For each combination of the 3 parameters (and the given amount of cathode active material),
   * the OCV curve is simulated and the error with the measured OCV curve is calculated.
   * The best fit (i.e. three parameters with the lowest error on the OCV curve) is returned
   *
   * IN
   * hierarchy 		level of the hierarchy in which we are
   * ap 				index in the search space for the cathode active material
   * AMp 				amount of cathode active material [m3]
   * cmaxp 			maximum lithium concentration in the cathode [mol m-3]
   * cmaxn			maximum lithium concentration in the anode [mol m-3]
   *
   * OUT
   * err 				lowest error for this amount of cathode active material
   * par 				parameters giving the lowest error (AMp, AMn, sp, sn)
   */

  //!< *********************************************************** 1 variables ***********************************************************************

  //!< variables
  double errmin = 1e10; //!< lowest error encountered so far

  //!< *********************************************************** 2 loop through the search space ***********************************************************************
  for (const auto AMn : AMn_space)   //!< loop for the search space of AMn
    for (const auto sp : sp_space)   //!< loop through the search space of sp / avoid that the starting points go out of range (the lithium fraction has to be between 0 and 1)
      for (const auto sn : sn_space) //!< loop through the search space of sn / avoid that the starting points go out of range (the lithium fraction has to be between 0 and 1)
      {
        //!< Simulate the OCV curve with these parameters and calculate the error between the simulated and measured OCV curve
        const auto erri = cost_OCV(OCVp, OCVn, AMp, AMn, sp, sn, cmaxp, cmaxn, OCVcell); //!< calculates total squared error.

        //!< Store the minimum error & parameters leading to this error
        if (erri < errmin) { //!< check if the error of this combination is lower than the best fit so far
          par = { AMp, AMn, sp, sn };
          errmin = erri;
        }
      }

  //!< Return the minimum error
  *err = std::sqrt(errmin / OCVcell.size()); //!< RMSE error.
}

auto hierarchicalOCVfit(int hmax, slide::FixedData<double> AMp_space, slide::FixedData<double> AMn_space, slide::FixedData<double> sp_space,
                        slide::FixedData<double> sn_space, std::string namepos, std::string nameneg, std::string namecell, double cmaxp, double cmaxn)
{
  /*
   * Hierarchical search algorithm to converge on the best fit. For a convex problem, the optimal point is found.
   * For a nonconvex problem, the algorithm might find a local minimum only.
   *
   * It iteratively 'zooms in' on the optimal point in the search space.
   * I.e. you first call it with a large range and a large step size for each parameter
   * 		it then finds the optimal combination (with this large step)
   * 		the search space is then refined to a region around the optimal point, with a smaller range and a smaller step size
   * 		the second iteration then finds the optimal point in this new (smaller) region (with the smaller step)
   * 		the search space is again refined to the region around this 2nd optimal point, with an even smaller range, and even smaller steps
   * 		etc.
   * The region is always refined to the points before and after the optimal point for each parameter.
   * E.g. if the optimal value is P and we were scanning with a step size of dp, then the new range in the next level is P-dp to P+dp
   * And the number of steps remains the same, so the new step size is ( (P+dp) - (P-dp) ) / (number of steps - 1) = 2dp/(step-1)
   *
   * So if the function is called with a step number nstep and an initial step size of dp,
   * then in hierarchical level n, the step size is dp* (2/(nstep-1))^(n-1)
   * So after 'hmax' level, the accuracy is dp* (2/(nstep-1))^(hmax-1)
   *
   * IN
   * hmax  		number of hierarchical steps we should take.
   * 				the higher the value, the higher the accuracy of the fit will be but the longer the calculation will take
   * namepos 		ame of the CSV file with the cathode OCV curve
   * nameneg 		name of the CSV file with the anode OCV curve
   * namecell		name of the CSV file with the cell's OCV curve
   * cmaxp 		maximum lithium concentration in the cathode [mol m-3]
   * cmaxn		maximum lithium concentration in the anode [mol m-3]
   *
   * OUT
   * err 			lowest error for this amount of cathode active material
   * par 			parameters giving the lowest error (AMp, AMn, sp, sn)
   */

  std::vector<std::array<double, 4>> par_arr(AMp_space.size()); //!< parameters [AMp AMn sp sn] giving the lowest error for that amount of cathode active material
  std::vector<double> err_arr(AMp_space.size());

  std::array<double, 4> par; //!< Optimal parameters.
  double err;                //!< lowest error.

  slide::XYdata_vv OCVp(100), OCVn(100), OCVcell(100); //!< Temp vectors
  readOCVinput(namepos, nameneg, namecell, OCVp, OCVn, OCVcell);

  //!< Loop for each level in the search
  for (int h = 0; h < hmax; h++) {

    //!< print the search space of this level
    std::cout << "Start hierarchy level " << h << " with the following search spaces: \n"
              << "AMp: from " << AMp_space.front() << " to " << AMp_space.back() << " in " << AMp_space.size() << " steps with magnitude " << AMp_space.dstep() << '\n'
              << "AMn: from " << AMn_space.front() << " to " << AMn_space.back() << " in " << AMn_space.size() << " steps with magnitude " << AMn_space.dstep() << '\n'
              << "sp: from " << sp_space.front() << " to " << sp_space.back() << " in " << sp_space.size() << " steps with magnitude " << sp_space.dstep() << '\n'
              << "sn: from " << sn_space.front() << " to " << sn_space.back() << " in " << sn_space.size() << " steps with magnitude " << sn_space.dstep() << '\n';

    auto task_indv = [&](int i) {
      const double AMp1 = AMp_space[i]; //!< amount of cathode active material for each
      fitAMnAndStartingPoints(h, i, AMp1, AMn_space, sp_space, sn_space, cmaxp, cmaxn, &err_arr[i], par_arr[i], OCVp, OCVn, OCVcell);
    };

    slide::run(task_indv, AMp_space.size()); //!< Loop through the search space of AMp

    const auto minIndex = std::min_element(err_arr.begin(), err_arr.end()) - err_arr.begin();
    par = par_arr[minIndex]; //!< Make the output parameters
    err = err_arr[minIndex];

    writeOCVParam(h, par); //!< Print the best fit, and write in a CSV file

    //!< Calculate the best fit in this level
    //!< Update the search space

    const auto [AMp, AMn, sp, sn] = par; //!< the best fit parameters in this level of the hierarchy

    AMp_space = slide::linspace_fix(AMp_space.prev(AMp), AMp_space.next(AMp), AMp_space.size());
    AMn_space = slide::linspace_fix(AMn_space.prev(AMn), AMn_space.next(AMn), AMn_space.size());
    sp_space = slide::linspace_fix(std::max(sp_space.prev(sp), 0.0), std::min(sp_space.next(sp), 1.0), sp_space.size()); //!< min to limit it to 1<.
    sn_space = slide::linspace_fix(std::max(sn_space.prev(sn), 0.0), std::min(sn_space.next(sn), 1.0), sn_space.size()); //!< max to limit it to >0
  }

  return std::pair(par, err); //!< Return parameters and error.
}

void estimateOCVparameters() // #TODO this function is slow and hand-tuned. Change with determining sn and AMp from boundary.
{
  /*
   * Function which will find the parameters which best fit the OCV curve of the user.
   *
   * As input, the user has to supply the electrode OCV curves.
   * These must come in the form of two csv files, each with two columns.
   * The first column of each file must give the lithium fraction in increasing order (i.e. going from 0 to 1)
   * The second column of each file must give the voltage vs li/li+ of the electrode at the corresponding li fraction
   *
   * Additionally, the user has to supply the measured cell OCV curve.
   * This must be one csv file with 2 columns.
   * The first column gives the discharged charge [Ah], starting from 0 (i.e. cell is fully charged) to c (i.e. the cell is fully discharged) [c is the capacity which can be discharged]
   * The second column gives the OCV of the cell [V].
   *
   * The function finds optimal values for 4 parameters:
   * 		cinip 		the initial lithium fraction for the cathode
   * 		cinin 		the initial lithium fraction for the anode
   * 		thickp		the thickness of the cathode
   * 		thickn		the thickness of the anode
   *
   * The values for the other geometric parameters are assumed known
   * 		elec_surf 	the geometric surface area of the electrode (assumed to be 0.982 m2, which is representative for an 18650 cell)
   * 		ep			the volume fraction of active material on the cathode (assumed to be 50%)
   * 		en 			the volume fraction of active material on the anode (assumed to be 50%)
   * If you have different values, you have to change them in the code below (lines 791-793)
   *
   * The values of the parameters which give the best fit are written in a CSV file.
   * Also, the simulated OCV curve with these parameters is written to a csv file so it can be compared with the measured one.
   * 		The first column gives the simulated discharged charge
   * 		The second column gives the simulated cell OCV
   * 		The third column gives the simulated anode lithium fraction
   * 		The fourth column gives the simulated anode OCV (vs Li/Li+)
   * 		The fifth column gives the simulated cathode lithium fraction
   * 		The sixth column gives the simulated cathode OCV (vs Li/Li+)
   */

  //!< *********************************************************** 1 USER INPUT ***********************************************************************
  using namespace PhyConst;

  //!< input parameters
  std::string namepos = "OCVfit_cathode.csv"; //!< name of the file with the OCV curve of the cathode
  const int np = 49;                          //!< number of data points on the cathode OCV curve
  std::string nameneg = "OCVfit_anode.csv";   //!< name of the file with the OCV curve of the anode
  const int nn = 63;                          //!< number of data points on the anode OCV curve
  std::string namecell = "OCVfit_cell.csv";   //!< name of the file with the OCV curve of the cell
  const int ncell = 74;                       //!< number of data points on the cell OCV curve
  double cmaxp = 51385;                       //!< maximum li-concentration in the cathode [mol m-3]
  double cmaxn = 30555;                       //!< maximum li-concentration in the anode [mol m-3]

  //!< names of the output files
  std::string nameparam = "OCVfit_parameters.csv"; //!< name of the output file in which the parameters of the best fit will be written
  std::string nameOCV = "OCVfit_sim.csv";          //!< name of the output file in which the OCV curve with the best parameters will be written

  //!< Read the OCV curves
  slide::XYdata_vv OCVp(np), OCVn(nn), OCVcell(ncell);
  try {
    readOCVinput(namepos, nameneg, namecell, OCVp, OCVn, OCVcell);
  } catch (int e) {
    //!< std::cout << "Throw test: " << 81 << '\n';
    std::cerr << "Error in determineOCV::estimateOCVparameters, the input files have the wrong format.\n";
    return; //!< stop calculating because we can't do anything
  }

  //!< ****************************************** 2 define the search space for fitting parameters ***********************************************************************

  //!< We can play with 4 parameters:
  //!< AMp 	the amount of active material on the positive electrode
  //!< AMn 	the amount of active material on the negative electrode
  //!< sp	the starting point on the positive electrode OCV curve
  //!< sn 	the starting point on the negative electrode OCV curve

  //!< estimate the amount of active material needed to reach the cell capacity:
  constexpr double n = 1.0;                                //!< number of electrons involved in the reaction
  const double cap = OCVcell.x.back();                     //!< capacity of the cell [Ah] (the discharge Ah at the last point on the OCV curve)
  const double AMp_guess = cap * 3600.0 / (n * F * cmaxp); //!< the total charge in an electrode in Ah is given by; n * F * cmax * AM / 3600
  const double AMn_guess = cap * 3600.0 / (n * F * cmaxn);

  //!< //!< Define the search space for the amount of active material on each electrode
  constexpr double step = 0.02; //!< take steps of 5% of the guessed active material
  constexpr double AMmax = 2.5; //!< the maximum amount of active material is 5 times the guessed amount
  constexpr double AMmin = 1;   //!< the minimum amount of active material is 0, actually should not be zero since it is a divisor.

  //!< Define the search space for the initial lithium fractions at each electrode
  auto sp_space = slide::linspace_fix(0.38, 0.4, 51); //!< From 0 to 1 lithium fraction 5% increase of the lithium fraction
  auto sn_space = slide::linspace_fix(0.55, 0.58, 51);

  //!< Define the search space for the amount of active material on each electrode
  auto AMp_space = slide::linspace_fix(AMmin * AMp_guess, AMmax * AMp_guess, 100); //!< From 0 to 5x of guessed active material with 10% steps.
  auto AMn_space = slide::linspace_fix(AMmin * AMn_guess, AMmax * AMn_guess, 100);


  //!< ***************************************************** 3 Fit the parameters ***********************************************************************

  //!< Call the hierarchical search algorithm, which does the fitting
  constexpr int hmax = 3;                                                                                                               //!< number of levels in the hierarchy to consider.
  const auto [par, err] = hierarchicalOCVfit(hmax, AMp_space, AMn_space, sp_space, sn_space, namepos, nameneg, namecell, cmaxp, cmaxn); //!< parameters of the best fit and lowest error.

  //!< ***************************************************** 4 write outputs ***********************************************************************

  //!< Print the best fit
  std::cout << "The best fit is: AMp = " << par[0] << ", AMn = " << par[1] << ", sp = " << par[2] << ", sn = " << par[3] << ".\n";
  std::ofstream output;

  //!< write the parameters in a csv file
  output.open(PathVar::results / nameparam, std::ios_base::out);
  output << "AMp" << ',' << par[0] << '\n'
         << "AMn" << ',' << par[1] << '\n'
         << "start pos" << ',' << par[2] << '\n'
         << "start neg" << ',' << par[3] << '\n';
  output.close();

  //!< Simulate the best-fit OCV curve
  const double Vend = OCVcell.y.back(); //!< minimum voltage of the OCV curve

  slide::XYdata_vv OCVsim, OCVnsim, OCVpsim;                       //!< array for the simulated OCV curve and voltage of each electrode
  constexpr int nin = 200;                                         //!< size hint for the array to store the simulated OCV curve
  OCVsim.reserve(nin), OCVnsim.reserve(nin), OCVpsim.reserve(nin); //!< reserve some size.

  double fp[3], fn[3]; //!< arrays to store the lithium fractions at 100%, 50% and 0% SOC
  discharge(OCVp, OCVn, cap, par[0], par[1], cmaxp, cmaxn, par[2], par[3], Vend, OCVsim, OCVnsim, OCVpsim, fp, fn);
  //!< Write this best-fit OCV curve in a csv file so the user can check it using the matlab script readEstimateOCV.m
  output.open(PathVar::results / nameOCV);
  for (size_t i = 0; i < OCVsim.x.size(); i++)
    output << OCVsim.x[i] << ',' << OCVsim.y[i] << ',' << OCVnsim.x[i] << ','
           << OCVnsim.y[i] << ',' << OCVpsim.x[i] << ',' << OCVpsim.y[i] << '\n';
  output.close();

  //!< From the 4 fitted values, we need to determine the following parameters, needed by the single particle model implemented in Cell
  //!< elec_surf 	the geometric surface area of the electrode
  //!< thickp		the thickness of the cathode
  //!< thickn		the thickness of the anode
  //!< ep			the volume fraction of active material on the cathode
  //!< en 			the volume fraction of active material on the anode
  //!< cinip 		the initial lithium fraction for the cathode
  //!< cinin 		the initial lithium fraction for the anode

  //!< Assume the volume fractions are 50% and the electrode surface is 0.0982 m2.
  //!< Then we can calculate the thickness of the electrodes to give the desired amount of active material, which is given by:
  //!< 		AM = elec_surf * thick * e
  constexpr double ep = 0.5;
  constexpr double en = 0.5;
  constexpr double elec_surf = 0.0982;
  const double thickp = par[0] / (elec_surf * ep);
  const double thickn = par[1] / (elec_surf * en);

  //!< Append these parameters in the csv file where we had written the fitted parameters
  output.open(PathVar::results / nameparam, std::ios_base::app);
  output << "cathode volume fraction ep" << ',' << ep << '\n';
  output << "anode volume fraction en" << ',' << en << '\n';
  output << "electrode surface ele_surf" << ',' << elec_surf << '\n';
  output << "cathode thickness thickp" << ',' << thickp << '\n';
  output << "anode thickness thickn" << ',' << thickn << '\n';
  output << "cathode lithium fraction at 100% SOC" << ',' << fp[0] << '\n';
  output << "anode lithium fraction at 100% SOC" << ',' << fn[0] << '\n';
  output << "cathode lithium fraction at 50% SOC" << ',' << fp[1] << '\n';
  output << "anode lithium fraction at 50% SOC" << ',' << fn[1] << '\n';
  output << "cathode lithium fraction at 0% SOC" << ',' << fp[2] << '\n';
  output << "anode lithium fraction at 0% SOC" << ',' << fn[2] << '\n';
  output << "capacity of the cell in Ah" << ',' << cap << '\n';
  output << "maximum voltage of the cell" << ',' << OCVcell.y[0] << '\n';
  output << "minimum voltage of the cell" << ',' << OCVcell.y.back() << '\n';
  output << "The error on the OCV curve with this fit is" << ',' << err << '\n';

  //!< then note down the settings used to get this fit
  output << '\n';
  output << "Below are the settings which produced this result\n";
  output << "maximum li-concentration in the cathode" << ',' << cmaxp << '\n';
  output << "maximum li-concentration in the anode" << ',' << cmaxn << '\n';
  output << "name of the file with the cathode OCV curve" << ',' << namepos << '\n';
  output << "name of the file with the anode OCV curve" << ',' << nameneg << '\n';
  output << "name of the file with the cell's OCV curve" << ',' << namecell << '\n';
  output << "capacity of the cell" << ',' << cap << '\n';

  //!< the search space
  output << '\n';
  output << "Below are the settings of the initial search space\n";
  output << "relative step size in the search for active material" << ',' << step << '\n';
  output << "relative minimum amount of active material" << ',' << AMmin << '\n';
  output << "relative maximum amount of active material" << ',' << AMmax << '\n';
  output << "step size in the search for the starting li-fraction" << ',' << sp_space.dstep() << '\n';
  output << "minimum li-fraction" << ',' << 0.0 << '\n';
  output << "maximum li-fraction" << ',' << 1.0 << '\n';
  output << "number of levels in the search hierarchy" << ',' << hmax << '\n';
  output.close();
}

void writeOCVParam(int h, const std::array<double, 4> &par)
{
  //!< Print the best fit, and write in a CSV file
  std::cout << "The best fit in hierarchy " << h << " is: AMp = " << par[0]
            << ", AMn = " << par[1] << ", sp = " << par[2] << ", sn = " << par[3] << ".\n";

  std::ofstream output; //!< write the parameters
  const auto na = "OCVFit_" + std::to_string(h) + "_param.csv";
  output.open(PathVar::results / na, std::ios_base::out);
  output << "AMp" << ',' << par[0] << '\n';
  output << "AMn" << ',' << par[1] << '\n';
  output << "start pos" << ',' << par[2] << '\n';
  output << "start neg" << ',' << par[3] << '\n';
  output.close();
}

} // namespace slide