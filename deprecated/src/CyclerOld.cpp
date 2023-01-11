/*
 * cycler.cpp
 *
 * Class to implement a 'cycler', which extends a 'basic cycler'.
 * A CyclerOld consists of the programs one would write on battery cycles.
 * I.e. it defines the procedures for loading a cell in a degradation experiment.
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "CyclerOld.hpp"
#include "../settings/settings.hpp"
#include "../utility/utility.hpp"

#include <vector>
#include <array>
#include <numeric>

void CyclerOld::getOCV(slide::FixedData<double> &Ah_fixed, std::vector<double> &OCVp, std::vector<double> &OCVn)
{
  /*
   * function to calculate the half cell OCV curves.
   * Each curve is fully simulated (i.e. from a lithium fraction of 0 to a lithium fraction of 1).
   * Both curves are aligned as they are in the actual cell, and the common x-axis gives the discharged charge [Ah].
   * The discharged charge is '0' when the cell is charged to its maximum voltage, and negative for all points with a higher cell voltage.
   * When one electrode has reached the full lithium concentration (but the other electrode can still further (dis)charge, its potential is kept constant.
   *
   * Degradation is never accounted for while recording the OCV curves because the code does not support this
   * After the OCV curves are recorded, the original states are restored
   *
   * IN	 *
   * OUT
   * Ah 		discharged charge, 0 when the cell OCV is at its maximum voltage, negative when the cell OCV is above its maximum [Ah]
   * OCVp		open circuit potential of the positive electrode, constant after we have reached the extreme points of the cathode OCV curve [V]
   * OCVn		open circuit potential of the negative electrode, constant after we have reached the extreme points of the cathode OCV curve [V]
   *
   * THROWS
   * 1002		the arrays provided for feedback are too short (ninocv < nout)
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::getOCV is starting\n";

  //!<*********************************************************** 1 variables & settings ***********************************************************************

  //!< variables
  slide::State sini, sCharged;                    //!< initial & fully charged battery state
  double Iini, ICharged;                          //!< initial & fully charged battery current
  double Tenvi, Trefi;                            //!< initial environmental & reference temperatures
  constexpr double Crate = 0.04;                  //!< Crate for the CC phase.
  double Ccutcha = 0.005;                         //!< C rate of the cut-off current during CV charge
  const double Icha = -c.getNominalCap() * Crate; //!< charging current [A]
  const double Idis = c.getNominalCap() * Crate;  //!< discharging current [A]
  constexpr double dt = 2;                        //!< time step for time integration
  bool blockDegradation = true;                   //!< don't account for degradation while we measure the OCV curves
  double ah, wh, tt;                              //!< unneeded feedback variables
  double ah_neg_p;                                //!< charge on the cathode OCV curve with a negative x-value (i.e. the points with a higher OCVp than the OCVp at the fully charged condition)
  double ah_neg_n;                                //!< charge on the anode OCV curve with a negative x-value (i.e. the points with a lower OCVn than the OCVn at the fully charged condition)
  bool bound = false;                             //!< in linear interpolation return the closest point if you are out of range
                                                  //!<	this ensures the electrode OCV remains constant after it has reached the full lithium concentration

  //!< Data collection: we don't want to store a point every time step because that leads to very large arrays.
  constexpr int nStore = 750;                                       //!< store a value every 750 time steps
  constexpr int nin2 = static_cast<int>(1 / Crate * 3600 / dt * 2); //!< expected duration 1/Crate [hours] * 3600/dt [time steps], *2 such that the arrays are for sure long enough //!<90k elements

  static thread_local std::vector<double> OCVni_vec, OCVpi_vec; //!<~50k and ~79k elements respectively.

  OCVni_vec.reserve(nin2), OCVpi_vec.reserve(nin2); //!< Reserve estimated number of elements.

  OCVni_vec.clear(); //!< clear elements of vectors.
  OCVpi_vec.clear();

  //!< Store the initial battery states so we can restore them later
  c.getStates(sini, &Iini);

  //!< Set the cell and environmental temperature to the reference temperature
  c.getTemperatures(&Tenvi, &Trefi);
  c.setTenv(Trefi);
  c.setT(Trefi);

  //!<*********************************************************** 2 get the half-cell curves ***********************************************************************

  try {
    //!< fully charge battery to its specified maximum voltage (Cell::getVmax())
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::getOCV is fully charging the cell to the maximum cell voltage.\n";

    CC_V_CV_I(Crate, c.getVmax(), Ccutcha, dt, blockDegradation, &ah, &wh, &tt); //!< do a CC charge until the maximum voltage, followed by a CV at this voltage

    //!< Store this battery state
    c.getStates(sCharged, &ICharged);

    //!< We can cycle each electrode separately using the function CC_halfCell_full.
    //!< This function returns the electrode OCV (the resistive voltage drop is just ignored)
    //!< and the electrode is cycled until the lithium is exhausted (i.e. li-concentration is 0 or the maximum concentration)

    //!< fully cycle the cathode until it reaches its extreme lithium concentrations
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::getOCV is fully charging the cathode.\n";

    CC_halfCell_full(Icha, dt, true, OCVpi_vec, &ah_neg_p, false); //!< ah_neg_p contains the capacity we charged additional (starting from the 'fully charged' state) on the positive electrode

    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::getOCV is fully discharging the cathode.\n";

    CC_halfCell_full(Idis, dt, true, OCVpi_vec, &ah, true); //!< now OCVpi contains the cathode OCV for the full lithium concentration range (from 0 to the maximum concentration)

    //!< restore the fully charged battery state (undo the cycling of the cathode)
    c.setStates(sCharged, ICharged);

    //!< fully cycle the anode until it reaches its extreme lithium concentrations
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::getOCV is fully charging the anode.\n";

    CC_halfCell_full(Icha, dt, false, OCVni_vec, &ah_neg_n, false); //!< ah_neg_n contains the capacity we charged additional (starting from the 'fully charged' state) on the negative electrode

    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::getOCV is fully discharging the anode.\n";

    CC_halfCell_full(Idis, dt, false, OCVni_vec, &ah, true); //!< now OCVni contains the anode OCV for the full lithium concentration range (from maximum concentration to 0)
  } catch (int e) {
    //!< std::cout << "Throw test: " << 43 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in a subfunction of CyclerOld::getOCV when half-cycling the cells" << e << ". Throwing it on.\n";
    throw e;
  }

  //!<*********************************************************** 3 make a common x-axis ***********************************************************************
  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::getOCV is making the x-axis\n";
  //!< the x-axis we want gives discharged charge [Ah]
  //!< We want the discharge charge to be 0 when the cell OCV is 4.2V (OCVp - OCVn = 4.2) [where 4.2 is the maximum cell voltage]
  //!< The points 'to the left', i.e. with a higher value for OCVp - OCVn must have a negative discharged charge
  //!< The points 'to the right', i.e. with a lower value for OCVp - OCVn must have a positive discharged charge

  //!< make the x-axis for the cathode
  //!<		the variable 'ah_neg_p' has the charge we charged additionally starting from the 'fully charged' position (i.e. how far we went 'to the left')
  //!<		so the 'leftmost' point of the cathode (with the highest OCVp) must have an x-value of ah_neg_p (which will be negative)
  //!<		then the point on the OCVp curve where the cell voltage is 4.2V will have an x-value of 0
  //!< From that leftmost point, the x-axis gives the discharged charge in Ah.
  //!<		the current during the discharge was constant (the variable Idis), so we just need to multiply the current by the time until now
  //!<		every time step takes 'dt' seconds, so the time [h] in step i is given by 'i*dt/3600'

  const auto np = OCVpi_vec.size(), nn = OCVni_vec.size(); //!< actual number of data points in the full discharge of each electrode

  slide::FixedData<double> Ahpi_vec_fixed(ah_neg_p, (dt / 3600.0 * Idis), np);

  //!< make the x-axis for the anode, exactly the same as for the cathode
  slide::FixedData<double> Ahni_vec_fixed(ah_neg_n, (dt / 3600.0 * Idis), nn);

  //!< Then we need to align both half-cell curves on a common x-axis
  //!<		the most negative point on the x-axis must be the leftmost point of both curves: xmin = min (Ahpi[0], Ahni[0])
  //!<		the most positive point on the x-axis must be the rightmost point of both curves: xmax = max (Ahpi[np-1], Ahni[nn-1])
  //!<			where np is the number of points on the cathode, i.e. Ahpi[np-1] is the 'rightmost' point of the cathode OCV. Similar for the anode
  const double Ahmin = std::min(Ahni_vec_fixed[0], Ahpi_vec_fixed[0]);         //!< the most negative point of the half-cell curves
  const double Ahmax = std::max(Ahni_vec_fixed.back(), Ahpi_vec_fixed.back()); //!< the most positive point of the half-cell curves

  //!< between these two extreme points, we make a data point every 'nstore' points.
  //!<		i.e. the second point on the x-axis is Ahmin + nstore*dt/3600*Idis
  //!<			 the thirds point on the x-axis is Ahmin + 2* nstore*dt/3600*Idis
  //!<			 etc.
  //!<		and we stop when the discharged charge is larger than Ahmax
  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::getOCV is aligning the half-cell curves\n";

  const double dAh_fixed = nStore * dt / 3600.0 * Idis;
  const int n_elem = static_cast<int>((Ahmax - Ahmin) / dAh_fixed) + 2; //!<+2 because while (Ah.back() <= Ahmax) so we pass Ahmax. stop if we've reached the most positive point of the OCV curve

  Ah_fixed = slide::FixedData<double>(Ahmin, dAh_fixed, n_elem);

  for (int i = 0; i < Ah_fixed.size(); ++i) {
    //!<*********************************************************** 4 align the half-cell curves on the common x-axis **********************************************************************
    //!< Now we have points of the electrodes each with their own x-axis
    //!< We can use linear interpolation to get the points on the common x-axis
    //!< electrode OCVs for this data point on the common x-axis
    constexpr bool is_Ahpi_fixed{ true }, is_Ahni_fixed{ true };                                                        //!< Always true probably. check_is_fixed(Ahpi, np);
    const double ocvp = linInt(printBool::printCrit, bound, Ahpi_vec_fixed, OCVpi_vec, np, Ah_fixed[i], is_Ahpi_fixed); //!< interpolate the negative half cell curve at the x-point we're at
    const double ocvn = linInt(printBool::printCrit, bound, Ahni_vec_fixed, OCVni_vec, nn, Ah_fixed[i], is_Ahni_fixed); //!< interpolate the positive half-cell curve at the x-point we're at

    OCVp.push_back(ocvp);
    OCVn.push_back(ocvn);
  }

  //!<*********************************************************** 5 output variables ***********************************************************************
  //!< restore the original battery state
  c.setStates(sini, Iini); //!< reset the states and current
  c.setTenv(Tenvi);        //!< reset the environmental temperature

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::getOCV terminating.\n";
}

double CyclerOld::checkUp_batteryStates(bool blockDegradation, bool checkCap, int cumCycle, double cumTime, double cumAh, double cumWh)
{
  /*
   * Function to do a capacity check and write the battery states as part of a check-up.
   * It will add one row of data in the csv file with the results (DegradationData_batteryState.csv in the subfolder of this CyclerOld)
   * the row has the following entries:
   * 		number of cycles until now
   * 		time the cell has been cycled until now [h]
   * 		cumulative Ah throughput up to now [Ah]
   * 		cumulative Wh throughput up to now [Wh]
   * 		the remaining capacity of the cell [Ah]
   * 		the transformed li concentration at the positive inner nodes of the positive particle (nch values)
   * 		the transformed li concentration at the positive inner nodes of the negative particle (nch values)
   * 		the cell temperature [K]
   * 		the thickness of the SEI layer [m]
   * 		the lost lithium [As]
   * 		the thickness of the cathode [m]
   * 		the thickness of the anode [m]
   * 		the volume fraction of active material in the cathode [-]
   * 		the volume fraction of active material in the anode [-]
   * 		the effective surface area of the cathode [m2 m-3]
   * 		the effective surface area of the anode [m2 m-3]
   * 		the surface area of the cracks at the surface of the negative particle [m2]
   * 		the diffusion constant at reference temperature of the cathode [m s-1]
   * 		the diffusion constant at reference temperature of the anode [m s-1]
   * 		the specific resistance of the combined electrodes [Ohm m2]
   * 		the thickness of the plated lithium layer [m]
   * 		the DC resistance of the cell [Ohm]
   * 		the active surface area of the anode [m2]
   *
   * IN
   * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * checkCap 			if true, the capacity is measured.
   * 						if false, the capacity is set to 0 and only the battery states are written
   * cumCycle				number of cycles up to now [-]
   * cumTime				time this cell has been cycled up to now [hour]
   * cumAh				cumulative Ah throughput up to now [Ah]
   * cumWh				cumulative Wh throughput up to now [Wh]
   *
   * OUT
   * cap 					the remaining capacity of the cell [Ah]
   *
   * THROWS
   * 1001 				the file in which to write the results couldn't be opened
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_batteryStates is starting.\n";

  //!< Get the capacity and the states
  const auto &states = c.getStates(); //!< array with the battery states
  double cap;                         //!< the remaining capacity [Ah]
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_batteryStates is getting the capacity.\n";
    if (checkCap)
      cap = getCapacity(blockDegradation);
    else
      cap = 0; //!< capacity is not measured
  } catch (int e) {
    //!< std::cout << "Throw test: " << 44 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkUp_batteryStates when measuring the cell capacity: " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< Write the results to the csv file
  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::checkUp_batteryStates is writing the capacity to the csv file.\n";

  const auto fol = PathVar::results / ID; //!< we want to write the file in a subfolder, so append the name of the subfolder before the name of the csv file
  std::ofstream output;

  const auto w_mode = !fileStatus.is_DegradationData_batteryState_created ? std::ios_base::out : std::ios_base::app; //!< Check if created earlier, if not then create, if created then append.
  output.open(fol / "DegradationData_batteryState.csv", w_mode);

  if (!output.is_open()) {
    //!< we couldn't open the file
    if constexpr (printBool::printCrit)
      std::cerr << "ERROR in CyclerOld::checkUp_batteryStates. File " << fol + "DegradationData_batteryState.csv"
                << " could not be opened. Throwing an error.\n";

    throw 1001;
  }

  fileStatus.is_DegradationData_batteryState_created = true; //!< Set it to true, because we just opened it.

  output << cumCycle << ',' << cumTime << ',' << cumAh << ',' << cumWh; //!< write the data points for where the cell is in it's life
  output << ',' << cap;                                                 //!< write the cell capacity

  for (const auto s : states)
    output << ',' << s; //!< write the state variables of the cell

  output << ',' << c.getR();            //!< write the total cell resistance (DC resistance [Ohm])
  output << ',' << c.getAnodeSurface(); //!< write the active anode surface area an*thickn*elec_surf excluding cracks [m2]
  output << '\n';                       //!< write an end-line (we have written everything we want)
  output.close();

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_batteryStates terminating.\n";

  return cap;
}

void CyclerOld::checkUp_OCVcurves(bool blockDegradation, double ocvpini, double ocvnini)
{
  /*
   * Function to measure and write the half-cell OCV curves as part of a check-up.
   * This function will add 4 rows to the csv file with the results (DegradationData_OCV.csv in the subfolder of this CyclerOld)
   * 		the first row contains the 'x-axis' of the half cell curves with the discharged charge [Ah]. It is 0 when the OCV of the cell is the maximum voltage, negative for points where the cell OCV is larger and positive for points where the cell OCV is smaller
   *		the second row contains the cathode potential at the corresponding points on the x-axis
   *		the third row contains the anode potential at the corresponding point on the x-axis
   *		the fourth row is an empty row to indicate the previous 3 rows are one 'data set'
   * At the end of the third row, also the 'operating point' of the cell is written, i.e. the potentials of the electrodes when the check-up procedure was started
   *		this is especially interesting when simulating calendar ageing because it allows to check the anode potential at which the cell is actually resting
   *
   * IN
   * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * ocvpini				cathode potential at the operating point when the check-up is called [V]
   * ocvnini				anode potential at the operating point when the check-up is called [V]
   *
   * THROWS
   * 1001 				the file in which to write the results couldn't be opened
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_OCVcurves is starting.\n";

  //!< variables
  constexpr int Ninocv = 25.0 * 3600.0 / 2.0 / 750.0 * 5.0; //!< length of the arrays with the OCV curves (0.04C -> 25 hours, 2s time steps, store one in every 750 points, *5 to ensure they are long enough)
  std::vector<double> ocvp, ocvn;                           //!< arrays to contain the OCV curves
  slide::FixedData<double> ocvAh;

  ocvp.reserve(Ninocv); //!< Reserve for the estimated length. #CHECK
  ocvn.reserve(Ninocv);

  //!< record the OCV cycles
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_OCVcurves is getting the OCV curves.\n";
    getOCV(ocvAh, ocvp, ocvn); //!< measure the OCV curves
  } catch (int e) {
    //!< std::cout << "Throw test: " << 45 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkUp_OCVcurves when measuring the OCV curves: " << e << ". Throwing the error on.\n";
    throw e;
  }

  //!< Write the curves to the csv file
  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::checkUp_OCVcurves is writing the OCV curves to a csv file.\n";

  const auto fol = PathVar::results / ID;                                                                   //!< we want to write the file in a subfolder, so append the name of the subfolder before the name of the csv file
  const auto w_mode = !fileStatus.is_DegradationData_OCV_created ? std::ios_base::out : std::ios_base::app; //!< Check if created earlier, if not then create, if created then append.

  std::ofstream output(fol + "DegradationData_OCV.csv", w_mode);
  if (!output.is_open()) {
    //!< we couldn't open the file
    if constexpr (printBool::printCrit)
      std::cerr << "ERROR on CyclerOld::checkUp_OCVcurves. File " << fol + "DegradationData_OCV.csv"
                << " could not be opened. Throwing an error.\n";
    throw 1001;
  }

  fileStatus.is_DegradationData_OCV_created = true; //!< Set it to true, because we just opened it.

  for (int i = 0; i < ocvAh.size(); i++) //!< write the x-axis, discharged Ah
    output << ocvAh[i] << ',';
  output << '\n';

  for (const auto ocvp_i : ocvp) //!< write the cathode OCV
    output << ocvp_i << ',';
  output << '\n';

  for (const auto ocvn_i : ocvn) //!< write the anode OCV
    output << ocvn_i << ',';

  //!< write the voltages of the operating point of the cell next to the anode OCV
  output << ",," << ocvpini << ',' << ocvnini << ',' << '\n';

  //!< leave one row blank to indicate these three lines were one measurement
  output << '\n';
  output.close();

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_OCVcurves terminating\n";
}

void CyclerOld::checkUp_CCCV(bool blockDegradation, int nCycles, double Crates[], double Ccut_cha, double Ccut_dis, bool includeCycleData)
{
  /*
   * Function to simulate and write some CCCV cycles as part of a check-up
   * The voltage and temperature while cycling are written in a csv file in the subfolder of this CyclerOld.
   * The data from each check-up is written in a new csv file, with as name DegradationData_CheckupCycle_x.csv
   * 		where 'x' is 1 (for the first check-up), 2 (for the second check-up), etc.
   * In the file, a data entry is recorded for every 2 seconds. Each data entry is written on one row of the csv file.
   * Each row has 15 columns (with the values for that entry):
   * 	 	total_time 		the total (cumulative) time since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	Ah_throughput 	the total (cumulative) charge throughput since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	Wh_throughput 	the total (cumulative) energy throughput since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	I 				the battery current at this point in time [A], > 0 for discharge, < 0 for charge
   * 	 	V  				the battery voltage at this point in time [V]
   * 	 	OCVp 			the cathode potential at this point in time [V]
   * 	 	OCVn  			the anode potential at this point in time [V]
   * 	 	Temperature 	the battery temperature at this point in time [K]
   * 	 	charge_time  	the cumulative time spent on charging since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	charge_Ah  		the cumulative charged charge since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	charge_Wh  		the cumulative charged energy since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	discharge_time 	the cumulative time spent on discharging since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	discharge_Ah 	the cumulative discharged charge since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	discharge_Wh 	the cumulative discharged energy since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	rest_time 		the cumulative time spent on resting since the start of this function, i.e. it is 0 for the first row [s]
   *
   * IN
   * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * nCycles				number of different cycles to be simulated (i.e. the length of the array crates)
   * Crates				array with the Crates of the CC phase of the different cycles to be checked, must be all positive
   * Ccut_cha				C rate of the cutoff current for the CV phase in the charges, must be positive
   * Ccut_dis				C rate of the cutoff current for the CV phase in the discharges, must be positive
   * includeCycleData	 	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
   * 							if false, the cycling data from the CCCV curves is only stored in a separate data file (and not part of the 'regular cycling data files': 'cyclingData_x.csv')
   * 							if true, the cycling data from the CCCV curves is stored both in a separate document and in the regular cycling data files ('cyclingData_x.csv')
   *
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_CCCV starting\n";

  //!< variables
  std::string nameCCCV = "DegradationData_CheckupCycle_" + std::to_string(indexdegr) + ".csv"; //!< name of the csv file in which the cycling data will be written
  double ahi, whi, timei;                                                                      //!< unneeded feedback variables (charge, energy and time spent in underlying functions)
  double dt = 2;                                                                               //!< use time steps of 2 seconds
  int feedb_old = CyclingDataTimeInterval;                                                     //!< the original data collection time interval
  int dataTimeInterval = 2;                                                                    //!< time resolution at which we want to store the cycling data from the CCCV curves
  CyclingDataTimeInterval = dataTimeInterval;                                                  //!< update the time resolution at which cycling data is stored for the CCCV cycles

  //!< fully charge the cell
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_CCCV is charging the cell before doing the CCCV cycles.\n";
    CC_V_CV_I(1.0, c.getVmax(), Ccut_cha, dt, blockDegradation, &ahi, &whi, &timei); //!< 1C charge to the specified cutoff Crate
  } catch (int e) {
    //!< std::cout << "Throw test: " << 46 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkup_CCCV when initially charging the cell: " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< std::flush the cycling data of the cell which was still stored in the arrays so we can start with an empty data buffer to store the results of the check-up
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_CCCV is flushing the previously stored cycling data.\n";
    writeCyclingData();
  } catch (int e) {
    //!< std::cout << "Throw test: " << 47 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkup_CCCV when writing the original cycling data: " << e << ".\n";
    throw e;
  }

  //!< do the cycles
  for (int i = 0; i < nCycles; i++) {
    try {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::checkUp_CCCV is simulating a CCCV cycle at Crate = " << Crates[i] << ".\n";
      CC_V_CV_I(Crates[i], c.getVmin(), Ccut_dis, dt, blockDegradation, &ahi, &whi, &timei); //!< discharge to the minimum cell voltage (CC at the given C rate and CV to the given cutoff C rate)
      CC_V_CV_I(Crates[i], c.getVmax(), Ccut_cha, dt, blockDegradation, &ahi, &whi, &timei); //!< charge to the maximum cell voltage (CC at the given C rate and CV to the given cutoff C rate)
    } catch (int e) {
      //!< std::cout << "Throw test: " << 48 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in CyclerOld::checkup_CCCV when following CCCV cycle at Crate: " << Crates[i]
                  << ". error " << e << ". Throwing it on.\n";

      throw e;
    }
  }

  //!< write the cycling data from the CCCV cycles
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_CCCV is writing the cycling data from the CCCV cycles.\n";
    writeCyclingData(nameCCCV, !includeCycleData); //!< if we don't want to include the cycling data, clear the buffer
  } catch (int e) {
    //!< std::cout << "Throw test: " << 49 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkup_CCCV when writing the cycling data from the CCCV cycles: " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< restore the old data-collection setting
  CyclingDataTimeInterval = feedb_old;

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_CCCV terminating.\n";
}

void CyclerOld::checkUp_pulse(bool blockDegradation, const std::vector<double> &I, const std::vector<double> &T, int profileLength, bool includeCycleData)
{
  /*
   * Function to simulate and write a pulse test as part of a check-up.
   * The voltage and temperature while cycling are written in a csv file in the subfolder of this CyclerOld.
   * The data from each check-up is written in a new csv file, with as name DegradationData_CheckupPulse_x.csv
   * 		where 'x' is 1 (for the first check-up), 2 (for the second check-up), etc.
   * In the file, a data entry is recorded for every 2 seconds. Each data entry is written on one row of the csv file.
   * Each row has 15 columns (with the values for that entry):
   * 	 	total_time 		the total (cumulative) time since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	Ah_throughput 	the total (cumulative) charge throughput since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	Wh_throughput 	the total (cumulative) energy throughput since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	I 				the battery current at this point in time [A], > 0 for discharge, < 0 for charge
   * 	 	V  				the battery voltage at this point in time [V]
   * 	 	OCVp 			the cathode potential at this point in time [V]
   * 	 	OCVn  			the anode potential at this point in time [V]
   * 	 	Temperature 	the battery temperature at this point in time [K]
   * 	 	charge_time  	the cumulative time spent on charging since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	charge_Ah  		the cumulative charged charge since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	charge_Wh  		the cumulative charged energy since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	discharge_time 	the cumulative time spent on discharging since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	discharge_Ah 	the cumulative discharged charge since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	discharge_Wh 	the cumulative discharged energy since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	rest_time 		the cumulative time spent on resting since the start of this function, i.e. it is 0 for the first row [s]
   *
   * IN
   * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * profileName 			name of the csv file which contains the current profile for the pulse test
   * 							the first column contains the current in [A] (positive for discharge, negative for charge)
   * 							the second column contains the time in [sec] the current should be maintained
   * 							the profile must be a net discharge, i.e. sum (I*dt) > 0
   * profileLength		length of the current profiles (number of rows in the csv file)
   * includeCycleData	 	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
   * 							if false, the cycling data from the pulse discharge is only stored in a separate data file (and not part of the 'regular cycling data files': 'cyclingData_x.csv')
   * 							if true, the cycling data from the pulse discharge is stored both in a separate document and in the regular cycling data files ('cyclingData_x.csv')
   *
   * THROWS
   * 1013 				the pulse profile is not a net discharge
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_pulse starting.\n";

  //!< variables
  std::string namePulse = "DegradationData_CheckupPulse_" + std::to_string(indexdegr) + ".csv"; //!< name of the csv file in which the cycling data will be written
  double dt = 2;                                                                                //!< take time steps of 2 seconds
  int limit = 1;                                                                                //!< if a voltage limit is reached during the pulse discharge profile, do a CV at that voltage for the rest of the time in that step
  int feedb_old = CyclingDataTimeInterval;                                                      //!< the old data collection time interval
  int dataTimeInterval = 2;                                                                     //!< time resolution at which we want to store the cycling data from the pulse discharge
  CyclingDataTimeInterval = dataTimeInterval;                                                   //!< update the time resolution of the cycling data collection to the new value for the pulse discharge
  double ahi, whi, timei;                                                                       //!< unneeded feedback variables
  double Ccut = 0.05;                                                                           //!< C rate of the cutoff current for the CV phase

  //!< check that the pulse profile is a discharge
  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::checkUp_pulse is checking if the current profile is a discharge.\n";

  const double aht = std::inner_product(I.begin(), I.end(), T.begin(), 0.0); //!< total charge throughput of the profile [As]

  if (aht <= 0) { //!< if the total charge throughput is negative, the profile is a net charge
    if constexpr (printBool::printCrit)
      std::cerr << "ERROR in CyclerOld::checkUp_pulse: the pulse profile is not a net discharge. "
                   "The total charge throughput was "
                << aht / 3600.0 << "Ah, where a negative number indicates a net charge.\n";

    //!< The total charge throughput MUST BE positive, i.e. there must be more 'discharge' than 'charge' in the profile.
    //!< Else, the user has to change the function 'checkUp_pulse' to deal with a net charge
    throw 1013;
  }

  //!<*********************************************************** 2 charge the cell ***********************************************************************
  //!< fully charge the cell
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_pulse is charging the cell.\n";
    CC_V_CV_I(1.0, c.getVmax(), Ccut, dt, blockDegradation, &ahi, &whi, &timei); //!< fully charge the cell to the maximum voltage (1C CC, CV until 0.05C cut off current)
  } catch (int e) {
    //!< std::cout << "Throw test: " << 50 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkup_pulse when initially charging the cell: " << e << ".\n";
  }

  //!< std::flush the cycling data of the cell which was still stored in the arrays so we can start with an empty data buffer to store the results of the check-up
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_pulse is flushing out the previously stored cycling data.\n";
    writeCyclingData();
  } catch (int e) {
    //!< std::cout << "Throw test: " << 51 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkup_pulse when writing the original cycling data: " << e << "\n.\n";
    throw e;
  }

  //!<*********************************************************** 3 do the pulse discharge until the lower voltage is reached ***********************************************************************
  //!< The current pulses must be in a file with the name given by 'CheckPulseProfile'
  //!<	and length given by 'length_checkPulseProfile'
  //!<	the first column gives the current [A]
  //!<	the second column gives the time each current should be maintained [s]
  //!< the cell starts the pulse tests at the maximum cell voltage (SOC = 100%). Keep doing the pulses until the lower voltage is reached
  bool Vmin = false;               //!< has the lower voltage limit been reached?
  const double Vupp = c.getVmax(); //!< the upper voltage limit is the maximum voltage of the cell
  const double Vlow = c.getVmin(); //!< the lower voltage limit is the minimum voltage of the cell
  int Vlim;                        //!< return-integer of followI indicating which voltage limit was hit
  try {
    //!< a loop to repeatedly apply the pulse profile until the lower voltage is reached
    while (!Vmin) {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::checkUp_pulse is applying the pulse profile and has so far discharged " << ahi << "Ah.\n";
      Vlim = followI(profileLength, I, T, blockDegradation, limit, Vupp, Vlow, &ahi, &whi, &timei); //!< follow the current pulses, the return-integer indicates which voltage limit was hit while following the profile
      Vmin = (Vlim == -1) || (Vlim == 10);                                                          //!< a return-integer of -1 or 10 means the lower voltage was hit, stop when it is
    }
  } catch (int e) {
    //!< std::cout << "Throw test: " << 52 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkup_pulse when following the pulses: " << e << ".\n";
  }

  //!< write the cycling data from the pulse cycles
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_pulse is writing the cycling data from the pulse profile.\n";
    writeCyclingData(namePulse, !includeCycleData); //!< if we don't want to include the cycling data, clear the buffer
  } catch (int e) {
    //!< std::cout << "Throw test: " << 53 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkup_pulse when writing the cycling data from the pulse discharge: " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< restore the old data-collection setting
  CyclingDataTimeInterval = feedb_old;

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp_pulse terminating.\n";
}

void CyclerOld::checkUp_pulse(bool blockDegradation, const std::string &profileName, int profileLength, bool includeCycleData)
{
  /*
   * Function to simulate and write a pulse test as part of a check-up.
   * The voltage and temperature while cycling are written in a csv file in the subfolder of this CyclerOld.
   * The data from each check-up is written in a new csv file, with as name DegradationData_CheckupPulse_x.csv
   * 		where 'x' is 1 (for the first check-up), 2 (for the second check-up), etc.
   * In the file, a data entry is recorded for every 2 seconds. Each data entry is written on one row of the csv file.
   * Each row has 15 columns (with the values for that entry):
   * 	 	total_time 		the total (cumulative) time since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	Ah_throughput 	the total (cumulative) charge throughput since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	Wh_throughput 	the total (cumulative) energy throughput since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	I 				the battery current at this point in time [A], > 0 for discharge, < 0 for charge
   * 	 	V  				the battery voltage at this point in time [V]
   * 	 	OCVp 			the cathode potential at this point in time [V]
   * 	 	OCVn  			the anode potential at this point in time [V]
   * 	 	Temperature 	the battery temperature at this point in time [K]
   * 	 	charge_time  	the cumulative time spent on charging since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	charge_Ah  		the cumulative charged charge since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	charge_Wh  		the cumulative charged energy since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	discharge_time 	the cumulative time spent on discharging since the start of this function, i.e. it is 0 for the first row [s]
   * 	 	discharge_Ah 	the cumulative discharged charge since the start of this function, i.e. it is 0 for the first row [Ah]
   * 	 	discharge_Wh 	the cumulative discharged energy since the start of this function, i.e. it is 0 for the first row [Wh]
   * 	 	rest_time 		the cumulative time spent on resting since the start of this function, i.e. it is 0 for the first row [s]
   *
   * IN
   * blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 							and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 						if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 							and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * profileName 			name of the csv file which contains the current profile for the pulse test
   * 							the first column contains the current in [A] (positive for discharge, negative for charge)
   * 							the second column contains the time in [sec] the current should be maintained
   * 							the profile must be a net discharge, i.e. sum (I*dt) > 0
   * profileLength		length of the current profiles (number of rows in the csv file)
   * includeCycleData	 	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
   * 							if false, the cycling data from the pulse discharge is only stored in a separate data file (and not part of the 'regular cycling data files': 'cyclingData_x.csv')
   * 							if true, the cycling data from the pulse discharge is stored both in a separate document and in the regular cycling data files ('cyclingData_x.csv')
   *
   * THROWS
   * 1013 				the pulse profile is not a net discharge
   */

  //!<*********************************************************** 1 get the pulse profile ***********************************************************************
  //!< read the pulse profile
  static thread_local std::vector<double> I(profileLength), T(profileLength); //!< array for the profile and array for the duration of each step as integers % integer did not make sense so I converted it to double.
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp_pulse is reading the current profile.\n";
    slide::loadCSV_2col(PathVar::data / profileName, I, T); //!< read the file
  } catch (int e) {
    //!< std::cout << "Throw test: " << 54 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkUp_pulse when reading the pulse profile: " << e << ". Throwing it on.\n";
    throw e;
  }

  checkUp_pulse(blockDegradation, I, T, profileLength, includeCycleData);
}

double CyclerOld::checkUp(struct checkUpProcedure &proc, int cumCycle, double cumTime, double cumAh, double cumWh)
{
  /*
   * Function to do a check-up during a degradation experiment
   * a check-up can consists of:
   * 		a capacity measurement
   * 		a half-cell OCV measurement
   * 		some CC CV cycles (CC discharge, CV discharge, CC charge, CV charge) at 0.5, 1 and 2 C
   * 		a pulse discharge test
   *
   * The results are written in csv files, in the subfolder of this CyclerOld.
   * 		DegradationData_batteryState.csv				one new row of data appended at the end of the existing csv file
   * 		DegradationData_OCV.csv							4 new rows of data appended at the end of the existing csv file
   * 		DegradationData_CheckupCycle_x.csv 				a new file for the data of this check-up, x is the index of the check-up (1 for the first check-up, 2 for the second, etc.)
   * 		DegradationData_CheckupPulse_x.csv				a new file for the data of this check-up, x is the index of the check-up (1 for the first check-up, 2 for the second, etc.)
   * See the individual functions (checkUp_yyy) for an exact description of what is in each file.
   *
   * IN
   * proc 				structure with the parameters of the check-up procedure with the following fields:
   * 		blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 								and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 								if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 								and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * 		capCheck 			boolean indicating if the capacity should be checked
   * 		OCVCheck			boolean indicating if the half-cell OCV curves should be checked
   * 		CCCVCheck			boolean indicating if some CCCV cycles should be done as part of the check-up procedure
   * 		pulseCheck			boolean indicating if a pulse discharge test should be done as part of the check-up procedure
   * 		includeCycleData	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
   * 		nCycles				number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
   * 		Crates				array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
   * 		Ccut_cha			C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
   * 		Ccut_dis			C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
   * 		profileName 		name of the csv file which contains the current profile for the pulse test
   * 								the first column contains the current in [A] (positive for discharge, negative for charge)
   * 								the second column contains the time in [sec] the current should be maintained
   * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
   * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
   * cumCycle				number of cycles up to now [-]
   * cumTime				time this cell has been cycled up to now [hour]
   * cumAh				cumulative Ah throughput up to now [Ah]
   * cumWh				cumulative Wh throughput up to now [Wh]
   *
   * OUT
   * cap 					the remaining capacity of the cell [Ah]
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp starting\n";

  //!< variables
  const auto sini = c.getStates();            //!< initial cell state
  const auto Iini = c.getI();                 //!< initial cell current [A]
  double v;                                   //!< initial cell voltage
  double tcellini, Tenvini, Trefi;            //!< initial cell, environmental and reference temperature
  double dt = 2;                              //!< time step when cycling and resting [s]
  double cap = 0;                             //!< capacity of the cell [Ah]
  double ocvpini, ocvnini, etap, etan, rdrop; //!< unneeded feedback variables
  double ahi, whi, timei;                     //!< unneeded feedback variables

  //!< Store the initial OCV and temperature from before the check-up
  c.getVoltage(printBool::printCrit, &v, &ocvpini, &ocvnini, &etap, &etan, &rdrop, &tcellini); //!< initial cell voltage and temperature
  c.getTemperatures(&Tenvini, &Trefi);                                                         //!< initial environmental and reference temperature
  int feedbackini = CyclingDataTimeInterval;                                                   //!< initial data collection time interval

  //!< if we don't want to include the cycling data from the check-up in the cycling data from the cell
  //!< push the cycling data and set the data collection to 0 to stop recording
  if (!proc.includeCycleData) {
    try {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::checkUp is flushing the previously stored cycling data.\n";
      writeCyclingData();          //!< write the cycling data which was still stored
      CyclingDataTimeInterval = 0; //!< don't record cycling data from the check-up
    } catch (int e) {
      //!< std::cout << "Throw test: " << 55 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in CyclerOld::checkUp when writing the previously stored cycling data: " << e << ".\n";
      throw e;
    }
  }

  //!<*********************************************************** 1 precharge the cell & set the temperature ***********************************************************************

  //!< Get the cell to 0.2V below the maximum voltage.
  //!< This is needed because later we are changing the temperature, which might affect the voltage a little bit.
  //!< If the cell is fully charged and the voltage changes a little bit, you can get an illegal voltage (e.g. Vmax + 0.001V), which results in an error
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::checkUp is bringing the cell to just below its maximum voltage.\n";
    CC_V_CV_I(1.0, c.getVmax() - 0.2, 0.05, dt, proc.blockDegradation, &ahi, &whi, &timei); //!< 1C current, 0.05C cutoff current
  } catch (int e) {
    //!< std::cout << "Throw test: " << 56 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::checkUp when setting the cell voltage to Vmax-0.2V: " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< do the check-up at reference temperature
  c.setTenv(Trefi); //!< set the environmental temperature to the reference temperature
  c.setT(Trefi);    //!< set the cell temperature to the reference temperature

  //!<*********************************************************** 2 do the different check-up tests ***********************************************************************
  //!< Measure the capacity
  if (proc.capCheck) {
    try {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::checkUp is starting a capacity check.\n";
      cap = checkUp_batteryStates(proc.blockDegradation, true, cumCycle, cumTime, cumAh, cumWh);
    } catch (int e) {
      //!< std::cout << "Throw test: " << 57 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in CyclerOld::checkUp when measuring the cell capacity: " << e << ". Skip the capacity measurement.\n";
    }
  }

  //!< measure the OCV curves
  if (proc.OCVCheck) {
    try {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::checkUp is starting an OCV check.\n";
      checkUp_OCVcurves(proc.blockDegradation, ocvpini, ocvnini);
    } catch (int e) {
      //!< std::cout << "Throw test: " << 58 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in CyclerOld::checkUp when measuring the OCV curves: " << e << ". Skip the OCV measurements.\n";
      throw e;
    }
  }

  //!< measure the CCCV cycles
  if (proc.CCCVCheck) {
    try {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::checkUp is starting a CCCV cycle check.\n";
      checkUp_CCCV(proc.blockDegradation, proc.nCycles, proc.Crates, proc.Ccut_cha, proc.Ccut_dis, proc.includeCycleData);
    } catch (int e) {
      //!< std::cout << "Throw test: " << 58 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in CyclerOld::checkUp when measuring the CCCV cycles: " << e << ". Skip the CCCV measurement.\n";
    }
  }

  //!< measure the discharge pulses
  if (proc.pulseCheck) {
    try {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::checkUp is starting a pulse profile check.\n";
      checkUp_pulse(proc.blockDegradation, proc.I, proc.T, proc.profileLength, proc.includeCycleData);
    } catch (int e) {
      //!< std::cout << "Throw test: " << 59 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in CyclerOld::checkUp when measuring the pulse test: " << e << ". Skip the CCCV measurement.\n";
    }
  }

  //!< increase the counter of the number of check-ups we have done
  indexdegr++;

  //!<*********************************************************** 3 reset the original battery state ***********************************************************************

  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::checkUp is restoring the initial battery state.\n";

  c.setTenv(Tenvini);                    //!< restore the environmental temperature
  CyclingDataTimeInterval = feedbackini; //!< restore the initial data-collection setting
  if (proc.blockDegradation)
    c.setStates(sini, Iini); //!< restore the exact initial battery state
  else {
    double Trest = 3600;                                                                  //!< rest time is 1 hour
    double Ccut = 0.05;                                                                   //!< C rate of cutoff current from the CV phase
    CC_V_CV_I(1, ocvpini - ocvnini, Ccut, dt, proc.blockDegradation, &ahi, &whi, &timei); //!< bring the cell back to the original OCV (OCV_cell = OCV_p - OCV_n)
    CC_t(0.0, dt, proc.blockDegradation, Trest, &ahi, &whi, &timei);                      //!< rest such that the cell temperatures goes back to the environmental temperature
  }

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::checkUp terminating with capacity " << cap << "Ah.\n";

  //!< return the remaining capacity
  return cap;
}

void CyclerOld::cycleAgeing(double dt, double Vma, double Vmi, double Ccha, bool CVcha, double Ccutcha,
                            double Cdis, bool CVdis, double Ccutdis, double Ti, int nrCycles, int nrCap, struct checkUpProcedure &proc)
{
  /*
   * Function to simulate a cycle degradation experiment where a cell is continuously cycled at constant current and/or constant voltage.
   * The parameters of the cycling regime are set by the inputs.
   * After a set number of cycles, a check-up is done where the cell capacity, OCV curves, etc. are measured.
   * These results are written to csv files.
   *
   * IN
   * dt 		the time step to be used in the cycling, small enough to ensure stability and accuracy (1-5 seconds) [s]
   * Vma 		maximum voltage of this cycle, must be below the maximum voltage of the cell [V]
   * Vmi 		minimum voltage of this cycle, must be above the minimum voltage of the cell [V]
   * Ccha 	C-rate at which the battery should be charged during the CC phase[-]
   * CVcha 	boolean indicating if a CV charge should be done after the CC charge (true) or not (false)
   * 			if true, charging is CC CV
   * 			if false, charging is CC only
   * Ccutcha 	C rate of the cutoff current for the CV charge [-], > 0
   * Cdis 	C rate at which the battery should be discharged [-]
   * CVdis 	boolean indicating if a CV discharge should be done after the CC discharge (true) or not (false)
   * 			if true, discharging is CC CV
   * 			if false, discharging is CC only
   * Ccutdis 	C rate of the cutoff current for the CV discharge [-], > 0
   * Ti 		environmental temperature at which the test should be performed [K]
   * nrCycles	the number of cycles that has to be simulated [-]
   * nrCap 	the number of cycles between consecutive check-ups [-]
   * proc 	structure with the parameters of the check-up procedure.
   * 			The struct is defined in CyclerOld.hpp and has the following fields:
   * 		blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 								and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 								if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 								and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * 		capCheck 			boolean indicating if the capacity should be checked
   * 		OCVCheck			boolean indicating if the half-cell OCV curves should be checked
   * 		CCCVCheck			boolean indicating if some CCCV cycles should be done as part of the check-up procedure
   * 		pulseCheck			boolean indicating if a pulse discharge test should be done as part of the check-up procedure
   * 		includeCycleData	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
   * 		nCycles				number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
   * 		Crates				array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
   * 		Ccut				C rate of the cutoff current for the CV phase for the CCCV check-up, must be positive
   * 		profileName 		name of the csv file which contains the current profile for the pulse test
   * 								the first column contains the current in [A] (positive for discharge, negative for charge)
   * 								the second column contains the time in [sec] the current should be maintained for
   * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
   * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
   *
   *
   * throws
   * 1014		the input parameters describing the cycling regime are invalid
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::cycleAgeing starting.\n";

  slide::util::error::checkInputParam_CycAge(c, Vma, Vmi, Ccha, Ccutcha, Cdis, Ccutdis, Ti, nrCycles, nrCap); //!< Check the input parameters

  //!<*********************************************************** 1 variables & settings ***********************************************************************

  //!< variables
  double ahi;                        //!< discharged charge in this charge or discharge [Ah]
  double whi;                        //!< discharged energy in this charge or discharge [Wh]
  double ti;                         //!< time spent in this charge or discharge [sec]
  double Ahtot = 0;                  //!< charge throughput from the cycling regime until now [Ah]
  double Whtot = 0;                  //!< energy throughput from the cycling regime until now [Wh]
  double timetot = 0;                //!< time spent while following the cycling regime until now [hour]
  double cap;                        //!< capacity of the cell at this point in time [Ah]
  bool blockDegradation = false;     //!< account for degradation while we cycle
  bool final = true;                 //!< boolean to indicate if a check-up at the end of the cycling regime is needed
  double capnom = c.getNominalCap(); //!< nominal cell capacity [Ah] to convert Crate to Amperes

  //!<*********************************************************** 2 cell initialisation ***********************************************************************

  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::cycleAgeing is initialising the cell\n";

  //!< set the temperature
  c.setT(Ti);    //!< set the cell temperature
  c.setTenv(Ti); //!< set the environmental temperature

  //!< Get the battery to Vma so the cycling can start with a discharge
  try {
    if (CVcha)
      CC_V_CV_I(Ccha, Vma, Ccutcha, dt, blockDegradation, &ahi, &whi, &ti); //!< CC and CV charge
    else
      CC_V(-Ccha * capnom, dt, blockDegradation, Vma, &ahi, &whi, &ti); //!< CC charge (must have current as input, not C rate)
  } catch (int e) {
    //!< std::cout << "Throw test: " << 60 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in a subfunction of CyclerOld::cycleAgeing when getting the cell to Vma, error " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< do an initial check up
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::cycleAgeing is doing an initial check-up.\n";
    cap = checkUp(proc, 0, timetot, Ahtot, Whtot);
  } catch (int e) {
    //!< std::cout << "Throw test: " << 61 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in the initial check-up CyclerOld::cycleAgeing " << e << ". Throwing it on.\n";
    throw e;
  }

  //!<*********************************************************** 3 cycle age the cell ***********************************************************************

  for (int i = 0; i < nrCycles; i++) { //!< loop for the specified number of cycles
    try {

      //!< discharge
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::cycleAgeing is discharging the cell in cycle number " << i << ".\n";
      if (CVdis)
        CC_V_CV_I(Cdis, Vmi, Ccutdis, dt, blockDegradation, &ahi, &whi, &ti); //!< CC and CV discharge
      else
        CC_V(Cdis * capnom, dt, blockDegradation, Vmi, &ahi, &whi, &ti); //!< CC discharge (must have current as input, not C rate)
      Ahtot += std::abs(ahi);                                            //!< increase the charge throughput with the throughput of this discharge
      Whtot += std::abs(whi);                                            //!< increase the energy throughput with the throughput of this discharge
      timetot += (ti / 3600);                                            //!< increase the total time the time of this discharge

      //!< charge
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::cycleAgeing is charging the cell in cycle number " << i << ".\n";
      if (CVcha)
        CC_V_CV_I(Ccha, Vma, Ccutcha, dt, blockDegradation, &ahi, &whi, &ti); //!< CC and CV charge
      else
        CC_V(-Ccha * capnom, dt, blockDegradation, Vma, &ahi, &whi, &ti); //!< CC charge (must have current as input, not C rate)
      Ahtot += std::abs(ahi);                                             //!< increase the charge throughput with the throughput of this charge
      Whtot += std::abs(whi);                                             //!< increase the energy throughput with the throughput of this charge
      timetot += (ti / 3600);                                             //!< increase the total time the time of this charge

      //!< do a check-up every nrCap cycles
      //!<	i is the cycle number, so when it is a multiple of nrCap we need to do a check-up
      //!<  do i+1 to avoid doing a check-up in the first cycle
      if (std::fmod(i + 1, nrCap) == 0) { //!< remainder is 0 -> it is a multiple of nrCap
        if constexpr (printBool::printCyclerHighLevel)
          std::cout << "CyclerOld::cycleAgeing is doing a check-up in cycle number " << i << ".\n";
        cap = checkUp(proc, i + 1, timetot, Ahtot, Whtot); //!< do the check-up procedure

        //!< End the experiment if the cell capacity has decreased too much
        if (cap < capnom / 2.0) { //!< end simulation if the cell has only 50% of capacity left
          std::cout << "CyclerOld::cycleAgeing has finished cycling regime " << ID << " early because the cell has already lost 50% of its capacity.";
          std::cout << " We have done " << i << " cycles instead of " << nrCycles << " and the remaining capacity now is " << cap << " [Ah].\n";
          final = false; //!< skip the final check-up because we just did one
          break;         //!< stop cycling
        }
      }
    } //!< end try block

    //!< Catch an error which occurred while cycling the cell (or during the check-up procedure)
    catch (int e) {
      //!< std::cout << "Throw test: " << 62 << '\n';
      //!< print an error message (this is a critical error)
      if constexpr (printBool::printCrit) {
        std::cout << "Error in CyclerOld::cycleAgeing while cycling the cell according to  cycling regime " << ID << ". Error encountered is " << e << ". Stop cycling now.";
        std::cout << " We have done " << i << " cycles instead of " << nrCycles << " and the capacity last measured is " << cap << " [Ah].\n";
      }

      //!< we probably cannot do a full check-up procedure because the cell is in an illegal state.
      //!< Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
      checkUp_batteryStates(proc.blockDegradation, false, i + 1, timetot, Ahtot, Whtot);

      final = false; //!< skip the final check-up
      break;         //!< stop cycling
    }                //!< end try-catch block
  }                  //!< end loop to cycle the cell

  //!<*********************************************************** 4 final check-up ***********************************************************************

  //!< Don't do the check-up if the simulation ended due to an error (or shortage in capacity)
  if (final) {
    try {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::cycleAgeing is doing a final check-up.\n";
      checkUp(proc, nrCycles, timetot, Ahtot, Whtot);
    } catch (int e) {
      //!< std::cout << "Throw test: " << 63 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in a in the final check-up of CyclerOld::cycleAgeing " << e << ". Throwing it on.\n";
      throw e;
    }
  }

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::cycleAgeing terminating\n";
}

void CyclerOld::calendarAgeing(double dt, double V, double Ti, int Time, int timeCheck, int mode, struct checkUpProcedure &proc)
{
  /*
   * Function to simulate a calendar degradation experiment where a cell is rested at a given voltage.
   * The parameters of the calendar regime are set by the inputs.
   * After a given amount of time, a check-up is done where the cell capacity, OCV curves, etc. are measured.
   * These results are written to csv files.
   *
   * IN
   * dt 		the time step to be used [sec]
   * V 		the voltage at which the battery has to rest [V]
   * Ti 		the temperature at which the battery has to rest [K]
   * Time 	the time for which the battery has to rest [days]
   * 			must be a multiple of timeCheck
   * timeCheck the time after which a check-up has to be done [days]
   * mode 	integer deciding how often to recharge to the specified voltage.
   * 			When a cell rests, the voltage will slip a bit due to degradation, so the question is how often to recharge.
   * 			0 	recharge only after a check-up
   * 			1 	recharge every day
   * 			2 	float the cell at the specified voltage, instead of letting it rest
   * 				this option is not recommended, it takes ages to calculate
   * 			Any other value will produce an error
   * proc 	structure with the parameters of the check-up procedure with the following fields:
   * 		blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 								and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 								if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 								and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * 		capCheck 			boolean indicating if the capacity should be checked
   * 		OCVCheck			boolean indicating if the half-cell OCV curves should be checked
   * 		CCCVCheck			boolean indicating if some CCCV cycles should be done as part of the check-up procedure
   * 		pulseCheck			boolean indicating if a pulse discharge test should be done as part of the check-up procedure
   * 		includeCycleData	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
   * 		nCycles				number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
   * 		Crates				array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
   * 		Ccut				C rate of the cutoff current for the CV phase for the CCCV check-up, must be positive
   * 		profileName 		name of the csv file which contains the current profile for the pulse test
   * 								the first column contains the current in [A] (positive for discharge, negative for charge)
   * 								the second column contains the time in [sec] the current should be maintained
   * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
   * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
   *
   * throws
   * 1014		the input parameters describing the calendar regime are invalid
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::calendarAgeing starting\n";

  slide::util::error::checkInputParam_CalAge(c, V, Ti, Time, timeCheck, mode); //!< Check the input parameters

  //!<*********************************************************** 1 variables & settings ***********************************************************************

  //!< Variables
  double ahi, ahi2;                  //!< discharged charge in this charge or discharge [Ah]
  double whi, whi2;                  //!< discharged energy in this charge or discharge [Wh]
  double ti, ti2;                    //!< time spent in this charge or discharge [sec]
  double Ahtot = 0;                  //!< cumulative charge throughput until now [Ah]
  double Whtot = 0;                  //!< cumulative energy throughput until now [Wh]
  double timetot = 0;                //!< cumulative time until now [hour]
  double cap;                        //!< cell capacity [Ah]
  double capnom = c.getNominalCap(); //!< nominal cell capacity
  double Ccut = 0.005;               //!< Crate of the cutoff current for CV phases [-]
  bool blockDegradation = false;     //!< do account for degradation during calendar

  //!< time at which a check-up has to be done
  int nrdt = Time / timeCheck;                //!< number of check-ups to be done
  double trest = timeCheck * (24.0 * 3600.0); //!< resting time between consecutive check-ups in seconds

  //!<*********************************************************** 2 cell initialisation ***********************************************************************

  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::calendarAgeing is initialising the cell\n";

  //!< set temperature
  c.setT(Ti);    //!< set the Cell temperature
  c.setTenv(Ti); //!< set the environmental temperature

  //!< Get the battery to the resting voltage
  try {
    CC_V_CV_I(1.0, V, Ccut, dt, blockDegradation, &ahi, &whi, &ti); //!< do a 1C CC charge until the voltage, followed by a CV at this voltage
  } catch (int e) {
    //!< std::cout << "Throw test: " << 64 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in a subfunction of CyclerOld::CalendarAgeing when getting the cell to the resting voltage "
                << e << ". Throwing it on.\n";
    throw e;
  }

  //!< initial check-up
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::calendarAgeing is doing an initial check-up.\n";
    cap = checkUp(proc, 0, timetot, Ahtot, Whtot);
  } catch (int e) {
    //!< std::cout << "Throw test: " << 65 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in a subfunction of CyclerOld::CalendarAgeing in the initial check up " << e << ". Throwing it on.\n";
    throw e;
  }

  //!<*********************************************************** 3 calendar age the cell ***********************************************************************

  //!< loop for each resting period between two check-ups
  for (int i = 0; i < nrdt; i++) {
    try {
      //!< Calendar ageing type depending on 'mode'
      //!< rest the cell for the full period between consecutive check-ups and recharge at the end
      if (mode == 0) {
        if constexpr (printBool::printCyclerHighLevel)
          std::cout << "CyclerOld::calendarAgeing is resting the cell in period " << i << ".\n";
        CC_t(0, dt, blockDegradation, trest, &ahi, &whi, &ti); //!< rest (= CC charge at 0A)

        if constexpr (printBool::printCyclerHighLevel)
          std::cout << "CyclerOld::calendarAgeing is recharging the cell in period " << i << ".\n";
        CV_I(V, dt, blockDegradation, Ccut, &ahi2, &whi2, &ti2); //!< recharge to the specified voltage
        timetot += (ti + ti2) / 3600.0;                          //!< number of hours we have rested additionally
      }
      //!< recharge every day
      else if (mode == 1) {
        for (int j = 0; j < timeCheck; j++) //!< #CHECK
        {                                   //!< a loop for every day
          if constexpr (printBool::printCyclerHighLevel)
            std::cout << "CyclerOld::calendarAgeing is resting the cell in day " << j << " of period " << i << ".\n";
          CC_t(0, dt, blockDegradation, 24.0 * 3600.0, &ahi, &whi, &ti); //!< rest for one day
          if constexpr (printBool::printCyclerHighLevel)
            std::cout << "CyclerOld::calendarAgeing is recharging the cell in day " << j << " of period " << i << ".\n";
          CV_I(V, dt, blockDegradation, Ccut, &ahi2, &whi2, &ti2); //!< recharge to the specified voltage
          timetot += (ti + ti2) / 3600.0;                          //!< number of hours we have rested additionally
        }
      }
      //!< float at the set voltage (takes very long to calculate)
      else if (mode == 2) {
        if constexpr (printBool::printCyclerHighLevel)
          std::cout << "CyclerOld::calendarAgeing is floating the cell in period " << i << ".\n";
        CV_t(V, dt, blockDegradation, trest, &ahi, &whi, &ti); //!< do a CV (dis)charge for the entire trest period. Takes ages to compute
        timetot += (ti) / 3600.0;                              //!< number of hours we have rested additionally
      } else
        assert(false); //!< not allowed, we have checked at the start that this can't happen

      //!< do a check-up
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::calendarAgeing is doing a check-up in period " << i << ".\n";
      cap = checkUp(proc, 0, timetot, Ahtot, Whtot);

      //!< End the experiment if the cell capacity has decreased too much
      if (cap < capnom / 2.0) { //!< end simulation if the cell has only 50% of capacity left
        std::cout << "CyclerOld::CalendarAgeing has finished calendar regime " << ID << " early do to too little remaining capacity.";
        std::cout << " We have rested " << i * timeCheck << " days instead of " << Time << " and the remaining capacity now is "
                  << cap << " [Ah].\n";
        break;
      }
    } //!< end try block
    catch (int e) {
      //!< std::cout << "Throw test: " << 66 << '\n';
      std::cout << "Error in CyclerOld::CalendarAgeing while resting the cell according to calendar regime " << ID << ". Error encountered is " << e << ". Stop resting now.";
      std::cout << " We have rested " << i * timeCheck << " days instead of " << Time << " and the capacity last measured is " << cap << ".\n";
      //!< we probably cannot do a full check-up because the cell is in an illegal state.
      //!< Therefore, only write the BatteryStates
      checkUp_batteryStates(proc.blockDegradation, false, 0, timetot, Ahtot, Whtot);
      break;
    } //!< end try-catch block

  } //!< end loop to rest and check-up

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::calendarAgeing terminating.\n";
}

void CyclerOld::profileAgeing(const std::string &nameI, int limit, double Vma, double Vmi, double Ti, int nrProfiles, int nrCap, struct checkUpProcedure &proc, size_t length)
{
  /*
   * Function to simulate degradation by continuously cycling a cell with a certain current profile.
   * the profile is repeated until the cell is fully charged/discharged, then the cell is discharged/charged at 1C CC CV, and the profile is repeated again.
   * The parameters of the cycling regime are set by the inputs.
   * After a set number of profile repetitions, a check-up is done where the cell capacity, OCV curves, etc. are measured.
   * These results are written to csv files.
   *
   * IN
   * nameI 	name of the CSV-file with the current profile
   * 				the first column contains the current in [A], positive for discharge, negative for charge
   * 				the second column contains the time in [sec] the current should be maintained
   * limit 	integer describing what to do if the current can't be maintained because a voltage limit is reached
   * 				0 	immediately go to the next current step of the profile
   * 				1 	keep the voltage constant for the rest of this step of the profile
   * Vma 		maximum voltage of this degradation experiment, must be below the maximum voltage of the cell [V]
   * Vmi 		minimum voltage of this degradation experiment, must be above the minimum voltage of the cell [V]
   * Ti 		temperature at which the test should be performed [K]
   * nrProfiles the number of times the the current profile should be repeated
   * nrCap 	the number of times the current profile is repeated approximately between consecutive check-ups
   * 				a check-up is always done after a voltage limit was hit and we have re(dis)charged the cell
   * 				e.g. you want to do a check-up every 20 cycles. Suppose you can do 15 cycles in the given voltage window (before you hit a voltage limit and you have to re(dis)charge the cell to the other voltage limit)
   * 				then the check-up will be done after 30 cycles.
   * proc 	structure with the parameters of the check-up procedure with the following fields:
   * 		blockDegradation 	if true [RECOMMENDED], degradation is not accounted for during this check-up,
   * 								and at the end of the check-up, the exact battery state from the start of the check-up is restored
   * 								if false [NOT RECOMMENDED], degradation is accounted for during this check-up,
   * 								and at the end of the check-up, the cell's voltage and temperature are brought back to the levels from before the check-up
   * 		capCheck 			boolean indicating if the capacity should be checked
   * 		OCVCheck			boolean indicating if the half-cell OCV curves should be checked
   * 		CCCVCheck			boolean indicating if some CCCV cycles should be done as part of the check-up procedure
   * 		pulseCheck			boolean indicating if a pulse discharge test should be done as part of the check-up procedure
   * 		includeCycleData	boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
   * 		nCycles				number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
   * 		Crates				array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
   * 		Ccut				C rate of the cutoff current for the CV phase for the CCCV check-up, must be positive
   * 		profileName 		name of the csv file which contains the current profile for the pulse test
   * 								the first column contains the current in [A] (positive for discharge, negative for charge)
   * 								the second column contains the time in [sec] the current should be maintained
   * 								the profile must be a net discharge, i.e. sum (I*dt) > 0
   * 		profileLength		length of the current profiles for the pulse test (number of rows in the csv file)
   *
   * length = 1000 in default.
   *
   *
   * throws
   * 1014 		the input parameters are invalid
   * 1015 		the profile is invalid or the voltage limits are too small: the full voltage window was 'traversed' by just one repetition of the current profile
   * 				i.e.: 	starting from a cell at the lower voltage, you still hit the upper voltage while following the profile
   * 						or starting from a cell at the upper voltage, you still hit the lower voltage while following the profile
   * 				reduce the absolute value of the currents in the current profile (or reduce the durations) to produce a valid current profile
   */

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::profileAgeing starting.\n";

  //!< Check the input parameters
  bool vmax = Vma > c.getVmax(); //!< check if the maximum voltage is below the cell maximum voltage
  if (vmax)
    std::cerr << "Error in CyclerOld::profileAgeing. The maximum voltage " << Vma
              << " is too high. The maximum value is " << c.getVmax() << ".\n";
  bool vmin = Vmi < c.getVmin(); //!< check if the minimum voltage is above the cell minimum voltage
  if (vmin)
    std::cerr << "Error in CyclerOld::profileAgeing. The minimum voltage " << Vmi
              << " is too low. The minimum value is " << c.getVmin() << ".\n";
  bool Temin = Ti < settings::Tmin_Cell_K; //!< check the temperature is above 0 degrees, TMIN is defined in State.hpp
  if (Temin)
    std::cerr << "Error in CyclerOld::profileAgeing. The temperature " << Ti << "K is too low. The minimum value is 273.\n";
  bool Temax = Ti > settings::Tmax_Cell_K; //!< check the temperature is below 60 degrees, TMIN is defined in State.hpp
  if (Temax)
    std::cerr << "Error in CyclerOld::profileAgeing. The temperature " << Ti << " is too high. The maximum value is (273+60).\n";
  bool cycles = nrProfiles <= nrCap; //!< check the number of cycles between consecutive check-ups is lower than the total number of cycles
  if (cycles)
    std::cerr << "Error in CyclerOld::profileAgeing. The number of cycles between two check ups " << nrCap
              << " is lower than the total number of cycles " << nrProfiles << ".\n";
  if (vmax || vmin || Temin || Temax || cycles)
    throw 1014;

  //!<*********************************************************** 1 variables & settings ***********************************************************************

  //!< Variables
  double capi;                   //!< capacity at this step [Ah]
  double timei;                  //!< time spent in this step [sec]
  double ahi;                    //!< charge discharged in this step [Ah]
  double whi;                    //!< energy discharged in this step [Wh]
  double Ahtot = 0;              //!< cumulative charge throughput until now [Ah]
  double Whtot = 0;              //!< cumulative energy throughput until now [Wh]
  double timetot = 0;            //!< cumulative time until now [hour]
  int nrep = 0;                  //!< number of times the profile was repeated before a voltage limit was hit
  int nreptot = 0;               //!< total number of times the profile was repeated
  int vlim;                      //!< integer to indicate which voltage limit was hit while we were following the current profile
  bool Vlimhit = false;          //!< boolean to indicate if we need to re(dis)charge the cell to continue following the profile
  bool check = false;            //!< boolean to indicate if we need to do a check-up after re(dis)charging the cell
  bool blockDegradation = false; //!< account for degradation while following the current profile
  bool final = true;             //!< boolean to check if a final check-up is needed
  double Ccc = 1;                //!< C rate for the CC phase of the recharge or redischarge between profiles [-]
  double Ccut = 0.05;            //!< C rate of the cutoff current for the CV phases  of the recharge or redischarge between profiles [-]

  //!< Read the current profile
  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::profileAgeing is reading the current profile.\n";

  static thread_local std::vector<double> I(length), T(length); //!< arrays to store the current profile as doubles
  try {
    slide::loadCSV_2col(PathVar::data / nameI, I, T); //!< read the file
  } catch (int e) {
    //!< std::cout << "Throw test: " << 67 << '\n';
    std::cout << "error in CyclerOld::profileAgeing when reading the file with the current profile called "
              << nameI << ", error " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< Determine if the profile is a net charge or a net discharge
  const double aht = std::inner_product(I.begin(), I.end(), T.begin(), 0.0); //!< charge throughput of the profile
  const int sign = (aht > 0) ? -1 : 1;                                       //!< the profile is a net discharge (-1) or net charge (1)

  //!<*********************************************************** 2 cell initialisation ***********************************************************************

  if constexpr (printBool::printCyclerHighLevel)
    std::cout << "CyclerOld::profileAgeing is initialising the cell.\n";

  //!< set the temperature
  c.setT(Ti);    //!< set the cell temperature
  c.setTenv(Ti); //!< set the environmental temperature

  //!<(dis)charge the cell
  try {
    if (sign == -1)                                                       //!< we need to recharge to Vma because the profile is a discharge
      CC_V_CV_I(Ccc, Vma, Ccut, 2, blockDegradation, &ahi, &whi, &timei); //!< CC charge, followed by CV at a time step of 2 seconds
    else                                                                  //!< we need to discharge to Vmi because the profile is a charge
      CC_V_CV_I(Ccc, Vmi, Ccut, 2, blockDegradation, &ahi, &whi, &timei); //!< CC discharge, followed by CV at a time step of 2 seconds
  } catch (int e) {
    //!< std::cout << "Throw test: " << 68 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in CyclerOld::profileAgeing when bringing the cell to the initial voltage " << e << ". Throwing it on.\n";
    throw e;
  }

  //!< initial check up
  try {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::profileAgeing is doing an initial check-up.\n";
    checkUp(proc, 0, timetot, Ahtot, Whtot);
  } catch (int e) {
    //!< std::cout << "Throw test: " << 69 << '\n';
    if constexpr (printBool::printCrit)
      std::cout << "Error in a subfunction of CyclerOld::profileAgeing in the initial check-up " << e << ". Throwing it on.\n";
    throw e;
  }

  //!<*********************************************************** 3 age the cell by continuously following the profile and re(dis)charging ***********************************************************************

  //!< Cycle the battery by following the profile by repeating the 3 steps
  //!<		keep applying the profile until you hit a voltage limit
  //!<		re(dis)charge the cell
  //!<		do a check-up if needed
  while (nreptot < nrProfiles) { //!< loop to follow the profile the set number of times

    //!<*********************************************************** 3A  follow the profile until you hit a voltage limit ***********************************************************************
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::profileAgeing has applied the profile " << nreptot
                << " times and is starting the next series of profile repetitions.\n";
    try {
      Vlimhit = false;   //!< boolean to check when we've hit the voltage limit
      nrep = 0;          //!< counter for how often we can follow the profile before hitting the voltage limit
      check = false;     //!< boolean to check if we need to do a check-up
      while (!Vlimhit) { //!< loop to keep applying the profile until you hit a voltage limit

        //!< follow the profile
        vlim = followI(length, I, T, blockDegradation, limit, Vma, Vmi, &ahi, &whi, &timei);

        //!< update the throughput
        Ahtot += std::abs(ahi);
        Whtot += std::abs(whi);
        timetot += timei / 3600.0;
        nrep++;                             //!< increase the counter of repetitions before hitting the voltage limit
        nreptot++;                          //!< increase the total counter
        if (std::fmod(nreptot, nrCap) == 0) //!< the total number of repetitions is a multiple of the nr profiles between a check, so store that we need to do a check-up
          check = true;

        //!< check if we have hit the voltage limit
        if (sign == -1)                           //!< the profile is a net discharge, so stop if you hit the lower voltage limit
          Vlimhit = (vlim == -1) || (vlim == 10); //!< a return-integer of -1 or 10 means the lower voltage was hit
        else
          Vlimhit = (vlim == 1) || (vlim == 10); //!< a return-integer of 1 or 10 means the upper voltage was hit
      }
    } catch (int e) {
      //!< std::cout << "Throw test: " << 70 << '\n';
      if constexpr (printBool::printCrit) {
        std::cout << "Error in CyclerOld::profileAgeing when cycling the cell with the current profile according to profile regime " << ID << ". Error encountered is " << e << ". Stop cycling now.";
        std::cout << " We have done " << nreptot << " repetitions instead of " << nrProfiles
                  << " and the capacity last measured is " << capi << ".\n";
      }

      //!< we probably cannot do a full check-up procedure because the cell is in an illegal state.
      //!< Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
      checkUp_batteryStates(proc.blockDegradation, false, nreptot, timetot, Ahtot, Whtot);

      final = false; //!< skip the final check-up
      break;         //!< stop cycling
    }

    //!< ensure we could repeat the profile multiple times before hitting the voltage limit
    if (nrep <= 1) {
      //!< this means that you can never fully follow the current profile in the specified voltage window.
      //!< this function does not allow this sort of behaviour, and an error is thrown.
      //!< the user has to decrease the currents in the profile, or enlarge the voltage window (if possible)
      std::cerr << "ERROR in CyclerOld::profileAgeing: the profile " << ID
                << " could only be repeated once before a voltage limit was hit.\n";
      //!< throw the error
      throw 1015;
    }

    //!<*********************************************************** 3B re(dis)charge ***********************************************************************
    //!< re(dis)charge the cell to the other voltage limit from the one you have hit
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::profileAgeing has applied the profile " << nreptot << " times and is going to re(dis)charge the cell.\n";
    try {
      if (sign == -1)                                                       //!< we need to recharge to Vma
        CC_V_CV_I(Ccc, Vma, Ccut, 2, blockDegradation, &ahi, &whi, &timei); //!< CC charge, followed by CV at a time step of 2 seconds
      else                                                                  //!< we need to discharge to Vmi
        CC_V_CV_I(Ccc, Vmi, Ccut, 2, blockDegradation, &ahi, &whi, &timei); //!< CC discharge, followed by CV at a time step of 2 seconds
      Ahtot += std::abs(ahi);
      Whtot += std::abs(whi);
      timetot += timei / 3600.0;
    } catch (int e) {
      //!< std::cout << "Throw test: " << 71 << '\n';
      if constexpr (printBool::printCrit) {
        std::cout << "Error in CyclerOld::profileAgeing when re(dis)charging the cell according to profile regime " << ID << ". Error encountered is " << e << ". Stop cycling now."
                  << " We have done " << nreptot << " repetitions instead of " << nrProfiles << " and the capacity last measured is " << capi << ".\n";
      }

      //!< we probably cannot do a full check-up procedure because the cell is in an illegal state.
      //!< Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
      checkUp_batteryStates(proc.blockDegradation, false, nreptot, timetot, Ahtot, Whtot);

      final = false; //!< skip the final check-up
      break;         //!< stop cycling
    }

    //!<*********************************************************** 3C check-up ***********************************************************************
    //!< do a check-up if needed
    if (check) {
      if constexpr (printBool::printCyclerHighLevel)
        std::cout << "CyclerOld::profileAgeing has applied the profile " << nreptot << " times and is going to do a check-up.\n";
      try {
        capi = checkUp(proc, nreptot, timetot, Ahtot, Whtot);
        check = false;
      } catch (int e) {
        //!< std::cout << "Throw test: " << 72 << '\n';
        if constexpr (printBool::printCrit) {
          std::cout << "Error in CyclerOld::profileAgeing when doing a check-up according to profile regime " << ID << ". Error encountered is " << e << ". Stop cycling now."
                    << " We have done " << nreptot << " repetitions instead of " << nrProfiles << " and the capacity last measured is " << capi << ".\n";
        }

        //!< we probably cannot do a full check-up procedure because the cell is in an illegal state.
        //!< Therefore, only write the BatteryStates with a capacity of 0 to indicate something went wrong
        checkUp_batteryStates(proc.blockDegradation, false, nreptot, timetot, Ahtot, Whtot);

        final = false; //!< skip the final check-up
        break;         //!< stop cycling
      }

      //!< End the experiment if the cell capacity has decreased too much
      if (capi < c.getNominalCap() / 2.0) { //!< end simulation if the cell has only 50% of capacity left
        std::cout << "CyclerOld::ProfileAgeing has finished cycling regime " << ID << " early do to too little remaining capacity."
                  << " We have done " << nreptot << " repetitions of the profile instead of " << nrProfiles << " and the remaining capacity now is " << capi << ".\n";
        final = false; //!< skip the final check-up because we just did one
        break;         //!< stop cycling
      }
    }

    //!< keep repeating these 3 steps (repeat profile until voltage limit, re(dis)charge, check-up) until you have done enough repetitions

  } //!< end loop of profile ageing

  //!<*********************************************************** 4 final check-up ***********************************************************************
  //!< Don't do the check-up if the simulation ended due to an error (or shortage in capacity)
  if (final) {
    if constexpr (printBool::printCyclerHighLevel)
      std::cout << "CyclerOld::profileAgeing has applied the profile " << nreptot << " times and is going to do a final check-up.\n";
    try {
      checkUp(proc, nreptot, timetot, Ahtot, Whtot);
    } catch (int e) {
      //!< std::cout << "Throw test: " << 73 << '\n';
      if constexpr (printBool::printCrit)
        std::cout << "Error in a subfunction of CyclerOld::profileAgeing in the final check-up " << e << ". Throwing it on.\n";
      throw e;
    }
  }

  if constexpr (printBool::printCyclerFunctions)
    std::cout << "CyclerOld::profileAgeing terminating.\n";
}
