/*
 * determine_characterisation.cpp
 *
 * The functions below can be used to fit the 'rate parameters' or 'characterisation parameters' of a cell at one temperature (currently 25 degrees):
 * 		diffusion constants (Dp and Dn)
 * 		rate constants of the main li-insertion reaction (kp and kn)
 * 		DC resistance of the cell (r)
 *
 * The user has to supply voltage measurements for a number of CC CV cycles (i.e. charge and discharge, both with a CC and CV phase).
 * These have to come in csv files with two columns
 * 		The first column has to have the charge throughput in Ah, starting at 0 and increasing as the cell is charged or discharged.
 * 		The second column has to have the cell voltage at the corresponding point.
 * In the top-level function ('fitCharacterisationAtReferenceT), the user has to define which cycles this are exactly (i.e. what C rates for the CC phase and which current threshold for the CV)
 *
 * The script then uses a hierarchical search algorithm to find the parameters which minimise the error between the simulated and measured voltage curves.
 * The results are written in a csv file. If the user wants to use these parameters, they have to be copied to the constructor of the respective cell-class.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include "../cells/cells.hpp"
#include "determine_OCV.hpp"
#include "CyclerOld.hpp"
#include "Cycler.hpp"
#include "determine_characterisation.hpp"
#include "../utility/utility.hpp"

#include <array>
#include <thread>
#include <algorithm>

namespace slide {
bool CCCV_fit(Cell_SPM c1, double Crate, double Ccut, double Tref, double Dp, double Dn, double kp, double kn, double R, const struct OCVparam &ocvfit, const struct slide::Model_SPM &M,
              slide::XYdata_vv &Vsim, slide::XYdata_vv &Tsim)
{
  /*
   * Function which simulates a full CC CV (dis)charge with the given rate parameters at reference temperature.
   * The reference temperature is defined in Cell (Tref), and is 25 degrees.
   * If you have data at a different reference temperature, you have to change the value of Tref in the constructor of Cell_SPM.
   *
   * IN
   * Crate 	the C rate to be used for the CC phase
   * 				> 0 for discharge
   * 				< 0 for charge
   * Ccut 	Crate of the cutoff current below which the CV phase should be stopped, > 0, [-]
   * Tref 	temperature at which the characterisation is done [K]
   * Dp		diffusion constant of the positive electrode at reference temperature
   * Dn		diffusion constant of the negative electrode at reference temperature
   * kp		rate constant of the lithium insertion at the positive electrode at reference temperature
   * kn		rate constant of the lithium insertion at the negative electrode at reference temperature
   * R		DC resistance of the cell
   * ocvfit 	structure with the parameters determined by the functions in determineOCV.cpp (the struct is defined in determineCharacterisation.h)
   * M 		structure with the matrices of the spatial discretisation of the solid diffusion PDE
   * 				defined in Cell.hpp
   *
   * OUT
   * Vsim		Matrix with the cell voltage [V] in the second column and the charge throughput [Ah] in the first one
   * Tsim		Matrix with the cell temperature [K] in the second column and the charge throughput [Ah] in the first one
   */

  //!< *********************************************************** 1 variables ***********************************************************************
  int verbose = 0; //!< integer deciding how verbose the simulation should be
                   //!< The higher the number, the more output there is.
                   //!< Recommended value is 1, only use higher value for debugging
                   //!< From 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
                   //!< 	0 	almost no messages are printed, only in case of critical errors related to illegal parameters
                   //!< 	1 	error messages are printed in case of critical errors which might crash the simulation
                   //!< 	2 	all error messages are printed, whether the simulation can recover from the errors or not
                   //!< 	3 	on top of the output from 2, a message is printed every time a function in the CyclerOld and BasicCycler is started and terminated
                   //!< 	4 	on top of the output from 3, the high-level flow of the program in the CyclerOld is printed (e.g. 'we are going to discharge the cell')
                   //!< 	5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
                   //!< 	6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
                   //!< 	7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated

  //!< Check all parameters are positive
  if (Dp <= 0 || Dn <= 0 || kp <= 0 || kn <= 0 || R < 0)
    return false;

  c1.setCharacterisationParam(Dp, Dn, kp, kn, R);

  //!< time steps
  double dt = 2.0;      //!< time step for cycling [s]
  double Istep = 0.1;   //!< current step for ramping, indicating how fast the current can change per 'ramp time step', [A s-1]
  double tstep = 0.001; //!< time step for ramping [s]
                        //!< the current can change at Istep/tstep, so currently 0.1A per 1 ms.

  //!< variables
  const auto s = c1.getStates(); //!< initial state of the cell, used to recover after an error
  std::span<const double> s_span{ s };

  double Ccut2 = 0.05;                              //!< Crate for the cutoff current when bringing the cell to the initial state (i.e. charge the cell first before you simulate the CCCV discharge)
  double Vset;                                      //!< voltage at which the CCCV cycle should end (i.e. the minimum voltage if you are simulating a discharge)
  bool blockDegradation = true;                     //!< don't account for degradation while doing the cycles
  std::string ID = "CharacterisationFit";           //!< identification std::string for the CyclerOld
  int timeCycleData = -1;                           //!< time interval at which cycling data has to be recorded [s]
                                                    //!< <0 means no folder is created and no data is stored
                                                    //!< 	0 means no data is recorded but a folder is still created for later use
                                                    //!<  >0 means data is recorded approximately every so many seconds
  CyclerOld cycler(c1, ID, verbose, timeCycleData); //!< Make the CyclerOld
  double ahi, whi, timei;                           //!< feedback variables we don't need

  //!< *********************************************************** 2 (dis)charge ***********************************************************************

  //!< We are going to start trying to simulate the (dis)charge with a normal time step (for time integration and ramping)
  //!< 	But smaller diffusion constants lead to decreased numerical stability and other discretiation errors
  //!< 	Therefore, these normal time steps might leads to errors in the code
  //!< So instead, there are loops which iteratively decrease the time steps to check if that solves the error.

  bool finished = false; //!< boolean to indicate whether we have completed the discharge

  //!< loop to decrease the integration time step (dt)
  for (int i = 0; i < 2 && !finished; i++) {
    //!< Set the ramping time steps to their normal value
    Istep = 0.1;   //!< current step for ramping
    tstep = 0.001; //!< time step for ramping

    //!< loop to decrease the time steps for ramping
    for (int j = 0; j < 2 && !finished; j++) {

      //!< Try to simulate the (dis)charge
      try {

        //!< restore the original battery state in case an error occurred earlier in the loop
        c1.setStates(s_span, false, false); //!< Does not throw.
        //!< c1.setRamping(Istep, tstep);			//!< Does not throw.
        cycler.setCyclingDataTimeResolution(0); //!< don't collect cycling data during the charging //!< Does not throw. (anymore)

        //!< Bring the cell to the correct soc before simulating the (dis)charge
        if (Crate > 0) { //!< simulate a discharge
          //!< first fully charge the cell to the maximum voltage at 1C
          cycler.CC_V_CV_I(1, ocvfit.Vmax, Ccut2, dt, blockDegradation, &ahi, &whi, &timei);
          Vset = ocvfit.Vmin;
        } else { //!< simulate a charge
          //!< first fully discharge the cell to the minimum voltage at 1C
          cycler.CC_V_CV_I(1, ocvfit.Vmin, Ccut2, dt, blockDegradation, &ahi, &whi, &timei);
          Vset = ocvfit.Vmax;
        }

        //!< simulate the CC CV (dis)charge
        cycler.setCyclingDataTimeResolution(dt); //!< collect cycling data of every time step
        cycler.CC_V_CV_I(std::abs(Crate), Vset, Ccut, dt, blockDegradation, &ahi, &whi, &timei);

        //!< If we get here, no errors were thrown and we can leave the loop
        finished = true; //!< indicate we have finished the simulation
        break;           //!< leave the loops
      } catch (int e) {
        //!< std::cout << "Throw test: " << 78 << '\n';
        //!< An error occurred while simulating the (dis)charge #TODO -> we throw here very much.
        //		Istep = Istep / 10.0; //!< reduce the ramping parameters by a factor of 10
        //!< now we are still in the loop which decreases the ramping time steps ('final' is still false), so you will try again
        //!< and the original battery state will be restored at the start of the loop, so the illegal battery state will be 'forgotten'
      }
    } //!< end loop to decrease the ramping time step

    //!< if we haven't finished the cycle, try again with a smaller time step
    if (!finished) //!< now we are still in the loop which decreases the time step, so you will try again
      dt /= 10.0;  //!< and the original battery state will be restored at the start of the loop, so the illegal battery state will be 'forgotten'
  }                //!< end loop to decrease the integration time step

  //!< *********************************************************** 3 output ***********************************************************************

  if (!finished || ahi == 0)
    return false; //!< the simulation wasn't successful. This happens if 'finished' is still false or if ahi == 0 i.e. no charge could be discharged

  //!< Get the cell voltage from the simulated (dis)charge from the CyclerOld
  cycler.returnCyclingData(Vsim.x, Vsim.y, Tsim.y);
  Tsim.x = Vsim.x;

  return true; //!< Simulation was successful.
}

void CCCV(double Crate, double Ccut, double Tref, double Dp, double Dn, double kp, double kn, double R, const struct OCVparam &ocvfit, const struct Model_SPM *M,
          slide::XYdata_vv &Vsim, slide::XYdata_vv &Tsim)
{
  /*
   * Function which simulates a full CC CV (dis)charge with the given rate parameters at reference temperature.
   * The reference temperature is defined in Cell (Tref), and is 25 degrees.
   * If you have data at a different reference temperature, you have to change the value of Tref in the constructor of Cell_SPM.
   *
   * IN
   * Crate 	the C rate to be used for the CC phase
   * 				> 0 for discharge
   * 				< 0 for charge
   * Ccut 	Crate of the cutoff current below which the CV phase should be stopped, > 0, [-]
   * Tref 	temperature at which the characterisation is done [K]
   * Dp		diffusion constant of the positive electrode at reference temperature
   * Dn		diffusion constant of the negative electrode at reference temperature
   * kp		rate constant of the lithium insertion at the positive electrode at reference temperature
   * kn		rate constant of the lithium insertion at the negative electrode at reference temperature
   * R		DC resistance of the cell
   * ocvfit 	structure with the parameters determined by the functions in determineOCV.cpp (the struct is defined in determineCharacterisation.h)
   * M 		structure with the matrices of the spatial discretisation of the solid diffusion PDE
   * 				defined in Cell.hpp
   *
   * OUT
   * Vsim		Matrix with the cell voltage [V] in the second column and the charge throughput [Ah] in the first one
   * Tsim		Matrix with the cell temperature [K] in the second column and the charge throughput [Ah] in the first one
   *
   * THROWS
   * 10002 	an error occurred while cycling the cell, and reducing the time step didn't help
   * 10003	one of the parameters is negative, which is not allowed;
   */

  //!< *********************************************************** 1 variables ***********************************************************************
  int verbose = 0; //!< integer deciding how verbose the simulation should be
                   //!< The higher the number, the more output there is.
                   //!< Recommended value is 1, only use higher value for debugging
                   //!< From 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
                   //!< 	0 	almost no messages are printed, only in case of critical errors related to illegal parameters
                   //!< 	1 	error messages are printed in case of critical errors which might crash the simulation
                   //!< 	2 	all error messages are printed, whether the simulation can recover from the errors or not
                   //!< 	3 	on top of the output from 2, a message is printed every time a function in the CyclerOld and BasicCycler is started and terminated
                   //!< 	4 	on top of the output from 3, the high-level flow of the program in the CyclerOld is printed (e.g. 'we are going to discharge the cell')
                   //!< 	5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
                   //!< 	6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
                   //!< 	7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated

  //!< Make a cell of the sub-class Cell_SPM. This sub-class defines some extra functions to change its parameters
  Cell_SPM c1(M, verbose);

  //!< Check all parameters are positive
  if (Dp <= 0 || Dn <= 0 || kp <= 0 || kn <= 0 || R < 0)
    throw 10003;

  //!< Set the characterisation parameters of the cell to the ones given as input to this function
  c1.setOCVcurve(ocvfit.namepos, ocvfit.nameneg);
  c1.setInitialConcentration(ocvfit.cmaxp, ocvfit.cmaxn, ocvfit.lifracpini, ocvfit.lifracnini);
  c1.setGeometricParameters(ocvfit.cap, ocvfit.elec_surf, ocvfit.ep, ocvfit.en, ocvfit.thickp, ocvfit.thickn);
  c1.setCharacterisationParam(Dp, Dn, kp, kn, R);
  //!< c1.setVlimits(ocvfit.Vmax, ocvfit.Vmin); #TODO
  c1.setT(Tref);    //!< set the temperature of the cell to the given value
  c1.setTenv(Tref); //!< set the environmental temperature to the given value

  //!< time steps
  double dt = 2.0;      //!< time step for cycling [s]
  double Istep = 0.1;   //!< current step for ramping, indicating how fast the current can change per 'ramp time step', [A s-1]
  double tstep = 0.001; //!< time step for ramping [s]
                        //!< the current can change at Istep/tstep, so currently 0.1A per 1 ms.

  //!< variables
  const auto s = c1.getStates();          //!< initial state of the cell, used to recover after an error
  const auto Iini = c1.getI();            //!< initial current of the cell, used to recover after an error
  double Ccut2 = 0.05;                    //!< Crate for the cutoff current when bringing the cell to the initial state (i.e. charge the cell first before you simulate the CCCV discharge)
  double Vset;                            //!< voltage at which the CCCV cycle should end (i.e. the minimum voltage if you are simulating a discharge)
  bool blockDegradation = true;           //!< don't account for degradation while doing the cycles
  std::string ID = "CharacterisationFit"; //!< identification std::string for the CyclerOld
  int timeCycleData = -1;                 //!< time interval at which cycling data has to be recorded [s]
                                          //!< <0 means no folder is created and no data is stored
                                          //!< 	0 means no data is recorded but a folder is still created for later use
                                          //!<  >0 means data is recorded approximately every so many seconds

  CyclerOld cycler(c1, ID, verbose, timeCycleData); //!< Make the CyclerOld
  double ahi, whi, timei;                           //!< feedback variables we don't need

  //!< *********************************************************** 2 (dis)charge ***********************************************************************

  //!< We are going to start trying to simulate the (dis)charge with a normal time step (for time integration and ramping)
  //!< 	But smaller diffusion constants lead to decreased numerical stability and other discretiation errors
  //!< 	Therefore, these normal time steps might leads to errors in the code
  //!< So instead, there are loops which iteratively decrease the time steps to check if that solves the error.

  bool finished = false; //!< boolean to indicate whether we have completed the discharge

  //!< loop to decrease the integration time step (dt)
  for (int i = 0; i < 2 && !finished; i++) {

    //!< Set the ramping time steps to their normal value
    Istep = 0.1;   //!< current step for ramping
    tstep = 0.001; //!< time step for ramping

    //!< loop to decrease the time steps for ramping
    for (int j = 0; j < 2 && !finished; j++) {

      //!< Try to simulate the (dis)charge
      try {

        //!< restore the original battery state in case an error occurred earlier in the loop
        c1.setStates(s, Iini);
        c1.setRamping(Istep, tstep);
        cycler.setCyclingDataTimeResolution(0); //!< don't collect cycling data during the charging

        //!< Bring the cell to the correct soc before simulating the (dis)charge
        if (Crate > 0) { //!< simulate a discharge
          //!< first fully charge the cell to the maximum voltage at 1C
          cycler.CC_V_CV_I(1, ocvfit.Vmax, Ccut2, dt, blockDegradation, &ahi, &whi, &timei);
          Vset = ocvfit.Vmin;
        } else { //!< simulate a charge
          //!< first fully discharge the cell to the minimum voltage at 1C
          cycler.CC_V_CV_I(1, ocvfit.Vmin, Ccut2, dt, blockDegradation, &ahi, &whi, &timei);
          Vset = ocvfit.Vmax;
        }

        //!< simulate the CC CV (dis)charge
        cycler.setCyclingDataTimeResolution(dt); //!< collect cycling data of every time step
        cycler.CC_V_CV_I(std::abs(Crate), Vset, Ccut, dt, blockDegradation, &ahi, &whi, &timei);

        //!< If we get here, no errors were thrown and we can leave the loop
        finished = true; //!< indicate we have finished the simulation
        break;           //!< leave the loops
      } catch (int e) {
        //!< std::cout << "Throw test: " << 79 << '\n';

        //!< An error occurred while simulating the (dis)charge
        Istep = Istep / 10.0; //!< reduce the ramping parameters by a factor of 10
                              //!< now we are still in the loop which decreases the ramping time steps ('final' is still false), so you will try again
                              //!< and the original battery state will be restored at the start of the loop, so the illegal battery state will be 'forgotten'
      }
    } //!< end loop to decrease the ramping time step

    //!< if we haven't finished the cycle, try again with a smaller time step
    if (!finished) {
      dt = dt / 10.0;
      //!< now we are still in the loop which decreases the time step, so you will try again
      //!< and the original battery state will be restored at the start of the loop, so the illegal battery state will be 'forgotten'
    }
  } //!< end loop to decrease the integration time step

  //!< *********************************************************** 3 output ***********************************************************************

  //!< Throw an error if the simulation wasn't successful. This happens if 'finished' is still false or if ahi == 0 i.e. no charge could be discharged
  if (!finished || ahi == 0)
    throw 10002;

  //!< Get the cell voltage from the simulated (dis)charge from the CyclerOld

  cycler.returnCyclingData(Vsim.x, Vsim.y, Tsim.y);
  Tsim.x = Vsim.x;
}

void fitDiffusionAndRate(int hierarchy, int ir, double R, slide::FixedData<double> Dp_space, slide::FixedData<double> Dn_space,
                         slide::FixedData<double> kp_space, slide::FixedData<double> kn_space,
                         std::vector<slide::XYdata_vv> &Vdata_all, double weights[],
                         double Crates[], double Ccuts[], double Tref, const struct OCVparam &ocvfit,
                         double *err, std::array<double, 5> &par)
{
  /*
   * Function which goes trough the specified search space for Dp, Dn, kp and kn, for a constant value of the DC resistance R.
   * The steps can be taken logarithmically (e.g. e-10, e-11, etc.) or linearly (e.g. 1e-10, 2e-10, etc.)
   * For each combination of the 4 parameters, the voltage of the CCCV cycles is simulated and the error with the measured voltages is calculated.
   * The best fit (i.e. 4 parameters with the lowest error on the OCV curve) is returned
   *
   * IN
   * hierarchy level of the hierarchy in which we are
   * ir 		index of where in the search space for r we are
   * R 		DC resistance of the cell to be used in the search
   * r_space  search space for r (DC resistance of the cell)
   * Dp_space search space for Dp (cathode diffusion constant)
   * Dn_space search space for Dn (anode diffusion constant)
   * kp_space search space for kp (cathode rate constant)
   * kn_space search space for kn (anode rate constant)
   *
   * nCCCV 	number of CCCV cycles
   * weights 	array with the weight that should be given to each CCCV curve when calculating the overall error
   * 			the sum of all weights should be one
   * names 	array with the names of the csv files with the measurements for the cell
   * Crates 	array with the C rates of the CC phases of each experiment, >0 for discharge, <0 for charge
   * Ccuts 	array with the C rates of the current threshold for the CV phase of each experiment, >0 . (set to a very large value if you don't want a CV phase)
   * Tref 	temperature at which the characterisation is done [K]
   * ocvfit 	structure with the values of the OCV parameters determined by determineOCV::estimateOCVparam
   *
   * OUT
   * err 		error of the best fit
   * par 		values of r, Dp, Dn, kp and kn which achieved the best fit
   *
   */

  //!< *********************************************************** 1 variables ***********************************************************************

  //!< Variables
  //!< auto M = Model_SPM::makeModel(); //!< structure with the matrices for the spatial discretisation for the solid diffusion PDE
  //!< constexpr double dt = 2;	  //!< time step to be used for the simulation
  slide::XYdata_vv Vsim, Tsim; //!< arrays to store the simulation results #TODO -> static thread_local
  double errmin = 10000000000; //!< lowest error encountered so far

  //!< ************************************************ Create an initial cell. ************************************************
  const int verbose = 0;
  Cell_SPM cell_init{};

  //!< Set the characterisation parameters of the cell to the ones given as input to this function
  cell_init.setOCVcurve(ocvfit.namepos, ocvfit.nameneg);
  cell_init.setInitialConcentration(ocvfit.cmaxp, ocvfit.cmaxn, ocvfit.lifracpini, ocvfit.lifracnini);
  cell_init.setGeometricParameters(ocvfit.cap, ocvfit.elec_surf, ocvfit.ep, ocvfit.en, ocvfit.thickp, ocvfit.thickn);
  //!< cell_init.setVlimits(ocvfit.Vmax, ocvfit.Vmin); #TODO set Vlimits?
  cell_init.setT(Tref);    //!< set the temperature of the cell to the given value
  cell_init.setTenv(Tref); //!< set the environmental temperature to the given value

  std::string ID = "CharacterisationFit"; //!< identification std::string for the CyclerOld
  int timeCycleData = -1;                 //!< time interval at which cycling data has to be recorded [s]
  //!< CyclerOld cycler(cell_init, ID, verbose, timeCycleData); //!< Make the CyclerOld

  //!< *********************************************************** 2 loop through the search space ***********************************************************************

  for (const auto Dp : Dp_space)       //!< scan the search space for Dp
    for (const auto Dn : Dn_space)     //!< scan the search space for Dn
      for (const auto kp : kp_space)   //!< scan the search space for kp
        for (const auto kn : kn_space) //!< scan the search space for kn
        {
          //!< Calculate the error for this set of parameters
          double errcomb = 0;                             //!< initialise the combined error of all CCCV experiments for this combination of Dp, Dn, kp and kn to 0
          for (size_t i = 0; i < Vdata_all.size(); i++) { //!< loop through all CCCV cycles

            //!< Simulate this CCCV experiment
            Vsim.clear(), Tsim.clear();
            auto flag = CCCV_fit(cell_init, Crates[i], Ccuts[i], Tref, Dp, Dn, kp, kn, R, ocvfit, M, Vsim, Tsim);

            if (flag) {
              const double erri = calculateError(false, Vdata_all[i], Vsim); //!< calculate the error of this CCCV cycle with the given parameters
              errcomb += std::abs(erri) * weights[i];                        //!< calculate the total (weighted) error
            } else                                                           //!< if flag is false, an error occured while simulating. This means the parameters were infeasible. High cost.
            {
              errcomb = 10000000000;
              break;
            }

          } //!< loop for CCCV experiments

          //!< Store the minimum error
          if (errcomb < errmin) { //!< check if the error of this combination is better than the best fit so far
            par = { R, Dp, Dn, kp, kn };
            errmin = errcomb;
          }
        } //!< loop for kn

  *err = errmin; //!< return the lowest error
}

void hierarchicalCharacterisationFit(int hmax, slide::FixedData<double> r_space, slide::FixedData<double> Dp_space,
                                     slide::FixedData<double> Dn_space, slide::FixedData<double> kp_space,
                                     slide::FixedData<double> kn_space, std::vector<slide::XYdata_vv> &Vdata_all,
                                     double weights[], double Crates[], double Ccuts[], double Tref,
                                     const struct OCVparam &ocvfit, double *err, std::array<double, 5> &par)
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
   * then in hierarchical level n, the step size is dp* (2/(nstep-1)^(n-1)
   * So after 'hmax' level, the accuracy is dp* (2/(nstep-1)^(hmax-1)
   *
   * IN
   * hmax  	number of hierarchical steps we should take.
   * 			the higher the value, the higher the accuracy of the fit will be but the longer the calculation will take
   *
   * r_space  search space for r (DC resistance of the cell)
   * Dp_space search space for Dp (cathode diffusion constant)
   * Dn_space search space for Dn (anode diffusion constant)
   * kp_space search space for kp (cathode rate constant)
   * kn_space search space for kn (anode rate constant)
   *
   * weights 	array with the weight that should be given to each CCCV curve when calculating the overall error
   * 			the sum of all weights should be one
   * Crates 	array with the C rates of the CC phases of each experiment, >0 for discharge, <0 for charge
   * Ccuts 	array with the Crate of the current threshold for the CV phase of each experiment, >0 . (set to a very large value if you don't want a CV phase)
   * Tref 	temperature at which the characterisation is done [K]
   * ocvfit 	structure with the values of the OCV parameters determined by determineOCV::estimateOCVparam
   *
   * OUT
   * err 		error of the best fit
   * par 		values of r, Dp, Dn, kp and kn which achieved the best fit
   */

  //!< variables
  std::vector<std::array<double, 5>> par_arr(r_space.size()); //!< array of parameters [R Dp Dn kp kn] giving the lowest error for that resistance
  std::vector<double> err_arr(r_space.size());

  //!< Loop for each level in the search
  for (int h = 0; h < hmax; h++) {

    //!< print the search space of this level
    std::cout << "Start hierarchy level " << h << " with the following search spaces: \n"
              << "R: from " << r_space.front() << " to " << r_space.back() << " in " << r_space.size() << " steps with magnitude " << r_space.dstep() << '\n'
              << "Dp: from " << Dp_space.front() << " to " << Dp_space.back() << " in " << Dp_space.size() << " steps with magnitude " << Dp_space.dstep() << '\n'
              << "Dn: from " << Dn_space.front() << " to " << Dn_space.back() << " in " << Dn_space.size() << " steps with magnitude " << Dn_space.dstep() << '\n'
              << "kp: from " << kp_space.front() << " to " << kp_space.back() << " in " << kp_space.size() << " steps with magnitude " << kp_space.dstep() << '\n'
              << "kn: from " << kn_space.front() << " to " << kn_space.back() << " in " << kn_space.size() << " steps with magnitude " << kn_space.dstep() << '\n';

    //!< Calculate the best fit in this level

    auto task_indv = [&](int i) {
      fitDiffusionAndRate(h, i, r_space[i], Dp_space, Dn_space, kp_space, kn_space, Vdata_all, weights, Crates, Ccuts, Tref, ocvfit, &err_arr[i], par_arr[i]);
    };

    slide::run(task_indv, r_space.size());

    const auto minIndex = std::min_element(err_arr.begin(), err_arr.end()) - err_arr.begin();

    writeCharacterisationParam(h, par_arr[minIndex], *err); //!< Print the best fit, and write in a CSV file

    //!< Update the search space
    //!< suppose the optimal value for a parameter p was in the search space at index i
    //!< then the new search space has as minimum value the value of p at i-1 and as maximum the value of p at i+1

    //!< only the first search level is logarithmic, afterwards the search has linear steps
    const auto [r, Dp, Dn, kp, kn] = par_arr[minIndex];

    r_space = slide::linspace_fix(r_space.prev(r), r_space.next(r), r_space.size());
    Dp_space = slide::linspace_fix(Dp_space.prev(Dp), Dp_space.next(Dp), Dp_space.size());
    Dn_space = slide::linspace_fix(Dn_space.prev(Dn), Dn_space.next(Dn), Dn_space.size());
    kp_space = slide::linspace_fix(kp_space.prev(kp), kp_space.next(kp), kp_space.size());
    kn_space = slide::linspace_fix(kn_space.prev(kn), kn_space.next(kn), kn_space.size());

    //!< Make the output parameters
    *err = err_arr[minIndex];
    par = par_arr[minIndex];
  }
}
void estimateCharacterisation()
{
  /*
   * Function which will find the parameters for the diffusion constants, rate constants and DC resistance
   * which best fit the measured voltage curve of the user.
   *
   * As input, the user has to supply csv files with voltage measurements of some full CCCV (dis)charges.
   * E.g. starting from the maximum voltage, do a CC discharge at 1C followed by a CV discharge until the current is below some threshold.
   * The first column has to have the charge throughput in Ah, starting at 0 and increasing as the cell is charged or discharged.
   * The second column has to have the cell voltage at the corresponding point.
   * The cycles have to be 'full', i.e. if it is a charge, than the first voltage must be Vmax and the last voltage must be Vmin (i.e. from fully charged to fully discharged) and vice versa.
   *
   * The user can specify in the code below which C rates were used for the CC phases, and what the cutoff current threshold is.
   * Additionally, the user has to specify the OCV-related parameters for the cell, as calculated by determineOCV.
   *
   * The function finds optimal values for 5 parameters:
   * 		diffusion constants (Dp and Dn)
   * 		rate constants of the main li-insertion reaction (kp and kn)
   * 		DC resistance of the cell (r)
   *
   * The values of the parameters which give the best fit are written in a CSV file.
   * Also, the simulated voltage curves with these parameters are written to csv files so it can be compared with the measured ones.
   */

  //!< *********************************************************** 1 USER INPUT ***********************************************************************

  //!< Specify the characterisation tests for which data is available
  double Tref = PhyConst::Kelvin + 25; //!< Temperature at which the characterisation should be done [K]

  std::string names[] = { "Characterisation_0.2C_CC_discharge.csv",
                          "Characterisation_0.5C_CC_discharge.csv",
                          "Characterisation_1C_CC_discharge.csv",
                          "Characterisation_2C_CC_discharge.csv",
                          "Characterisation_3C_CC_discharge.csv" }; //!< Name of the files with the voltage curve from the cell

  double Crates[] = { 0.2, 0.5, 1, 2, 3 };        //!< C rates of the CC phases for each voltage curve, <0 for charge, >0 for discharge
  double Ccuts[] = { 100, 100, 100, 100, 100 };   //!< C rates of the current threshold for the CV phases [A], > 0
                                                  //!< 	a very high value (above the Crate of the CC) will avoid there is a CV phase at all, so there is only a CC phase
  double weights[] = { 0.2, 0.2, 0.2, 0.2, 0.2 }; //!< array with the weight that should be attributed to each individual CCCV curve when calculating the overall error
                                                  //!< 	no restrictions are imposed, but it is recommended that all weights sum up to 1
                                                  //!< 	if a weight is set to 0, that curve is ignored. If a weight is negative, the error for that curve will be maximised so this is not recommended.

  const size_t nCCCV = std::size(names); //!< number of data sets for CCCV (dis)charges

  //!< names of output files
  const std::string nameparam = "characterisationFit_parameters.csv"; //!< name of the output csv file in which the optimal parameters will be written
  const std::string nameCCCVfit = "characterisationFit_";             //!< prefix appended before the name of the output data files with the simulations for the best fit
                                                                      //!< each file will have 3 columns: charge throughput, voltage and temperature

  //!< Specify the OCV parameters (calculated by determineOCV::estimateOCVparameters)
  OCVparam ocvfit;
  ocvfit.elec_surf = 0.0982;  //!< electrode surface
  ocvfit.ep = 0.5;            //!< volume fraction of active material in the cathode
  ocvfit.en = 0.5;            //!< volume fraction of active material in the anode
  ocvfit.thickp = 70e-6;      //!< thickness of the cathode
  ocvfit.thickn = 73.5e-6;    //!< thickness of the anode
  ocvfit.lifracpini = 0.6862; //!< lithium fraction in the cathode at 50% soC
  ocvfit.lifracnini = 0.4843; //!< lithium fraction in the anode at 50% SOC
  ocvfit.cmaxp = 51385;       //!< maximum lithium concentration in the cathode [mol m-3]
  ocvfit.cmaxn = 30555;       //!< maximum lithium concentration in the anode [mol m-3]
  ocvfit.cap = 2.7;           //!< the capacity of the cell [Ah]
  ocvfit.Vmax = 4.2;          //!< maximum voltage of the cell [V]
  ocvfit.Vmin = 2.7;          //!< minimum voltage of the cell [V]

  ocvfit.namepos = "OCVfit_cathode.csv"; //!< name of the CSV file with the cathode OCV curve
  ocvfit.nameneg = "OCVfit_anode.csv";   //!< name of the CSV file with the anode OCV curve
  ocvfit.np = 49;                        //!< number of points in the cathode OCV curve
  ocvfit.nn = 63;                        //!< number of points in the anode OCV curve

  //!< ****************************************** 2 define the search space for fitting parameters ***********************************************************************

  //!< define the search space for the characterisation parameters at reference temperature
  //	 Dp		diffusion constant of the positive electrode at reference temperature
  //	 Dn		diffusion constant of the negative electrode at reference temperature
  //	 kp		rate constant of the lithium insertion at the positive electrode at reference temperature
  //	 kn		rate constant of the lithium insertion at the negative electrode at reference temperature
  //	 Rdc	DC resistance of the cell

  bool logDstep = true; //!< if true, the steps in the first level are logarithmically, i.e. D = Dmin * Dstep^i
                        //!< if false, the steps are linearly, i.e. D = Dmin + i*Dstep
                        //!< steps in the later search levels are always linearly
  bool logkstep = true; //!< if true, the steps in the first level are logarithmically, i.e. k = kmin * kstep^i
                        //!< if false, the steps are linearly, i.e. k = kmin + i*kstep
                        //!< steps in the later search levels are always linearly

  auto Dp_space = slide::logstep_fix(1e-18, 13, 10); //!< lowest value/step size/Nsteps in the search space for Dp [m s-1]
  auto Dn_space = slide::logstep_fix(1e-18, 13, 10); //!< lowest value/step size/Nsteps in the search space for Dn [m s-1]

  auto kp_space = slide::logstep_fix(1e-18, 13, 10); //!< lowest value/step size/Nsteps in the search space for kp [m s-1]
  auto kn_space = slide::logstep_fix(1e-18, 13, 10); //!< lowest value/step size/Nsteps in the search space for kn [m s-1]

  auto r_space = slide::linstep_fix(1e-6, 5e-3, 9); //!< lowest value/step size/Nsteps in the search space for Rdc [Ohm]

  //!< Read the measured voltage profiles
  bool checkRange = false; //!< the first column has to start at 0, but not end at 1

  std::vector<slide::XYdata_vv> Vdata_all(nCCCV);

  for (auto &Vdata_i : Vdata_all)
    Vdata_i.reserve(100);

  for (size_t i = 0; i < nCCCV; i++) { //!< loop to read the voltage profile of each cycle

    slide::loadCSV_2col(PathVar::data / names[i], Vdata_all[i].x, Vdata_all[i].y); //!< read the csv file with the voltage profile

    //!< check the data is in the correct format
    const bool val = validOCV(checkRange, Vdata_all[i]); //!< boolean indicating if the data is in the correct format
    if (!val) {
      std::cerr << "ERROR in determineCharacterisation::fitDiffusionAndRate. Input file "
                << names[i] << " has the wrong format. throwing an error.\n";
      throw 10000;
    }
  }

  //!< ***************************************************** 3 Fit the parameters ***********************************************************************

  //!< Call the hierarchical search algorithm, which does the fitting
  int hmax = 3;              //!< number of hierarchical levels to use. Increasing this number will improve the accuracy, but take longer to calculate
  double err;                //!< error in the best fit
  std::array<double, 5> par; //!< parameters giving the lowest error [R Dp Dn kp kn]
  hierarchicalCharacterisationFit(hmax, r_space, Dp_space, Dn_space, kp_space, kn_space, Vdata_all, weights, Crates, Ccuts, Tref, ocvfit, &err, par);

  //!< ***************************************************** 4 write outputs ***********************************************************************

  //!< Print the best fit, and write in a CSV file
  auto [Rdc, Dp, Dn, kp, kn] = par;
  std::cout << "The best fit is: Rdc = " << Rdc << ", Dp = " << Dp << ", Dn = " << Dn << ", kp = " << kp << ", kn = " << kn << ".\n";

  std::ofstream output; //!< write the parameters in a csv file
  output.open(PathVar::results / nameparam, std::ios_base::out);
  output << "Rdc" << ',' << Rdc << '\n';
  output << "Dp" << ',' << Dp << '\n';
  output << "Dn" << ',' << Dn << '\n';
  output << "kp" << ',' << kp << '\n';
  output << "kn" << ',' << kn << '\n';
  output << "temperature" << ',' << Tref << '\n';
  output.close();

  //!< Simulate the voltage curves at the best fit
  auto M_ptr = Model_SPM::makeModel(); //!< structure with the matrices for the spatial discretisation for the solid diffusion PDE #TODO possible duplication somewhere of below code.

  constexpr double dt = 2;                      //!< time step to be used for the simulation
  constexpr int nin = 1 / 0.5 * 3600 / dt * 10; //!< length of the arrays to store the simulation results. (for a 0.5 C rate CC phase, *10 for the CV phase and safety margin)

  slide::XYdata_vv Vsim, Tsim; //!< arrays to store the simulation results #TODO
  Vsim.reserve(nin), Tsim.reserve(nin);

  for (int i = 0; i < nCCCV; i++) { //!< loop through all CCCV experiments
    Vsim.clear(), Tsim.clear();
    CCCV(Crates[i], Ccuts[i], Tref, Dp, Dn, kp, kn, Rdc, ocvfit, M_ptr, Vsim, Tsim);

    //!< write the simulated voltages in a csv file
    output.open(PathVar::results / (nameCCCVfit + names[i]));
    for (size_t i = 0; i < Vsim.size(); i++)
      output << Vsim.x[i] << ',' << Vsim.y[i] << ',' << Tsim.y[i] << '\n';
    output.close();
  }

  //!< Calculate the error of the best fit
  double errcomb = 0;                                            //!< the combined error of all CCCV experiments for this combination of Dp, Dn, kp and kn
  std::vector<double> erri(nCCCV);                               //!< the error for this CCCV cycle
  output.open(PathVar::results / nameparam, std::ios_base::app); //!< append the errors in the file with the parameter values

  for (int i = 0; i < nCCCV; i++) { //!< loop through all CCCV experiments
    //!< Simulate this CCCV experiment
    Vsim.clear(), Tsim.clear();
    CCCV(Crates[i], Ccuts[i], Tref, Dp, Dn, kp, kn, Rdc, ocvfit, M, Vsim, Tsim);

    //!< calculate the error (function defined in determineOCV.cpp)
    erri[i] = calculateError(false, Vdata_all[i], Vsim);
    errcomb += std::abs(erri[i]) * weights[i];
  }
  output << "combined RMSE" << ',' << errcomb << '\n';

  //!< Write on which cycles this fit is based
  output << "\n this was for the following cycles\n"
         << "name of the data file" << ',' << "length of the data file" << ','
         << "C rate of the CC phase" << ',' << "C rate of the limit current for the CV phase" << ','
         << "RMSE" << ',' << "weight\n";

  for (int i = 0; i < nCCCV; i++) {
    output << names[i] << ',' << Vdata_all[i].size() << ',' << Crates[i] << ','
           << Ccuts[i] << ',' << erri[i] << ',' << weights[i] << '\n';
  }
  output << '\n';

  //!< then note down the settings used to get this fit
  output << '\n'
         << "Below are the settings which produced this result\n"
         << "temperature" << ',' << Tref << '\n'
         << "electrode surface" << ',' << ocvfit.elec_surf << '\n'
         << "cathode active volume fraction" << ',' << ocvfit.ep << '\n'
         << "anode active volume fraction" << ',' << ocvfit.en << '\n'
         << "cathode thickness" << ',' << ocvfit.thickp << '\n'
         << "anode thickness" << ',' << ocvfit.thickn << '\n'
         << "cathode initial li-fraction at 50% SOC" << ',' << ocvfit.lifracpini << '\n'
         << "anode initial li-fraction at 50% SOC" << ',' << ocvfit.lifracnini << '\n'
         << "maximum voltage" << ',' << ocvfit.Vmax << '\n'
         << "minimum voltage" << ',' << ocvfit.Vmin << '\n'
         << "maximum li-concentration in the cathode" << ',' << ocvfit.cmaxp << '\n'
         << "maximum li-concentration in the anode" << ',' << ocvfit.cmaxn << '\n'
         << "name of the file with the cathode OCV curve" << ',' << ocvfit.namepos << '\n'
         << "name of the file with the anode OCV curve" << ',' << ocvfit.nameneg << '\n';

  //!< the search space
  output << '\n'
         << "Below are the settings of the initial search space" << '\n'
         << "number of steps in the search space for the diffusion constant" << ',' << Dp_space.size() << '\n'
         << "minimum cathode diffusion constant" << ',' << Dp_space.front() << '\n'
         << "minimum anode diffusion constant" << ',' << Dn_space.front() << '\n'
         << "initial step size for the cathode diffusion constant" << ',' << Dp_space.dstep() << '\n'
         << "initial step size for the anode diffusion constant" << ',' << Dn_space.dstep() << '\n'
         << "logarithmic steps for the diffusion constant?" << ',' << logDstep << '\n'
         << "number of steps in the search space for the rate constant" << ',' << kp_space.size() << '\n'
         << "minimum cathode rate constant" << ',' << kp_space.front() << '\n'
         << "minimum anode rate constant" << ',' << kn_space.front() << '\n'
         << "initial step size for the cathode rate constant" << ',' << kp_space.dstep() << '\n'
         << "initial step size for the anode rate constant" << ',' << kn_space.dstep() << '\n'
         << "logarithmic steps for the rate constant?" << ',' << logkstep << '\n'
         << "number of steps in the search space for the DC resistance" << ',' << r_space.size() << '\n'
         << "minimum DC resistance" << ',' << r_space.front() << '\n'
         << "initial step size in the search for the DC constant" << ',' << r_space.dstep() << '\n'
         << "number of levels in the search hierarchy" << ',' << hmax << '\n';

  output.close();
}

void writeCharacterisationParam(int h, const std::array<double, 5> &par, double err)
{
  //!< Print the best fit, and write in a CSV file
  const auto [Rdc, Dp, Dn, kp, kn] = par;

  std::cout << "The best fit in hierarchy " << h << " is: Rdc = " << Rdc << ", Dp = "
            << Dp << ", Dn = " << Dn << ", kp = " << kp << ", kn = " << kn << ".\n";
  std::ofstream output;
  //!< write the parameters and the magnitude of the error
  const auto na = PathVar::results / ("characterisationFit_" + std::to_string(h) + "_param.csv");
  output.open(na, std::ios_base::out);
  output << "Rdc" << ',' << Rdc << '\n'
         << "Dp" << ',' << Dp << '\n'
         << "Dn" << ',' << Dn << '\n'
         << "kp" << ',' << kp << '\n'
         << "kn" << ',' << kn << '\n'
         << "total RMSE" << ',' << err << '\n';
  output.close();
}
} // namespace slide