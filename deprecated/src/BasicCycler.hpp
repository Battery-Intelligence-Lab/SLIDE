/*
 * BasicCycler.hpp
 *
 * Header for the class implementing a basic cycler.
 * A basic cycler simulates a battery tester (and can be programmed similarly)
 * It offers functions to load a cell with a CC and/or CV.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <memory>
#include <vector>
#include <array>
#include <filesystem>
#include <string>

#include "../cells/cells.hpp"

namespace slide {

namespace fs = std::filesystem;

struct FileStatus
{
  bool is_DegradationData_batteryState_created{ false };
  bool is_DegradationData_OCV_created{ false };
};

struct CyclerData
{
  double I;       //!< current of the cell at every step [A]
  double V;       //!< voltage of the cell at every step [V]
  double OCVp;    //!< cathode potential of the cell at every step [V]
  double OCVn;    //!< anode potential of the cell at every step [V]
  double T;       //!< temperature of the cell at every step [K]
  double timeCha; //!< cumulative time spent on charging since the start at every step [s]
  double AhCha;   //!< cumulative charged throughput since the start at every step [A]
  double WhCha;   //!< cumulative charged energy throughput since the start at every step [Wh]
  double timeDis; //!< cumulative time spent on discharging since the start at every step [s]
  double AhDis;   //!< cumulative discharged charge throughput since the start at every step [A]
  double WhDis;   //!< cumulative discharged energy throughput since the start at every step [Wh]
  double timeRes; //!< cumulative time spent on rest since the start at every step [s]
};

class BasicCycler
{

protected:
  Cell c;         //!< cell connected to the basicCycler
  std::string ID; //!< identification string for this basicCycler, will be included in the name of the subfolder in which data files are stored.
  int verbose{ settings::verbose };
  //!< integer deciding how verbose the simulation should be
  //!< The higher the number, the more output there is.
  //!< Recommended level is 1, only use higher levels for debugging
  //!< From level 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
  //!<	0 	almost no messages are printed, only in case of critical errors related to illegal parameters
  //!<	1 	error messages are printed in case of critical errors which might crash the simulation
  //!<	2 	all error messages are printed, whether the simulation can recover from the errors or not
  //!<	3 	on top of the output from 2, a message is printed every time a function in the Cycler and BasicCycler is started and terminated
  //!<	4 	on top of the output from 3, the high-level flow of the program in the Cycler is printed (e.g. 'we are going to discharge the cell')
  //!<	5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
  //!<	6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
  //!<	7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated

  //!< data of the cell which is being cycled
  int CyclingDataTimeInterval;                     //!< time resolution at which the cell data is stored [s], if 0 no data is stored
                                                   //!< the higher the time resolution, the more accurate the stored data is
  int fileIndex{ 0 };                              //!< number of data files written so far
  size_t maxLength{ 100000 };                      //!< length of the arrays in which the results are stored before they are written to a csv file, default = 100000
  double timeCha{ 0 }, timeDis{ 0 }, timeRes{ 0 }; //!< cumulative time spent on charging/discharging/rest since the start of this data collection [s]
  double AhCha{ 0 }, AhDis{ 0 };                   //!< cumulative charged/discharged throughput since the start of this data collection [A]
  double WhCha{ 0 }, WhDis{ 0 };                   //!< cumulative charged/discharged energy throughput since the start of this data collection [Wh]

  std::vector<CyclerData> outVec; //!< Output vector.

  FileStatus fileStatus;

  void storeResults(double I, double v, double ocvp, double ocvn, double tem); //!< store the cycling data of a cell
  int setCurrent(double I, double Vupp, double Vlow);                          //!< auxiliary function of CC_t_V to set the current
  void findCVcurrent_recursive(double Imin, double Imax, int sgn, double Vset, double dt, bool blockDegradation, double *Il, double *Vl);
  void findCVcurrent_recursive_bin(double Imin, double Imax, unsigned short Nrecursion, double Vset, double dt, bool blockDegradation, const slide::State &s_ini, const double I_ini, double *Il, double *Vl); //!< binary search version

  //!< auxiliary function to solve the nonlinear equation to keep the voltage constant

public:
  BasicCycler(Cell &ci, std::string IDi, int verbose, int CyclingDataTimeIntervali);

  Cell &getCell() { return c; }                          //!< returns (a reference to) the cell of the basicCycler
  void setCyclingDataTimeResolution(int timeResolution); //!< change the time resolution of the data collection

  void clearData(); //!< Clear the recorded data.
  void reset();

  //!< Functions to write the cycling data of the cell
  void openFolder(int CyclingDataTimeIntervali);

  void writeCyclingData();                                                                         //!< writes the data of the cell to a csv file with the standard name
  void writeCyclingData(const std::string &name, bool clear);                                      //!< writes the data of the cell to the specified csv file
  void returnCyclingData(std::vector<double> &Ah, std::vector<double> &V, std::vector<double> &T); //!< return the voltage of the cell instead of writing it to the file

  //!< cycle battery at constant current
  int CC_t_V(double I, double dt, bool blockDegradation, double time, double Vupp, double Vlow, double *ahi, double *whi, double *timei); //!< CC cycle with a time and two voltage end-condition
  int CC_t(double I, double dt, bool blockDegradation, double time, double *ahi, double *whi, double *timei);                             //!< CC cycle for a fixed amount of time
  int CC_V(double I, double dt, bool blockDegradation, double Vset, double *ahi, double *whi, double *timei);                             //!< CC cycle until a given voltage is reached
  void CC_halfCell_full(double I, double dt, bool pos, std::vector<double> &OCVi, double *ahi, bool isWritten);                           //!< CC cycle only one electrode

  //!< cycle battery at constant voltage
  void findCVcurrent(double Vset, double dt, bool blockDegradation, double *Il, double *Vl);     //!< find the current needed to keep the voltage constant at the specified value
  void findCVcurrent_bin(double Vset, double dt, bool blockDegradation, double *Il, double *Vl); //!< findCVcurrent with binary search

  int CV_t_I(double V, double dt, bool blockDegradation, double time, double Icut, double *ahi, double *whi, double *timei); //!< CV cycle with both a time and current limit
  void CV_t(double V, double dt, bool blockDegradation, double time, double *ahi, double *whi, double *timei);               //!< CV cycle for a fixed amount of time
  void CV_I(double V, double dt, bool blockDegradation, double Icut, double *ahi, double *whi, double *timei);               //!< CV cycle until the current is below a threshold value

  //!< hybrid CC CV
  int CC_t_CV_t(double I, double dt, bool blockDegradation, double time, double Vupp, double Vlow, double *ahi, double *whi, double *timei); //!< load the cell for a given time, with a CC and switch automatically to a CV if a voltage limit is reached
  void CC_V_CV_I(double Crate, double Vset, double Icut, double dt, bool blockDegradation, double *ahi, double *whi, double *timei);         //!< bring the cell to a specified voltage, starting with a CC phase, followed by a CV phase

  //!< current profile
  int followI(int nI, const std::string &nameI, bool blockDegradation, int limit, double Vupp, double Vlow, double *ahi, double *whi, double *timei);                                   //!< follow a predefined current pattern (reads CSV)
  int followI(int nI, const std::vector<double> &I, const std::vector<double> &T, bool blockDegradation, int limit, double Vupp, double Vlow, double *ahi, double *whi, double *timei); //!< follow a predefined current pattern (does not read CSV)
};
} // namespace slide