/*
 * CyclerOld.hpp
 *
 * Header for the class implementing a cycler (check-up and degradation procedures)
 * A CyclerOld consists of the programs one would write on battery cycles.
 * I.e. it defines the procedures for loading a cell in a degradation experiment.
 *
 * It also defines a struct which describes the settings for the check-up procedure
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

#include "BasicCycler.hpp"
#include "../utility/utility.hpp"

//!< Define a structure which outlines the check-up procedure.
//!< A check-up can consist of 4 things:
//!<		capacity measurement
//!<		half-cell OCV measurement
//!<		CCCV (or CC) cycles
//!<		a pulse discharge
struct checkUpProcedure
{
  bool blockDegradation;   //!< boolean indicating if degradation is accounted for during the check-up, [RECOMMENDED: TRUE]
  bool capCheck;           //!< boolean indicating if the capacity should be checked
  bool OCVCheck;           //!< boolean indicating if the half-cell OCV curves should be checked
  bool CCCVCheck;          //!< boolean indicating if some CCCV cycles should be done as part of the check-up procedure
  bool pulseCheck;         //!< boolean indicating if a pulse discharge test should be done as part of the check-up procedure
  bool includeCycleData;   //!< boolean indicating if the cycling data from the check-up should be included in the cycling data of the cell or not
  int nCycles;             //!< number of different cycles to be simulated for the CCCV check-up (i.e. the length of the array crates)
  double Crates[100];      //!< array with the Crates of the different cycles to be checked for the CCCV check-up, must be all positive
  double Ccut_cha;         //!< C rate of the cutoff current for the CV phase for the charges in the CCCV check-up, must be positive
  double Ccut_dis;         //!< C rate of the cutoff current for the CV phase for the discharges in the CCCV check-up, must be positive
  std::string profileName; //!< name of the csv file which contains the current profile for the pulse test
                           //	the first column contains the current in [A] (positive for discharge, negative for charge)
                           //	the second column contains the time in [sec] the current should be maintained
                           //	the profile must be a net discharge, i.e. sum (I*dt) > 0
  int profileLength;       //!< length of the current profiles for the pulse test (number of rows in the csv file)

  std::vector<double> I, T; //!< profile data;

  void set_profileName(const std::string &_profileName)
  {
    profileName = _profileName;

    //!< read the pulse profile
    try {
      if constexpr (settings::printBool::printCyclerHighLevel)
        std::cout << "checkUpProcedure::set_profileName is reading the current profile.\n";
      slide::loadCSV_2col(PathVar::data / profileName, I, T); //!< read the file
    } catch (int e) {
      //!< std::cout << "Throw test: " << 74 << '\n';
      if constexpr (settings::printBool::printCrit)
        std::cout << "Error in checkUpProcedure::set_profileName when reading the pulse profile: " << e << ".\n";
      throw e;
    }
  }
};

class CyclerOld : public BasicCycler
{
private:
  int indexdegr; //!< index number of the check-up (how many check-ups have we done so far)

  //!< functions for a check-up
  double getCapacity(bool blockDegradation);                                                                                                       //!< measure the remaining cell capacity
  void getOCV(slide::FixedData<double> &Ah, std::vector<double> &OCVp, std::vector<double> &OCVn);                                                 //!< measure the half-cell OCV curves
  double checkUp_batteryStates(bool blockDegradation, bool checkCap, int cumCycle, double cumTime, double cumAh, double cumWh);                    //!< measure the capacity and battery state & write to a file
  void checkUp_OCVcurves(bool blockDegradation, double ocvpini, double ocvnini);                                                                   //!< measure the half-cell OCV curves & write them to a file
  void checkUp_CCCV(bool blockDegradation, int nCycles, double Crates[], double Ccut_cha, double Ccut_dis, bool includeCycleData);                 //!< measure the voltage and temperature during some CCCV cycles & write to a file
  void checkUp_pulse(bool blockDegradation, const std::string &profileName, int profileLength, bool includeCycleData);                             //!< measure the voltage and temperature during a pulse discharge & write to a file
  void checkUp_pulse(bool blockDegradation, const std::vector<double> &I, const std::vector<double> &T, int profileLength, bool includeCycleData); //!< measure the voltage and temperature during a pulse discharge & write to a file

  double checkUp(struct checkUpProcedure &proc, int cumCycle, double cumTime, double cumAh, double cumWh); //!< function to do a check-up of a cell

public:
  CyclerOld(Cell &ci, std::string IDi, int verbosei, int feedbacki) : BasicCycler(ci, IDi, verbosei, feedbacki), indexdegr(0) {} //!< constructor

  //!< Degradation procedures
  void cycleAgeing(double dt, double Vma, double Vmi, double Ccha, bool CVcha, double Icutcha, //!< cycle ageing by continuously repeating the same cycle
                   double Cdis, bool CVdis, double Icutdis, double Ti, int nrCycles, int nrCap, struct checkUpProcedure &proc);
  void calendarAgeing(double dt, double V, double Ti, int Time, //!< calendar ageing by resting a cell
                      int timeCheck, int mode, struct checkUpProcedure &proc);
  void profileAgeing(const std::string &nameI, int limit, //!< profile ageing by repeating the same current profile
                     double Vma, double Vmi, double Ti, int nrProfiles, int nrCap, struct checkUpProcedure &proc, size_t length = 1000);
};
