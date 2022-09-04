/*
 * settings.hpp
 *
 * Author : Volkan Kumtepeli
 *
 * Defines constants to be used in the program.
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <string>

#include "enum_definitions.hpp"
#include "constants.hpp"
#include "slide_paths.hpp"
#include "tolerances.hpp"

namespace slide::settings::cool //!< Cooling System Settings.
{
constexpr double flowrate_perCell{ 0.0005 }; //!< flow rate m3/s, per cell value
}

namespace slide::settings {
constexpr bool isParallel{ true };                //!< Parallelises the code if possible.
constexpr unsigned int numMaxParallelWorkers = 8; //!< Maximum number of threads to use if isParallel true.

//!< if this assertion fails, the user has changed something in the code at some point, without accounting for this change somewhere else.
//!< e.g. if you add an extra state-variable, you have to increase the value of 'ns' (defined in Constants.hpp), and add it in all functions in State.
constexpr size_t nch{ 5 }; ////!< number of points in the spatial discretisation of the solid diffusion equation
//!< number of points in the spatial discretisation of the solid diffusion equation *** DON'T CHANGE THE VALUE ***
//!< this is the number of positive inner Chebyshev nodes
//!< 		the full Chebyshev interval is from x = -1 to x = 1
//!< 		the positive points go from x = 0 to x = 1
//!< 		the inner positive points are the positive points excluding the point at x=0 and at x=1
//!< 		so nch is the number of Chebyshev points with 0 < x < 1
//!< do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied MATLAB scripts.
//!< See the word document '2 overview of the code', section 'MATLAB setup before running the C++ code'

constexpr double Tmin_Cell_C{ 0 };  //!< the minimum temperature allowed in the simulation [oC]
constexpr double Tmax_Cell_C{ 60 }; //!< the maximum temperature allowed in the simulation [oC]

constexpr double Tmin_Cell_K{ PhyConst::Kelvin + Tmin_Cell_C }; //!< the minimum temperature allowed in the simulation [K]
constexpr double Tmax_Cell_K{ PhyConst::Kelvin + Tmax_Cell_C }; //!< the maximum temperature allowed in the simulation [K]

constexpr bool overwrite_data = true; //!< if this is false then folder overwriting is forbidden so you need to delete folders in results.

//!< Data storage  //!< from slidepack
constexpr int DATASTORE_NHIST = 100; //!< length of the arrays with the histograms (if 1)

constexpr auto DATASTORE_CELL = cellDataStorageLevel::storeCumulativeData; //!< if 0, no cell-level data is stored
                                                                           //!< if 1, statistics about I, V and T are stored, as well as overall utilisation (throughput)
                                                                           //!< if 2, current, voltage, temperature, soc is stored at every time step, as well as overall utilisation (throughput)

//!< constexpr int DATASTORE_MODULE = 0; //!< if 0, no module-level data is stored
//!< if 2, current, voltage, temperature, soc is stored at every time step, as well as overall utilisation (throughput)
//!< constexpr int DATASTORE_BATT = 0; //!< if 0, no module-level data is stored
//!< if 2, current, voltage, temperature, soc is stored at every time step, as well as overall utilisation (throughput)
constexpr int DATASTORE_COOL = 0; //!< if 0, no data is stored
                                  //!< if 1, statistics about the cooling system is stored
                                  //!< if 2, operating power etc is stored every time step

constexpr int CELL_NDEG = 10; //!< maximum number of any type of degradation model

constexpr int MODULE_NSUs_MAX = 20; //!< maximum number of cells in a base module
                                    //!< note: CELL_NSTATE_MAX * MODULE_NCELL_MAX <= StorageUnit_NSTATES_MAX

constexpr double T_ENV = 273.0 + 21.0; //!< environmental temperature

constexpr size_t CYCLER_NDATA_MAX{ 10000 }; //!< length of the arrays which hold the cycling data (if 2)
                                            //!< Large battery ~ 3000 SPM cells * (7+4) arrays * 8 Byte per double * N doubles
                                            //!< = 3000 * 88 * 10,000 = 2.64 GB
                                            //!< 1 CC cycle ~ 2 hours. So if you store data every 20s, 1 cycle gives 360 data points -> 10k is about 30 cycles

constexpr size_t CELL_NDATA_HIST_MAX = 10'000'000; //!< If histogram then write data every 10 millionth data.
constexpr size_t CELL_NDATA_INST_MAX = 10'000;     //!< If it is storing every data then we should have much less storage.

constexpr size_t CELL_NDATA_MAX = DATASTORE_CELL <= cellDataStorageLevel::storeHistogramData ? CELL_NDATA_HIST_MAX : CELL_NDATA_INST_MAX;
//!< length of arrays in which we store cycling data at every time step
constexpr size_t MODULE_NDATA_MAX{ CYCLER_NDATA_MAX }; //!< if we store cell-level data, make the array as long as the one in Cycler

constexpr size_t CELL_NSTATE_MAX{ 30 }; //!< maximum number of states of all types of cells //!< used to prepare arrays which are long enough for getStates()

//!< Data storage
//#define DATASTORE_CELL 0   //!< if 0, no cell-level data is stored
//!< if 1, statistics about I, V and T are stored, as well as overall utilisation (throughput)
//!< if 2, current, voltage, temperature, soc is stored at every time step, as well as overall utilisation (throughput)
#define DATASTORE_MODULE 0 //!< if 0, no module-level data is stored
                           //!< if 2, current, voltage, temperature, soc is stored at every time step, as well as overall utilisation (throughput)
#define DATASTORE_BATT 0   //!< if 0, no module-level data is stored
                           //!< if 2, current, voltage, temperature, soc is stored at every time step, as well as overall utilisation (throughput)

//!< timing of various functions
#define TIMING false

//!< temperature
constexpr int T_MODEL{ 0 }; //!< which thermal model to use
                            //!< 	0 no thermal model
                            //!< 	1 individual cell bulk thermal model
                            //!< 	2 coupled cell thermal model with cooling from modules
                            //!< if 1, statistics about the cooling system is stored
                            //!< if 2, operating power etc is stored every time step

//!< Cell types
//!< #define CELLTYPE_CELL 1     //!< regular cells
//!< #define CELLTYPE_CELL_ECM 2 //!< ECM cells

//!< Choose how much messages should be printed to the terminal
constexpr int verbose{ printLevel::printCrit };
//!< integer deciding how verbose the simulation should be
//!< The higher the number, the more output there is.
//!< Recommended value is 1, only use higher value for debugging
//!< From 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for
//!< debugging to find where the error is and why it is happening
//!< 	0 	almost no messages are printed, only in case of critical errors related to illegal parameters
//!< 	1 	error messages are printed in case of critical errors which might crash the simulation
//!< 	2 	all error messages are printed, whether the simulation can recover from the errors or not
//!< 	3 	on top of the output from 2, a message is printed every time a function in the Cycler and BasicCycler is started and terminated
//!< 	4 	on top of the output from 3, the high-level flow of the program in the Cycler is printed (e.g. 'we are going to discharge the cell')
//!< 	5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
//!< 	6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
//!< 	7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated

} // namespace slide::settings

//!< Non-user related settings, please do not change!
namespace slide::settings {
constexpr auto CVcurrentFindingMethod = CVcurrentAlgorithm::falsePosition;
//!< printLevel_bool::
} // namespace slide::settings


namespace slide::settings::data {
constexpr bool storeCumulativeData = (DATASTORE_CELL >= cellDataStorageLevel::storeCumulativeData);
constexpr size_t N_CumulativeData = storeCumulativeData ? 3 : 0;
} // namespace slide::settings::data


namespace slide::settings::printBool {
constexpr auto printCrit = settings::verbose >= printLevel::printCrit;                               //!< threshold of verbose of when to print error messages for critical errors
constexpr auto printNonCrit = settings::verbose >= printLevel::printNonCrit;                         //!< threshold of verbose of when to print error messages for noncritical errors
constexpr auto printCyclerFunctions = settings::verbose >= printLevel::printCyclerFunctions;         //!< threshold of verbose of when to print the start and end of functions of the BasicCycler
constexpr auto printCyclerHighLevel = settings::verbose >= printLevel::printCyclerHighLevel;         //!< threshold of verbose of when to print the high-level flow of the program in the BasicCycler
constexpr auto printCyclerDetail = settings::verbose >= printLevel::printCyclerDetail;               //!< threshold of verbose of when to print the low-level detailed flow of the program in the BasicCycler
constexpr auto printfindCVcurrentDetail = settings::verbose >= printLevel::printfindCVcurrentDetail; //!< threshold of verbose of when to print the details of how the current for the CV phase is found
constexpr auto printCellFunctions = settings::verbose >= printLevel::printCellFunctions;             //!< threshold of verbose of when to print messages at the start and end of functions of the Cell
};                                                                                                   // namespace slide::settings::printBool
