/*
 * derived_settings.hpp
 *
 * These settings are derived settings based on the user settings;
 * Please do not change this file.
 *
 *  Created on: 11 Oct 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

namespace slide::settings {
constexpr size_t CELL_NDATA_MAX = DATASTORE_CELL <= cellDataStorageLevel::storeHistogramData ? CELL_NDATA_HIST_MAX : CELL_NDATA_INST_MAX;
//!< length of arrays in which we store cycling data at every time step
} // namespace slide::settings

namespace slide::settings::data {
constexpr bool storeCumulativeData = (DATASTORE_CELL >= cellDataStorageLevel::storeCumulativeData);
constexpr bool writeCumulativeData = storeCumulativeData;

constexpr size_t N_CumulativeCell = storeCumulativeData ? 3 : 0;
constexpr size_t N_CumulativeModule = (DATASTORE_MODULE >= moduleDataStorageLevel::storeCumulativeData) ? 3 : 0;


} // namespace slide::settings::data


namespace slide::settings::printBool {
constexpr auto printCrit = settings::verbose >= printLevel::printCrit;                               //!< threshold of verbose of when to print error messages for critical errors
constexpr auto printNonCrit = settings::verbose >= printLevel::printNonCrit;                         //!< threshold of verbose of when to print error messages for noncritical errors
constexpr auto printCyclerFunctions = settings::verbose >= printLevel::printCyclerFunctions;         //!< threshold of verbose of when to print the start and end of functions of the BasicCycler
constexpr auto printCyclerHighLevel = settings::verbose >= printLevel::printCyclerHighLevel;         //!< threshold of verbose of when to print the high-level flow of the program in the BasicCycler
constexpr auto printCyclerDetail = settings::verbose >= printLevel::printCyclerDetail;               //!< threshold of verbose of when to print the low-level detailed flow of the program in the BasicCycler
constexpr auto printfindCVcurrentDetail = settings::verbose >= printLevel::printfindCVcurrentDetail; //!< threshold of verbose of when to print the details of how the current for the CV phase is found
} // namespace slide::settings::printBool
