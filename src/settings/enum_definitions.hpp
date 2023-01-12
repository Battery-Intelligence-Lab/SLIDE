/*
 * enum_definitions.hpp
 *
 *  Created on: 10 Apr 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

namespace slide::settings {
enum printLevel {
  printCrit = 1,           //!< threshold of verbose of when to print error messages for critical errors
  printNonCrit,            //!< threshold of verbose of when to print error messages for noncritical errors
  printCyclerFunctions,    //!< threshold of verbose of when to print the start and end of functions of the BasicCycler
  printCyclerHighLevel,    //!< threshold of verbose of when to print the high-level flow of the program in the BasicCycler
  printCyclerDetail,       //!< threshold of verbose of when to print the low-level detailed flow of the program in the BasicCycler
  printfindCVcurrentDetail //!< threshold of verbose of when to print the details of how the current for the CV phase is found
};

enum class cellDataStorageLevel {
  noStorage = 0,
  storeCumulativeData,
  storeHistogramData,
  storeTimeData
};

enum class moduleDataStorageLevel {
  noStorage = 0,
  storeCumulativeData,
  storeTimeData
};

enum CVcurrentAlgorithm //!< Current finding method;
{
  linearSearch = 0,
  falsePosition = 1

};

} // namespace slide::settings

namespace slide {
enum cellType {
  //!< which cell to use for the simulation.
  //!< 0 	high power Kokam NMC cell (18650)
  //!< 1 	high energy LG Chem NMC cell (18650)
  //!< 2 	user cell (template class provided for where the user can define his own parameters)
  KokamNMC = 0,
  LGChemNMC = 1,
  UserCell = 2
};

}