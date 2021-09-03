/*
 * Constants.hpp
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
#include <filesystem>

inline std::filesystem::path operator+(const std::filesystem::path &lhs, const std::string &rhs) // To make path type compatible with strings.
{
    const std::filesystem::path temp{rhs};
    return lhs / temp;
}

namespace PathVar
{
    namespace fs = std::filesystem;

#ifdef SLIDE_ROOT_DIR //SLIDE_CMAKE_MACROS
    // If path macros are defined in CMake, then use them.
    const auto root_folder = fs::path(SLIDE_ROOT_DIR).make_preferred();
#else
    const std::string root_folder{".."};
#endif

    const fs::path results_folder = "results";
    const fs::path data_folder = "data";

    const fs::path results = root_folder / results_folder;
    const fs::path data = root_folder / data_folder;
} // namespace PathVar

namespace PhyConst
{
    constexpr double Kelvin = 273;
    constexpr double F = 96487;  // Faraday's constant
    constexpr double Rg = 8.314; // ideal gas constant
}

namespace slide
{

    enum cellType
    {
        // which cell to use for the simulation.
        // 0 	high power Kokam NMC cell (18650)
        // 1 	high energy LG Chem NMC cell (18650)
        // 2 	user cell (template class provided for where the user can define his own parameters)

        KokamNMC = 0,
        LGChemNMC = 1,
        UserCell = 2
    };

}

namespace settings
{
    constexpr bool isParallel{true};                  // Parallelises the code if possible.
    constexpr unsigned int numMaxParallelWorkers = 8; // Maximum number of threads to use if isParallel true.

    // if this assertion fails, the user has changed something in the code at some point, without accounting for this change somewhere else.
    // e.g. if you add an extra state-variable, you have to increase the value of 'ns' (defined in Constants.hpp), and add it in all functions in State.
    constexpr int nch{5};
    // number of points in the spatial discretisation of the solid diffusion equation *** DON'T CHANGE THE VALUE ***
    // this is the number of positive inner Chebyshev nodes
    // 		the full Chebyshev interval is from x = -1 to x = 1
    // 		the positive points go from x = 0 to x = 1
    // 		the inner positive points are the positive points excluding the point at x=0 and at x=1
    // 		so nch is the number of Chebyshev points with 0 < x < 1
    // do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied Matlab scripts.
    // See the word document '2 overview of the code', section 'Matlab setup before running the C++ code'
    constexpr int ns{2 * nch + 14};

    constexpr double Tmin_C{0};  // the minimum temperature allowed in the simulation [oC]
    constexpr double Tmax_C{60}; // the maximum temperature allowed in the simulation [oC]

    constexpr double Tmin_K{PhyConst::Kelvin + Tmin_C}; // the minimum temperature allowed in the simulation [K]
    constexpr double Tmax_K{PhyConst::Kelvin + Tmax_C}; // the maximum temperature allowed in the simulation [K]

    constexpr bool overwrite_data = true; // if this is false then folder overwriting is forbidden so you need to delete folders in results.

    // Choose how much messages should be printed to the terminal
    constexpr int verbose{0}; // integer deciding how verbose the simulation should be
                              // The higher the number, the more output there is.
                              // Recommended value is 1, only use higher value for debugging
                              // From 4 (and above) there will be too much info printed to follow what is going on, but this might be useful for debugging to find where the error is and why it is happening
                              // 	0 	almost no messages are printed, only in case of critical errors related to illegal parameters
                              // 	1 	error messages are printed in case of critical errors which might crash the simulation
                              // 	2 	all error messages are printed, whether the simulation can recover from the errors or not
                              // 	3 	on top of the output from 2, a message is printed every time a function in the Cycler and BasicCycler is started and terminated
                              // 	4 	on top of the output from 3, the high-level flow of the program in the Cycler is printed (e.g. 'we are going to discharge the cell')
                              // 	5 	on top of the output from 4, the low-level flow of the program in the BasicCycler is printed (e.g. 'in time step 101, the voltage is 3.65V')
                              // 	6 	on top of the output from 5, we also print details of the nonlinear search for the current needed to do a CV phase
                              // 	7 	on top of the output from 6, a message is printed every time a function in the Cell is started and terminated

    namespace path::Kokam
    {
        const std::string namepos{"Kokam_OCV_NMC.csv"};
        const std::string nameneg{"Kokam_OCV_C.csv"};
        const std::string nameentropicC{"Kokam_entropic_C.csv"};
        const std::string nameentropicCell{"Kokam_entropic_cell.csv"};

    } // namespace path::Kokam

} // namespace settings

// Non-user related settings, please do not change!
namespace settings
{
    //printLevel_bool::
}

enum printLevel
{
    printCrit = 1,            // threshold of verbose of when to print error messages for critical errors
    printNonCrit,             // threshold of verbose of when to print error messages for noncritical errors
    printCyclerFunctions,     // threshold of verbose of when to print the start and end of functions of the BasicCycler
    printCyclerHighLevel,     // threshold of verbose of when to print the high-level flow of the program in the BasicCycler
    printCyclerDetail,        // threshold of verbose of when to print the low-level detailed flow of the program in the BasicCycler
    printfindCVcurrentDetail, // threshold of verbose of when to print the details of how the current for the CV phase is found
    printCellFunctions        // threshold of verbose of when to print messages at the start and end of functions of the Cell

};
