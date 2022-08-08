/*
 * slide_paths.hpp
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

inline std::filesystem::path operator+(const std::filesystem::path &lhs, const std::string &rhs)
{ // To make path type compatible with strings.
    const std::filesystem::path temp{rhs};
    return lhs / temp;
}

namespace PathVar
{
    namespace fs = std::filesystem;

#ifdef SLIDE_ROOT_DIR // SLIDE_CMAKE_MACROS
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

namespace settings::path::Kokam
{
    const std::string namepos{"Kokam_OCV_NMC.csv"};
    const std::string nameneg{"Kokam_OCV_C.csv"};
    const std::string nameentropicC{"Kokam_entropic_C.csv"};
    const std::string nameentropicCell{"Kokam_entropic_cell.csv"};

} // namespace path::Kokam