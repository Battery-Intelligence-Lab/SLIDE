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
  const std::filesystem::path temp{ rhs };
  return lhs / temp;
}

namespace PathVar {
namespace fs = std::filesystem;

#ifdef SLIDE_ROOT_DIR // SLIDE_CMAKE_MACROS
// If path macros are defined in CMake, then use them.
const auto root_folder = fs::path(SLIDE_ROOT_DIR).make_preferred();
#else
const static fs::path root_folder{ "../.." };
#endif

const static fs::path results_folder = "results";
const static fs::path data_folder = "data";

static fs::path results = root_folder / results_folder;
static fs::path data = root_folder / data_folder;
} // namespace PathVar

namespace slide::settings::path::Kokam {
const static std::string namepos{ "Kokam_OCV_NMC.csv" };
const static std::string nameneg{ "Kokam_OCV_C.csv" };
const static std::string nameentropicC{ "Kokam_entropic_C.csv" };
const static std::string nameentropicCell{ "Kokam_entropic_cell.csv" };

} // namespace slide::settings::path::Kokam