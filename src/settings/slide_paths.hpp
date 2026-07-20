/**
 * @file slide_paths.hpp
 * @brief Defines paths to be used in the program.
 * @author Volkan Kumtepeli
 * @date 3 Mar 2020
 */

#pragma once

#include <string>
#include <filesystem>

inline std::filesystem::path operator+(const std::filesystem::path &lhs, const std::string &rhs)
{ //!< To make path type compatible with strings.
  const std::filesystem::path temp{ rhs };
  return lhs / temp;
}

namespace PathVar {
namespace fs = std::filesystem;

//!< `inline` (not `static`): one object for the whole program. As `static` these
//!< were a separate dynamically-initialised copy in every translation unit that
//!< included this header, and `results`/`data` were additionally non-const, so a
//!< write would have silently changed only the writing TU's view. Nothing assigns
//!< to them, so `inline const` is the honest declaration.
#ifdef SLIDE_ROOT_DIR //!< SLIDE_CMAKE_MACROS
//!< If path macros are defined in CMake, then use them.
inline const auto root_folder = fs::path(SLIDE_ROOT_DIR).make_preferred();
#else
//!< Fallback for builds that do not define the macro: resolved against the
//!< *current working directory*, so a binary run from anywhere else cannot find
//!< `data/`. The build system defines SLIDE_ROOT_DIR precisely to avoid this.
inline const fs::path root_folder{ "../.." };
#endif

inline const fs::path results_folder = "results";
inline const fs::path data_folder = "data";

inline const fs::path results = root_folder / results_folder;
inline const fs::path data = root_folder / data_folder;
} // namespace PathVar

namespace slide::settings::path::Kokam {
const static std::string namepos{ "Kokam_OCV_NMC.csv" };
const static std::string nameneg{ "Kokam_OCV_C.csv" };
const static std::string nameentropicC{ "Kokam_entropic_C.csv" };
const static std::string nameentropicCell{ "Kokam_entropic_cell.csv" };

} // namespace slide::settings::path::Kokam