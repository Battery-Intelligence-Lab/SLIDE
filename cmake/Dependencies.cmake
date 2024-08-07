# This cmake file is to add external dependency projects.
# Adapted from https://github.com/cpp-best-practices/cmake_template/tree/main
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")

include(cmake/recipes/CPM.cmake)

# Range-v3 library:
CPMAddPackage(
  NAME range-v3
  URL "https://github.com/ericniebler/range-v3/archive/refs/tags/0.12.0.tar.gz"
  DOWNLOAD_ONLY YES
)

include(eigen)
include(boost)

add_library(range-v3 INTERFACE)
target_include_directories(range-v3 SYSTEM INTERFACE ${range-v3_SOURCE_DIR}/include)

# Catch2 library:
if(NOT TARGET Catch2::Catch2WithMain)
  CPMAddPackage(
    NAME Catch2
    URL "https://github.com/catchorg/Catch2/archive/refs/tags/v3.6.0.tar.gz"
  )
endif()

# fmt library:
CPMAddPackage(
  NAME fmt
  URL "https://github.com/fmtlib/fmt/archive/refs/tags/11.0.2.tar.gz"
)