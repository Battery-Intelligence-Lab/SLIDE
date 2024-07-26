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

add_library(range-v3 INTERFACE)
target_include_directories(range-v3 SYSTEM INTERFACE ${range-v3_SOURCE_DIR}/include)

# Catch2 library:
if(NOT TARGET Catch2::Catch2WithMain)
  CPMAddPackage(
    NAME Catch2
    URL "https://github.com/catchorg/Catch2/archive/refs/tags/v3.3.2.tar.gz"
  )
endif()


# CPMAddPackage("gh:fmtlib/fmt#10.1.1") #fmt library
# fmt library:
CPMAddPackage(
  NAME fmt
  URL "https://github.com/fmtlib/fmt/archive/refs/tags/10.1.1.tar.gz"
)


# # Glaze library:
# CPMAddPackage(
#   NAME glaze
#   URL "https://github.com/stephenberry/glaze/archive/refs/tags/v1.3.4.tar.gz"
# )