# This cmake file is to add external dependency projects.
# Adapted from https://github.com/cpp-best-practices/cmake_template/tree/main

include(cmake/CPM.cmake)

# Range-v3 library:
CPMAddPackage(
  NAME range-v3
  URL "https://github.com/ericniebler/range-v3/archive/refs/tags/0.12.0.tar.gz"
  DOWNLOAD_ONLY YES
)

add_library(range-v3 INTERFACE)
target_include_directories(range-v3 SYSTEM INTERFACE ${range-v3_SOURCE_DIR}/include)

# Eigen library:
CPMAddPackage(
  NAME eigen
  URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
  DOWNLOAD_ONLY YES
)

add_library(eigen INTERFACE)
target_include_directories(eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})

# Catch2 library:
if(NOT TARGET Catch2::Catch2WithMain)
  CPMAddPackage(
    NAME Catch2
    URL "https://github.com/catchorg/Catch2/archive/refs/tags/v3.3.2.tar.gz"
  )
endif()

# Arrow library:
CPMAddPackage(
  NAME Arrow
  URL "https://www.apache.org/dyn/closer.lua?action=download&filename=arrow/arrow-12.0.0/apache-arrow-12.0.0.tar.gz"
  SOURCE_SUBDIR "cpp"
  OPTIONS "xsimd_SOURCE BUNDLED"
)