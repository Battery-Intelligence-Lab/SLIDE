cmake_minimum_required(VERSION 3.31)

# include(cmake/ProjectSettings.cmake)
project(SLIDE VERSION "3.0.0"
  DESCRIPTION "A simulator for lithium-ion battery degradation"
  HOMEPAGE_URL https://github.com/Battery-Intelligence-Lab/SLIDE/
  LANGUAGES CXX)

# The project contains no C++20 modules. Disabling dependency scanning globally
# also keeps CMake try_compile probes (for example FindThreads) independent of
# the separately packaged clang-scan-deps executable.
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

option(SLIDE_CORE_ONLY "Configure only the dependency-light v4 core" OFF)
option(SLIDE_WITH_CUDA "Enable the optional CUDA SPM batch backend" OFF)
option(SLIDE_WITH_MATLAB "Build the optional MATLAB MEX and +slide package" OFF)
option(SLIDE_BUILD_CORE_TESTS "Build the dependency-light v4 core smoke tests" OFF)
option(SLIDE_BUILD_FUZZERS "Build Clang libFuzzer parser targets" OFF)
if(SLIDE_WITH_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "SLIDE_WITH_CUDA=ON but no CUDA compiler is available")
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)
include(cmake/Dependencies.cmake) # Include external projects

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_20)

# Sanitizer lanes are homogeneous whole-project builds. Apply their directory
# options before creating any SLIDE library or executable so static core,
# legacy objects, test translation units, and final links are all instrumented.
include(cmake/Sanitizers.cmake)
enable_sanitizers()

if(SLIDE_CORE_ONLY)
  add_subdirectory(src/core bin/core)
else()
  enable_testing()
  add_subdirectory(src bin)
endif()

if(SLIDE_BUILD_CORE_TESTS)
  enable_testing()
  add_subdirectory(tests/core_only)
endif()

if(SLIDE_BUILD_FUZZERS)
  add_subdirectory(tests/fuzz)
endif()

option(SLIDE_BUILD_PYTHON "Build the nanobind Python extension" OFF)
if(SLIDE_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(SLIDE_WITH_MATLAB)
  find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM MX_LIBRARY MEX_COMPILER)
  add_subdirectory(matlab)
endif()

if(NOT SLIDE_CORE_ONLY)
  include(cmake/Coverage.cmake)
  add_subdirectory(tests)
  add_subdirectory(benchmark)

  target_link_libraries(src PUBLIC Eigen3::Eigen range-v3 fmt Boost::boost) # Link dependencies to the src library.
  # glaze::glaze  -> maybe in future.

  add_executable(slide
    src/main.cpp
  )

  target_link_libraries(slide
    PRIVATE
    src
    project_warnings
    project_options
    nlopt
  )

  target_include_directories(slide
    PUBLIC
    data
  )
endif()

# message(STATUS "Project: ${PROJECT_NAME} version ${PROJECT_HOMEPAGE_URL}")

# cmake --build ./build/ --parallel $(($(nproc)-1))
# cmake --build ./build/ --clean-first
# cmake --install build --strip
# cmake --install build --config Debug
# ccmake -G "Unix Makefiles" -S . -B ./build
if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)

  if(ENABLE_BUILD_WITH_TIME_TRACE)
    target_compile_options(project_options INTERFACE -ftime-trace)
  endif()
endif()

# enable cache system
# include(cmake/Cache.cmake)

# standard compiler warnings
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# enable doxygen
# include(cmake/Doxygen.cmake)
# enable_doxygen(src doxygen-docs)

# allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

# set(CMAKE_INCLUDE_CURRENT_DIR ON)
# set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

# include(GenerateExportHeader)

# generate_export_header(${PROJECT_NAME})
if(NOT SLIDE_CORE_ONLY)
  set(LIB_NAME slide-lib)

  add_library(${LIB_NAME} INTERFACE)

  target_link_libraries(${LIB_NAME}
    INTERFACE
    src
  )

  # # install header
  include(GNUInstallDirs)

  install(TARGETS ${LIB_NAME} EXPORT ${LIB_NAME}Targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
  )
endif()

# Prints out the licenses of third party packages in binary folder.
if(NOT SLIDE_CORE_ONLY)
  cpm_licenses_create_disclaimer_target(write-licenses "${CMAKE_CURRENT_SOURCE_DIR}/bin/third_party.txt" "${CPM_PACKAGES}")
endif()
