cmake_minimum_required(VERSION 3.17)

# include(cmake/ProjectSettings.cmake)
project(SLIDE VERSION "3.0.0"
  DESCRIPTION "A simulator for lithium-ion battery degradation"
  HOMEPAGE_URL https://github.com/Battery-Intelligence-Lab/SLIDE/
  LANGUAGES CXX)

 
include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)

enable_testing()
add_subdirectory(src bin)

include(cmake/Coverage.cmake)
add_subdirectory(tests)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_20)

add_executable(slide
  src/main.cpp
)

target_link_libraries(slide
  PRIVATE
  src
  project_warnings
  project_options
)

target_include_directories(slide
  PUBLIC
  data
)



# message(STATUS "Project: ${PROJECT_NAME} version ${PROJECT_HOMEPAGE_URL}")

# cmake --build ./build/ --parallel $(($(nproc)-1))
# cmake --build ./build/ --clean-first
# cmake --install build --strip
# cmake --install build --config Debug
# ccmake -G "Unix Makefiles" -S . -B ./build
if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)

  if(ENABLE_BUILD_WITH_TIME_TRACE)
    add_compile_definitions(project_options INTERFACE -ftime-trace)
  endif()
endif()

# # enable cache system
# include(cmake/Cache.cmake)

# standard compiler warnings
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# # sanitizer options if supported by compiler
# include(cmake/Sanitizers.cmake)
# enable_sanitizers(project_options)

# enable doxygen
# include(cmake/Doxygen.cmake)
# enable_doxygen(src doxygen-docs)

# # allow for static analysis options
# include(cmake/StaticAnalyzers.cmake)

# set(CMAKE_INCLUDE_CURRENT_DIR ON)
# set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

# include(GenerateExportHeader)

# generate_export_header(${PROJECT_NAME})
set(LIB_NAME slide-lib)

add_library(${LIB_NAME} INTERFACE)

target_link_libraries(${LIB_NAME}
  INTERFACE
  src
)

# # install header
# install(FILES "src/slide.hpp" DESTINATION include)
include(GNUInstallDirs)

install(TARGETS ${LIB_NAME} EXPORT ${LIB_NAME}Targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

# install(EXPORT MyProject-targets DESTINATION lib/cmake/MyProject)

# # # install the library
# # install(TARGETS MyProject
# # 	EXPORT MyProject-targets
# # 	LIBRARY DESTINATION lib
# # 	ARCHIVE DESTINATION lib
# # )

# # # install export target and config for find_package
# # install(EXPORT ${LIB_NAME}Targets DESTINATION lib/cmake/${LIB_NAME})
# # include(CMakePackageConfigHelpers)
# # configure_package_config_file(
# # 	"MyProjectConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/MyProjectConfig.cmake"
# # 	INSTALL_DESTINATION "lib/cmake/${LIB_NAME}"
# # )
# # install(FILES "${CMAKE_CURRENT_BINARY_DIR}/MyProjectConfig.cmake" DESTINATION "lib/cmake/${LIB_NAME}")

# # defines targets and sources CMAKE_BUILD_PARALLEL_LEVEL -> Try this to see if it gets better than -j
# add_subdirectory(src)

# install(TARGETS slide_v3
# CONFIGURATIONS Debug
# RUNTIME DESTINATION Debug/bin)
# install(TARGETS slide_v3
# CONFIGURATIONS Release
# RUNTIME DESTINATION Release/bin)