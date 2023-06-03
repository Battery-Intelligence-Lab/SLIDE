cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

#include(GNUInstallDirs)
#set(PACKAGE_NAME matplotlib_cpp)

# Library target
add_library(matplotlib_cpp INTERFACE)


# target_compile_features(matplotlib_cpp INTERFACE
#   cxx_std_11
# )

# TODO: Use `Development.Embed` component when requiring cmake >= 3.18
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Python3 executable: ${Python3_EXECUTABLE}")
message(STATUS "Python3 include dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
message(STATUS "Python3 NumPy include dir: ${Python3_NumPy_INCLUDE_DIRS}")
target_link_libraries(matplotlib_cpp INTERFACE
  Python3::Python
  Python3::Module
  ${Python3_LIBRARIES}
)

target_include_directories(matplotlib_cpp
  INTERFACE
  .
  ${Python3_INCLUDE_DIRS}
)

find_package(Python3 COMPONENTS NumPy)
if(Python3_NumPy_FOUND)
  target_link_libraries(matplotlib_cpp INTERFACE
    Python3::NumPy
  )
else()
  target_compile_definitions(matplotlib_cpp INTERFACE WITHOUT_NUMPY)
endif()

# if(Python3_NumPy_FOUND)
#   add_executable(surface examples/surface.cpp)
#   target_link_libraries(surface PRIVATE matplotlib_cpp)
#   set_target_properties(surface PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

#   add_executable(colorbar examples/colorbar.cpp)
#   target_link_libraries(colorbar PRIVATE matplotlib_cpp)
#   set_target_properties(colorbar PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
#   add_executable(contour examples/contour.cpp)
#   target_link_libraries(contour PRIVATE matplotlib_cpp)
#   set_target_properties(contour PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

#   add_executable(spy examples/spy.cpp)
#   target_link_libraries(spy PRIVATE matplotlib_cpp)
#   set_target_properties(spy PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
# endif()


# Install headers
# install(FILES
#   "${PROJECT_SOURCE_DIR}/matplotlibcpp.h"
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


