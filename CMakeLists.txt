cmake_minimum_required(VERSION 3.16)

project (SLIDE_V3 VERSION 3.0)

# defines targets and sources CMAKE_BUILD_PARALLEL_LEVEL -> Try this to see if it gets better than -j 
add_subdirectory(src)

set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON) #export compiler flags for code completion engines

# include(CheckIPOSupported)
# check_ipo_supported(RESULT ipo_supported)
# if(ipo_supported)
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION True)
# endif()

add_library(slidelib INTERFACE)

target_link_libraries(slidelib
  INTERFACE
    SLIDE
  )

add_executable(slide_v3
    src/main.cpp
    )

target_link_libraries(slide_v3
  PRIVATE
    SLIDE
  ) # "-Wl,--stack,8000000000" -> No need anymore. % pthread

# 1F15CC8FAF2E004105282ADDDC78DC21098DFF10