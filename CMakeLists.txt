cmake_minimum_required (VERSION 3.6)

project (SLIDE VERSION 2.0)

set (CMAKE_DEBUG_POSTFIX d)

# ensure C++17
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

#export compiler flags for code completion engines
set (CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/Debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_CURRENT_SOURCE_DIR}/RelWithDebInfo)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${CMAKE_CURRENT_SOURCE_DIR}/MinSizeRel)


# Set global compiler warnings
if (MSVC)
    set(common_options /W3 /permissive-) #/WX
    set(debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)
    set(not_debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)

    add_compile_options(
        "$<$<CONFIG:Debug>:${debug_options}>$<$<NOT:$<CONFIG:Debug>>:${not_debug_options}>"
        )
else ()
	add_compile_options ("-Wall" 
                       "$<$<CONFIG:Release>:-Ofast;-march=native>") #-Ofast -march=native -g -fno-omit-frame-pointer -gdwarf-2 (flto not good)
endif ()

message(STATUS "Host system: ${CMAKE_HOST_SYSTEM}")
add_definitions(-DSLIDE_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")


set(slide_headers
  src/cell_KokamNMC.hpp   
  src/determine_OCV.h
  src/read_CSVfiles.h   
  src/basic_cycler.hpp
  src/cycling.h 
  src/cell_LGChemNMC.hpp
  src/interpolation.h 
  src/state.hpp        
  src/degradation.h
  src/cell_fit.hpp
  src/cell_user.hpp                           
  src/cell.hpp
  src/determine_characterisation.h
  src/model.h
  src/cycler.hpp
  src/constants.hpp
  src/util_error.hpp
  src/util.hpp
  src/slide_aux.hpp
  )

set (slide_source
  src/basic_cycler.cpp
  src/cell_LGChemNMC.cpp
  src/degradation.cpp
  src/cell.cpp
  src/cell_user.cpp
  src/determine_characterisation.cpp
  src/model.cpp
  src/cell_fit.cpp
  src/cycler.cpp 
  src/determine_OCV.cpp
  src/read_CSVfiles.cpp
  src/cell_KokamNMC.cpp
  src/cycling.cpp
  src/interpolation.cpp
  src/state.cpp
  src/util_error.cpp
  src/util.cpp
  )


add_executable (slide src/main.cpp
    ${slide_source}
    ${slide_headers}
    )

target_include_directories (slide PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(slide) # "-Wl,--stack,8000000000" -> No need anymore. % pthread

# 1F15CC8FAF2E004105282ADDDC78DC21098DFF10