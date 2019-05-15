cmake_minimum_required (VERSION 3.6)

project (SLIDE)

set (CMAKE_DEBUG_POSTFIX d)

# ensure C++11
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

#export compiler flags for code completion engines
set (CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Set global compiler warnings
if (MSVC)
    set(common_options /W3 /WX)
    set(debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)
    set(not_debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)

    add_compile_options(
        "$<$<CONFIG:Debug>:${debug_options}>$<$<NOT:$<CONFIG:Debug>>:${not_debug_options}>"
        )
else ()
	add_compile_options (-Wall)
endif ()


set (slide_headers
  src/Cell_KokamNMC.hpp   
  src/determineOCV.h
  src/ReadCSVfiles.h   
  src/BasicCycler.hpp
  src/Cycling.h 
  src/Cell_LGChemNMC.hpp
  src/Interpolation.h 
  src/State.hpp        
  src/Degradation.h
  src/Cell_Fit.hpp
  src/Cell_user.hpp                           
  src/Cell.hpp
  src/determineCharacterisation.h
  src/Model.h
  src/Cycler.hpp      
  )

set (slide_source
  src/BasicCycler.cpp
  src/Cell_LGChemNMC.cpp
  src/Degradation.cpp
  src/Cell.cpp
  src/Cell_user.cpp
  src/determineCharacterisation.cpp
  src/Model.cpp
  src/Cell_Fit.cpp
  src/Cycler.cpp 
  src/determineOCV.cpp
  src/ReadCSVfiles.cpp
  src/Cell_KokamNMC.cpp
  src/Cycling.cpp
  src/Interpolation.cpp
  src/State.cpp
  )


add_executable (slide src/main.cpp
    ${slide_source}
    ${slide_headers}
    )

target_include_directories (slide PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
TARGET_LINK_LIBRARIES(slide pthread)

