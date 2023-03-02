option (SLIDE_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" OFF)

# Setup coverage testing for GCC or Clang
if (SLIDE_ENABLE_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        message (STATUS "Configuring with coverage")
        target_compile_options (slide PUBLIC --coverage -O0)
        target_link_libraries (slide PUBLIC --coverage)
    else ()
        message (FATAL_ERROR "GCC or Clang required with SLIDE_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()