option (SLIDE_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" ON)

# Setup coverage testing for GCC or Clang
if (SLIDE_ENABLE_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        message (STATUS "Configuring with coverage")
        target_compile_options (unit_test_Cell_Bucket PUBLIC --coverage -O0)
        target_link_libraries (unit_test_Cell_Bucket PUBLIC --coverage)

        target_compile_options (unit_test_Cell_ECM PUBLIC --coverage -O0)
        target_link_libraries (unit_test_Cell_ECM PUBLIC --coverage)

        target_compile_options (unit_test_Cell_SPM PUBLIC --coverage -O0)
        target_link_libraries (unit_test_Cell_SPM PUBLIC --coverage)

        target_compile_options (unit_test_Converter PUBLIC --coverage -O0)
        target_link_libraries (unit_test_Converter PUBLIC --coverage)

    else ()
        message (FATAL_ERROR "GCC or Clang required with SLIDE_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()