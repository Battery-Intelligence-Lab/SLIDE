# Phase 1 P1-G1: full single-cell legacy/core trajectory parity.
add_executable_with_coverage_and_test(parity_P1G1_spm P1G1_spm_test.cpp)
target_link_libraries(parity_P1G1_spm PRIVATE slide_core)
