# Phase 1 P1-G1: full single-cell legacy/core trajectory parity.
add_executable_with_coverage_and_test(parity_P1G1_spm P1G1_spm_test.cpp)
target_link_libraries(parity_P1G1_spm PRIVATE slide_core)

# Phase 2 P2-G1: compiled 3s2p pack vs legacy nested Module_s/Module_p
add_executable_with_coverage_and_test(parity_P2G1_pack P2G1_pack_test.cpp)
target_link_libraries(parity_P2G1_pack PRIVATE slide_core)
