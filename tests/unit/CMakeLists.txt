# Add executables and tests with coverage flags applied
add_executable_with_coverage_and_test(unit_test_Cell_Bucket Cell_Bucket_test.cpp)
add_executable_with_coverage_and_test(unit_test_Cell_ECM Cell_ECM_test.cpp)
add_executable_with_coverage_and_test(unit_test_Cell_SPM Cell_SPM_test.cpp)
add_executable_with_coverage_and_test(unit_test_Chebyshev Chebyshev_test.cpp)

add_executable_with_coverage_and_test(unit_test_Converter Converter_test.cpp)
add_executable_with_coverage_and_test(unit_test_Module_p Module_p_test.cpp)
add_executable_with_coverage_and_test(unit_test_Module_p_phase0 Module_p_phase0_test.cpp)

# Phase 0 (agent B) regression test for bug B1
add_executable_with_coverage_and_test(unit_test_CellDataStorage CellDataStorage_test.cpp)
# Phase 0 (agent B) regression test for bug B2
add_executable_with_coverage_and_test(unit_test_Histogram Histogram_test.cpp)
# Phase 0 (agent B) regression test for bug B3
add_executable_with_coverage_and_test(unit_test_Cycler_energy Cycler_energy_test.cpp)

# Phase 1 (v4 core): StateArena/BatchBuilder (PLAN.md §3.1) — links slide_core, not legacy physics
add_executable_with_coverage_and_test(unit_test_core_StateArena core_StateArena_test.cpp)
target_link_libraries(unit_test_core_StateArena PRIVATE slide_core)

# Phase 1: cold physical description hierarchy + canonical v4 domain mapping (PLAN.md §3.3)
add_executable_with_coverage_and_test(unit_test_core_CellDesign core_CellDesign_test.cpp)
target_link_libraries(unit_test_core_CellDesign PRIVATE slide_core)

# Phase 1 P1-G0 pilot (Q8): legacy vs v4 core diffusion parity-drift — needs BOTH legacy src and slide_core
add_executable_with_coverage_and_test(unit_test_core_P1G0_pilot core_P1G0_pilot_test.cpp)
target_link_libraries(unit_test_core_P1G0_pilot PRIVATE slide_core)

# Phase 1: production SpectralDiffusion<NCH> vs legacy-shaped oracle (PLAN.md §3.2/§3.5) — needs legacy src (Model_SPM) + slide_core
add_executable_with_coverage_and_test(unit_test_core_SpectralDiffusion core_SpectralDiffusion_test.cpp)
target_link_libraries(unit_test_core_SpectralDiffusion PRIVATE slide_core)

# Phase 1 P1-G3 (part a): analytic Chebyshev eigenvalue oracle (roots of tan(mu)=mu) — needs Model_SPM
add_executable_with_coverage_and_test(unit_test_core_ChebyshevEigenvalues core_ChebyshevEigenvalues_test.cpp)

# Phase 1 P1-G3 (part b): constant-flux spherical transient incl. centre-node output oracle
add_executable_with_coverage_and_test(unit_test_core_ChebyshevTransient core_ChebyshevTransient_test.cpp)
target_link_libraries(unit_test_core_ChebyshevTransient PRIVATE slide_core)

# Phase 1: observable-reconstruction layer — surface concentration vs legacy Cell_SPM + analytic round-trip
add_executable_with_coverage_and_test(unit_test_core_SpmObservables core_SpmObservables_test.cpp)
target_link_libraries(unit_test_core_SpmObservables PRIVATE slide_core)

# add_executable_with_coverage_and_test(unit_test_Module_s Module_s_test.cpp)
# add_executable_with_coverage_and_test(unit_test_Cycler Cycler_test.cpp)
# add_executable_with_coverage_and_test(unit_test_Procedure Procedure_test.cpp)
