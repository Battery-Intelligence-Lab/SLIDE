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

# Phase 1: parameter-function canonicalisation to exact indexed curves / uniform LUTs (§3.11)
add_executable_with_coverage_and_test(unit_test_core_CompiledCurve core_CompiledCurve_test.cpp)
target_link_libraries(unit_test_core_CompiledCurve PRIVATE slide_core)

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

# Phase 1: shared kinetics/OCV/resistance/voltage/heat observable stage vs legacy Cell_SPM
add_executable_with_coverage_and_test(unit_test_core_SpmElectrical core_SpmElectrical_test.cpp)
target_link_libraries(unit_test_core_SpmElectrical PRIVATE slide_core)

# Phase 1: allocation-free lumped thermal RHS over shared heat observables + q_ext seam
add_executable_with_coverage_and_test(unit_test_core_ThermalLumped core_ThermalLumped_test.cpp)
target_link_libraries(unit_test_core_ThermalLumped PRIVATE slide_core)

# Phase 1 ageing: SEI mechanisms 1--4, porosity coupling, and additive arena RHS
add_executable_with_coverage_and_test(unit_test_core_Sei core_Sei_test.cpp)
target_link_libraries(unit_test_core_Sei PRIVATE slide_core)

# Phase 1 ageing: shared Dai/Laresgoiti observables + restart-safe stress history
add_executable_with_coverage_and_test(unit_test_core_SpmStress core_SpmStress_test.cpp)
target_link_libraries(unit_test_core_SpmStress PRIVATE slide_core)

# Phase 1 ageing: surface cracking 1--5, extra SEI, and diffusivity coupling
add_executable_with_coverage_and_test(unit_test_core_SurfaceCrack core_SurfaceCrack_test.cpp)
target_link_libraries(unit_test_core_SurfaceCrack PRIVATE slide_core)

# Phase 1 ageing: loss-of-active-material mechanisms 1--4
add_executable_with_coverage_and_test(unit_test_core_Lam core_Lam_test.cpp)
target_link_libraries(unit_test_core_Lam PRIVATE slide_core)

# Phase 1 ageing: Yang lithium-plating side reaction and state coupling
add_executable_with_coverage_and_test(unit_test_core_LithiumPlating core_LithiumPlating_test.cpp)
target_link_libraries(unit_test_core_LithiumPlating PRIVATE slide_core)

# Phase 1: compile-time fixed SPM RHS pipeline, mandatory zeroing, and rebind-safe evaluation
add_executable_with_coverage_and_test(unit_test_core_SpmPipeline core_SpmPipeline_test.cpp)
target_link_libraries(unit_test_core_SpmPipeline PRIVATE slide_core)

# Phase 1 factory prerequisite: validated per-batch spectral model build (C1-C3)
add_executable_with_coverage_and_test(unit_test_core_SpectralModel core_SpectralModel_test.cpp)
target_link_libraries(unit_test_core_SpectralModel PRIVATE slide_core)

# Phase 1 D-02: runtime registry selects one concrete pipeline per homogeneous batch
add_executable_with_coverage_and_test(unit_test_core_SpmFactory core_SpmFactory_test.cpp)
target_link_libraries(unit_test_core_SpmFactory PRIVATE slide_core)

# Phase 1: allocation-free EulerLegacy stepping and minimal single-batch Simulation façade
add_executable_with_coverage_and_test(unit_test_core_Simulation core_Simulation_test.cpp)
target_link_libraries(unit_test_core_Simulation PRIVATE slide_core)

# Phase 1 gate P1-G2: exact allocation counter over one accepted 10,000-lane step
add_executable_with_coverage_and_test(unit_test_core_P1G2_allocation core_P1G2_allocation_test.cpp)
target_link_libraries(unit_test_core_P1G2_allocation PRIVATE slide_core)

# Phase 1 gate P1-G4: arena snapshot into a rebuilt pipeline resumes bit-for-bit
add_executable_with_coverage_and_test(unit_test_core_P1G4_restart core_P1G4_restart_test.cpp)
target_link_libraries(unit_test_core_P1G4_restart PRIVATE slide_core)

# Phase 2 foundation: value-semantic combinators, flat netlist, and D-21 thermal graph
add_executable_with_coverage_and_test(unit_test_core_PackTopology core_PackTopology_test.cpp)
target_link_libraries(unit_test_core_PackTopology PRIVATE slide_core)

# Phase 2: batch Thevenin boundary, persistent sparse workspace, and ladder Mode B
add_executable_with_coverage_and_test(unit_test_core_PackSolver core_PackSolver_test.cpp)
target_link_libraries(unit_test_core_PackSolver PRIVATE slide_core)

# Phase 2: transactional electrical/thermal pack stepping and restore invalidation
add_executable_with_coverage_and_test(unit_test_core_PackStepper core_PackStepper_test.cpp)
target_link_libraries(unit_test_core_PackStepper PRIVATE slide_core)

# Phase 2 P2-G5: Mode-B admission against Mode A at 16p/256p
add_executable_with_coverage_and_test(unit_test_core_P2G5_ModeB core_P2G5_ModeB_test.cpp)
target_link_libraries(unit_test_core_P2G5_ModeB PRIVATE slide_core)

# Phase 2 D-21/P2-G1: allocation counter over one accepted coupled pack step
add_executable_with_coverage_and_test(unit_test_core_P2G1_allocation core_P2G1_allocation_test.cpp)
target_link_libraries(unit_test_core_P2G1_allocation PRIVATE slide_core)

# Phase 3: exact modal propagator, conservation, and nch=12 stability gates
add_executable_with_coverage_and_test(unit_test_core_ExponentialModal core_ExponentialModal_test.cpp)
target_link_libraries(unit_test_core_ExponentialModal PRIVATE slide_core)

# Phase 4: waveform-relaxation/Baumgarte and 100k-cell scale gates
add_executable_with_coverage_and_test(unit_test_core_ModeC core_ModeC_test.cpp)
target_link_libraries(unit_test_core_ModeC PRIVATE slide_core)

# Phase 5: PyBaMM-style grammar, controls, exact events, and legacy Cycler parity
add_executable_with_coverage_and_test(unit_test_core_Experiment core_Experiment_test.cpp)
target_link_libraries(unit_test_core_Experiment PRIVATE slide_core)

# Phase 6: lazy snapshot recording, CSV, and hardened mmap binary I/O
add_executable_with_coverage_and_test(unit_test_core_Recorder core_Recorder_test.cpp)
target_link_libraries(unit_test_core_Recorder PRIVATE slide_core)

add_executable_with_coverage_and_test(unit_test_core_RecorderAllocation core_RecorderAllocation_test.cpp)
target_link_libraries(unit_test_core_RecorderAllocation PRIVATE slide_core)

# Phase 7: PyBaMM-named Chen2020 and BPX 1.x parameter absorption
add_executable_with_coverage_and_test(unit_test_core_ParameterSet core_ParameterSet_test.cpp)
target_link_libraries(unit_test_core_ParameterSet PRIVATE slide_core)

# Phase 7 P7-G3: hermetic behaviour parity against committed PyBaMM traces
add_executable_with_coverage_and_test(unit_test_core_P7G3_PyBaMM core_P7G3_PyBaMM_test.cpp)
target_link_libraries(unit_test_core_P7G3_PyBaMM PRIVATE slide_core)
target_compile_definitions(unit_test_core_P7G3_PyBaMM PRIVATE SLIDE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# Phase 7 P7-G2: true forward sensitivities against centered finite differences
add_executable_with_coverage_and_test(unit_test_core_ForwardSensitivity core_ForwardSensitivity_test.cpp)
target_link_libraries(unit_test_core_ForwardSensitivity PRIVATE slide_core)

# Phase 8 P8-G4: persistent dependency-free pool and measured diagnostic
add_executable_with_coverage_and_test(unit_test_core_ThreadPool core_ThreadPool_test.cpp)
target_link_libraries(unit_test_core_ThreadPool PRIVATE slide_core)

# Phase 8 P8-G3: asynchronous shuffled/compressed recording and producer allocations
add_executable_with_coverage_and_test(unit_test_core_AsyncRecorder core_AsyncRecorder_test.cpp)
target_link_libraries(unit_test_core_AsyncRecorder PRIVATE slide_core)

add_executable_with_coverage_and_test(unit_test_core_AsyncRecorderAllocation core_AsyncRecorderAllocation_test.cpp)
target_link_libraries(unit_test_core_AsyncRecorderAllocation PRIVATE slide_core)

# add_executable_with_coverage_and_test(unit_test_Module_s Module_s_test.cpp)
# add_executable_with_coverage_and_test(unit_test_Cycler Cycler_test.cpp)
# add_executable_with_coverage_and_test(unit_test_Procedure Procedure_test.cpp)
