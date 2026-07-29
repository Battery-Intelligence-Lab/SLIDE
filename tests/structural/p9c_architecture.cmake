# M0.6--M0.9 architecture aggregate. Keep one CTest slot while retaining the
# independently runnable 9C-2, 9C-3, 9C-4, and 9C-5 structural proofs.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/p9c2_ageing_kernel.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/p9c3_cold_file_split.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/p9c4_shared_test_harness.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/p9c5_public_surface.cmake")

# MQ.2 O2: finite-math builds may assume that an ordinary range comparison never
# sees NaN or infinity.  Pin the fast-math-safe classifier to the exact
# z_surface rejection slice rather than accepting an unrelated use elsewhere.
load_compact("src/core/SpmObservables.hpp" mq2_spm_observables)
require_tokens("MQ.2 fast-math surface validity include" mq2_spm_observables
  "#include\"Numeric.hpp\"")

set(mq2_surface_slice_begin
  "constRealz_surface=spm_scalar::surfaceStoichiometry(")
set(mq2_surface_slice_end "constRealreaction_rate=")
string(FIND "${mq2_spm_observables}" "${mq2_surface_slice_begin}"
  mq2_surface_slice_begin_position)
string(FIND "${mq2_spm_observables}" "${mq2_surface_slice_end}"
  mq2_surface_slice_end_position)
if(mq2_surface_slice_begin_position EQUAL -1
   OR mq2_surface_slice_end_position EQUAL -1
   OR mq2_surface_slice_end_position LESS_EQUAL mq2_surface_slice_begin_position)
  message(FATAL_ERROR
    "MQ.2 fast-math surface validity: z_surface rejection slice is missing")
endif()
math(EXPR mq2_surface_slice_length
  "${mq2_surface_slice_end_position} - ${mq2_surface_slice_begin_position}")
string(SUBSTRING "${mq2_spm_observables}"
  ${mq2_surface_slice_begin_position}
  ${mq2_surface_slice_length}
  mq2_surface_rejection_slice)
require_token_count("MQ.2 fast-math surface validity predicate"
  mq2_surface_rejection_slice
  "if(!is_finite_primal(z_surface)||!(primal_value(z_surface)>0.0&&primal_value(z_surface)<1.0))"
  1)

# MQ.2 O1: two fixed-field scratch owners share one checked cold extent,
# observable storage is consumed exactly, and cache ownership must not absorb
# any of the scalar transport physics.
load_compact("src/core/detail/CheckedLaneExtent.hpp" mq2_checked_lane_extent)
load_compact("src/core/AgeingKernel.hpp" mq2_ageing_kernel)
load_compact("src/core/Sei.hpp" mq2_sei)
require_token_count("MQ.2 checked lane extent owner"
  mq2_checked_lane_extent "checked_lane_extent(" 1)
require_token_count("MQ.2 checked ageing extent consumer"
  mq2_ageing_kernel "checked_lane_extent(" 1)
require_token_count("MQ.2 checked observable extent consumer"
  mq2_spm_observables "checked_lane_extent(" 1)
require_token_count("MQ.2 observable BatchView addressing"
  mq2_spm_observables "state.raw()" 0)
require_token_count("MQ.2 transport cache load boundary"
  mq2_spm_observables "try_load(" 2)
require_token_count("MQ.2 transport cache store boundary"
  mq2_spm_observables "store(" 2)
require_tokens("MQ.2 named observable scratch fields" mq2_spm_observables
  "staticconstexprstd::size_tper_domain_lane_fields=6;"
  "staticconstexprstd::size_tshared_lane_fields=9;"
  "assert(cursor<=storage_.size()&&count<=storage_.size()-cursor);"
  "assert(cursor==static_cast<std::size_t>(active_lanes)*per_lane);")
require_token_count("MQ.2 observable scratch bare field count"
  mq2_spm_observables "+21" 0)
require_token_count("MQ.2 SEI lane index conversion"
  mq2_sei "static_cast<std::size_t>(lane)" 3)

set(mq2_transport_cache_slice_begin "structSpmTransportCache")
set(mq2_transport_cache_slice_end
  "template<intNCH,classReal>voidcomputeSpmTransportLane(")
string(FIND "${mq2_spm_observables}" "${mq2_transport_cache_slice_begin}"
  mq2_transport_cache_slice_begin_position)
string(FIND "${mq2_spm_observables}" "${mq2_transport_cache_slice_end}"
  mq2_transport_cache_slice_end_position)
if(mq2_transport_cache_slice_begin_position EQUAL -1
   OR mq2_transport_cache_slice_end_position EQUAL -1
   OR mq2_transport_cache_slice_end_position
      LESS_EQUAL mq2_transport_cache_slice_begin_position)
  message(FATAL_ERROR "MQ.2 O1 transport-cache slice is missing")
endif()
math(EXPR mq2_transport_cache_slice_length
  "${mq2_transport_cache_slice_end_position} - ${mq2_transport_cache_slice_begin_position}")
string(SUBSTRING "${mq2_spm_observables}"
  ${mq2_transport_cache_slice_begin_position}
  ${mq2_transport_cache_slice_length}
  mq2_transport_cache_slice)
forbid_tokens("MQ.2 transport cache owns no scalar physics"
  mq2_transport_cache_slice "spm_scalar::")

set(mq2_transport_slice_begin "voidcomputeSpmTransportLane(")
set(mq2_transport_slice_end
  "template<intNCH,classReal>voidcomputeSpmTransport(")
string(FIND "${mq2_spm_observables}" "${mq2_transport_slice_begin}"
  mq2_transport_slice_begin_position)
string(FIND "${mq2_spm_observables}" "${mq2_transport_slice_end}"
  mq2_transport_slice_end_position)
if(mq2_transport_slice_begin_position EQUAL -1
   OR mq2_transport_slice_end_position EQUAL -1
   OR mq2_transport_slice_end_position LESS_EQUAL mq2_transport_slice_begin_position)
  message(FATAL_ERROR "MQ.2 O1 transport-kernel slice is missing")
endif()
math(EXPR mq2_transport_slice_length
  "${mq2_transport_slice_end_position} - ${mq2_transport_slice_begin_position}")
string(SUBSTRING "${mq2_spm_observables}"
  ${mq2_transport_slice_begin_position}
  ${mq2_transport_slice_length}
  mq2_transport_slice)
require_token_count("MQ.2 transport arrhenius ownership"
  mq2_transport_slice "spm_scalar::arrheniusFactor(" 2)
require_token_count("MQ.2 transport activation ownership"
  mq2_transport_slice "spm_scalar::activatedValue(" 2)
require_token_count("MQ.2 transport denominator ownership"
  mq2_transport_slice "spm_scalar::fluxDenominator(" 2)
require_token_count("MQ.2 transport flux ownership"
  mq2_transport_slice "spm_scalar::molarFlux(" 3)

# MQ.2 P1: the affine pack algebra has one owner per ordered operation.  Exact
# per-file counts prevent a duplicated caller from being hidden by a deleted
# caller elsewhere, while the body slices keep both sparse reconstruction paths
# independently pinned.
load_compact("src/core/PackSolverInternal.hpp" mq2_pack_internal)
load_compact("src/core/PackSolver.hpp" mq2_pack_header)
load_compact("src/core/PackSolver.cpp" mq2_pack_sparse)
load_compact("src/core/PackSolverIterative.cpp" mq2_pack_iterative)

require_token_count("MQ.2 P1 cell-current owner"
  mq2_pack_internal "cellCurrentFromDrop(" 1)
require_token_count("MQ.2 P1 cell-current sparse consumers"
  mq2_pack_sparse "cellCurrentFromDrop(" 2)
require_token_count("MQ.2 P1 cell-current iterative consumers"
  mq2_pack_iterative "cellCurrentFromDrop(" 2)
require_token_count("MQ.2 P1 branch-drop owner"
  mq2_pack_internal "branchDrop(" 1)
require_token_count("MQ.2 P1 branch-drop sparse consumers"
  mq2_pack_sparse "branchDrop(" 4)
require_token_count("MQ.2 P1 branch-drop iterative consumers"
  mq2_pack_iterative "branchDrop(" 2)
require_token_count("MQ.2 P1 branch-affine owner"
  mq2_pack_internal "branchAffine(" 1)
require_token_count("MQ.2 P1 branch-affine sparse consumers"
  mq2_pack_sparse "branchAffine(" 2)
require_token_count("MQ.2 P1 branch-affine iterative consumers"
  mq2_pack_iterative "branchAffine(" 2)
require_token_count("MQ.2 P1 branch-numerator owner"
  mq2_pack_internal "branchCurrentNumerator(" 1)
require_token_count("MQ.2 P1 branch-numerator sparse consumer"
  mq2_pack_sparse "branchCurrentNumerator(" 1)
require_token_count("MQ.2 P1 branch-numerator iterative consumer"
  mq2_pack_iterative "branchCurrentNumerator(" 1)
require_token_count("MQ.2 P1 branch-current owner"
  mq2_pack_internal "branchCurrentOut(" 1)
require_token_count("MQ.2 P1 branch-current sparse consumer"
  mq2_pack_sparse "branchCurrentOut(" 1)
require_token_count("MQ.2 P1 branch-current iterative consumer"
  mq2_pack_iterative "branchCurrentOut(" 1)
require_token_count("MQ.2 P1 exact branch-drop owner"
  mq2_pack_internal
  "[[nodiscard]]inlinereal_tbranchDrop(constCompiledElectricalBranch&branch,std::span<constreal_t>node_voltage)noexcept{returnnode_voltage[branch.node_positive]-node_voltage[branch.node_negative];}"
  1)
require_token_count("MQ.2 P1 exact branch-affine record"
  mq2_pack_internal
  "structBranchAffine{real_tresistance{};real_tsource{};};"
  1)
require_token_count("MQ.2 P1 exact branch-affine owner"
  mq2_pack_internal
  "[[nodiscard]]inlineBranchAffinebranchAffine(constCompiledElectricalBranch&branch,std::span<constreal_t>cell_ocv,std::span<constreal_t>cell_resistance)noexcept{if(branch.kind==ElectricalBranchKind::cell)return{.resistance=cell_resistance[branch.cell],.source=cell_ocv[branch.cell]};return{.resistance=branch.resistance,.source=real_t{}};}"
  1)
require_token_count("MQ.2 P1 exact branch-numerator owner"
  mq2_pack_internal
  "[[nodiscard]]inlinereal_tbranchCurrentNumerator(constreal_t&drop,constreal_t&source)noexcept{returndrop-source;}"
  1)
require_token_count("MQ.2 P1 exact branch-current owner"
  mq2_pack_internal
  "[[nodiscard]]inlinereal_tbranchCurrentOut(constreal_t&numerator,constreal_t&resistance)noexcept{returnnumerator/resistance;}"
  1)
require_token_count("MQ.2 P1 exact cell-current owner"
  mq2_pack_internal
  "[[nodiscard]]inlineboolcellCurrentFromDrop(constreal_t&drop,constreal_t&ocv,constreal_t&resistance,real_t&current)noexcept{constreal_tnumerator=ocv-drop;if(!is_finite(drop)||!is_finite(numerator))returnfalse;current=numerator/resistance;returnis_finite(current);}"
  1)
require_token_count("MQ.2 P1 sparse affine argument order"
  mq2_pack_sparse "branchAffine(branch,ocv_,resistance_)" 2)
require_token_count("MQ.2 P1 iterative affine argument order"
  mq2_pack_iterative "branchAffine(branch,ocv_,resistance_)" 2)
require_token_count("MQ.2 P1 sparse cell-current argument order"
  mq2_pack_sparse
  "cellCurrentFromDrop(voltage,ocv_[branch.cell],resistance_[branch.cell],candidate_current_[branch.cell])"
  2)
require_token_count("MQ.2 P1 iterative node-current argument order"
  mq2_pack_iterative
  "cellCurrentFromDrop(voltage,ocv_[branch.cell],resistance_[branch.cell],candidate_current_[branch.cell])"
  1)
require_token_count("MQ.2 P1 iterative ladder-current argument order"
  mq2_pack_iterative
  "cellCurrentFromDrop(layer_voltage_[layer],ocv_[cell],resistance_[cell],candidate_current_[cell])"
  1)

foreach(mq2_pack_source IN ITEMS mq2_pack_sparse mq2_pack_iterative)
  require_token_count("MQ.2 P1 raw branch drop removed"
    ${mq2_pack_source}
    "candidate_node_voltage_[branch.node_positive]-candidate_node_voltage_[branch.node_negative]"
    0)
  require_token_count("MQ.2 P1 raw branch current removed"
    ${mq2_pack_source}
    "candidate_current_[branch.cell]=numerator/resistance_[branch.cell]"
    0)
  require_token_count("MQ.2 P1 all raw branch divisions removed"
    ${mq2_pack_source}
    "/resistance_[branch.cell]"
    0)
  require_token_count("MQ.2 P1 all raw branch-current assignments removed"
    ${mq2_pack_source}
    "candidate_current_[branch.cell]="
    0)
  require_token_count("MQ.2 P1 raw affine resistance ternary removed"
    ${mq2_pack_source}
    "branch.kind==ElectricalBranchKind::cell?resistance_[branch.cell]:branch.resistance"
    0)
  require_token_count("MQ.2 P1 aliased affine resistance ternary removed"
    ${mq2_pack_source}
    "cell?resistance_[branch.cell]:branch.resistance"
    0)
  require_token_count("MQ.2 P1 raw affine source ternary removed"
    ${mq2_pack_source}
    "branch.kind==ElectricalBranchKind::cell?ocv_[branch.cell]:0.0"
    0)
  require_token_count("MQ.2 P1 aliased affine source ternary removed"
    ${mq2_pack_source}
    "cell?ocv_[branch.cell]:0.0"
    0)
endforeach()
require_token_count("MQ.2 P1 raw ladder reconstruction removed"
  mq2_pack_iterative
  "candidate_current_[cell]=current_numerator/resistance_[cell]"
  0)
require_token_count("MQ.2 P1 all raw ladder-current assignments removed"
  mq2_pack_iterative
  "candidate_current_[cell]="
  0)
require_token_count("MQ.2 P1 sparse KCL policy and order"
  mq2_pack_sparse
  "constreal_tvoltage=branchDrop(branch,candidate_node_voltage_);assert(is_finite(voltage));constautoaffine=branchAffine(branch,ocv_,resistance_);constreal_tnumerator=branchCurrentNumerator(voltage,affine.source);if(!is_finite(numerator))returnslide::Status::Invalid_states;constreal_tbranch_current=branchCurrentOut(numerator,affine.resistance);"
  1)
require_token_count("MQ.2 P1 relaxation KCL policy and order"
  mq2_pack_iterative
  "constreal_tvoltage=branchDrop(branch,candidate_node_voltage_);if(!is_finite(voltage))returnslide::Status::Invalid_states;constautoaffine=branchAffine(branch,ocv_,resistance_);constreal_tnumerator=branchCurrentNumerator(voltage,affine.source);assert(is_finite(numerator));constreal_tbranch_current=branchCurrentOut(numerator,affine.resistance);"
  1)

set(mq2_sparse_undamped_begin "real_tdamping=1.0;")
set(mq2_sparse_undamped_end
  "assert(is_finite(damping)&&damping>=0.0&&damping<=1.0);")
require_token_count("MQ.2 P1 sparse undamped slice begin"
  mq2_pack_sparse "${mq2_sparse_undamped_begin}" 1)
require_token_count("MQ.2 P1 sparse undamped slice end"
  mq2_pack_sparse "${mq2_sparse_undamped_end}" 1)
string(FIND "${mq2_pack_sparse}" "${mq2_sparse_undamped_begin}"
  mq2_sparse_undamped_begin_position)
string(FIND "${mq2_pack_sparse}" "${mq2_sparse_undamped_end}"
  mq2_sparse_undamped_end_position)
if(mq2_sparse_undamped_begin_position EQUAL -1
   OR mq2_sparse_undamped_end_position EQUAL -1
   OR mq2_sparse_undamped_end_position LESS_EQUAL mq2_sparse_undamped_begin_position)
  message(FATAL_ERROR "MQ.2 P1 sparse undamped reconstruction slice is missing")
endif()
math(EXPR mq2_sparse_undamped_length
  "${mq2_sparse_undamped_end_position} - ${mq2_sparse_undamped_begin_position}")
string(SUBSTRING "${mq2_pack_sparse}"
  ${mq2_sparse_undamped_begin_position}
  ${mq2_sparse_undamped_length}
  mq2_sparse_undamped_slice)
require_token_count("MQ.2 P1 sparse undamped reconstruction"
  mq2_sparse_undamped_slice "cellCurrentFromDrop(" 1)

set(mq2_sparse_damped_begin "if(damping<1.0){")
set(mq2_sparse_damped_end
  "for(constauto&branch:topology_.electrical.branches)if(branch.kind==ElectricalBranchKind::resistor)")
require_token_count("MQ.2 P1 sparse damped slice begin"
  mq2_pack_sparse "${mq2_sparse_damped_begin}" 1)
require_token_count("MQ.2 P1 sparse damped slice end"
  mq2_pack_sparse "${mq2_sparse_damped_end}" 1)
string(FIND "${mq2_pack_sparse}" "${mq2_sparse_damped_begin}"
  mq2_sparse_damped_begin_position)
string(FIND "${mq2_pack_sparse}" "${mq2_sparse_damped_end}"
  mq2_sparse_damped_end_position)
if(mq2_sparse_damped_begin_position EQUAL -1
   OR mq2_sparse_damped_end_position EQUAL -1
   OR mq2_sparse_damped_end_position LESS_EQUAL mq2_sparse_damped_begin_position)
  message(FATAL_ERROR "MQ.2 P1 sparse damped reconstruction slice is missing")
endif()
math(EXPR mq2_sparse_damped_length
  "${mq2_sparse_damped_end_position} - ${mq2_sparse_damped_begin_position}")
string(SUBSTRING "${mq2_pack_sparse}"
  ${mq2_sparse_damped_begin_position}
  ${mq2_sparse_damped_length}
  mq2_sparse_damped_slice)
require_token_count("MQ.2 P1 sparse damped reconstruction"
  mq2_sparse_damped_slice "cellCurrentFromDrop(" 1)

forbid_tokens("MQ.2 P1 dead sparse declarations"
  mq2_pack_sparse
  "#include<cstring>"
  "usingdetail::addCompensatedFinite;")
forbid_tokens("MQ.2 P1 dead iterative declarations"
  mq2_pack_iterative
  "usingdetail::conservativePackRoundoffBound;"
  "usingdetail::finiteCandidate;")

file(READ "${SLIDE_SOURCE_DIR}/src/core/PackSolverInternal.hpp"
  mq2_pack_internal_with_comments)
string(REGEX REPLACE "[ \t\r\n]" ""
  mq2_pack_internal_with_comments "${mq2_pack_internal_with_comments}")
require_tokens("MQ.2 P1 compensated-sum rationale adjacency"
  mq2_pack_internal_with_comments
  "//Keepthesestrict-FPpragmaswiththecompensatedfinitesum:theKahan//evaluationorderisthecontract,andfastmathwouldreassociateit.#ifdefined(_MSC_VER)")
forbid_tokens("MQ.2 P1 false compensated-sum consumers"
  mq2_pack_internal_with_comments
  "directsparsesolve"
  "matrix-freestrategies")

# MQ.2 P2: Mode-C relaxation storage has seven named roles rather than
# reusing target/compensation for coefficient and KCL state. All allocation
# remains in configure(), publication is one no-throw aggregate move, and the
# solve retains exactly its six cold-sized fills.
set(mq2_pack_private_begin
  "private:[[nodiscard]]slide::StatussolveImpl(")
set(mq2_pack_private_end "friendclassPackStepper;};")
require_token_count("MQ.2 P2 PackSolver private-slice begin"
  mq2_pack_header "${mq2_pack_private_begin}" 1)
require_token_count("MQ.2 P2 PackSolver private-slice end"
  mq2_pack_header "${mq2_pack_private_end}" 1)
string(FIND "${mq2_pack_header}" "${mq2_pack_private_begin}"
  mq2_pack_private_begin_position)
string(FIND "${mq2_pack_header}" "${mq2_pack_private_end}"
  mq2_pack_private_end_position)
if(mq2_pack_private_begin_position EQUAL -1
   OR mq2_pack_private_end_position EQUAL -1
   OR mq2_pack_private_end_position LESS_EQUAL mq2_pack_private_begin_position)
  message(FATAL_ERROR "MQ.2 P2 PackSolver private slice is missing")
endif()
string(LENGTH "${mq2_pack_private_end}" mq2_pack_private_end_length)
math(EXPR mq2_pack_private_length
  "${mq2_pack_private_end_position} - ${mq2_pack_private_begin_position} + ${mq2_pack_private_end_length}")
string(SUBSTRING "${mq2_pack_header}"
  ${mq2_pack_private_begin_position}
  ${mq2_pack_private_length}
  mq2_pack_private)
require_token_count("MQ.2 P2 exact relaxation scratch owner"
  mq2_pack_private
  "structRelaxationScratch{std::vector<real_t>diagonal{};std::vector<real_t>diagonal_compensation{};std::vector<real_t>rhs{};std::vector<real_t>rhs_compensation{};std::vector<real_t>target{};std::vector<real_t>residual{};std::vector<real_t>residual_compensation{};};RelaxationScratchrelaxation_{};"
  1)

set(mq2_relaxation_candidate
  "RelaxationScratchrelaxation{.diagonal=std::vector<real_t>(nodes,0.0),.diagonal_compensation=std::vector<real_t>(nodes,0.0),.rhs=std::vector<real_t>(nodes,0.0),.rhs_compensation=std::vector<real_t>(nodes,0.0),.target=std::vector<real_t>(nodes,0.0),.residual=std::vector<real_t>(nodes,0.0),.residual_compensation=std::vector<real_t>(nodes,0.0)};")
set(mq2_relaxation_nothrow
  "static_assert(std::is_nothrow_move_assignable_v<RelaxationScratch>);")
set(mq2_relaxation_publication
  "topology_=std::move(candidate_topology);thevenin_=std::move(thevenin);batch_executor_=std::move(batch_executor);workspace_=std::move(workspace);solution_=std::move(solution);current_guess_=std::move(current_guess);candidate_current_=std::move(candidate_current);ocv_=std::move(ocv);resistance_=std::move(resistance);candidate_node_voltage_=std::move(candidate_node_voltage);layer_voltage_=std::move(layer_voltage);rollback_cell_current_=std::move(rollback_cell_current);rollback_node_voltage_=std::move(rollback_node_voltage);relaxation_=std::move(relaxation);")
require_token_count("MQ.2 P2 exact relaxation scratch candidate"
  mq2_pack_sparse "${mq2_relaxation_candidate}" 1)
require_token_count("MQ.2 P2 no-throw scratch publication"
  mq2_pack_sparse "${mq2_relaxation_nothrow}" 1)
require_token_count("MQ.2 P2 contiguous no-throw publication block"
  mq2_pack_sparse "${mq2_relaxation_publication}" 1)
require_token_count("MQ.2 P2 single scratch publication"
  mq2_pack_sparse "relaxation_=std::move(relaxation);" 1)
require_token_count("MQ.2 P2 no other scratch publication"
  mq2_pack_sparse "relaxation_=" 1)

string(FIND "${mq2_pack_sparse}" "${mq2_relaxation_candidate}"
  mq2_relaxation_candidate_position)
string(FIND "${mq2_pack_sparse}" "${mq2_relaxation_nothrow}"
  mq2_relaxation_nothrow_position)
string(FIND "${mq2_pack_sparse}" "topology_=std::move(candidate_topology);"
  mq2_first_publication_position)
if(mq2_relaxation_candidate_position EQUAL -1
   OR mq2_relaxation_nothrow_position EQUAL -1
   OR mq2_first_publication_position EQUAL -1
   OR NOT mq2_relaxation_candidate_position LESS mq2_relaxation_nothrow_position
   OR NOT mq2_relaxation_nothrow_position LESS mq2_first_publication_position)
  message(FATAL_ERROR
    "MQ.2 P2 all scratch allocation and its no-throw proof must precede member publication")
endif()

file(GLOB_RECURSE mq2_core_production_sources
  LIST_DIRECTORIES false
  "${SLIDE_SOURCE_DIR}/src/core/*.cpp"
  "${SLIDE_SOURCE_DIR}/src/core/*.cu"
  "${SLIDE_SOURCE_DIR}/src/core/*.hpp")
foreach(mq2_old_relaxation_member IN ITEMS
        relaxation_diagonal_
        relaxation_rhs_
        relaxation_target_
        relaxation_compensation_)
  foreach(mq2_core_production_source IN LISTS mq2_core_production_sources)
    file(READ "${mq2_core_production_source}" mq2_core_production_content)
    string(FIND "${mq2_core_production_content}"
      "${mq2_old_relaxation_member}" mq2_old_relaxation_member_position)
    if(NOT mq2_old_relaxation_member_position EQUAL -1)
      file(RELATIVE_PATH mq2_core_production_relative
        "${SLIDE_SOURCE_DIR}" "${mq2_core_production_source}")
      message(FATAL_ERROR
        "MQ.2 P2 old relaxation member ${mq2_old_relaxation_member} remains in production file ${mq2_core_production_relative}")
    endif()
  endforeach()
endforeach()

set(mq2_relaxation_solve_begin
  "slide::StatusPackSolver::solveRelaxation(real_tapplied_current)")
set(mq2_relaxation_solve_end
  "diagnostics_.constraint_drift=drift;returnslide::Status::Success;}")
require_token_count("MQ.2 P2 relaxation solve owner"
  mq2_pack_iterative "${mq2_relaxation_solve_begin}" 1)
require_token_count("MQ.2 P2 relaxation solve tail"
  mq2_pack_iterative "${mq2_relaxation_solve_end}" 1)
string(FIND "${mq2_pack_iterative}" "${mq2_relaxation_solve_begin}"
  mq2_relaxation_solve_begin_position)
string(FIND "${mq2_pack_iterative}" "${mq2_relaxation_solve_end}"
  mq2_relaxation_solve_end_position)
if(mq2_relaxation_solve_begin_position EQUAL -1
   OR mq2_relaxation_solve_end_position EQUAL -1
   OR mq2_relaxation_solve_end_position LESS_EQUAL mq2_relaxation_solve_begin_position)
  message(FATAL_ERROR "MQ.2 P2 relaxation solve slice is missing")
endif()
string(LENGTH "${mq2_relaxation_solve_end}"
  mq2_relaxation_solve_end_length)
math(EXPR mq2_relaxation_solve_length
  "${mq2_relaxation_solve_end_position} - ${mq2_relaxation_solve_begin_position} + ${mq2_relaxation_solve_end_length}")
string(SUBSTRING "${mq2_pack_iterative}"
  ${mq2_relaxation_solve_begin_position}
  ${mq2_relaxation_solve_length}
  mq2_relaxation_solve)

require_token_count("MQ.2 P2 relaxation diagonal uses"
  mq2_relaxation_solve "relaxation_.diagonal[" 4)
require_token_count("MQ.2 P2 relaxation diagonal-compensation uses"
  mq2_relaxation_solve "relaxation_.diagonal_compensation[" 2)
require_token_count("MQ.2 P2 relaxation rhs uses"
  mq2_relaxation_solve "relaxation_.rhs[" 6)
require_token_count("MQ.2 P2 relaxation rhs-compensation uses"
  mq2_relaxation_solve "relaxation_.rhs_compensation[" 4)
require_token_count("MQ.2 P2 relaxation target uses"
  mq2_relaxation_solve "relaxation_.target[" 4)
require_token_count("MQ.2 P2 relaxation residual uses"
  mq2_relaxation_solve "relaxation_.residual[" 5)
require_token_count("MQ.2 P2 relaxation residual-compensation uses"
  mq2_relaxation_solve "relaxation_.residual_compensation[" 4)
require_token_count("MQ.2 P2 six relaxation fills"
  mq2_relaxation_solve "std::fill(" 6)
require_token_count("MQ.2 P2 target is assigned, never filled"
  mq2_relaxation_solve "std::fill(relaxation_.target.begin()" 0)
require_token_count("MQ.2 P2 solve constructs no vector"
  mq2_relaxation_solve "std::vector<real_t>" 0)
foreach(mq2_forbidden_relaxation_mutator IN ITEMS
        "std::ranges::fill"
        "std::fill_n"
        ".assign("
        ".resize("
        ".clear("
        ".reserve("
        ".push_back("
        ".emplace_back("
        ".swap("
        "std::swap("
        "std::exchange("
        "std::copy"
        "memset(")
  require_token_count("MQ.2 P2 alternate scratch mutation is forbidden"
    mq2_relaxation_solve "${mq2_forbidden_relaxation_mutator}" 0)
endforeach()
require_token_count("MQ.2 P2 exact coefficient reset"
  mq2_relaxation_solve
  "std::fill(relaxation_.diagonal.begin(),relaxation_.diagonal.end(),0.0);std::fill(relaxation_.diagonal_compensation.begin(),relaxation_.diagonal_compensation.end(),0.0);std::fill(relaxation_.rhs.begin(),relaxation_.rhs.end(),0.0);std::fill(relaxation_.rhs_compensation.begin(),relaxation_.rhs_compensation.end(),0.0);autostamp="
  1)
require_token_count("MQ.2 P2 positive diagonal association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.diagonal[p],relaxation_.diagonal_compensation[p],conductance)"
  1)
require_token_count("MQ.2 P2 negative diagonal association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.diagonal[n],relaxation_.diagonal_compensation[n],conductance)"
  1)
require_token_count("MQ.2 P2 positive rhs association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.rhs[p],relaxation_.rhs_compensation[p],positive_rhs)"
  1)
require_token_count("MQ.2 P2 negative rhs association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.rhs[n],relaxation_.rhs_compensation[n],negative_rhs)"
  1)
require_token_count("MQ.2 P2 positive terminal rhs association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.rhs[netlist.terminal_positive],relaxation_.rhs_compensation[netlist.terminal_positive],-applied_current)"
  1)
require_token_count("MQ.2 P2 negative terminal rhs association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.rhs[netlist.terminal_negative],relaxation_.rhs_compensation[netlist.terminal_negative],applied_current)"
  1)
require_token_count("MQ.2 P2 every node receives one target"
  mq2_relaxation_solve
  "for(std::uint32_tnode=0;node<netlist.node_count;++node){if(node==netlist.terminal_negative){relaxation_.target[node]=real_t{};continue;}if(!(is_strictly_positive_finite(relaxation_.diagonal[node])&&is_finite(relaxation_.rhs[node])))returnslide::Status::Numerical_failure;relaxation_.target[node]=relaxation_.rhs[node]/relaxation_.diagonal[node];if(!is_finite(relaxation_.target[node]))returnslide::Status::Invalid_states;}"
  1)
require_token_count("MQ.2 P2 exact KCL reset immediately before assembly"
  mq2_relaxation_solve
  "std::fill(relaxation_.residual.begin(),relaxation_.residual.end(),0.0);std::fill(relaxation_.residual_compensation.begin(),relaxation_.residual_compensation.end(),0.0);real_troundoff_operation_scale="
  1)
require_token_count("MQ.2 P2 positive branch residual association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.residual[branch.node_positive],relaxation_.residual_compensation[branch.node_positive],branch_current)"
  1)
require_token_count("MQ.2 P2 negative branch residual association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.residual[branch.node_negative],relaxation_.residual_compensation[branch.node_negative],-branch_current)"
  1)
require_token_count("MQ.2 P2 positive terminal residual association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.residual[netlist.terminal_positive],relaxation_.residual_compensation[netlist.terminal_positive],applied_current)"
  1)
require_token_count("MQ.2 P2 negative terminal residual association"
  mq2_relaxation_solve
  "addCompensatedFinite(relaxation_.residual[netlist.terminal_negative],relaxation_.residual_compensation[netlist.terminal_negative],-applied_current)"
  1)
require_token_count("MQ.2 P2 residual diagnostic reads residual owner"
  mq2_relaxation_solve
  "constreal_tresidual=relaxation_.residual[node];"
  1)

# MQ.2 S1: one cold prefix-uniqueness owner, two exact arena-copy
# directions, and the public N*dt/frozen-solve semantics are executable.
load_compact("src/core/PackTopologyInternal.hpp" mq2_pack_topology_internal)
load_compact("src/core/PackTopology.hpp" mq2_pack_topology_header)
load_compact("src/core/PackTopology.cpp" mq2_pack_topology)
load_compact("src/core/PackStepper.hpp" mq2_pack_stepper_header)
load_compact("src/core/PackStepper.cpp" mq2_pack_stepper)
load_compact("src/core/PackSolverOwnership.cpp" mq2_pack_ownership)
load_compact("cmake/SlideCoreTarget.cmake" mq2_core_target)

require_token_count("MQ.2 S1.1 PackSolver explicit move contract"
  mq2_pack_header
  "PackSolver()=default;PackSolver(constPackSolver&)=delete;PackSolver&operator=(constPackSolver&)=delete;PackSolver(PackSolver&&other)noexcept;PackSolver&operator=(PackSolver&&other)noexcept;"
  1)
require_token_count("MQ.2 S1.1 PackStepper explicit move contract"
  mq2_pack_stepper_header
  "PackStepper()=default;PackStepper(constPackStepper&)=delete;PackStepper&operator=(constPackStepper&)=delete;PackStepper(PackStepper&&other)noexcept;PackStepper&operator=(PackStepper&&other)noexcept;"
  1)
require_token_count("MQ.2 S1.1 PackSolver direct utility include"
  mq2_pack_sparse "#include<utility>" 1)
require_token_count("MQ.2 S1.1 PackSolver ownership utility include"
  mq2_pack_ownership "#include<utility>" 1)
require_token_count("MQ.2 S1.1 PackSolver reset temporaries stay no-throw"
  mq2_pack_ownership
  "static_assert(std::is_nothrow_default_constructible_v<PackSolution>&&std::is_nothrow_default_constructible_v<PackSolveDiagnostics>);"
  1)
require_token_count("MQ.2 S1.1 PackStepper direct utility include"
  mq2_pack_stepper "#include<utility>" 1)
require_token_count("MQ.2 S1.1 PackSolver ownership source"
  mq2_core_target
  "\"\${PROJECT_SOURCE_DIR}/src/core/PackSolverOwnership.cpp\""
  1)
require_token_count("MQ.2 S1.1 exact SolverWorkspace moves"
  mq2_pack_sparse
  "SolverWorkspace::SolverWorkspace(SolverWorkspace&&other)noexcept:impl_{std::move(other.impl_)},factorized_resistance_{std::move(other.factorized_resistance_)},valid_{std::exchange(other.valid_,false)},age_{std::exchange(other.age_,0)},numeric_factorizations_{std::exchange(other.numeric_factorizations_,0)},symbolic_factorizations_{std::exchange(other.symbolic_factorizations_,0)}{static_assert(detail::nothrow_movable<std::unique_ptr<Impl>,std::vector<real_t>,bool,int>);}SolverWorkspace&SolverWorkspace::operator=(SolverWorkspace&&other)noexcept{if(this!=&other){impl_=std::move(other.impl_);factorized_resistance_=std::move(other.factorized_resistance_);valid_=std::exchange(other.valid_,false);age_=std::exchange(other.age_,0);numeric_factorizations_=std::exchange(other.numeric_factorizations_,0);symbolic_factorizations_=std::exchange(other.symbolic_factorizations_,0);}return*this;}"
  1)
require_token_count("MQ.2 S1.1 exact PackSolver moves"
  mq2_pack_ownership
  "PackSolver::PackSolver(PackSolver&&other)noexcept:topology_{std::move(other.topology_)},thevenin_{std::move(other.thevenin_)},batch_executor_{std::move(other.batch_executor_)},workspace_{std::move(other.workspace_)},solution_{std::exchange(other.solution_,PackSolution{})},diagnostics_{std::exchange(other.diagnostics_,PackSolveDiagnostics{})},current_guess_{std::move(other.current_guess_)},candidate_current_{std::move(other.candidate_current_)},ocv_{std::move(other.ocv_)},resistance_{std::move(other.resistance_)},candidate_node_voltage_{std::move(other.candidate_node_voltage_)},layer_voltage_{std::move(other.layer_voltage_)},rollback_cell_current_{std::move(other.rollback_cell_current_)},rollback_node_voltage_{std::move(other.rollback_node_voltage_)},relaxation_{std::move(other.relaxation_)},candidate_terminal_voltage_{std::exchange(other.candidate_terminal_voltage_,real_t{})},residual_norm_{std::exchange(other.residual_norm_,real_t{})},relaxation_alpha_{std::exchange(other.relaxation_alpha_,real_t{2.0/3.0})},configured_{std::exchange(other.configured_,false)},has_solution_{std::exchange(other.has_solution_,false)}{static_assert(nothrow_movable<CompiledPackTopology,PackTheveninSystem,BatchExecutor,SolverWorkspace,PackSolution,PackSolveDiagnostics,std::vector<real_t>,RelaxationScratch,real_t,bool>);}PackSolver&PackSolver::operator=(PackSolver&&other)noexcept{if(this!=&other){topology_=std::move(other.topology_);thevenin_=std::move(other.thevenin_);batch_executor_=std::move(other.batch_executor_);workspace_=std::move(other.workspace_);solution_=std::exchange(other.solution_,PackSolution{});diagnostics_=std::exchange(other.diagnostics_,PackSolveDiagnostics{});current_guess_=std::move(other.current_guess_);candidate_current_=std::move(other.candidate_current_);ocv_=std::move(other.ocv_);resistance_=std::move(other.resistance_);candidate_node_voltage_=std::move(other.candidate_node_voltage_);layer_voltage_=std::move(other.layer_voltage_);rollback_cell_current_=std::move(other.rollback_cell_current_);rollback_node_voltage_=std::move(other.rollback_node_voltage_);relaxation_=std::move(other.relaxation_);candidate_terminal_voltage_=std::exchange(other.candidate_terminal_voltage_,real_t{});residual_norm_=std::exchange(other.residual_norm_,real_t{});relaxation_alpha_=std::exchange(other.relaxation_alpha_,real_t{2.0/3.0});configured_=std::exchange(other.configured_,false);has_solution_=std::exchange(other.has_solution_,false);}return*this;}"
  1)
require_token_count("MQ.2 S1.1 exact PackStepper moves"
  mq2_pack_stepper
  "PackStepper::PackStepper(PackStepper&&other)noexcept:topology_{std::move(other.topology_)},batches_{std::move(other.batches_)},steppers_{std::move(other.steppers_)},exponential_steppers_{std::move(other.exponential_steppers_)},solver_{std::move(other.solver_)},current_density_{std::move(other.current_density_)},checkpoint_offsets_{std::move(other.checkpoint_offsets_)},checkpoint_{std::move(other.checkpoint_)},cell_temperature_{std::move(other.cell_temperature_)},cell_external_heat_{std::move(other.cell_external_heat_)},boundary_heat_{std::move(other.boundary_heat_)},solver_checkpoint_solution_{std::move(other.solver_checkpoint_solution_)},solver_checkpoint_diagnostics_{std::move(other.solver_checkpoint_diagnostics_)},cell_external_heat_checkpoint_{std::move(other.cell_external_heat_checkpoint_)},boundary_heat_checkpoint_{std::move(other.boundary_heat_checkpoint_)},solver_checkpoint_has_solution_{std::exchange(other.solver_checkpoint_has_solution_,false)},configured_{std::exchange(other.configured_,false)}{other.checkpoint_.clear();other.cell_external_heat_.clear();other.boundary_heat_.clear();}PackStepper&PackStepper::operator=(PackStepper&&other)noexcept{if(this!=&other){topology_=std::move(other.topology_);batches_=std::move(other.batches_);steppers_=std::move(other.steppers_);exponential_steppers_=std::move(other.exponential_steppers_);solver_=std::move(other.solver_);current_density_=std::move(other.current_density_);checkpoint_offsets_=std::move(other.checkpoint_offsets_);checkpoint_=std::move(other.checkpoint_);cell_temperature_=std::move(other.cell_temperature_);cell_external_heat_=std::move(other.cell_external_heat_);boundary_heat_=std::move(other.boundary_heat_);solver_checkpoint_solution_=std::move(other.solver_checkpoint_solution_);solver_checkpoint_diagnostics_=std::move(other.solver_checkpoint_diagnostics_);cell_external_heat_checkpoint_=std::move(other.cell_external_heat_checkpoint_);boundary_heat_checkpoint_=std::move(other.boundary_heat_checkpoint_);solver_checkpoint_has_solution_=std::exchange(other.solver_checkpoint_has_solution_,false);configured_=std::exchange(other.configured_,false);other.checkpoint_.clear();other.cell_external_heat_.clear();other.boundary_heat_.clear();}return*this;}"
  1)

require_token_count("MQ.2 S1 exact prefix-uniqueness owner"
  mq2_pack_topology_internal
  "template<classRange,classProjection=std::identity>requiresstd::ranges::random_access_range<constRange>&&std::ranges::sized_range<constRange>[[nodiscard]]constexprboolfirstOccurrence(constRange&range,std::size_tindex,Projectionprojection={}){if(index>=std::ranges::size(range))returnfalse;constautofirst=std::ranges::begin(range);constautocurrent=first+static_cast<std::ranges::range_difference_t<constRange>>(index);returnstd::ranges::find(first,current,std::invoke(projection,*current),projection)==current;}"
  1)
foreach(mq2_prefix_owner_include IN ITEMS
        "#include<algorithm>"
        "#include<functional>"
        "#include<ranges>")
  require_token_count("MQ.2 S1 direct prefix-owner include"
    mq2_pack_topology_internal "${mq2_prefix_owner_include}" 1)
endforeach()
require_token_count("MQ.2 S1 prefix owner count"
  mq2_pack_topology_internal "firstOccurrence(" 1)
require_token_count("MQ.2 S1 prefix PackSolver consumers"
  mq2_pack_sparse "firstOccurrence(" 2)
require_token_count("MQ.2 S1 prefix PackStepper consumer"
  mq2_pack_stepper "firstOccurrence(" 1)
require_token_count("MQ.2 S1 archetype-prefix rejection"
  mq2_pack_sparse
  "batch_archetypes[batch].empty()||!detail::firstOccurrence(batch_archetypes,batch)"
  1)
require_token_count("MQ.2 S1 identity-prefix rejection"
  mq2_pack_sparse
  "!detail::firstOccurrence(batches,batch,&TheveninBatchView::identity)"
  1)
require_token_count("MQ.2 S1 pointer-prefix rejection"
  mq2_pack_stepper
  "!detail::firstOccurrence(batches,batch)"
  1)
require_token_count("MQ.2 S1 old archetype prefix search removed"
  mq2_pack_sparse "std::find(batch_archetypes.begin()" 0)
require_token_count("MQ.2 S1 old identity prefix loop removed"
  mq2_pack_sparse "for(std::size_tprior=0;prior<batch;++prior)" 0)
require_token_count("MQ.2 S1 old pointer prefix search removed"
  mq2_pack_stepper "std::find(batches.begin()" 0)

set(mq2_pack_stepper_private_begin
  "private:[[nodiscard]]slide::StatusstepImpl(")
set(mq2_pack_stepper_private_end "boolconfigured_{};};")
require_token_count("MQ.2 S1 PackStepper private-slice begin"
  mq2_pack_stepper_header "${mq2_pack_stepper_private_begin}" 1)
require_token_count("MQ.2 S1 PackStepper private-slice end"
  mq2_pack_stepper_header "${mq2_pack_stepper_private_end}" 1)
string(FIND "${mq2_pack_stepper_header}" "${mq2_pack_stepper_private_begin}"
  mq2_pack_stepper_private_begin_position)
string(FIND "${mq2_pack_stepper_header}" "${mq2_pack_stepper_private_end}"
  mq2_pack_stepper_private_end_position)
if(mq2_pack_stepper_private_begin_position EQUAL -1
   OR mq2_pack_stepper_private_end_position EQUAL -1
   OR mq2_pack_stepper_private_end_position
      LESS_EQUAL mq2_pack_stepper_private_begin_position)
  message(FATAL_ERROR "MQ.2 S1 PackStepper private slice is missing")
endif()
string(LENGTH "${mq2_pack_stepper_private_end}"
  mq2_pack_stepper_private_end_length)
math(EXPR mq2_pack_stepper_private_length
  "${mq2_pack_stepper_private_end_position} - ${mq2_pack_stepper_private_begin_position} + ${mq2_pack_stepper_private_end_length}")
string(SUBSTRING "${mq2_pack_stepper_header}"
  ${mq2_pack_stepper_private_begin_position}
  ${mq2_pack_stepper_private_length}
  mq2_pack_stepper_private)
require_token_count("MQ.2 S1 private gather declaration"
  mq2_pack_stepper_private
  "voidgatherStates(std::span<real_t>destination)const;"
  1)
require_token_count("MQ.2 S1 private scatter declaration"
  mq2_pack_stepper_private
  "voidscatterStates(std::span<constreal_t>source);"
  1)
require_token_count("MQ.2 S1 exact gather owner"
  mq2_pack_stepper
  "voidPackStepper::gatherStates(std::span<real_t>destination)const{assert(destination.size()==checkpoint_.size());for(std::size_tbatch=0;batch<batches_.size();++batch){constautostate=batches_[batch]->state().raw();std::memcpy(destination.data()+checkpoint_offsets_[batch],state.data(),state.size_bytes());}}"
  1)
require_token_count("MQ.2 S1 exact scatter owner"
  mq2_pack_stepper
  "voidPackStepper::scatterStates(std::span<constreal_t>source){assert(source.size()==checkpoint_.size());for(std::size_tbatch=0;batch<batches_.size();++batch){autostate=batches_[batch]->state().raw();std::memcpy(state.data(),source.data()+checkpoint_offsets_[batch],state.size_bytes());}}"
  1)
require_token_count("MQ.2 S1 gather owner and callers"
  mq2_pack_stepper "gatherStates(" 3)
require_token_count("MQ.2 S1 scatter owner and callers"
  mq2_pack_stepper "scatterStates(" 3)
require_token_count("MQ.2 S1 exact arena-copy directions"
  mq2_pack_stepper "std::memcpy(" 2)
require_token_count("MQ.2 S1 internal gather call"
  mq2_pack_stepper "gatherStates(checkpoint_);" 1)
require_token_count("MQ.2 S1 public gather call"
  mq2_pack_stepper "gatherStates(destination);" 1)
require_token_count("MQ.2 S1 internal scatter call"
  mq2_pack_stepper "scatterStates(checkpoint_);" 1)
require_token_count("MQ.2 S1 public scatter call"
  mq2_pack_stepper "scatterStates(source);" 1)
require_token_count("MQ.2 S1 internal gather placement"
  mq2_pack_stepper
  "voidPackStepper::saveCheckpoint(){gatherStates(checkpoint_);solver_checkpoint_solution_.cell_current="
  1)
require_token_count("MQ.2 S1 internal scatter placement"
  mq2_pack_stepper
  "voidPackStepper::restoreCheckpoint(){scatterStates(checkpoint_);solver_.solution_.cell_current="
  1)
require_token_count("MQ.2 S1 exact public gather transaction"
  mq2_pack_stepper
  "slide::StatusPackStepper::checkpoint(std::span<real_t>destination)const{if(!configured_||destination.size()!=checkpoint_.size())returnslide::Status::Invalid_parameters;gatherStates(destination);returnslide::Status::Success;}"
  1)
require_token_count("MQ.2 S1 exact public scatter transaction"
  mq2_pack_stepper
  "slide::StatusPackStepper::restore(std::span<constreal_t>source){if(!configured_||source.size()!=checkpoint_.size())returnslide::Status::Invalid_parameters;scatterStates(source);solver_.invalidate();returnslide::Status::Success;}"
  1)

file(READ "${SLIDE_SOURCE_DIR}/src/core/PackStepper.hpp"
  mq2_pack_stepper_header_with_comments)
string(REGEX REPLACE "[ \t\r\n]" ""
  mq2_pack_stepper_header_with_comments
  "${mq2_pack_stepper_header_with_comments}")
require_token_count("MQ.2 S1 public N*dt contract"
  mq2_pack_stepper_header_with_comments
  "Advancesthepackby`substeps*dt`,NOTby`dt`."
  1)
require_token_count("MQ.2 S1 public frozen-current contract"
  mq2_pack_stepper_header_with_comments
  "electricalcurrentisre-solvedonce,thenheldfrozenacrossthem"
  1)

set(mq2_pack_step_impl_begin
  "slide::StatusPackStepper::stepImpl(")
set(mq2_pack_step_impl_end
  "if(status!=slide::Status::Success)returnfail(status);}returnslide::Status::Success;}")
require_token_count("MQ.2 S1 PackStepper stepImpl owner"
  mq2_pack_stepper "${mq2_pack_step_impl_begin}" 1)
require_token_count("MQ.2 S1 PackStepper stepImpl tail"
  mq2_pack_stepper "${mq2_pack_step_impl_end}" 1)
string(FIND "${mq2_pack_stepper}" "${mq2_pack_step_impl_begin}"
  mq2_pack_step_impl_begin_position)
string(FIND "${mq2_pack_stepper}" "${mq2_pack_step_impl_end}"
  mq2_pack_step_impl_end_position)
if(mq2_pack_step_impl_begin_position EQUAL -1
   OR mq2_pack_step_impl_end_position EQUAL -1
   OR mq2_pack_step_impl_end_position LESS_EQUAL mq2_pack_step_impl_begin_position)
  message(FATAL_ERROR "MQ.2 S1 PackStepper stepImpl slice is missing")
endif()
string(LENGTH "${mq2_pack_step_impl_end}" mq2_pack_step_impl_end_length)
math(EXPR mq2_pack_step_impl_length
  "${mq2_pack_step_impl_end_position} - ${mq2_pack_step_impl_begin_position} + ${mq2_pack_step_impl_end_length}")
string(SUBSTRING "${mq2_pack_stepper}"
  ${mq2_pack_step_impl_begin_position}
  ${mq2_pack_step_impl_length}
  mq2_pack_step_impl)

set(mq2_substeps_loop_begin
  "for(intsubstep=0;substep<substeps;++substep){")
set(mq2_substeps_loop_end
  "if(status!=slide::Status::Success)returnfail(status);}")
require_token_count("MQ.2 S1 exact substeps loop"
  mq2_pack_step_impl "${mq2_substeps_loop_begin}" 1)
require_token_count("MQ.2 S1 exact substeps loop tail"
  mq2_pack_step_impl "${mq2_substeps_loop_end}" 1)
string(FIND "${mq2_pack_step_impl}" "${mq2_substeps_loop_begin}"
  mq2_substeps_loop_begin_position)
string(FIND "${mq2_pack_step_impl}" "${mq2_substeps_loop_end}"
  mq2_substeps_loop_end_position)
if(mq2_substeps_loop_begin_position EQUAL -1
   OR mq2_substeps_loop_end_position EQUAL -1
   OR mq2_substeps_loop_end_position LESS_EQUAL mq2_substeps_loop_begin_position)
  message(FATAL_ERROR "MQ.2 S1 PackStepper substeps-loop slice is missing")
endif()
string(LENGTH "${mq2_substeps_loop_end}" mq2_substeps_loop_end_length)
math(EXPR mq2_substeps_loop_length
  "${mq2_substeps_loop_end_position} - ${mq2_substeps_loop_begin_position} + ${mq2_substeps_loop_end_length}")
string(SUBSTRING "${mq2_pack_step_impl}"
  ${mq2_substeps_loop_begin_position}
  ${mq2_substeps_loop_length}
  mq2_substeps_loop)
string(SUBSTRING "${mq2_pack_step_impl}"
  0 ${mq2_substeps_loop_begin_position} mq2_pack_step_preloop)

require_token_count("MQ.2 S1 one frozen electrical solve"
  mq2_pack_step_impl "solver_.solve(" 1)
require_token_count("MQ.2 S1 electrical solve precedes substeps"
  mq2_pack_step_preloop "solver_.solve(" 1)
require_token_count("MQ.2 S1 no electrical solve inside substeps"
  mq2_substeps_loop "solver_.solve(" 0)
require_token_count("MQ.2 S1 one frozen thermal assembly"
  mq2_pack_step_impl "topology_.thermal.assemble(" 1)
require_token_count("MQ.2 S1 thermal assembly precedes substeps"
  mq2_pack_step_preloop "topology_.thermal.assemble(" 1)
require_token_count("MQ.2 S1 no thermal assembly inside substeps"
  mq2_substeps_loop "topology_.thermal.assemble(" 0)
require_token_count("MQ.2 S1 full-dt substep calls"
  mq2_substeps_loop
  "time+static_cast<real_t>(substep)*dt,dt)"
  2)

# MQ.2 T1: one cold owner derives branch adjacency/sparsity for both generated
# and imported topologies, and one sorted slot record owns all archetype batch
# metadata. The validator ordering keeps the P9-B18 allocation bound and every
# endpoint guard ahead of graph construction.
set(mq2_batch_locations_begin
  "boolassignBatchLocations(CompiledPackTopology&pack){")
set(mq2_batch_locations_end
  "voidcompileElectricalMetadata(CompiledPackTopology&pack,std::uint32_tnode_count){")
require_token_count("MQ.2 T1 batch-location slice begin"
  mq2_pack_topology "${mq2_batch_locations_begin}" 1)
require_token_count("MQ.2 T1 batch-location slice end"
  mq2_pack_topology "${mq2_batch_locations_end}" 1)
string(FIND "${mq2_pack_topology}" "${mq2_batch_locations_begin}"
  mq2_batch_locations_begin_position)
string(FIND "${mq2_pack_topology}" "${mq2_batch_locations_end}"
  mq2_batch_locations_end_position)
if(mq2_batch_locations_begin_position EQUAL -1
   OR mq2_batch_locations_end_position EQUAL -1
   OR mq2_batch_locations_end_position LESS_EQUAL mq2_batch_locations_begin_position)
  message(FATAL_ERROR "MQ.2 T1 batch-location slice is missing")
endif()
math(EXPR mq2_batch_locations_length
  "${mq2_batch_locations_end_position} - ${mq2_batch_locations_begin_position}")
string(SUBSTRING "${mq2_pack_topology}"
  ${mq2_batch_locations_begin_position}
  ${mq2_batch_locations_length}
  mq2_batch_locations)

require_token_count("MQ.2 T1 one BatchSlot record"
  mq2_batch_locations
  "structBatchSlot{std::uint32_tbatch{};std::uint32_tnext_lane{};boolthermal{};};"
  1)
require_token_count("MQ.2 T1 one archetype map"
  mq2_batch_locations "std::map<" 1)
require_token_count("MQ.2 T1 exact slot map"
  mq2_batch_locations "std::map<std::string,BatchSlot>slots;" 1)
require_token_count("MQ.2 T1 one slot insertion"
  mq2_batch_locations
  "slots.try_emplace(cell.archetype,BatchSlot{.thermal=cell.thermal})"
  1)
require_token_count("MQ.2 T1 one checked slot lookup"
  mq2_batch_locations "slots.at(cell.archetype)" 1)
require_token_count("MQ.2 T1 one post-increment lane publication"
  mq2_batch_locations
  "cell.location={.batch=slot.batch,.lane=slot.next_lane++};"
  1)
require_token_count("MQ.2 T1 old integer maps removed"
  mq2_batch_locations "std::map<std::string,std::uint32_t>" 0)
require_token_count("MQ.2 T1 old thermal map removed"
  mq2_batch_locations "std::map<std::string,bool>" 0)
require_ordered_tokens("MQ.2 T1 validate before batch publication"
  mq2_batch_locations
  "slots.try_emplace(cell.archetype,BatchSlot{.thermal=cell.thermal})"
  "if(!inserted&&it->second.thermal!=cell.thermal)"
  "for(auto&[name,slot]:slots)"
  "slots.at(cell.archetype)"
  "cell.location={.batch=slot.batch,.lane=slot.next_lane++};")

require_token_count("MQ.2 T1 BranchGraph owner"
  mq2_pack_topology "BranchGraphbuildBranchGraph(" 1)
require_token_count("MQ.2 T1 BranchGraph owner and consumers"
  mq2_pack_topology "buildBranchGraph(" 3)
require_token_count("MQ.2 T1 connectivity owner"
  mq2_pack_topology "boolisConnectedFrom(" 1)
require_token_count("MQ.2 T1 connectivity owner and consumers"
  mq2_pack_topology "isConnectedFrom(" 3)
require_token_count("MQ.2 T1 exact connectivity consumers"
  mq2_pack_topology
  "isConnectedFrom(graph.adjacency,netlist.terminal_positive)"
  2)
require_token_count("MQ.2 T1 positive adjacency owner"
  mq2_pack_topology
  "graph.adjacency[branch.node_positive].push_back(branch.node_negative);"
  1)
require_token_count("MQ.2 T1 negative adjacency owner"
  mq2_pack_topology
  "graph.adjacency[branch.node_negative].push_back(branch.node_positive);"
  1)
require_token_count("MQ.2 T1 global positive adjacency census"
  mq2_pack_topology
  "adjacency[branch.node_positive].push_back(branch.node_negative);"
  1)
require_token_count("MQ.2 T1 global negative adjacency census"
  mq2_pack_topology
  "adjacency[branch.node_negative].push_back(branch.node_positive);"
  1)
require_token_count("MQ.2 T1 global sparsity-contribution census"
  mq2_pack_topology "sparsity.emplace_back(" 3)
require_token_count("MQ.2 T1 one graph traversal queue"
  mq2_pack_topology "std::queue<std::uint32_t>pending;" 1)
require_token_count("MQ.2 T1 one sparsity overflow proof"
  mq2_pack_topology "sizeof(CompiledElectricalBranch)>3" 1)
require_token_count("MQ.2 T1 one sparsity reserve"
  mq2_pack_topology "graph.sparsity.reserve(branches.size()*3);" 1)
require_token_count("MQ.2 T1 old validator sparsity owner removed"
  mq2_pack_topology "expected_sparsity" 0)

set(mq2_branch_graph_begin "BranchGraphbuildBranchGraph(")
set(mq2_branch_graph_end "boolisConnectedFrom(")
string(FIND "${mq2_pack_topology}" "${mq2_branch_graph_begin}"
  mq2_branch_graph_begin_position)
string(FIND "${mq2_pack_topology}" "${mq2_branch_graph_end}"
  mq2_branch_graph_end_position)
if(mq2_branch_graph_begin_position EQUAL -1
   OR mq2_branch_graph_end_position EQUAL -1
   OR mq2_branch_graph_end_position LESS_EQUAL mq2_branch_graph_begin_position)
  message(FATAL_ERROR "MQ.2 T1 BranchGraph owner slice is missing")
endif()
math(EXPR mq2_branch_graph_length
  "${mq2_branch_graph_end_position} - ${mq2_branch_graph_begin_position}")
string(SUBSTRING "${mq2_pack_topology}"
  ${mq2_branch_graph_begin_position}
  ${mq2_branch_graph_length}
  mq2_branch_graph)
require_token_count("MQ.2 T1 positive adjacency lives in graph owner"
  mq2_branch_graph
  "graph.adjacency[branch.node_positive].push_back(branch.node_negative);"
  1)
require_token_count("MQ.2 T1 negative adjacency lives in graph owner"
  mq2_branch_graph
  "graph.adjacency[branch.node_negative].push_back(branch.node_positive);"
  1)
require_token_count("MQ.2 T1 three sparsity contributions live in graph owner"
  mq2_branch_graph "graph.sparsity.emplace_back(" 3)

set(mq2_netlist_validator_begin
  "slide::Statusdetail::validateElectricalNetlist(")
set(mq2_netlist_validator_end
  "slide::Statusdetail::finalizeImportedPackTopology(")
string(FIND "${mq2_pack_topology}" "${mq2_netlist_validator_begin}"
  mq2_netlist_validator_begin_position)
string(FIND "${mq2_pack_topology}" "${mq2_netlist_validator_end}"
  mq2_netlist_validator_end_position)
if(mq2_netlist_validator_begin_position EQUAL -1
   OR mq2_netlist_validator_end_position EQUAL -1
   OR mq2_netlist_validator_end_position
      LESS_EQUAL mq2_netlist_validator_begin_position)
  message(FATAL_ERROR "MQ.2 T1 electrical-validator slice is missing")
endif()
math(EXPR mq2_netlist_validator_length
  "${mq2_netlist_validator_end_position} - ${mq2_netlist_validator_begin_position}")
string(SUBSTRING "${mq2_pack_topology}"
  ${mq2_netlist_validator_begin_position}
  ${mq2_netlist_validator_length}
  mq2_netlist_validator)
require_ordered_tokens("MQ.2 T1 guarded validator graph construction"
  mq2_netlist_validator
  "static_cast<std::size_t>(netlist.node_count-1)>netlist.branches.size()"
  "netlist.terminal_positive>=netlist.node_count"
  "for(constauto&branch:netlist.branches){"
  "if(branch.node_positive>=netlist.node_count||branch.node_negative>=netlist.node_count||branch.node_positive==branch.node_negative)"
  "if(cell_branches!=cell_count||std::any_of(seen_cell.begin(),seen_cell.end(),[](unsignedcharseen){returnseen==0;}))"
  "constautograph=buildBranchGraph(netlist.branches,netlist.node_count);"
  "if(netlist.nodal_sparsity!=graph.sparsity)"
  "if(!isConnectedFrom(graph.adjacency,netlist.terminal_positive))")

# MQ.2 T2: the per-step thermal assembly is explicitly no-throw and the
# leading PackTopology contract names its owner plus the hot/cold boundary.
require_token_count("MQ.2 T2 noexcept declaration"
  mq2_pack_topology_header
  "[[nodiscard]]slide::Statusassemble(std::span<constreal_t>cell_temperature,std::span<constreal_t>boundary_temperature,std::span<real_t>q_ext,std::span<real_t>boundary_heat)noexcept;"
  1)
require_token_count("MQ.2 T2 noexcept definition"
  mq2_pack_topology
  "slide::StatusCompiledThermalGraph::assemble(std::span<constreal_t>cell_temperature,std::span<constreal_t>boundary_temperature,std::span<real_t>q_ext,std::span<real_t>boundary_heat)noexcept{"
  1)

file(READ "${SLIDE_SOURCE_DIR}/src/core/PackTopology.hpp"
  mq2_pack_topology_header_with_comments)
string(REGEX REPLACE "[ \t\r\n]" ""
  mq2_pack_topology_header_with_comments
  "${mq2_pack_topology_header_with_comments}")
string(FIND "${mq2_pack_topology_header_with_comments}" "#pragmaonce"
  mq2_pack_topology_pragma_position)
if(mq2_pack_topology_pragma_position EQUAL -1)
  message(FATAL_ERROR "MQ.2 T2 PackTopology #pragma once is missing")
endif()
string(SUBSTRING "${mq2_pack_topology_header_with_comments}"
  0 ${mq2_pack_topology_pragma_position}
  mq2_pack_topology_leading_contract)
require_token_count("MQ.2 T2 substantive topology brief"
  mq2_pack_topology_header_with_comments
  "@briefValue-semanticpackcombinatorsandcompiledelectrical/thermaltopology."
  1)
require_token_count("MQ.2 T2 topology ownership"
  mq2_pack_topology_header_with_comments
  "OwnsD-19packauthoring/compilationpluscanonicalelectricalandD-21thermaloutputs."
  1)
require_token_count("MQ.2 T2 exact PLAN owner"
  mq2_pack_topology_header_with_comments
  "ImplementsPLAN.mdsection3.4(D-19andD-21)."
  1)
require_token_count("MQ.2 T2 cold boundary"
  mq2_pack_topology_header_with_comments
  "Cold:descriptionauthoringandcompilePackDescription()."
  1)
require_token_count("MQ.2 T2 hot boundary"
  mq2_pack_topology_header_with_comments
  "Hot:CompiledThermalGraph::assemble(),atmostonceperPackStepperstepattempt."
  1)
require_token_count("MQ.2 T2 stale all-cold brief removed"
  mq2_pack_topology_header_with_comments
  "@briefValue-semanticpackcombinatorsandcold-compiledelectrical/thermaltopology."
  0)
require_ordered_tokens("MQ.2 T2 leading contract placement"
  mq2_pack_topology_leading_contract
  "@briefValue-semanticpackcombinatorsandcompiledelectrical/thermaltopology."
  "OwnsD-19packauthoring/compilationpluscanonicalelectricalandD-21thermaloutputs."
  "ImplementsPLAN.mdsection3.4(D-19andD-21)."
  "Cold:descriptionauthoringandcompilePackDescription()."
  "Hot:CompiledThermalGraph::assemble(),atmostonceperPackStepperstepattempt.")

# MQ.2 C1: NetlistCsv retains one source for each limit diagnostic, one
# ASCII-only digit scanner, and one Status-returning semantic-failure owner.
# Exact cursor/status/grammar slices keep the behavior-preserving extraction
# from widening the accepted grammar or changing first-failure diagnostics.
load_compact("src/core/NetlistCsv.cpp" mq2_netlist_csv)
load_compact("src/core/NetlistCsv.hpp" mq2_netlist_csv_header)
file(READ "${SLIDE_SOURCE_DIR}/src/core/NetlistCsv.cpp"
  mq2_netlist_csv_with_comments)
string(REGEX REPLACE "[ \t\r\n]" ""
  mq2_netlist_csv_with_comments
  "${mq2_netlist_csv_with_comments}")
file(READ "${SLIDE_SOURCE_DIR}/src/core/NetlistCsv.hpp"
  mq2_netlist_csv_header_with_comments)
string(REGEX REPLACE "[ \t\r\n]" ""
  mq2_netlist_csv_header_with_comments
  "${mq2_netlist_csv_header_with_comments}")

require_token_count("MQ.2 C1 msg_too_large owner"
  mq2_netlist_csv
  "constexprstd::string_viewmsg_too_large=\"liionpackCSVexceeds4194304bytes\";"
  1)
require_token_count("MQ.2 C1 msg_too_wide owner"
  mq2_netlist_csv
  "constexprstd::string_viewmsg_too_wide=\"CSVrowexceeds32columns\";"
  1)
require_token_count("MQ.2 C1 msg_field_large owner"
  mq2_netlist_csv
  "constexprstd::string_viewmsg_field_large=\"CSVfieldexceeds65536bytes\";"
  1)
require_token_count("MQ.2 C1 msg_too_many_rows owner"
  mq2_netlist_csv
  "constexprstd::string_viewmsg_too_many_rows=\"liionpackCSVexceeds100000rows\";"
  1)
require_token_count("MQ.2 C1 numeric limit/message pin"
  mq2_netlist_csv
  "static_assert(max_csv_bytes==4194304U&&max_csv_rows==100000U&&max_csv_columns==32U&&max_csv_field_bytes==65536U,\"diagnosticmessagetextquotestheselimitsverbatim\");"
  1)

require_token_count("MQ.2 C1 too-large literal owner"
  mq2_netlist_csv "\"liionpackCSVexceeds4194304bytes\"" 1)
require_token_count("MQ.2 C1 too-wide literal owner"
  mq2_netlist_csv "\"CSVrowexceeds32columns\"" 1)
require_token_count("MQ.2 C1 field-large literal owner"
  mq2_netlist_csv "\"CSVfieldexceeds65536bytes\"" 1)
require_token_count("MQ.2 C1 too-many-rows literal owner"
  mq2_netlist_csv "\"liionpackCSVexceeds100000rows\"" 1)
require_token_count("MQ.2 C1 too-large owner and consumers"
  mq2_netlist_csv "msg_too_large" 3)
require_token_count("MQ.2 C1 too-wide owner and consumer"
  mq2_netlist_csv "msg_too_wide" 2)
require_token_count("MQ.2 C1 field-large owner and consumers"
  mq2_netlist_csv "msg_field_large" 3)
require_token_count("MQ.2 C1 too-many-rows owner and consumer"
  mq2_netlist_csv "msg_too_many_rows" 2)

require_token_count("MQ.2 C1 direct too-large consumer"
  mq2_netlist_csv "diagnostic.message=msg_too_large;" 1)
require_token_count("MQ.2 C1 load too-large consumer"
  mq2_netlist_csv
  "diagnostic={.message=std::string{msg_too_large}};"
  1)
require_token_count("MQ.2 C1 bounded-file byte limit"
  mq2_netlist_csv
  "detail::readBoundedFile(path,max_csv_bytes,contents)"
  1)
require_token_count("MQ.2 C1 too-wide consumer"
  mq2_netlist_csv "returnfail(msg_too_wide);" 1)
require_token_count("MQ.2 C1 quoted and unquoted field consumers"
  mq2_netlist_csv "returnfail(msg_field_large);" 2)
require_token_count("MQ.2 C1 too-many-rows consumer"
  mq2_netlist_csv "diagnostic.message=msg_too_many_rows;" 1)

require_token_count("MQ.2 C1 exact column guard"
  mq2_netlist_csv
  "if(fields.size()>=max_csv_columns)returnfail(msg_too_wide);"
  1)
require_token_count("MQ.2 C1 exact field guards"
  mq2_netlist_csv
  "if(field.size()>max_csv_field_bytes)returnfail(msg_field_large);"
  2)
require_token_count("MQ.2 C1 exact direct-size guard"
  mq2_netlist_csv
  "if(csv.size()>max_csv_bytes){diagnostic.message=msg_too_large;returnslide::Status::Invalid_parameters;}"
  1)
require_token_count("MQ.2 C1 exact data-row guard"
  mq2_netlist_csv
  "if(data_rows>max_csv_rows){diagnostic.row=reader.row()-1;diagnostic.message=msg_too_many_rows;returnslide::Status::Invalid_parameters;}"
  1)

set(mq2_netlist_parse_entry_begin
  "slide::StatusparseLiionpackNetlistCsv(")
set(mq2_netlist_parse_entry_end "CsvReaderreader{csv,diagnostic};")
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_parse_entry_begin}"
  mq2_netlist_parse_entry_begin_position)
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_parse_entry_end}"
  mq2_netlist_parse_entry_end_position)
if(mq2_netlist_parse_entry_begin_position EQUAL -1
   OR mq2_netlist_parse_entry_end_position EQUAL -1
   OR mq2_netlist_parse_entry_end_position
      LESS_EQUAL mq2_netlist_parse_entry_begin_position)
  message(FATAL_ERROR "MQ.2 C1 parse-entry slice is missing")
endif()
math(EXPR mq2_netlist_parse_entry_length
  "${mq2_netlist_parse_entry_end_position} - ${mq2_netlist_parse_entry_begin_position}")
string(SUBSTRING "${mq2_netlist_csv}"
  ${mq2_netlist_parse_entry_begin_position}
  ${mq2_netlist_parse_entry_length}
  mq2_netlist_parse_entry)
require_ordered_tokens("MQ.2 C1 diagnostic reset precedes direct limit"
  mq2_netlist_parse_entry
  "try{diagnostic={};"
  "if(csv.empty())"
  "if(csv.size()>max_csv_bytes)"
  "diagnostic.message=msg_too_large;")

require_token_count("MQ.2 C1 non-reader offset reset census"
  mq2_netlist_csv "diagnostic.offset" 1)
require_token_count("MQ.2 C1 reader offset publication census"
  mq2_netlist_csv "diagnostic_.offset" 1)
require_token_count("MQ.2 C1 exact allocation-failure diagnostic owner"
  mq2_netlist_csv
  "slide::StatusallocationFailure(NetlistCsvDiagnostic&diagnostic,std::string_viewmessage)noexcept{diagnostic.row=0;diagnostic.offset=0;assignDiagnosticNoThrow(diagnostic.message,message);returnslide::Status::Numerical_failure;}"
  1)
require_token_count("MQ.2 C1 exact reader failure owner"
  mq2_netlist_csv
  "boolfail(std::string_viewmessage){diagnostic_.row=row_;diagnostic_.offset=cursor_;diagnostic_.message.assign(message);returnfalse;}"
  1)

set(mq2_netlist_quoted_begin "if(source_[cursor_]=='\"'){")
set(mq2_netlist_unquoted_begin
  "}else{while(cursor_<source_.size()")
set(mq2_netlist_fields_publish "fields.push_back(std::move(field));")
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_quoted_begin}"
  mq2_netlist_quoted_begin_position)
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_unquoted_begin}"
  mq2_netlist_unquoted_begin_position)
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_fields_publish}"
  mq2_netlist_fields_publish_position)
if(mq2_netlist_quoted_begin_position EQUAL -1
   OR mq2_netlist_unquoted_begin_position EQUAL -1
   OR mq2_netlist_fields_publish_position EQUAL -1
   OR mq2_netlist_unquoted_begin_position
      LESS_EQUAL mq2_netlist_quoted_begin_position
   OR mq2_netlist_fields_publish_position
      LESS_EQUAL mq2_netlist_unquoted_begin_position)
  message(FATAL_ERROR "MQ.2 C1 quoted/unquoted reader slices are missing")
endif()
math(EXPR mq2_netlist_quoted_length
  "${mq2_netlist_unquoted_begin_position} - ${mq2_netlist_quoted_begin_position}")
math(EXPR mq2_netlist_unquoted_length
  "${mq2_netlist_fields_publish_position} - ${mq2_netlist_unquoted_begin_position}")
string(SUBSTRING "${mq2_netlist_csv}"
  ${mq2_netlist_quoted_begin_position}
  ${mq2_netlist_quoted_length}
  mq2_netlist_quoted)
string(SUBSTRING "${mq2_netlist_csv}"
  ${mq2_netlist_unquoted_begin_position}
  ${mq2_netlist_unquoted_length}
  mq2_netlist_unquoted)
require_token_count("MQ.2 C1 quoted field-limit consumer"
  mq2_netlist_quoted
  "if(field.size()>max_csv_field_bytes)returnfail(msg_field_large);"
  1)
require_token_count("MQ.2 C1 unquoted field-limit consumer"
  mq2_netlist_unquoted
  "if(field.size()>max_csv_field_bytes)returnfail(msg_field_large);"
  1)

require_token_count("MQ.2 C1 exact ASCII digit owner"
  mq2_netlist_csv
  "[[nodiscard]]constexprboolisAsciiDigit(charvalue)noexcept{returnvalue>='0'&&value<='9';}"
  1)
require_token_count("MQ.2 C1 exact digit scanner owner"
  mq2_netlist_csv
  "constexprstd::size_tscanDigits(std::string_viewtext,std::size_t&cursor)noexcept{conststd::size_tbegin=cursor;while(cursor<text.size()&&isAsciiDigit(text[cursor]))++cursor;returncursor-begin;}"
  1)
require_token_count("MQ.2 C1 digit predicate census"
  mq2_netlist_csv "isAsciiDigit(" 4)
require_token_count("MQ.2 C1 digit scanner census"
  mq2_netlist_csv "scanDigits(" 4)
require_token_count("MQ.2 C1 scanner return may be discarded"
  mq2_netlist_csv
  "[[nodiscard]]constexprstd::size_tscanDigits(" 0)
require_token_count("MQ.2 C1 no raw zero-to-nine cursor tests"
  mq2_netlist_csv
  "text[cursor]>='0'&&text[cursor]<='9'" 0)
require_token_count("MQ.2 C1 locale-independent digit rationale"
  mq2_netlist_csv_with_comments
  "//std::isdigitislocale-sensitiveandrequiresanunsigned-chardomain.[[nodiscard]]constexprboolisAsciiDigit"
  1)
require_token_count("MQ.2 C1 explicit descriptor digit conversion"
  mq2_netlist_csv "isAsciiDigit(static_cast<char>(value))" 1)

set(mq2_netlist_descriptor_begin "boolvalidDescriptor(")
set(mq2_netlist_descriptor_end "failSemantic(")
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_descriptor_begin}"
  mq2_netlist_descriptor_begin_position)
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_descriptor_end}"
  mq2_netlist_descriptor_end_position)
if(mq2_netlist_descriptor_begin_position EQUAL -1
   OR mq2_netlist_descriptor_end_position EQUAL -1
   OR mq2_netlist_descriptor_end_position
      LESS_EQUAL mq2_netlist_descriptor_begin_position)
  message(FATAL_ERROR "MQ.2 C1 validDescriptor slice is missing")
endif()
math(EXPR mq2_netlist_descriptor_length
  "${mq2_netlist_descriptor_end_position} - ${mq2_netlist_descriptor_begin_position}")
string(SUBSTRING "${mq2_netlist_csv}"
  ${mq2_netlist_descriptor_begin_position}
  ${mq2_netlist_descriptor_length}
  mq2_netlist_descriptor)
require_token_count("MQ.2 C1 descriptor digit consumer"
  mq2_netlist_descriptor
  "isAsciiDigit(static_cast<char>(value))" 1)
require_token_count("MQ.2 C1 descriptor raw digit test removed"
  mq2_netlist_descriptor
  "value>='0'&&value<='9'" 0)

set(mq2_netlist_value_begin "boolparseValue(")
set(mq2_netlist_value_end "boolvalidDescriptor(")
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_value_begin}"
  mq2_netlist_value_begin_position)
string(FIND "${mq2_netlist_csv}" "${mq2_netlist_value_end}"
  mq2_netlist_value_end_position)
if(mq2_netlist_value_begin_position EQUAL -1
   OR mq2_netlist_value_end_position EQUAL -1
   OR mq2_netlist_value_end_position LESS_EQUAL mq2_netlist_value_begin_position)
  message(FATAL_ERROR "MQ.2 C1 parseValue slice is missing")
endif()
math(EXPR mq2_netlist_value_length
  "${mq2_netlist_value_end_position} - ${mq2_netlist_value_begin_position}")
string(SUBSTRING "${mq2_netlist_csv}"
  ${mq2_netlist_value_begin_position}
  ${mq2_netlist_value_length}
  mq2_netlist_value)
require_token_count("MQ.2 C1 parseValue scanner consumers"
  mq2_netlist_value "scanDigits(" 3)
require_token_count("MQ.2 C1 parseValue direct digit consumer"
  mq2_netlist_value "isAsciiDigit(" 1)
require_token_count("MQ.2 C1 minus-only leading sign"
  mq2_netlist_value
  "if(text[cursor]=='-'&&++cursor==text.size())returnfalse;" 1)
require_token_count("MQ.2 C1 plus sign remains exponent-only"
  mq2_netlist_value "'+'" 1)
require_token_count("MQ.2 C1 zero integer branch"
  mq2_netlist_value "if(text[cursor]=='0'){" 1)
require_token_count("MQ.2 C1 no-leading-zero rejection"
  mq2_netlist_value
  "if(cursor<text.size()&&isAsciiDigit(text[cursor]))returnfalse;" 1)
require_token_count("MQ.2 C1 nonzero integer branch"
  mq2_netlist_value
  "elseif(text[cursor]>='1'&&text[cursor]<='9'){" 1)
require_token_count("MQ.2 C1 known-nonzero digit scan"
  mq2_netlist_value "scanDigits(text,cursor);" 1)
require_token_count("MQ.2 C1 nonempty fraction/exponent scans"
  mq2_netlist_value
  "if(scanDigits(text,cursor)==0)returnfalse;" 2)
require_token_count("MQ.2 C1 complete value consumption"
  mq2_netlist_value "if(cursor!=text.size())returnfalse;" 1)

require_token_count("MQ.2 C1 exact semantic Status owner"
  mq2_netlist_csv
  "[[nodiscard]]slide::StatusfailSemantic(NetlistCsvDiagnostic&diagnostic,std::size_trow,std::string_viewmessage){diagnostic.row=row;diagnostic.message.assign(message);returnslide::Status::Invalid_parameters;}"
  1)
require_token_count("MQ.2 C1 semantic owner and consumers"
  mq2_netlist_csv "failSemantic(" 9)
require_token_count("MQ.2 C1 direct semantic returns"
  mq2_netlist_csv "returnfailSemantic(" 8)
require_token_count("MQ.2 C1 row semantic returns"
  mq2_netlist_csv "returnfailSemantic(diagnostic,row," 7)
require_token_count("MQ.2 C1 element-row semantic return"
  mq2_netlist_csv "returnfailSemantic(diagnostic,element.row," 1)
require_token_count("MQ.2 C1 direct Invalid-parameters census"
  mq2_netlist_csv "returnslide::Status::Invalid_parameters;" 12)

require_token_count("MQ.2 C1 offset-member contract"
  mq2_netlist_csv_header_with_comments
  "std::size_toffset{};//!<sourcecursorat/immediatelyafterareader-detected//!<syntax/limitfailure;zerowhenunavailablestd::stringmessage{};"
  1)
require_token_count("MQ.2 C1 public limit contract"
  mq2_netlist_csv_header_with_comments
  "Inputisboundedto4MiB,100,000rows,32columns,and65,536bytes/field."
  1)
require_token_count("MQ.2 C1 stale unrestricted offset removed"
  mq2_netlist_csv_header_with_comments
  "sourcebyteatorimmediatelyafterthefailure" 0)

# MQ.2 R1: recording format facts, eager-reader indexing, schema order, and
# widened row addressing each have one explicit owner. Missing future headers
# are treated as empty only so old production reaches the registered first-red
# common-owner boundary instead of failing inside file(READ).
load_compact("src/core/Recorder.cpp" mq2_r1_recorder)
load_compact("src/core/RecordingFormat.cpp" mq2_r1_recording_format)
load_compact("src/core/AsyncRecorder.cpp" mq2_r1_async_recorder)
load_compact("src/core/AsyncRecordingCodec.cpp" mq2_r1_async_codec)
load_compact("src/core/detail/AsyncRecordingFormat.hpp"
  mq2_r1_async_format)

set(mq2_r1_common_path
  "${SLIDE_SOURCE_DIR}/src/core/detail/RecordingFormatCommon.hpp")
if(EXISTS "${mq2_r1_common_path}")
  load_compact("src/core/detail/RecordingFormatCommon.hpp" mq2_r1_common)
  file(READ "${mq2_r1_common_path}" mq2_r1_common_with_comments)
  string(REGEX REPLACE "[ \t\r\n]" ""
    mq2_r1_common_with_comments "${mq2_r1_common_with_comments}")
else()
  set(mq2_r1_common "")
  set(mq2_r1_common_with_comments "")
endif()

set(mq2_r1_snapshot_path
  "${SLIDE_SOURCE_DIR}/src/core/detail/SnapshotIndexing.hpp")
if(EXISTS "${mq2_r1_snapshot_path}")
  load_compact("src/core/detail/SnapshotIndexing.hpp" mq2_r1_snapshot)
  file(READ "${mq2_r1_snapshot_path}" mq2_r1_snapshot_with_comments)
  string(REGEX REPLACE "[ \t\r\n]" ""
    mq2_r1_snapshot_with_comments "${mq2_r1_snapshot_with_comments}")
else()
  set(mq2_r1_snapshot "")
  set(mq2_r1_snapshot_with_comments "")
endif()

# This is deliberately the first R1 assertion: old production must stop here
# at 0/1, after every preceding structural suite has passed.
require_token_count("MQ.2 R1 RecordingFormatCommon owner"
  mq2_r1_common
  "inlineconstexprstd::uint32_tendian_marker=0x01020304U;"
  1)
require_token_count("MQ.2 R1 exact CRC owner"
  mq2_r1_common
  "inlinestd::uint32_tcrc32(std::span<conststd::byte>bytes){std::uint32_tcrc=0xffffffffU;for(constautobyte:bytes){crc^=std::to_integer<std::uint8_t>(byte);for(intbit=0;bit<8;++bit)crc=(crc>>1U)^(0xedb88320U&(0U-(crc&1U)));}return~crc;}"
  1)
require_token_count("MQ.2 R1 generic header CRC owner"
  mq2_r1_common
  "template<classHeader>[[nodiscard]]inlinestd::uint32_theaderCrc(Headerheader){static_assert(std::is_trivially_copyable_v<Header>);header.header_crc32=0;returncrc32(std::as_bytes(std::span{&header,1}));}"
  1)
require_token_count("MQ.2 R1 allocation-failure owner"
  mq2_r1_common
  "inlineslide::StatusallocationFailureStatus()noexcept{returnslide::Status::Numerical_failure;}"
  1)
require_token_count("MQ.2 R1 common Status include"
  mq2_r1_common "#include\"../../types/Status.hpp\"" 1)
require_token_count("MQ.2 R1 common byte include"
  mq2_r1_common "#include<cstddef>" 1)
require_token_count("MQ.2 R1 common integer include"
  mq2_r1_common "#include<cstdint>" 1)
require_token_count("MQ.2 R1 common span include"
  mq2_r1_common "#include<span>" 1)
require_token_count("MQ.2 R1 common type-trait include"
  mq2_r1_common "#include<type_traits>" 1)
require_token_count("MQ.2 R1 common owns no major version"
  mq2_r1_common "format_major" 0)
require_token_count("MQ.2 R1 common owns no minor version"
  mq2_r1_common "format_minor" 0)
string(FIND "${mq2_r1_common_with_comments}"
  "/***@fileRecordingFormatCommon.hpp" mq2_r1_common_leading_position)
string(FIND "${mq2_r1_common_with_comments}"
  "#pragmaonce" mq2_r1_common_pragma_position)
if(NOT mq2_r1_common_leading_position EQUAL 0
   OR mq2_r1_common_pragma_position EQUAL -1)
  message(FATAL_ERROR
    "MQ.2 R1 RecordingFormatCommon contract is not the leading file contract")
endif()
string(SUBSTRING "${mq2_r1_common_with_comments}"
  0 ${mq2_r1_common_pragma_position} mq2_r1_common_leading_contract)
require_ordered_tokens("MQ.2 R1 common leading contract"
  mq2_r1_common_leading_contract
  "/***@fileRecordingFormatCommon.hpp"
  "@briefSharedchecksum,byte-order,andallocation-failurefactsforrecording."
  "Owns:`endian_marker`,`crc32`,`headerCrc`,and`allocationFailureStatus`."
  "ImplementsPLAN.mdsection3.7.Cold:file-formatchecksandallocationtranslationonly."
  "@surfaceinternal*/")

set(mq2_r1_format_sources
  "${mq2_r1_common}${mq2_r1_recording_format}${mq2_r1_recorder}${mq2_r1_async_format}${mq2_r1_async_recorder}${mq2_r1_async_codec}${mq2_r1_snapshot}")
require_token_count("MQ.2 R1 global crc32 census"
  mq2_r1_format_sources "crc32(" 6)
require_token_count("MQ.2 R1 global endian census"
  mq2_r1_format_sources "endian_marker" 8)
require_token_count("MQ.2 R1 global allocation mapper census"
  mq2_r1_format_sources "allocationFailureStatus(" 13)
require_token_count("MQ.2 R1 global header CRC census"
  mq2_r1_format_sources "headerCrc(" 7)
require_token_count("MQ.2 R1 file-header wrapper removed"
  mq2_r1_format_sources "fileHeaderCrc(" 0)
require_token_count("MQ.2 R1 block-header wrapper removed"
  mq2_r1_format_sources "blockHeaderCrc(" 0)
require_token_count("MQ.2 R1 direct common include census"
  mq2_r1_format_sources "RecordingFormatCommon.hpp" 3)
require_token_count("MQ.2 R1 direct snapshot include census"
  mq2_r1_format_sources "SnapshotIndexing.hpp" 2)

# Per-file counts prevent a globally balanced migration from moving a format
# fact to the wrong reader/writer or allocation boundary.
require_token_count("MQ.2 R1 common crc32 census"
  mq2_r1_common "crc32(" 2)
require_token_count("MQ.2 R1 common endian census"
  mq2_r1_common "endian_marker" 1)
require_token_count("MQ.2 R1 common allocation census"
  mq2_r1_common "allocationFailureStatus(" 1)
require_token_count("MQ.2 R1 common header CRC census"
  mq2_r1_common "headerCrc(" 1)

require_token_count("MQ.2 R1 synchronous crc32 census"
  mq2_r1_recording_format "crc32(" 0)
require_token_count("MQ.2 R1 synchronous endian census"
  mq2_r1_recording_format "endian_marker" 3)
require_token_count("MQ.2 R1 synchronous allocation census"
  mq2_r1_recording_format "allocationFailureStatus(" 2)
require_token_count("MQ.2 R1 synchronous header CRC census"
  mq2_r1_recording_format "headerCrc(" 2)

require_token_count("MQ.2 R1 Recorder crc32 census"
  mq2_r1_recorder "crc32(" 0)
require_token_count("MQ.2 R1 Recorder endian census"
  mq2_r1_recorder "endian_marker" 0)
require_token_count("MQ.2 R1 Recorder allocation census"
  mq2_r1_recorder "allocationFailureStatus(" 6)
require_token_count("MQ.2 R1 Recorder header CRC census"
  mq2_r1_recorder "headerCrc(" 0)

foreach(token IN ITEMS "crc32(" "endian_marker"
    "allocationFailureStatus(" "headerCrc(")
  require_token_count("MQ.2 R1 async format has no common fact"
    mq2_r1_async_format "${token}" 0)
endforeach()

require_token_count("MQ.2 R1 async writer crc32 census"
  mq2_r1_async_recorder "crc32(" 2)
require_token_count("MQ.2 R1 async writer endian census"
  mq2_r1_async_recorder "endian_marker" 2)
require_token_count("MQ.2 R1 async writer allocation census"
  mq2_r1_async_recorder "allocationFailureStatus(" 2)
require_token_count("MQ.2 R1 async writer header CRC census"
  mq2_r1_async_recorder "headerCrc(" 2)

require_token_count("MQ.2 R1 async reader crc32 census"
  mq2_r1_async_codec "crc32(" 2)
require_token_count("MQ.2 R1 async reader endian census"
  mq2_r1_async_codec "endian_marker" 2)
require_token_count("MQ.2 R1 async reader allocation census"
  mq2_r1_async_codec "allocationFailureStatus(" 2)
require_token_count("MQ.2 R1 async reader header CRC census"
  mq2_r1_async_codec "headerCrc(" 2)

foreach(token IN ITEMS "crc32(" "endian_marker"
    "allocationFailureStatus(" "headerCrc(")
  require_token_count("MQ.2 R1 snapshot header has no format fact"
    mq2_r1_snapshot "${token}" 0)
endforeach()

require_token_count("MQ.2 R1 Recorder direct common include"
  mq2_r1_recorder
  "#include\"detail/RecordingFormatCommon.hpp\"" 1)
require_token_count("MQ.2 R1 synchronous direct common include"
  mq2_r1_recording_format
  "#include\"detail/RecordingFormatCommon.hpp\"" 1)
require_token_count("MQ.2 R1 async format direct common include"
  mq2_r1_async_format
  "#include\"RecordingFormatCommon.hpp\"" 1)
require_token_count("MQ.2 R1 async writer no direct common include"
  mq2_r1_async_recorder "RecordingFormatCommon.hpp" 0)
require_token_count("MQ.2 R1 async reader no direct common include"
  mq2_r1_async_codec "RecordingFormatCommon.hpp" 0)

require_token_count("MQ.2 R1 synchronous major owner"
  mq2_r1_recording_format
  "constexprstd::uint16_tformat_major=1;" 1)
require_token_count("MQ.2 R1 synchronous minor owner"
  mq2_r1_recording_format
  "constexprstd::uint16_tformat_minor=0;" 1)
require_token_count("MQ.2 R1 async major owner"
  mq2_r1_async_format
  "constexprstd::uint16_tformat_major=1;" 1)
require_token_count("MQ.2 R1 async minor owner"
  mq2_r1_async_format
  "constexprstd::uint16_tformat_minor=0;" 1)
require_token_count("MQ.2 R1 synchronous minor comparison"
  mq2_r1_recording_format "header.minor>format_minor" 1)
require_token_count("MQ.2 R1 async minor comparison"
  mq2_r1_async_codec "header.minor!=format_minor" 1)
require_token_count("MQ.2 R2 legal zero CRC has no stored-value guard"
  mq2_r1_recording_format "header.header_crc32==0" 0)
forbid_tokens("MQ.2 R2 no alternate zero-valued CRC special case"
  mq2_r1_recording_format
  "0==header.header_crc32"
  "!header.header_crc32"
  "headerCrc(header)==0"
  "headerCrc(header)!=0")
require_token_count("MQ.2 R1 synchronous writer header CRC"
  mq2_r1_recording_format
  "header.header_crc32=headerCrc(header);" 1)
require_token_count("MQ.2 R1 synchronous reader header CRC"
  mq2_r1_recording_format
  "headerCrc(header)!=header.header_crc32" 1)
require_token_count("MQ.2 R1 synchronous exchange removed"
  mq2_r1_recording_format "std::exchange(" 0)
require_token_count("MQ.2 R1 async block writer header CRC"
  mq2_r1_async_recorder
  "header.header_crc32=headerCrc(header);" 2)
require_token_count("MQ.2 R1 async file reader header CRC"
  mq2_r1_async_codec
  "header.header_crc32!=headerCrc(header)" 1)
require_token_count("MQ.2 R1 async block reader header CRC"
  mq2_r1_async_codec
  "block.header_crc32!=headerCrc(block)" 1)

require_token_count("MQ.2 R1 raw and payload writer roles"
  mq2_r1_async_recorder
  ".raw_crc32=crc32(raw),.payload_crc32=crc32(payload)," 1)
require_ordered_tokens("MQ.2 R1 raw and payload reader roles"
  mq2_r1_async_codec
  "autopayload_view=std::span<std::byte>{payload}.first("
  "constautopayload_status=readExact(input,payload_view,remaining)"
  "if(payload_status!=slide::Status::Success)"
  "returnpayload_status"
  "if(block.payload_crc32!=crc32(payload_view))"
  "constautounshuffle_status=byteUnshuffle(shuffled,raw,sizeof(real_t))"
  "if(block.raw_crc32!=crc32(raw))")

set(mq2_r1_drain_begin
  "slide::StatusAsyncRecorder::drainSlot(Slot&slot){")
set(mq2_r1_drain_end "voidAsyncRecorder::drainLoop(){")
string(FIND "${mq2_r1_async_recorder}" "${mq2_r1_drain_begin}"
  mq2_r1_drain_begin_position)
string(FIND "${mq2_r1_async_recorder}" "${mq2_r1_drain_end}"
  mq2_r1_drain_end_position)
if(mq2_r1_drain_begin_position EQUAL -1
   OR mq2_r1_drain_end_position EQUAL -1
   OR mq2_r1_drain_end_position LESS_EQUAL mq2_r1_drain_begin_position)
  message(FATAL_ERROR "MQ.2 R1 async block-writer slice is missing")
endif()
math(EXPR mq2_r1_drain_length
  "${mq2_r1_drain_end_position} - ${mq2_r1_drain_begin_position}")
string(SUBSTRING "${mq2_r1_async_recorder}"
  ${mq2_r1_drain_begin_position} ${mq2_r1_drain_length}
  mq2_r1_drain_slice)
require_token_count("MQ.2 R1 block-writer header CRC role"
  mq2_r1_drain_slice "header.header_crc32=headerCrc(header);" 1)
require_token_count("MQ.2 R1 block-writer one header CRC"
  mq2_r1_drain_slice "headerCrc(" 1)

set(mq2_r1_finalize_begin
  "slide::StatusAsyncRecorder::finalizeFile(){")
set(mq2_r1_finalize_end "slide::StatusAsyncRecorder::finish(){")
string(FIND "${mq2_r1_async_recorder}" "${mq2_r1_finalize_begin}"
  mq2_r1_finalize_begin_position)
string(FIND "${mq2_r1_async_recorder}" "${mq2_r1_finalize_end}"
  mq2_r1_finalize_end_position)
if(mq2_r1_finalize_begin_position EQUAL -1
   OR mq2_r1_finalize_end_position EQUAL -1
   OR mq2_r1_finalize_end_position
      LESS_EQUAL mq2_r1_finalize_begin_position)
  message(FATAL_ERROR "MQ.2 R1 async file-writer slice is missing")
endif()
math(EXPR mq2_r1_finalize_length
  "${mq2_r1_finalize_end_position} - ${mq2_r1_finalize_begin_position}")
string(SUBSTRING "${mq2_r1_async_recorder}"
  ${mq2_r1_finalize_begin_position} ${mq2_r1_finalize_length}
  mq2_r1_finalize_slice)
require_token_count("MQ.2 R1 file-writer header CRC role"
  mq2_r1_finalize_slice "header.header_crc32=headerCrc(header);" 1)
require_token_count("MQ.2 R1 file-writer one header CRC"
  mq2_r1_finalize_slice "headerCrc(" 1)

# The eager SoA readers share one indexing body, retain their caller-specific
# assertions, and pass every same-typed argument in the frozen order.
require_token_count("MQ.2 R1 exact eager snapshot owner"
  mq2_r1_snapshot
  "[[nodiscard]]inlineSnapshotViewsnapshotView(std::size_tindex,std::size_tlanes,std::size_tstate_values,std::span<conststd::uint64_t>steps,std::span<constreal_t>times,std::span<constreal_t>currents,std::span<constreal_t>states)noexcept{return{.accepted_step=steps[index],.time=times[index],.current_density=currents.subspan(index*lanes,lanes),.state=states.subspan(index*state_values,state_values)};}"
  1)
require_token_count("MQ.2 R1 eager snapshot global census"
  mq2_r1_format_sources "snapshotView(" 3)
require_token_count("MQ.2 R1 Recorder snapshot include"
  mq2_r1_recorder "#include\"detail/SnapshotIndexing.hpp\"" 1)
require_token_count("MQ.2 R1 compressed snapshot include"
  mq2_r1_async_codec "#include\"detail/SnapshotIndexing.hpp\"" 1)
require_token_count("MQ.2 R1 synchronous mapped reader has no eager include"
  mq2_r1_recording_format "SnapshotIndexing.hpp" 0)
require_token_count("MQ.2 R1 eager owner includes public view"
  mq2_r1_snapshot "#include\"../Recorder.hpp\"" 1)
require_token_count("MQ.2 R1 eager owner cstddef include"
  mq2_r1_snapshot "#include<cstddef>" 1)
require_token_count("MQ.2 R1 eager owner cstdint include"
  mq2_r1_snapshot "#include<cstdint>" 1)
require_token_count("MQ.2 R1 eager owner span include"
  mq2_r1_snapshot "#include<span>" 1)
require_token_count("MQ.2 R1 Recorder assert include"
  mq2_r1_recorder "#include<cassert>" 1)
require_token_count("MQ.2 R1 Recorder caller assertion"
  mq2_r1_recorder "assert(index<count_);" 1)
require_token_count("MQ.2 R1 compressed caller assertion"
  mq2_r1_async_codec "assert(valid_&&index<steps_.size());" 1)
require_token_count("MQ.2 R1 Recorder ordered eager call"
  mq2_r1_recorder
  "returndetail::snapshotView(index,static_cast<std::size_t>(lanes_),state_values_,accepted_steps_,times_,current_density_,states_);"
  1)
require_token_count("MQ.2 R1 compressed ordered eager call"
  mq2_r1_async_codec
  "returndetail::snapshotView(index,static_cast<std::size_t>(lanes_),state_values_,steps_,times_,currents_,states_);"
  1)
require_token_count("MQ.2 R1 mapped snapshot owner retained"
  mq2_r1_recording_format
  "SnapshotViewBinaryRecording::snapshot(std::size_tindex)const" 1)
require_token_count("MQ.2 R1 mapped snapshot table retained"
  mq2_r1_recording_format
  "constauto*table=mapping_->data+table_offset_;" 1)
require_token_count("MQ.2 R1 mapped snapshot typed views retained"
  mq2_r1_recording_format
  "reinterpret_cast<constreal_t*>(cursor)" 2)

# The schema has one canonical order. CSV joins it directly; the optional
# Parquet branch advances one monotone cursor through the same five append
# families and proves that every name was consumed.
require_token_count("MQ.2 R1 helper anonymous namespace opening"
  mq2_r1_recorder
  "namespace{[[nodiscard]]std::vector<std::string>recorderColumnNames"
  1)
require_token_count("MQ.2 R1 helper anonymous namespace closing"
  mq2_r1_recorder
  "static_cast<std::size_t>(lanes));}}slide::StatusRecorder::configure("
  1)
require_token_count("MQ.2 R1 schema owner"
  mq2_r1_recorder
  "[[nodiscard]]std::vector<std::string>recorderColumnNames(introws,intlanes)"
  1)
require_token_count("MQ.2 R1 exact canonical schema owner"
  mq2_r1_recorder
  "[[nodiscard]]std::vector<std::string>recorderColumnNames(introws,intlanes){std::vector<std::string>names{\"accepted_step\",\"time_s\"};for(intlane=0;lane<lanes;++lane)names.push_back(\"current_density_lane\"+std::to_string(lane)+\"_A_m2\");for(introw=0;row<rows;++row)for(intlane=0;lane<lanes;++lane)names.push_back(\"state_r\"+std::to_string(row)+\"_lane\"+std::to_string(lane));for(intlane=0;lane<lanes;++lane)names.push_back(\"terminal_voltage_lane\"+std::to_string(lane)+\"_V\");returnnames;}"
  1)
require_token_count("MQ.2 R1 no noexcept default schema vector"
  mq2_r1_recorder "std::vector<std::string>names;" 0)
require_token_count("MQ.2 R1 schema owner and consumers"
  mq2_r1_recorder "recorderColumnNames(" 3)
require_token_count("MQ.2 R1 exact schema call arguments"
  mq2_r1_recorder "recorderColumnNames(rows_,lanes_)" 2)
foreach(token IN ITEMS
    "\"accepted_step\"" "\"time_s\"" "\"current_density_lane\""
    "\"_A_m2\"" "\"state_r\"" "\"_lane\""
    "\"terminal_voltage_lane\"" "\"_V\"")
  require_token_count("MQ.2 R1 once-only schema fragment"
    mq2_r1_recorder "${token}" 1)
endforeach()
require_ordered_tokens("MQ.2 R1 canonical schema order"
  mq2_r1_recorder
  "std::vector<std::string>names{\"accepted_step\",\"time_s\"}"
  "\"current_density_lane\"+std::to_string(lane)+\"_A_m2\""
  "\"state_r\"+std::to_string(row)+\"_lane\"+std::to_string(lane)"
  "\"terminal_voltage_lane\"+std::to_string(lane)+\"_V\"")
require_token_count("MQ.2 R1 exact CSV schema join"
  mq2_r1_recorder
  "constautonames=recorderColumnNames(rows_,lanes_);for(std::size_tindex=0;index<names.size();++index)output<<(index==0?\"\":\",\")<<names[index];"
  1)
require_token_count("MQ.2 R1 one Parquet column cursor"
  mq2_r1_recorder "std::size_tcolumn{};" 1)
require_token_count("MQ.2 R1 Parquet schema precedes cursor"
  mq2_r1_recorder
  "constautonames=recorderColumnNames(rows_,lanes_);std::size_tcolumn{};"
  1)
require_token_count("MQ.2 R1 Parquet name consumption"
  mq2_r1_recorder "names[column++]" 5)
require_token_count("MQ.2 R1 only five cursor advances"
  mq2_r1_recorder "column++" 5)
require_token_count("MQ.2 R1 Parquet uint64 name consumption"
  mq2_r1_recorder "append_uint64(names[column++]" 1)
require_token_count("MQ.2 R1 Parquet double name consumption"
  mq2_r1_recorder "append_double(names[column++]" 4)
require_token_count("MQ.2 R1 complete Parquet schema consumption"
  mq2_r1_recorder "assert(column==names.size());" 1)
require_ordered_tokens("MQ.2 R1 Parquet name-to-data associations"
  mq2_r1_recorder
  "std::size_tcolumn{}"
  "append_uint64(names[column++],[&](std::size_ti){returnaccepted_steps_[i]"
  "append_double(names[column++],[&](std::size_ti){returntimes_[i]"
  "append_double(names[column++],[&](std::size_ti){returnsnapshot(i).current_density[static_cast<std::size_t>(lane)]"
  "append_double(names[column++],[&](std::size_ti){returnsnapshotRow(snapshot(i),row,rows_,stride_,lanes_)[static_cast<std::size_t>(lane)]"
  "append_double(names[column++],[&](std::size_ti){returnvoltage[static_cast<std::size_t>(lane)][i]"
  "assert(column==names.size())"
  "constautotable=arrow::Table::Make(")
require_token_count("MQ.2 R1 Parquet step association"
  mq2_r1_recorder
  "append_uint64(names[column++],[&](std::size_ti){returnaccepted_steps_[i];})"
  1)
require_token_count("MQ.2 R1 Parquet time association"
  mq2_r1_recorder
  "append_double(names[column++],[&](std::size_ti){returntimes_[i];})"
  1)
require_token_count("MQ.2 R1 Parquet current loop and association"
  mq2_r1_recorder
  "for(intlane=0;lane<lanes_;++lane){if(!append_double(names[column++],[&](std::size_ti){returnsnapshot(i).current_density[static_cast<std::size_t>(lane)];}))returnslide::Status::Numerical_failure;}"
  1)
require_token_count("MQ.2 R1 Parquet state loops and association"
  mq2_r1_recorder
  "for(introw=0;row<rows_;++row)for(intlane=0;lane<lanes_;++lane)if(!append_double(names[column++],[&](std::size_ti){returnsnapshotRow(snapshot(i),row,rows_,stride_,lanes_)[static_cast<std::size_t>(lane)];}))returnslide::Status::Numerical_failure;"
  1)
require_token_count("MQ.2 R1 Parquet voltage loop and association"
  mq2_r1_recorder
  "for(intlane=0;lane<lanes_;++lane)if(!append_double(names[column++],[&](std::size_ti){returnvoltage[static_cast<std::size_t>(lane)][i];}))returnslide::Status::Numerical_failure;"
  1)
forbid_tokens("MQ.2 R1 no alternate Parquet cursor mutation"
  mq2_r1_recorder
  "++column" "--column" "column--"
  "column+=" "column-=" "column*=" "column/=" "column%="
  "column^=" "column|=" "column&=" "column<<=" "column>>=")
require_token_count("MQ.2 R1 only cursor equality spelling"
  mq2_r1_recorder "column=" 1)
require_token_count("MQ.2 R1 stale CSV schema literal removed"
  mq2_r1_recorder "output<<\"accepted_step,time_s\"" 0)
require_token_count("MQ.2 R1 stale Parquet schema literal removed"
  mq2_r1_recorder "append_uint64(\"accepted_step\"" 0)

# Cast both multiplicands before row addressing. The exact consumer calls pin
# rows/stride/lanes against same-typed argument swaps.
require_token_count("MQ.2 R1 exact widened row owner"
  mq2_r1_recorder
  "[[nodiscard]]std::span<constreal_t>snapshotRow(constSnapshotView&snapshot,introw,introws,intstride,intlanes)noexcept{assert(0<=row&&row<rows&&0<lanes&&lanes<=stride);returnsnapshot.state.subspan(static_cast<std::size_t>(row)*static_cast<std::size_t>(stride),static_cast<std::size_t>(lanes));}"
  1)
require_token_count("MQ.2 R1 row owner and consumers"
  mq2_r1_recorder "snapshotRow(" 3)
require_token_count("MQ.2 R1 CSV ordered row call"
  mq2_r1_recorder
  "snapshotRow(recorded,row,rows_,stride_,lanes_)" 1)
require_token_count("MQ.2 R1 Parquet ordered row call"
  mq2_r1_recorder
  "snapshotRow(snapshot(i),row,rows_,stride_,lanes_)" 1)
require_token_count("MQ.2 R1 CSV row result is consumed"
  mq2_r1_recorder
  "for(introw=0;row<rows_;++row)for(constreal_tvalue:snapshotRow(recorded,row,rows_,stride_,lanes_))output<<','<<value;"
  1)
require_token_count("MQ.2 R1 Parquet row result is consumed"
  mq2_r1_recorder
  "returnsnapshotRow(snapshot(i),row,rows_,stride_,lanes_)[static_cast<std::size_t>(lane)];"
  1)
require_token_count("MQ.2 R1 late row-stride cast removed"
  mq2_r1_recorder "static_cast<std::size_t>(row*stride_" 0)
require_token_count("MQ.2 R1 raw row-stride-lane arithmetic removed"
  mq2_r1_recorder "row*stride_+lane" 0)

# AsyncRecorder never consumed these checked-arithmetic names. The codec
# retains its live layout dependencies at their exact old-production counts.
require_token_count("MQ.2 R1 dead async checked include removed"
  mq2_r1_async_recorder "detail/CheckedArithmetic.hpp" 0)
require_token_count("MQ.2 R1 dead async checkedAdd removed"
  mq2_r1_async_recorder "checkedAdd" 0)
require_token_count("MQ.2 R1 dead async checkedMultiply removed"
  mq2_r1_async_recorder "checkedMultiply" 0)
require_token_count("MQ.2 R1 dead async compressionBound removed"
  mq2_r1_async_recorder "compressionBound" 0)
require_token_count("MQ.2 R1 live codec checkedAdd retained"
  mq2_r1_async_codec "checkedAdd" 10)
require_token_count("MQ.2 R1 live codec checkedMultiply retained"
  mq2_r1_async_codec "checkedMultiply" 14)
require_token_count("MQ.2 R1 live codec compressionBound retained"
  mq2_r1_async_codec "compressionBound" 2)
require_token_count("MQ.2 R1 live codec checked include retained"
  mq2_r1_async_codec "#include\"detail/CheckedArithmetic.hpp\"" 1)
require_token_count("MQ.2 R1 mapped move include retained"
  mq2_r1_recording_format "#include<utility>" 1)

require_token_count("MQ.2 R1 synchronous endian import"
  mq2_r1_recording_format "usingdetail::endian_marker;" 1)
require_token_count("MQ.2 R1 async writer endian import"
  mq2_r1_async_recorder "usingdetail::endian_marker;" 1)
require_token_count("MQ.2 R1 async reader endian import"
  mq2_r1_async_codec "usingdetail::endian_marker;" 1)
require_token_count("MQ.2 R1 synchronous writer endian consumer"
  mq2_r1_recording_format ".endian=endian_marker," 1)
require_token_count("MQ.2 R1 synchronous reader endian consumer"
  mq2_r1_recording_format "header.endian!=endian_marker" 1)
require_token_count("MQ.2 R1 async writer endian consumer"
  mq2_r1_async_recorder ".endian=endian_marker," 1)
require_token_count("MQ.2 R1 async reader endian consumer"
  mq2_r1_async_codec "header.endian!=endian_marker" 1)
require_token_count("MQ.2 R1 no synchronous private endian owner"
  mq2_r1_recording_format
  "constexprstd::uint32_tendian_marker=" 0)
require_token_count("MQ.2 R1 no async-writer private endian owner"
  mq2_r1_async_recorder
  "constexprstd::uint32_tendian_marker=" 0)
require_token_count("MQ.2 R1 no async-reader private endian owner"
  mq2_r1_async_codec
  "constexprstd::uint32_tendian_marker=" 0)
require_token_count("MQ.2 R1 exactly two major-version owners"
  mq2_r1_format_sources
  "constexprstd::uint16_tformat_major=" 2)
require_token_count("MQ.2 R1 exactly two minor-version owners"
  mq2_r1_format_sources
  "constexprstd::uint16_tformat_minor=" 2)

set(mq2_r1_allocation_catch_pair
  "catch(conststd::bad_alloc&){returnallocationFailureStatus();}catch(conststd::length_error&){returnallocationFailureStatus();}")
require_token_count("MQ.2 R1 Recorder allocation catch delegation"
  mq2_r1_recorder "${mq2_r1_allocation_catch_pair}" 3)
require_token_count("MQ.2 R1 synchronous allocation catch delegation"
  mq2_r1_recording_format "${mq2_r1_allocation_catch_pair}" 1)
require_token_count("MQ.2 R1 async writer allocation catch delegation"
  mq2_r1_async_recorder "${mq2_r1_allocation_catch_pair}" 1)
require_token_count("MQ.2 R1 async reader allocation catch delegation"
  mq2_r1_async_codec "${mq2_r1_allocation_catch_pair}" 1)

string(FIND "${mq2_r1_snapshot_with_comments}"
  "/***@fileSnapshotIndexing.hpp" mq2_r1_snapshot_leading_position)
string(FIND "${mq2_r1_snapshot_with_comments}"
  "#pragmaonce" mq2_r1_snapshot_pragma_position)
if(NOT mq2_r1_snapshot_leading_position EQUAL 0
   OR mq2_r1_snapshot_pragma_position EQUAL -1)
  message(FATAL_ERROR
    "MQ.2 R1 SnapshotIndexing contract is not the leading file contract")
endif()
string(SUBSTRING "${mq2_r1_snapshot_with_comments}"
  0 ${mq2_r1_snapshot_pragma_position} mq2_r1_snapshot_leading_contract)
require_ordered_tokens("MQ.2 R1 snapshot leading contract"
  mq2_r1_snapshot_leading_contract
  "/***@fileSnapshotIndexing.hpp"
  "@briefSharedSoA-to-SnapshotViewindexingforeagerrecordingreaders."
  "Owns:`snapshotView`.ImplementsPLAN.mdsection3.7."
  "Cold:allocation-freeeager-readerindexingonly."
  "@surfaceinternal*/")

# MQ.2 R2: the test-side oracle is itself pinned to three complete rows, an
# independent locale-free parser, the precomputed legal-zero vector, and a
# CRC-field-only corruption. Exact production roles remain in the R1 block.
load_compact("tests/unit/core_Recorder_test.cpp" mq2_r2_recorder_test)
require_token_count("MQ.2 R2 exact real-field count"
  mq2_r2_recorder_test "constexprstd::size_tcsv_real_fields=63;" 1)
require_token_count("MQ.2 R2 parser owner"
  mq2_r2_recorder_test
  "std::optional<std::array<CsvRow,3>>parseRecorderCsv(std::string_viewtext)"
  1)
require_token_count("MQ.2 R2 locale-free numeric parser calls"
  mq2_r2_recorder_test "std::from_chars(" 2)
require_token_count("MQ.2 R2 scientific real parser"
  mq2_r2_recorder_test "std::chars_format::scientific" 1)
forbid_tokens("MQ.2 R2 parser has no locale or throwing numeric fallback"
  mq2_r2_recorder_test "std::stod(" "std::stringstream")
require_token_count("MQ.2 R2 nondegenerate accepted steps"
  mq2_r2_recorder_test
  "constexprstd::array<std::uint64_t,3>accepted_steps{7,11,19};"
  1)
require_token_count("MQ.2 R2 independent live-row expectation"
  mq2_r2_recorder_test
  "batch.state().row(row)[static_cast<std::size_t>(lane)]" 1)
require_token_count("MQ.2 R2 independent voltage destination"
  mq2_r2_recorder_test
  "std::span{expected.values}.subspan(61,2)" 1)
require_token_count("MQ.2 R2 all live rows captured"
  mq2_r2_recorder_test
  "for(std::size_tindex=0;index<accepted_steps.size();++index)" 1)
require_token_count("MQ.2 R2 all expected rows compared"
  mq2_r2_recorder_test
  "for(std::size_tindex=0;index<expected_csv.size();++index)" 1)
require_token_count("MQ.2 R2 all real fields compared"
  mq2_r2_recorder_test
  "value<expected_csv[index].values.size()" 1)
require_token_count("MQ.2 R2 exact zero-CRC header owner"
  mq2_r2_recorder_test "zero_crc_header_hex=" 1)
require_token_count("MQ.2 R2 frozen zero-CRC header prefix"
  mq2_r2_recorder_test
  "534c4944455245430100000004030201400000000000000070000000c7000000"
  1)
require_token_count("MQ.2 R2 frozen zero-CRC header suffix"
  mq2_r2_recorder_test
  "87030000000000004000000000000000000e000000000000000e000000000000"
  1)
require_token_count("MQ.2 R2 CRC-field-only corruption"
  mq2_r2_recorder_test "corrupt[20]=std::byte{1};" 1)
require_token_count("MQ.2 R2 independent header CRC owner and witness"
  mq2_r2_recorder_test "testHeaderCrc(" 2)
require_token_count("MQ.2 R2 failed-open atomicity witness"
  mq2_r2_recorder_test
  "CHECK(recording.open(corrupt_path)==Status::Invalid_parameters);CHECK(recording.valid());"
  1)

message(STATUS "9C architecture aggregate structural gate passed")
