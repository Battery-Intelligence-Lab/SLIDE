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

message(STATUS "9C architecture aggregate structural gate passed")
