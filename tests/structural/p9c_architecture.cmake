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

message(STATUS "9C architecture aggregate structural gate passed")
