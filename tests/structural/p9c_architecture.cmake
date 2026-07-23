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

message(STATUS "9C architecture aggregate structural gate passed")
