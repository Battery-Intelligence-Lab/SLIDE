# M0.6 / 9C-2: four named ageing equation bodies, one production scaffold.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

function(load_compact relative_path output)
  set(path "${SLIDE_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "9C-2 source is missing: ${relative_path}")
  endif()
  file(READ "${path}" content)
  string(REGEX REPLACE "//[^\r\n]*" "" content "${content}")
  string(REGEX REPLACE "/\\*([^*]|\\*+[^*/])*\\*+/" "" content "${content}")
  string(REGEX REPLACE "[ \t\r\n]" "" content "${content}")
  set(${output} "${content}" PARENT_SCOPE)
endfunction()

function(require_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR "9C-2 ${label}: required token is absent: ${token}")
    endif()
  endforeach()
endfunction()

function(forbid_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR "9C-2 ${label}: private scaffold remains: ${token}")
    endif()
  endforeach()
endfunction()

function(require_token_count label content_variable token expected)
  string(LENGTH "${${content_variable}}" original_length)
  string(REPLACE "${token}" "" without_token "${${content_variable}}")
  string(LENGTH "${without_token}" stripped_length)
  string(LENGTH "${token}" token_length)
  math(EXPR count "(${original_length} - ${stripped_length}) / ${token_length}")
  if(NOT count EQUAL expected)
    message(
      FATAL_ERROR
      "9C-2 ${label}: expected ${expected} occurrences of ${token}, found ${count}")
  endif()
endfunction()

function(require_ordered_tokens label content_variable)
  set(search_offset 0)
  foreach(token IN LISTS ARGN)
    string(SUBSTRING "${${content_variable}}" ${search_offset} -1 remaining)
    string(FIND "${remaining}" "${token}" relative_position)
    if(relative_position EQUAL -1)
      message(
        FATAL_ERROR
        "9C-2 ${label}: required ordered token is absent after offset ${search_offset}: ${token}")
    endif()
    string(LENGTH "${token}" token_length)
    math(EXPR search_offset
      "${search_offset} + ${relative_position} + ${token_length}")
  endforeach()
endfunction()

load_compact("src/core/Sei.hpp" sei)
load_compact("src/core/AgeingKernel.hpp" scaffold)
load_compact("src/core/SurfaceCrack.hpp" crack)
load_compact("src/core/Lam.hpp" lam)
load_compact("src/core/LithiumPlating.hpp" plating)
# M0.9 / 9C-5 moved the cold parameter blocks and their shared model-mask vocabulary out of
# the mechanism headers so the public factory does not drag a kernel in. The 9C-2 invariant
# is unchanged -- one mask idiom, not four hand-rolled ones -- so it is now checked where the
# masks live, and each mechanism must consume its own parameter header.
load_compact("src/core/SeiParams.hpp" sei_params)
load_compact("src/core/SurfaceCrackParams.hpp" crack_params)
load_compact("src/core/LamParams.hpp" lam_params)
load_compact("src/core/LithiumPlatingParams.hpp" plating_params)
load_compact("src/core/SpmStress.hpp" stress)
load_compact("src/core/SpmPipeline.hpp" pipeline)
load_compact("src/core/SpmFactory.cpp" factory)
load_compact("src/cells/Cell_SPM/Cell_SPM_degradation.cpp" legacy)

foreach(consumer IN ITEMS sei crack lam plating stress pipeline)
  require_tokens("${consumer}" ${consumer} "#include\"AgeingKernel.hpp\"")
endforeach()

require_tokens("shared scaffold" scaffold
  "for(unsignedmodel=1;model<=ModelCount;++model)"
  "SLIDE_AGEING_FORCE_INLINEslide::Statusfor_each_enabled_ageing_model_lane("
  "SLIDE_AGEING_FORCE_INLINEslide::Statusevaluate_ageing_stage(")

require_tokens("SEI masks" sei_params
  "ageing_model_bit<4>(model)"
  "valid_ageing_model_mask<4>(p.model_mask)")
require_tokens("surface-crack masks" crack_params
  "ageing_model_bit<5>(model)"
  "valid_ageing_model_mask<5>(p.model_mask)")
require_tokens("LAM masks" lam_params
  "ageing_model_bit<4>(model)"
  "valid_ageing_model_mask<4>(p.model_mask)")
require_tokens("plating scales" plating_params
  "tryLithiumPlatingScales(")

# Include-presence alone would be weaker than what 9C-2 defended before M0.9: a mechanism could
# include its params header and still stop consuming the shared mask. So each mechanism must
# still hand its own `p.model_mask` into the shared traversal, which is the actual invariant.
require_tokens("SEI mask use" sei
  "detail::for_each_enabled_ageing_model_lane<4>(p.model_mask,")
require_tokens("surface-crack mask use" crack
  "detail::for_each_enabled_ageing_model_lane<5>(p.model_mask,")
require_tokens("LAM mask use" lam
  "detail::for_each_enabled_ageing_model_lane<4>(p.model_mask,")

require_tokens("SEI" sei
  "#include\"SeiParams.hpp\""
  "detail::AgeingScratchStorage<Real,2>"
  "detail::clear_ageing_fields<Real,2>("
  "detail::for_each_enabled_ageing_model_lane<4>("
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane("
  "return{storage_.field(0),storage_.field(1)}")

require_tokens("surface crack" crack
  "#include\"SurfaceCrackParams.hpp\""
  "detail::AgeingScratchStorage<Real,3>"
  "detail::clear_ageing_fields<Real,3>("
  "detail::for_each_enabled_ageing_model_lane<5>("
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane("
  "return{storage_.field(0),storage_.field(1),storage_.field(2)}")

require_tokens("LAM" lam
  "#include\"LamParams.hpp\""
  "detail::AgeingScratchStorage<Real,6>"
  "detail::clear_ageing_fields<Real,6>("
  "detail::for_each_enabled_ageing_model_lane<4>("
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane("
  "(*field)[domain_index(domain)]=storage_.field(cursor)"
  "++cursor")

require_tokens("lithium plating" plating
  "#include\"LithiumPlatingParams.hpp\""
  "BasicLithiumPlatingOutput"
  "detail::AgeingScratchStorage<Real,1>"
  "return{storage_.field(0)}"
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane(")

require_tokens("stress scratch" stress "detail::AgeingScratchStorage<Real,3>")
require_token_count("pipeline stages" pipeline "detail::evaluate_ageing_stage(" 4)
require_ordered_tokens("pipeline stage order" pipeline
  "detail::evaluate_ageing_stage(params_.enable_sei,"
  "detail::evaluate_ageing_stage(params_.enable_surface_crack,"
  "detail::evaluate_ageing_stage(params_.enable_lam,"
  "detail::evaluate_ageing_stage(params_.enable_lithium_plating,")
require_tokens("pipeline disabled-SEI dependency" pipeline
  "detail::clear_ageing_fields<real_t,2>(")
require_tokens("factory masks" factory
  "valid_optional_ageing_model_mask<4>(options.sei_model_mask)"
  "valid_optional_ageing_model_mask<5>(options.surface_crack_model_mask)"
  "valid_optional_ageing_model_mask<4>(options.lam_model_mask)")

foreach(mechanism IN ITEMS sei crack lam plating)
  forbid_tokens("${mechanism}" ${mechanism}
    "std::vector<Real>storage_"
    "std::fill("
    "for(unsignedmodel=1;model<="
    "for(intlane=0;lane<lanes;++lane)")
endforeach()

# The legacy Cell_SPM equations are the independent parity oracle. Sharing the
# production helper here would turn mechanism parity into a common-mode test.
forbid_tokens("legacy oracle" legacy
  "AgeingKernel.hpp"
  "AgeingScratchStorage"
  "for_each_enabled_ageing_model_lane"
  "evaluate_ageing_stage")

message(STATUS "9C-2 shared ageing-kernel structural gate passed")
