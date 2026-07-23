# M0.9 / 9C-5: the public surface is minimal, classified, and intentional (MC-5).
#
# Every src/core header carries exactly one `@surface api|support|internal` tag.
# `api` names what a binding may include; `support` is the vocabulary the API's signatures
# are written in; `internal` is implementation detail. The load-bearing rule is R3: no
# api/support header may include an internal one, so no user translation unit compiles a
# kernel. The api set and each api header's declaration count are pinned, so a new public
# entity cannot appear without a deliberate edit to this file.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

set(P9C5_API
  AsyncRecorder CellDesign CudaSpmBatch EulerLegacy Experiment ExponentialModal
  ForwardSensitivity NetlistCsv PackSolver PackStepper PackTopology ParameterSet
  Recorder Simulation SpmFactory ThreadPool)

set(P9C5_SUPPORT
  AgeingModelMask BatchBuilder BatchView CompiledCurve LamParams LithiumPlatingParams
  Numeric SeiParams SpmBatchLayout SpmScalarKernels StateArena SurfaceCrackParams)

# Column-0 declaration anchors, api AND support: newline followed by an identifier start or
# `[[nodiscard]]`. Adding or removing any namespace-scope entity changes the count.
set(P9C5_ANCHORS
  "AsyncRecorder=14" "CellDesign=24" "CudaSpmBatch=8" "EulerLegacy=4" "Experiment=16"
  "ExponentialModal=4" "ForwardSensitivity=8" "NetlistCsv=5" "PackSolver=16"
  "PackStepper=4" "PackTopology=23" "ParameterSet=6" "Recorder=11" "Simulation=6"
  "SpmFactory=9" "ThreadPool=12"
  "AgeingModelMask=9" "BatchBuilder=4" "BatchView=13" "CompiledCurve=8" "LamParams=4"
  "LithiumPlatingParams=5" "Numeric=9" "SeiParams=4" "SpmBatchLayout=5"
  "SpmScalarKernels=36" "StateArena=9" "SurfaceCrackParams=4")

# Member-level anchors (exactly two spaces of indent): class methods, enumerators, and struct
# fields. Without these, the entire member surface of every public class was unpinned.
set(P9C5_MEMBERS
  "AsyncRecorder=93" "CellDesign=65" "CudaSpmBatch=38" "EulerLegacy=14" "Experiment=60"
  "ExponentialModal=14" "ForwardSensitivity=36" "NetlistCsv=13" "PackSolver=85"
  "PackStepper=33" "PackTopology=55" "ParameterSet=18" "Recorder=69" "Simulation=20"
  "SpmFactory=90" "ThreadPool=54"
  "AgeingModelMask=6" "BatchBuilder=14" "BatchView=28" "CompiledCurve=27" "LamParams=21"
  "LithiumPlatingParams=27" "Numeric=18" "SeiParams=37" "SpmBatchLayout=25"
  "SpmScalarKernels=122" "StateArena=56" "SurfaceCrackParams=25")

# Public types each api header must still declare (removal is as much a surface change as
# addition, and an aggregate count alone cannot see a swap).
set(P9C5_TYPES_AsyncRecorder AsyncBackpressurePolicy CompressionCodec AsyncRecorderConfig AsyncRecorder CompressedRecording)
set(P9C5_TYPES_CellDesign Domain PerDomain OCVCurve Arrhenius ActiveMaterial StressParams AgingMechanismKind AgingMechanismSpec ElectrodeDesign SeparatorDesign ElectrolyteDesign ThermalDesign CellDesign ElectrodeParams)
set(P9C5_TYPES_CudaSpmBatch CudaAsyncRecorder CudaSpmBatch)
set(P9C5_TYPES_EulerLegacy EulerLegacy)
set(P9C5_TYPES_Experiment ControlMode Direction ExperimentVariables ExperimentFunction CustomTermination ExperimentSegment ParseDiagnostic Experiment TerminationReason ExperimentSolution DriveCycle CyclerIntegrator CyclerV2)
set(P9C5_TYPES_ExponentialModal ExponentialModal)
set(P9C5_TYPES_ForwardSensitivity SensitivityParameter ForwardSensitivitySolution)
set(P9C5_TYPES_NetlistCsv NetlistCsvOptions NetlistCsvDiagnostic)
set(P9C5_TYPES_PackSolver TheveninBatchView PackTheveninSystem PackSolveMode PackSolveDiagnostics PackSolution SolverWorkspace PackSolver)
set(P9C5_TYPES_PackStepper PackStepper)
set(P9C5_TYPES_PackTopology PackCellSpec PackLink PackNodeKind PackNode ThermalBoundarySpec ThermalLinkSpec PackDescription ElectricalBranchKind CompiledElectricalBranch BatchLaneLocation CompiledCell CompiledElectricalNetlist ThermalEdge ThermalIncident CompiledThermalGraph CompiledPackTopology)
set(P9C5_TYPES_ParameterSet ParameterValue ParameterDescription ParameterSet)
set(P9C5_TYPES_Recorder BackpressurePolicy RecorderConfig SnapshotView Recorder BinaryRecording)
set(P9C5_TYPES_Simulation ConstantCurrentExperiment SimulationSolution Simulation)
set(P9C5_TYPES_SpmFactory SpmComposition SpmModelOptions SpmFactoryInput SpmBatch)
set(P9C5_TYPES_ThreadPool ThreadPool BatchExecutor ParallelisationDiagnostic)

function(p9c5_strip_comments raw output)
  set(content "${raw}")
  string(REGEX REPLACE "/\\*([^*]|\\*+[^*/])*\\*+/" "" content "${content}")
  string(REGEX REPLACE "//[^\r\n]*" "" content "${content}")
  set(${output} "${content}" PARENT_SCOPE)
endfunction()

function(p9c5_compact raw output)
  p9c5_strip_comments("${raw}" content)
  string(REGEX REPLACE "[ \t\r\n]" "" content "${content}")
  set(${output} "${content}" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------------------
# R1 -- every core header is classified, exactly once, with a valid tier.
# ---------------------------------------------------------------------------------------
file(GLOB P9C5_HEADERS
  "${SLIDE_SOURCE_DIR}/src/core/*.hpp"
  "${SLIDE_SOURCE_DIR}/src/core/detail/*.hpp")
if(P9C5_HEADERS STREQUAL "")
  message(FATAL_ERROR "9C-5: no src/core headers found")
endif()

set(P9C5_TIER_MAP "")
foreach(header IN LISTS P9C5_HEADERS)
  get_filename_component(stem "${header}" NAME_WE)
  file(READ "${header}" raw)
  string(REGEX MATCHALL "@surface[ \t]+([A-Za-z]+)" tags "${raw}")
  list(LENGTH tags tag_count)
  if(NOT tag_count EQUAL 1)
    message(FATAL_ERROR
      "9C-5 R1: ${stem}.hpp must carry exactly one @surface tag, found ${tag_count}")
  endif()
  string(REGEX REPLACE "@surface[ \t]+" "" tier "${tags}")
  if(NOT tier MATCHES "^(api|support|internal)$")
    message(FATAL_ERROR
      "9C-5 R1: ${stem}.hpp has an unknown @surface tier: '${tier}'")
  endif()

  # R2 -- the tier a header claims must match the pinned classification.
  if(stem IN_LIST P9C5_API)
    set(expected api)
  elseif(stem IN_LIST P9C5_SUPPORT)
    set(expected support)
  else()
    set(expected internal)
  endif()
  if(NOT tier STREQUAL expected)
    message(FATAL_ERROR
      "9C-5 R2: ${stem}.hpp is tagged '${tier}' but this gate pins it as '${expected}'. "
      "Changing the public surface is a deliberate act: update tests/structural/p9c5_public_surface.cmake.")
  endif()
  list(APPEND P9C5_TIER_MAP "${stem}=${tier}")
endforeach()

# Every pinned api/support header must exist (a rename must not silently drop a rule).
foreach(stem IN LISTS P9C5_API P9C5_SUPPORT)
  if(NOT EXISTS "${SLIDE_SOURCE_DIR}/src/core/${stem}.hpp")
    message(FATAL_ERROR "9C-5 R2: pinned surface header is missing: ${stem}.hpp")
  endif()
endforeach()

function(p9c5_tier stem output)
  foreach(entry IN LISTS P9C5_TIER_MAP)
    if(entry MATCHES "^${stem}=(.*)$")
      set(${output} "${CMAKE_MATCH_1}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${output} "unknown" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------------------
# R3 -- no api/support header includes an internal header.
#
# This enforces exactly one thing: no header classified `internal` is reachable from the api
# set. That kills the ageing/pipeline kernel stack (SpmPipeline, SpmObservables, Sei, Lam,
# SurfaceCrack, LithiumPlating, ThermalLumped, SpmStress). It does NOT claim a user TU
# compiles no physics at all: `SpmScalarKernels.hpp` is `support` and IS reachable, because
# `CompiledCurve`'s inline evaluator calls `linearInterpolate` and duplicating it would break
# PC-10 (one physics source), which outranks MC-5. That residual is deliberate and recorded in
# `.claude/designs/m0-9-public-surface.md`; do not read this rule as more than it is.
#
# Both include spellings are checked, and a `../` prefix does not exempt a path that climbs
# back into src/core -- an adversarial review walked through both of those holes.
# ---------------------------------------------------------------------------------------
foreach(header IN LISTS P9C5_HEADERS)
  get_filename_component(stem "${header}" NAME_WE)
  p9c5_tier("${stem}" tier)
  if(tier STREQUAL "internal")
    continue()
  endif()
  file(READ "${header}" raw)
  p9c5_strip_comments("${raw}" content)
  string(REGEX MATCHALL "#include[ \t]+[\"<][A-Za-z_/.]+[\">]" includes "${content}")
  set(seen_include_targets)
  foreach(include IN LISTS includes)
    string(REGEX REPLACE "#include[ \t]+[\"<](.*)[\">]" "\\1" target "${include}")
    if(target IN_LIST seen_include_targets)
      message(FATAL_ERROR
        "9C-5 R3: ${tier} header ${stem}.hpp includes ${target} more than once")
    endif()
    list(APPEND seen_include_targets "${target}")
    get_filename_component(target_stem "${target}" NAME_WE)
    if(target MATCHES "^\\.\\./" AND NOT target MATCHES "core/")
      continue() # ../types/Status.hpp and friends genuinely live outside src/core
    endif()
    p9c5_tier("${target_stem}" target_tier)
    if(target_tier STREQUAL "internal")
      message(FATAL_ERROR
        "9C-5 R3: ${tier} header ${stem}.hpp includes internal header ${target}. "
        "Implementation detail must not reach a user translation unit (MC-5).")
    endif()
    if(target_tier STREQUAL "unknown")
      # Standard-library and external headers carry no tier and need none. A path that names
      # src/core, however, must be classified -- otherwise a new unclassified core header would
      # be a free pass through this rule.
      if(EXISTS "${SLIDE_SOURCE_DIR}/src/core/${target_stem}.hpp"
         OR EXISTS "${SLIDE_SOURCE_DIR}/src/core/detail/${target_stem}.hpp")
        message(FATAL_ERROR
          "9C-5 R3: ${stem}.hpp includes ${target}, which carries no surface classification")
      endif()
    endif()
  endforeach()
endforeach()

# ---------------------------------------------------------------------------------------
# R4 -- the public surface is pinned, at BOTH levels, across api AND support headers.
#
# The first version of this rule counted column-0 declarations only, so every class member,
# every enumerator, and every parameter field was unpinned: an adversarial review added a
# public method to `SpmBatch`, an enumerator to `SpmComposition`, and a field to `SeiParams`,
# and the gate stayed green all three times. Support headers were not pinned at all, even
# though the `*Params` structs ARE the public parameter surface -- a new field in `SeiParams`
# changes the layout of `SpmFactoryInput`, an api type.
#
# So both anchors are pinned now: namespace-scope declarations (column 0) and member-level
# declarations (exactly two spaces of indent, which is where this codebase puts class members
# and enumerators). This is a tripwire on the shape of these headers, not a parser: it fires
# on any added or removed declaration, and two offsetting edits could still net out.
# ---------------------------------------------------------------------------------------
function(p9c5_count_anchors content pattern output)
  # Count by length growth, not list(LENGTH): a match containing '[' does not survive CMake's
  # list re-parsing, and would silently undercount `[[nodiscard]]` declarations. Every match is
  # two or three characters, so measure the growth per match rather than assuming a width.
  string(REGEX REPLACE "${pattern}" "@@@@@@@@" marked "${content}")
  string(LENGTH "${content}" plain_length)
  string(LENGTH "${marked}" marked_length)
  math(EXPR delta "${marked_length} - ${plain_length}")
  set(${output} "${delta}" PARENT_SCOPE)
endfunction()

foreach(entry IN LISTS P9C5_ANCHORS)
  string(REGEX REPLACE "=.*$" "" stem "${entry}")
  string(REGEX REPLACE "^.*=" "" expected "${entry}")
  file(READ "${SLIDE_SOURCE_DIR}/src/core/${stem}.hpp" raw)
  p9c5_strip_comments("${raw}" content)
  # Each match is 2 chars ("\n" + 1) replaced by 8, so growth is 6 per match.
  p9c5_count_anchors("${content}" "\n([A-Za-z_]|\\[)" growth)
  math(EXPR count "${growth} / 6")
  if(NOT count EQUAL expected)
    message(FATAL_ERROR
      "9C-5 R4: ${stem}.hpp declares ${count} namespace-scope entities, the pinned surface "
      "has ${expected}. A public entity was added or removed -- update this gate on purpose.")
  endif()

  # Member level: "\n" + two spaces + 1 char = 4 chars replaced by 8, so growth is 4 per match.
  foreach(member_entry IN LISTS P9C5_MEMBERS)
    if(member_entry MATCHES "^${stem}=(.*)$")
      set(expected_members "${CMAKE_MATCH_1}")
      p9c5_count_anchors("${content}" "\n  ([A-Za-z_~]|\\[)" member_growth)
      math(EXPR member_count "${member_growth} / 4")
      if(NOT member_count EQUAL expected_members)
        message(FATAL_ERROR
          "9C-5 R4: ${stem}.hpp declares ${member_count} member-level entities, the pinned "
          "surface has ${expected_members}. A public method, enumerator, or field was added or "
          "removed -- update this gate on purpose.")
      endif()
    endif()
  endforeach()

  p9c5_compact("${raw}" compact)
  foreach(type IN LISTS P9C5_TYPES_${stem})
    # Match the declaration, not the substring: a plain FIND for `NetlistCsvDiagnostic`
    # also succeeds against `NetlistCsvDiagnosticRenamed`, so a rename of a public type
    # would slip through with the anchor count unchanged.
    if(NOT compact MATCHES "(class|struct|enumclass|using)${type}[{:=;<]")
      message(FATAL_ERROR
        "9C-5 R4: ${stem}.hpp no longer declares the pinned public type ${type}")
    endif()
  endforeach()
endforeach()

# ---------------------------------------------------------------------------------------
# R5 -- language bindings and docs see only the api surface.
# ---------------------------------------------------------------------------------------
set(P9C5_CONSUMERS
  "python/bindings.cpp"
  "matlab/slide_mex.cpp"
  "docs/v4/quickstart-cpp.md")
foreach(consumer IN LISTS P9C5_CONSUMERS)
  set(path "${SLIDE_SOURCE_DIR}/${consumer}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "9C-5 R5: registered public consumer is missing: ${consumer}")
  endif()
  file(READ "${path}" raw)
  # `core/[A-Za-z_]+\.hpp` could not match `core/detail/StrictJson.hpp`, so a binding could
  # include a detail header and pass. Match nested paths too.
  string(REGEX MATCHALL "core/[A-Za-z_/]+\\.hpp" includes "${raw}")
  foreach(include IN LISTS includes)
    get_filename_component(target_stem "${include}" NAME_WE)
    p9c5_tier("${target_stem}" target_tier)
    if(NOT target_tier STREQUAL "api")
      message(FATAL_ERROR
        "9C-5 R5: ${consumer} includes ${include}, which is '${target_tier}'. "
        "Bindings and docs may include api headers only.")
    endif()
  endforeach()
endforeach()

# ---------------------------------------------------------------------------------------
# R6 -- one name per concept (the M0.9 lexicon).
#   * lane/row counts are n_lanes()/n_rows(), never nLanes()/nRows()
#   * the device counters are ...Count / deviceArenaBytes
#   * the private lane check is checked_lane_count / checked_shape
#   * validate... returns slide::Status; valid.../is... returns bool
#   * no get-prefixed accessors
# ---------------------------------------------------------------------------------------
# Forbid the dead spelling itself, not one syntactic shape of it: `intnLanes(` missed
# `auto nLanes()`, `a->nLanes()`, and `std::size_t nLanes()`. No legitimate use of these
# spellings exists anywhere in core, so the bare token is the right prohibition.
set(P9C5_FORBIDDEN_TOKENS
  "nLanes("
  "nRows("
  "checkedLaneCount"
  "checkedShape"
  "deviceAllocations("
  "deviceWideSynchronizations("
  "deviceBytes("
  "SpmPipelineLayout")
set(P9C5_FORBIDDEN_REGEX
  "boolvalidate[A-Z]"      # validate... must return slide::Status
  "Statusvalid[A-Z]")      # valid... must be a bool predicate

# Checked against comment-stripped text that still has its whitespace: in the compacted form
# a leading `get` is glued to the return type (`inlineintgetWorkerCount`) and no word
# boundary survives to anchor against.
set(P9C5_FORBIDDEN_SPACED_REGEX
  "[ \t(*&,:~]get[A-Z]"    # no get-prefixed accessors in core
  "(^|[^_A-Za-z0-9])lanes[ \t\r\n]*\\(") # lane-count accessors are n_lanes()

file(GLOB P9C5_SOURCES
  "${SLIDE_SOURCE_DIR}/src/core/*.hpp"
  "${SLIDE_SOURCE_DIR}/src/core/*.cpp"
  "${SLIDE_SOURCE_DIR}/src/core/*.cu"
  "${SLIDE_SOURCE_DIR}/src/core/detail/*.hpp")
foreach(source IN LISTS P9C5_SOURCES)
  get_filename_component(name "${source}" NAME)
  file(READ "${source}" raw)
  p9c5_compact("${raw}" compact)
  foreach(token IN LISTS P9C5_FORBIDDEN_TOKENS)
    string(FIND "${compact}" "${token}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR
        "9C-5 R6: ${name} uses '${token}', which is not the canonical name for that concept")
    endif()
  endforeach()
  foreach(pattern IN LISTS P9C5_FORBIDDEN_REGEX)
    if(compact MATCHES "${pattern}")
      message(FATAL_ERROR
        "9C-5 R6: ${name} matches the forbidden naming pattern '${pattern}' (${CMAKE_MATCH_0})")
    endif()
  endforeach()
  p9c5_strip_comments("${raw}" spaced)
  foreach(pattern IN LISTS P9C5_FORBIDDEN_SPACED_REGEX)
    if(spaced MATCHES "${pattern}")
      message(FATAL_ERROR
        "9C-5 R6: ${name} matches the forbidden naming pattern '${pattern}' (${CMAKE_MATCH_0})")
    endif()
  endforeach()
endforeach()

# The canonical spellings must actually be there -- a gate of prohibitions alone would pass
# on an empty file.
set(P9C5_REQUIRED
  "src/core/StateArena.hpp=intn_lanes()"
  "src/core/StateArena.hpp=intn_rows()"
  "src/core/BatchView.hpp=intn_lanes()"
  "src/core/BatchView.hpp=intn_rows()"
  "src/core/Recorder.hpp=intn_lanes()"
  "src/core/Recorder.hpp=intn_rows()"
  "src/core/AsyncRecorder.hpp=intn_lanes()"
  "src/core/AsyncRecorder.hpp=intn_rows()"
  "src/core/CudaSpmBatch.hpp=intn_lanes()"
  "src/core/SpmFactory.hpp=intn_lanes()"
  "src/core/CudaSpmBatch.hpp=deviceAllocationCount()"
  "src/core/CudaSpmData.hpp=deviceAllocationCount("
  "src/core/SeiParams.hpp=StatusvalidateSeiParams("
  "src/core/LamParams.hpp=StatusvalidateLamParams("
  "src/core/SurfaceCrackParams.hpp=StatusvalidateSurfaceCrackParams("
  "src/core/LithiumPlatingParams.hpp=StatusvalidateLithiumPlatingParams(")
foreach(entry IN LISTS P9C5_REQUIRED)
  string(REGEX REPLACE "=.*$" "" relative "${entry}")
  string(REGEX REPLACE "^[^=]*=" "" token "${entry}")
  file(READ "${SLIDE_SOURCE_DIR}/${relative}" raw)
  p9c5_compact("${raw}" compact)
  string(FIND "${compact}" "${token}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR
      "9C-5 R6: ${relative} no longer spells the canonical name '${token}'")
  endif()
endforeach()

message(STATUS "9C-5 public-surface gate passed")
