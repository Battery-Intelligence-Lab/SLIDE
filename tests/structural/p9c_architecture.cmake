# M0.6--M0.7 architecture aggregate. Keep one CTest slot while retaining the
# independently runnable 9C-2 and 9C-3 structural proofs.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/p9c2_ageing_kernel.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/p9c3_cold_file_split.cmake")

message(STATUS "9C architecture aggregate structural gate passed")
