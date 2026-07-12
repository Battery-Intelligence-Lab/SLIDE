include_guard(GLOBAL)

option(ENABLE_SANITIZER_ADDRESS "Enable AddressSanitizer" OFF)
option(ENABLE_SANITIZER_LEAK "Enable LeakSanitizer" OFF)
option(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR "Enable UndefinedBehaviorSanitizer" OFF)
option(ENABLE_SANITIZER_THREAD "Enable ThreadSanitizer" OFF)
option(ENABLE_SANITIZER_MEMORY "Enable MemorySanitizer" OFF)

function(enable_sanitizers)
  set(sanitizers)
  if(ENABLE_SANITIZER_ADDRESS)
    list(APPEND sanitizers address)
  endif()
  if(ENABLE_SANITIZER_LEAK)
    list(APPEND sanitizers leak)
  endif()
  if(ENABLE_SANITIZER_UNDEFINED_BEHAVIOR)
    list(APPEND sanitizers undefined)
  endif()
  if(ENABLE_SANITIZER_THREAD)
    list(APPEND sanitizers thread)
  endif()
  if(ENABLE_SANITIZER_MEMORY)
    list(APPEND sanitizers memory)
  endif()

  if(NOT sanitizers)
    return()
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
     AND NOT CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    message(FATAL_ERROR
      "Sanitizers require GCC or Clang; found ${CMAKE_CXX_COMPILER_ID}")
  endif()
  if(ENABLE_IPO)
    message(FATAL_ERROR
      "Sanitizer builds require -DENABLE_IPO=OFF so instrumentation is inspectable")
  endif()
  if(ENABLE_SANITIZER_THREAD AND WIN32)
    message(FATAL_ERROR "ThreadSanitizer is unsupported on Windows; use Linux or WSL")
  endif()
  if(ENABLE_SANITIZER_THREAD AND (ENABLE_SANITIZER_ADDRESS
                                  OR ENABLE_SANITIZER_LEAK
                                  OR ENABLE_SANITIZER_UNDEFINED_BEHAVIOR
                                  OR ENABLE_SANITIZER_MEMORY))
    message(FATAL_ERROR "ThreadSanitizer must run in its own build tree")
  endif()
  if(ENABLE_SANITIZER_MEMORY AND (ENABLE_SANITIZER_ADDRESS
                                  OR ENABLE_SANITIZER_LEAK
                                  OR ENABLE_SANITIZER_THREAD
                                  OR ENABLE_SANITIZER_UNDEFINED_BEHAVIOR))
    message(FATAL_ERROR "MemorySanitizer must run in its own build tree")
  endif()
  if(ENABLE_SANITIZER_MEMORY AND NOT CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    message(FATAL_ERROR "MemorySanitizer requires Clang")
  endif()

  list(JOIN sanitizers "," sanitizer_list)
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=${sanitizer_list}>"
    "$<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer>"
    "$<$<COMPILE_LANGUAGE:CXX>:-fno-sanitize-recover=all>")
  add_link_options("$<$<LINK_LANGUAGE:CXX>:-fsanitize=${sanitizer_list}>")
  message(STATUS "SLIDE project sanitizer instrumentation: ${sanitizer_list}")
endfunction()
