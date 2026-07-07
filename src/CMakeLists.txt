add_library(src STATIC)


target_include_directories(src PUBLIC .)

add_subdirectory(core) # v4 slide::core — NOT linked into legacy `src` (strangler, PLAN.md §5)
add_subdirectory(types)
add_subdirectory(cells)
add_subdirectory(modules)
add_subdirectory(cooling)
add_subdirectory(power_conversion)
add_subdirectory(system)
add_subdirectory(procedures)
add_subdirectory(factories)
add_subdirectory(settings)
add_subdirectory(utility)

if (UNIX)
  if (NOT APPLE)
    set(PTHREADS "-pthread")
  endif()
endif()

target_link_libraries(src
  PUBLIC
  cells
  modules
  cooling
  power_conversion
  system
  procedures
  factories
  settings
  types
  utility
  ${PTHREADS} # -> for WSL
  Boost::boost # -> for odeint
)

add_executable(papercode paperCode.cpp)
target_link_libraries(papercode src)