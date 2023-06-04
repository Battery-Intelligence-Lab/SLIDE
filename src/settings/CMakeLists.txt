add_library(settings INTERFACE)
target_include_directories(settings INTERFACE .)  # #TODO try to make this PUBLIC so other libraries use headers without relative path.

if(MSVC)
add_compile_definitions(DETECTED_COMPILER_MSVC)
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
add_compile_definitions(DETECTED_COMPILER_CLANG)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
add_compile_definitions(DETECTED_COMPILER_GCC)
else()
message(AUTHOR_WARNING "No compiler warnings set for '${CMAKE_CXX_COMPILER_ID}' compiler.")
endif()
