cmake_minimum_required(VERSION 3.21)

add_library(utility INTERFACE)

target_sources(utility INTERFACE
    free_functions.hpp
)

target_include_directories(utility INTERFACE .)
target_include_directories(utility INTERFACE ./io)
target_include_directories(utility INTERFACE ${fmt_SOURCE_DIR}/include)

target_link_libraries(utility INTERFACE fmt)
