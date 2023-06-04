cmake_minimum_required(VERSION 3.17)
add_library (eigen INTERFACE)
target_include_directories(eigen SYSTEM INTERFACE .)