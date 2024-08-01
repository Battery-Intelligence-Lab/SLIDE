cmake_minimum_required(VERSION 3.21)

add_executable(parallel_model
  parallel_model_sln.cpp
)

target_link_libraries(parallel_model
  PRIVATE
  src
  project_warnings
  project_options
)