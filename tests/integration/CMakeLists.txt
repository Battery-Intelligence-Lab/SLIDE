cmake_minimum_required(VERSION 3.21)

add_executable(parallel_model
  parallel_model_sln.cpp
)

target_link_libraries(parallel_model
  PRIVATE
  src
  project_warnings
  project_options
  Boost::boost
)

add_executable(parallel_model_slide
  parallel_model_sln_slide.cpp
)

target_link_libraries(parallel_model_slide
  PRIVATE
  src
  project_warnings
  project_options
  Boost::boost
)