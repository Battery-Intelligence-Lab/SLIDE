
add_executable(benchmark_Cell_Bucket benchmark_Cell_Bucket.cpp)

target_link_libraries(benchmark_Cell_Bucket
  PRIVATE
  src
  project_warnings
  project_options
  Boost::boost
)

add_executable(benchmark_Cell_ECM benchmark_Cell_ECM.cpp)

target_link_libraries(benchmark_Cell_ECM
  PRIVATE
  src
  project_warnings
  project_options
  Boost::boost
)

add_executable(benchmark_Cell_SPM benchmark_Cell_SPM.cpp)

target_link_libraries(benchmark_Cell_SPM
  PRIVATE
  src
  project_warnings
  project_options
  Boost::boost
)

add_executable(benchmark_LP_cases benchmark_LP_cases.cpp)

target_link_libraries(benchmark_LP_cases
  PRIVATE
  src
  project_warnings
  project_options
  Boost::boost
)