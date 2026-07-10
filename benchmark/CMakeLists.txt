
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

add_executable(benchmark_PAY1_spm benchmark_PAY1_spm.cpp)
target_link_libraries(benchmark_PAY1_spm
  PRIVATE
  src
  slide_core
  project_warnings
  project_options
)

add_executable(benchmark_PAY2_pack benchmark_PAY2_pack.cpp)
target_link_libraries(benchmark_PAY2_pack
  PRIVATE
  src
  slide_core
  project_warnings
  project_options
)

add_executable(benchmark_PAY4_slide_pack benchmark_PAY4_slide_pack.cpp)
target_link_libraries(benchmark_PAY4_slide_pack
  PRIVATE
  slide_core
  project_warnings
  project_options
)

add_executable(benchmark_PAY4_slide_single benchmark_PAY4_slide_single.cpp)
target_link_libraries(benchmark_PAY4_slide_single
  PRIVATE
  slide_core
  project_warnings
  project_options
)

if(SLIDE_WITH_CUDA)
  add_executable(benchmark_PAY5_cuda benchmark_PAY5_cuda.cpp)
  target_link_libraries(benchmark_PAY5_cuda
    PRIVATE
    slide_core
    project_warnings
    project_options
  )
endif()
