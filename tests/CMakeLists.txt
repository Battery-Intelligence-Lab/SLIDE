cmake_minimum_required(VERSION 3.17)


add_subdirectory(unit)

add_executable(tests tests_main.cpp)


target_link_libraries(tests
  PRIVATE
  src
  unit
  )

add_test(slide_tests tests)