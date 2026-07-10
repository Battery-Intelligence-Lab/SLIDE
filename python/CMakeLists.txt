find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)

# nanobind recommends asking the active interpreter for its CMake package path;
# this works both in scikit-build isolation and in a developer virtualenv.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE nanobind_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_slide_core bindings.cpp)
target_link_libraries(_slide_core PRIVATE slide_core)
target_compile_features(_slide_core PRIVATE cxx_std_20)

install(TARGETS _slide_core LIBRARY DESTINATION slide RUNTIME DESTINATION slide)
