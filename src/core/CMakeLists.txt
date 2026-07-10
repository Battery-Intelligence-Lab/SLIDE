# v4 core (slide::core) — grows alongside legacy src/ (strangler pattern, PLAN.md §5).
# The registry is compiled once so downstream callers do not instantiate every common composition.
include(${PROJECT_SOURCE_DIR}/cmake/SlideCoreTarget.cmake)
slide_add_core_library(slide_core)

if(SLIDE_WITH_CUDA)
  # Keep CUDA language metadata and device-runtime linkage inside one private
  # shared boundary. This prevents nvcc's Windows device libraries leaking as
  # bare MSVC-style names into Clang-linked CPU consumers.
  add_library(slide_cuda_runtime SHARED CudaSpmRuntime.cu)
  target_include_directories(slide_cuda_runtime PRIVATE ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(slide_cuda_runtime PRIVATE CUDA::cudart)
  set_target_properties(slide_cuda_runtime PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
    CUDA_EXTENSIONS OFF
    CUDA_RUNTIME_LIBRARY Shared)
  # CUDA is an implementation detail behind CudaSpmBatch's CPU-only pimpl.
  # Do not leak backend metadata into unrelated consumers of slide_core.
  target_compile_definitions(slide_core PRIVATE SLIDE_WITH_CUDA=1)
  target_link_libraries(slide_core PRIVATE slide_cuda_runtime)
endif()

option(SLIDE_WITH_ARROW "Enable the optional Arrow/Parquet recorder sink" OFF)
if(SLIDE_WITH_ARROW)
  # Deliberately find-package-first: Arrow-via-CPM is not portable on Windows+Clang.
  find_package(Arrow CONFIG REQUIRED)
  find_package(Parquet CONFIG REQUIRED)
  if(TARGET Arrow::arrow_shared)
    set(SLIDE_ARROW_TARGET Arrow::arrow_shared)
  elseif(TARGET Arrow::arrow_static)
    set(SLIDE_ARROW_TARGET Arrow::arrow_static)
  else()
    message(FATAL_ERROR "Arrow package found without arrow_shared/arrow_static target")
  endif()
  if(TARGET Parquet::parquet_shared)
    set(SLIDE_PARQUET_TARGET Parquet::parquet_shared)
  elseif(TARGET Parquet::parquet_static)
    set(SLIDE_PARQUET_TARGET Parquet::parquet_static)
  else()
    message(FATAL_ERROR "Parquet package found without parquet_shared/parquet_static target")
  endif()
  target_compile_definitions(slide_core PUBLIC SLIDE_WITH_ARROW=1)
  target_link_libraries(slide_core PUBLIC ${SLIDE_ARROW_TARGET} ${SLIDE_PARQUET_TARGET})
endif()

option(SLIDE_WITH_ZSTD "Enable shuffled zstd async recording blocks" OFF)
if(SLIDE_WITH_ZSTD)
  find_package(zstd CONFIG QUIET)
  if(TARGET zstd::libzstd_static)
    set(SLIDE_ZSTD_TARGET zstd::libzstd_static)
  elseif(TARGET zstd::libzstd_shared)
    set(SLIDE_ZSTD_TARGET zstd::libzstd_shared)
  elseif(TARGET libzstd_static)
    set(SLIDE_ZSTD_TARGET libzstd_static)
  else()
    CPMAddPackage(
      NAME zstd
      GITHUB_REPOSITORY facebook/zstd
      GIT_TAG v1.5.7
      SOURCE_SUBDIR build/cmake
      OPTIONS
        "ZSTD_BUILD_PROGRAMS OFF"
        "ZSTD_BUILD_SHARED OFF"
        "ZSTD_BUILD_STATIC ON"
        "ZSTD_BUILD_TESTS OFF"
    )
    set(SLIDE_ZSTD_TARGET libzstd_static)
  endif()
  if(NOT TARGET ${SLIDE_ZSTD_TARGET})
    message(FATAL_ERROR "SLIDE_WITH_ZSTD=ON but no zstd library target is available")
  endif()
  target_compile_definitions(slide_core PUBLIC SLIDE_WITH_ZSTD=1)
  target_link_libraries(slide_core PRIVATE ${SLIDE_ZSTD_TARGET})
endif()
