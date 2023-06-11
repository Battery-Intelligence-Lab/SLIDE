# This cmake file is to add external projects. 

include(cmake/CPM.cmake)

CPMAddPackage(
  NAME range-v3
  URL "https://github.com/ericniebler/range-v3/archive/refs/tags/0.12.0.tar.gz"
  DOWNLOAD_ONLY YES
  )

CPMAddPackage(
  NAME eigen
  URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
  DOWNLOAD_ONLY YES
)

add_library(eigen INTERFACE)
target_include_directories(eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})

add_library(range-v3 INTERFACE)
target_include_directories(range-v3 SYSTEM INTERFACE ${range-v3_SOURCE_DIR}/include)