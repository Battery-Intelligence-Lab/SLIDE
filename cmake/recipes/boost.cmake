if(TARGET boost_headers)
    return()
endif()

message(STATUS "Third-party: creating targets 'boost_headers'...")

# Try to find system Boost first (e.g., from Homebrew)
find_package(Boost QUIET)

if(Boost_FOUND)
    message(STATUS "Found system Boost: ${Boost_VERSION}")
    # Create an interface target for consistency
    add_library(boost_headers INTERFACE)
    target_include_directories(boost_headers INTERFACE ${Boost_INCLUDE_DIRS})
    # Create alias for compatibility
    return()
endif()

# If system Boost not found, download header-only version
message(STATUS "System Boost not found, downloading header-only version...")

set(BOOST_URL "https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.gz" CACHE STRING "Boost download URL")
set(BOOST_URL_SHA256 "9de758db755e8330a01d995b0a24d09798048400ac25c03fc5ea9be364b13c93" CACHE STRING "Boost download URL SHA256 checksum")

include(CPM)
CPMAddPackage(
    NAME boost
    URL ${BOOST_URL}
    DOWNLOAD_ONLY ON
)

# Create interface target for header-only Boost
add_library(boost_headers INTERFACE)
target_include_directories(boost_headers INTERFACE ${boost_SOURCE_DIR})

# Create alias for compatibility
add_library(Boost::boost ALIAS boost_headers)

set(Boost_POPULATED ON)
