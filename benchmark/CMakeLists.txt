cmake_minimum_required(VERSION 3.14...3.22)

project(MPQsortBenchmarks LANGUAGES CXX)

# ---- Options ----

# option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)

CPMAddPackage(
    NAME benchmark
    GITHUB_REPOSITORY google/benchmark
    VERSION 1.7.0
    OPTIONS "BENCHMARK_ENABLE_TESTING Off"
)
# Check if TBB installed and if not, build from sources
# TBB is not currenty compatible with GCC 12 so it would be better to download binary and extract that to /usr/local and rename folder to tbb
# Tutorial: https://stackoverflow.com/questions/36534706/trouble-linking-intel-tbb-library
CPMFindPackage(
    NAME TBB
    GITHUB_REPOSITORY oneapi-src/oneTBB
    VERSION 2021.5.0
)

# Cant use CPM because we need the whole CUDA toolkit, which contains thrust and other needed headers
find_package(Thrust REQUIRED)
thrust_create_target(ThrustOMP HOST OMP)

CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")

CPMAddPackage(NAME MPQsort SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME} ${sources})

target_link_libraries(${PROJECT_NAME} Thrust::OMP::Host TBB::tbb benchmark::benchmark MPQsort::MPQsort)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
