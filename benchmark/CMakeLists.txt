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

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
    set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)
endif()

# Cant use CPM because we need the whole CUDA toolkit, which contains thrust and other needed
# headers find_package(Thrust REQUIRED) thrust_create_target(ThrustOMP HOST OMP)

CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")

CPMAddPackage(NAME MPQsort SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME} ${sources})

target_link_libraries(${PROJECT_NAME} benchmark::benchmark MPQsort::MPQsort)

# Link OpenMP dependencies
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
