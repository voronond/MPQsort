cmake_minimum_required(VERSION 3.14...3.22)

project(MPQsortBenchmarks LANGUAGES CXX)

# ---- Options ----

option(TEST_INSTALLED_VERSION "Test the version found by find_package" OFF)

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
CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")

if(TEST_INSTALLED_VERSION)
    find_package(MPQsort REQUIRED)
else()
    CPMAddPackage(NAME MPQsort SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
endif()

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} benchmark::benchmark MPQsort::MPQsort)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)

# enable compiler warnings
if(NOT TEST_INSTALLED_VERSION)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(MPQsort PUBLIC -Wall -Wpedantic -Wextra -Werror)
    elseif(MSVC)
        target_compile_options(MPQsort PUBLIC /W4 /WX)
    endif()
endif()
