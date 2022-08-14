cmake_minimum_required(VERSION 3.14...3.22)

project(MPQsortTests LANGUAGES CXX)

# ---- Options ----

option(ENABLE_TEST_COVERAGE "Enable test coverage" OFF)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)

CPMAddPackage("gh:catchorg/Catch2@2.13.9")
CPMAddPackage("gh:TheLartians/Format.cmake@1.7.3")
CPMAddPackage(NAME MPQsort SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME} ${sources})
target_link_libraries(${PROJECT_NAME} Catch2::Catch2 MPQsort::MPQsort)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)

# enable compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(MPQsort PUBLIC -Wall -Wpedantic -Wextra -Werror)
elseif(MSVC)
    target_compile_options(MPQsort PUBLIC /W4 /WX)
endif()

# ---- Add MPQsortTests ----

enable_testing()

# Note: doctest and similar testing frameworks can automatically configure CMake tests. For other
# testing frameworks add the tests target instead: add_test(NAME ${PROJECT_NAME} COMMAND
# ${PROJECT_NAME})

include(${Catch2_SOURCE_DIR}/contrib/Catch.cmake)
catch_discover_tests(${PROJECT_NAME})

# ---- code coverage ----

if(ENABLE_TEST_COVERAGE)
    target_compile_options(MPQsort PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
    target_link_options(MPQsort PUBLIC -fprofile-arcs -ftest-coverage)
endif()
