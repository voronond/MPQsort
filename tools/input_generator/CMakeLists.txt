cmake_minimum_required(VERSION 3.14...3.22)

project(MPQsortInputGenerator LANGUAGES CXX)

# --- Import tools ----

include(../../cmake/tools.cmake)

# ---- Dependencies ----

include(../../cmake/CPM.cmake)

CPMAddPackage(
    GITHUB_REPOSITORY jarro2783/cxxopts
    VERSION 3.0.0
    OPTIONS "CXXOPTS_BUILD_EXAMPLES NO" "CXXOPTS_BUILD_TESTS NO" "CXXOPTS_ENABLE_INSTALL YES"
)

# ---- Create binary ----

file(GLOB sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
add_executable(${PROJECT_NAME} ${sources})

target_link_libraries(${PROJECT_NAME} cxxopts)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
