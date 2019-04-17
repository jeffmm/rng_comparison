cmake_minimum_required(VERSION 3.12.4)
project(rng_comparison)

# Set binary output directory to build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Add executables
add_subdirectory("src")

# Reset binary output directory to default
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Add catch
include_directories(lib/Catch2/single_include/catch2)

# Configure tests
enable_testing()
add_subdirectory(tests)

