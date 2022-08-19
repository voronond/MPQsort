#!/bin/bash

# Exit when any command fails
set -e

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

echo "Building all..."
cmake -S all -B build
cmake --build build -j`nproc`
echo "Building done..."

echo "Running tests..."
CTEST_OUTPUT_ON_FAILURE=1 cmake --build build/test --target test
echo "Running tests done..."

echo "Formatting code..."
cmake --build build --target fix-format
echo "Formatting code done..."

echo "Building docs..."
cmake --build build --target GenerateDocs
echo "Building docs done"
