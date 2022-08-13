#!/bin/bash

# exit if any command fails
set -e

cmake -S test -B build/test
cmake --build build/test -j`nproc`
CTEST_OUTPUT_ON_FAILURE=1 cmake --build build/test --target test
