#!/bin/bash

# exit if any command fails
set -e

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

cmake -S tools/input_generator -B build/tools/input_generator
cmake --build build/tools/input_generator -j`nproc`
./build/tools/input_generator/MPQsortInputGenerator "$@"
