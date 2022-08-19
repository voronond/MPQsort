#!/bin/bash

# exit if any command fails
set -e

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

cmake -S test -B build/test
cmake --build build/test -j`nproc`
./build/test/MPQsortTests --use-colour yes "$@"
