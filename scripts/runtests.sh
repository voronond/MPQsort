#!/bin/bash

# exit if any command fails
set -e

cmake -S test -B build/test
cmake --build build/test -j`nproc`
./build/test/MPQsortTests --use-colour yes "$@"
