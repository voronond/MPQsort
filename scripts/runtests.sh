#!/bin/bash

# exit if any command fails
set -e

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

if [[ "$1" == "sanitizers" ]]
then
    # Remove sanitizers from arguments
    shift

    echo "Running with address sanitizer..."
    cmake -S test -B build/test/address -DUSE_SANITIZER=Address -DENABLE_COMPILER_WARNINGS=ON
    cmake --build build/test/address -j`nproc`
    ./build/test/address/MPQsortTests --use-colour yes "$@"

    exit 0
fi

cmake -S test -B build/test -DENABLE_COMPILER_WARNINGS=ON
cmake --build build/test -j`nproc`
# Break at first fail
./build/test/MPQsortTests -a --use-colour yes "$@"
