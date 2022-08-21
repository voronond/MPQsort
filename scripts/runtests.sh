#!/bin/bash

# exit if any command fails
set -e

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

if [[ "$1" == "sanitizers" ]]
then
    # Remove sanitizers from arguments
    shift

    echo "Running with address sanitizer..."
    cmake -S test -B build/test/address -DUSE_SANITIZER=Address
    cmake --build build/test/address -j`nproc`
    ./build/test/address/MPQsortTests --use-colour yes "$@"

    echo "Running with thread sanitizer..."
    cmake -S test -B build/test/thread -DUSE_SANITIZER=Thread
    cmake --build build/test/thread -j`nproc`
    ./build/test/thread/MPQsortTests --use-colour yes "$@"

    echo "Running with memory sanitizer..."
    cmake -S test -B build/test/memory -DUSE_SANITIZER=Memory
    cmake --build build/test/memory -j`nproc`
    ./build/test/memory/MPQsortTests --use-colour yes "$@"

    exit 0
fi

cmake -S test -B build/test
cmake --build build/test -j`nproc`
./build/test/MPQsortTests --use-colour yes "$@"
