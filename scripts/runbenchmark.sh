#!/bin/bash

# exit if any command fails
set -e

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

function build {
    cmake -S benchmark -B build/benchmark -DCMAKE_BUILD_TYPE=Release
    cmake --build build/benchmark --config Release -j`nproc`
}

if [[ "${1,,}" == "silent" ]]
then
    shift
    build &>/dev/null
else
    build
fi

./build/benchmark/MPQsortBenchmarks "$@"
