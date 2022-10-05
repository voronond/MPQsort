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

# Set deprecated OpenMP nested paralelism for GNU algoritms
export OMP_NESTED=true
./build/benchmark/MPQsortBenchmarks --benchmark_time_unit=ms "$@"
unset OMP_NESTED
