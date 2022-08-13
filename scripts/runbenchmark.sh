#!/bin/bash

# exit if any command fails
set -e

function build {
    cmake -S benchmark -B build/benchmark -DCMAKE_BUILD_TYPE=Release
    cmake --build build/benchmark --config Release -j`nproc`
}

if [[ "${1,,}" == "silent" ]]
then
    build &>/dev/null
else
    build
fi

./build/benchmark/MPQsortBenchmarks
