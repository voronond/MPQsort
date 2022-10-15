#!/bin/bash

# Check required python version
python -c 'import sys; exit(0) if (sys.version_info.major >= 3 and sys.version_info.minor >= 3) or sys.version_info.major > 3 else exit(1)'

if [ $? -ne 0 ]
then
    echo "Requires Python version 3.3 or newer."
    exit 1
fi

set -e

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

function build {
    cmake -S benchmark -B build/benchmark -DCMAKE_BUILD_TYPE=Release
    cmake --build build/benchmark --config Release -j`nproc`
}

build

# Set deprecated OpenMP nested paralelism for GNU algoritms
export OMP_NESTED=true
TOOLS_DIR=$HOME/.cache/CPM/benchmark/864e6ccb444b9a631c8b0950e9e6d1a05fb7abca/tools

# Create vnevn if not exists
if [[ ! -a build/venv ]]
then
    python -m venv build/venv
    # Activate
    source build/venv/bin/activate
    # Install requirements
    pip3 install -r ${TOOLS_DIR}/requirements.txt
    pip3 install pandas
else
    # Activate
    source build/venv/bin/activate
fi

# Run script
python ${TOOLS_DIR}/compare.py "$@" --benchmark_time_unit=ms

# Deactivate
deactivate
unset OMP_NESTED
