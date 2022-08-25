#!/bin/bash

# Check required python version
python -c 'import sys; exit(0) if (sys.version_info.major >= 3 and sys.version_info.minor >= 3) or sys.version_info.major > 3 else exit(1)'

if [ $? -ne 0 ]
then
    echo "Requires Python version 3.3 or newer."
    exit 1
fi

# Create vnevn if not exists
if [[ ! -a build/venv && -d build/venv ]]
then
    python -m venv build/venv
    # Activate
    source build/venv/bin/activate
    # Install requirements
    pip3 install -r build/_deps/benchmark-src/tools/requirements.txt
    pip3 install pandas
else
    # Activate
    source build/venv/bin/activate
fi

# Run script
python build/_deps/benchmark-src/tools/compare.py "$@"

# Deactivate
deactivate
