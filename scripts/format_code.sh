#!/bin/bash

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

cmake -S test -B build/format_code

# format code
cmake --build build/format_code --target fix-format
