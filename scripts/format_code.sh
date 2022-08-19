#!/bin/bash

export CPM_SOURCE_CACHE=$HOME/.cache/CPM

cmake -S test -B build/test

# format code
cmake --build build/test --target fix-format
