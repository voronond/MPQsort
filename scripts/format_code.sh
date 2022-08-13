#!/bin/bash

cmake -S test -B build/test

# format code
cmake --build build --target fix-format
