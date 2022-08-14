#!/bin/bash

cmake -S test -B build/test

# format code
cmake --build build/test --target fix-format
