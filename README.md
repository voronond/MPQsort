[![Actions Status](https://github.com/TheLartians/ModernCppStarter/workflows/MacOS/badge.svg)](https://github.com/voronond/MPQsort/actions)
[![Actions Status](https://github.com/TheLartians/ModernCppStarter/workflows/Windows/badge.svg)](https://github.com/voronond/MPQsort/actions)
[![Actions Status](https://github.com/TheLartians/ModernCppStarter/workflows/Ubuntu/badge.svg)](https://github.com/voronond/MPQsort/actions)
[![Actions Status](https://github.com/TheLartians/ModernCppStarter/workflows/Style/badge.svg)](https://github.com/voronond/MPQsort/actions)
[![Actions Status](https://github.com/TheLartians/ModernCppStarter/workflows/Install/badge.svg)](https://github.com/voronond/MPQsort/actions)

# Multiway parallel Quicksort algorithm

Master thesis which addresses a problem of multiway parallel Quicksort algorithm. The goal was to design and implement parallel partitioning with multiple pivots, because the parallel implementations until now use only one pivot (Hoare's partitioning) to do so. The parallelization is implemented using OpenMP standard.

## Usage

This library is implemented as header-only. So if you want to use the algorithm in your application, the only thing you need to do is including file in *include/mpqsort/mpqsort.h* in your application.
