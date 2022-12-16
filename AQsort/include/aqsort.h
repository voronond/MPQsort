#ifndef AQSORT_AQSORT_H
#define AQSORT_AQSORT_H

#include <cstddef>

// threshold for using insert sort instead of quick sort in sequential sorting
#ifndef AQSORT_INSERTION_SORT_THRESHOLD 
#define AQSORT_INSERTION_SORT_THRESHOLD 16
#endif

#include "impl/sequential_sort.h"

// parallel sorting stuff:
#ifdef _OPENMP

// blocks size used by paralell partitioning algorithm
#ifndef AQSORT_PARALLEL_PARTITION_BLOCK_SIZE
#define AQSORT_PARALLEL_PARTITION_BLOCK_SIZE 1024
#endif

#include "impl/parallel_sort.h"

#endif

namespace aqsort
{
#ifdef _OPENMP
    template <typename Comp, typename Swap>
    inline void parallel_sort(std::size_t length, Comp* const comp, Swap* const swap)
    {
        impl::parallel_sort(length, comp, swap);
    }
#endif

    template <typename Comp, typename Swap>
    inline void sequential_sort(std::size_t length, Comp* const comp, Swap* const swap)
    {
        impl::sequential_sort(length, comp, swap);
    }

    template <typename Comp, typename Swap>
    inline void sort(std::size_t length, Comp* const comp, Swap* const swap)
    {
#ifdef _OPENMP
        parallel_sort(length, comp, swap);
#else
        sequential_sort(length, comp, swap);
#endif
    }
}

#endif
