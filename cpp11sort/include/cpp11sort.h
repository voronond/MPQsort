#ifndef CPP11SORT_CPP11SORT_H
#define CPP11SORT_CPP11SORT_H

#include <cstddef>
#include <functional>
#include <iterator>
#include <thread>

// Threshold to switch to the sequential sort without further dividing to smaller tasks.
#ifndef CPP11SORT_SEQUENTIAL_THRESHOLD
#define CPP11SORT_SEQUENTIAL_THRESHOLD 50000
#endif

// Size of the block used in the parallel partitioning.
#ifndef CPP11SORT_PARTITION_BLOCK_SIZE
#define CPP11SORT_PARTITION_BLOCK_SIZE 6000
#endif

// Number of elements to pick the pivot from
#ifndef CPP11SORT_PIVOT_SAMPLES
#define CPP11SORT_PIVOT_SAMPLES 30
#endif

// need to be here after configuration macro parameters definition:
#include "impl/parallel_sort.h"
#include "impl/thread_manager.h"
#include "impl/thread_manager_pool.h"

namespace cpp11sort
{

template<typename RAIter, typename Comp>
void sort(RAIter first, RAIter last, Comp comp, unsigned int num_threads)
{
   impl::ThreadManagerPool thread_manager(num_threads);
   impl::ParallelSortQs(first, last, comp, num_threads, 0, thread_manager);
}

template<typename RAIter>
void sort(RAIter first, RAIter last)
{
   sort(first, last,
         std::less< typename std::iterator_traits<RAIter>::value_type >{}, 
         std::thread::hardware_concurrency() );
}

template<typename RAIter, typename Comp>
void sort(RAIter first, RAIter last, Comp comp)
{
   sort(first, last, comp, std::thread::hardware_concurrency());
}

template<typename RAIter>
void sort(RAIter first, RAIter last, unsigned int num_threads)
{
   sort(first, last,
         std::less< typename std::iterator_traits<RAIter>::value_type >{}, 
         num_threads);
}

}

#endif
