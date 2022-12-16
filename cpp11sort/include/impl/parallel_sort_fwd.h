#ifndef CPP11SORT_IMPL_PARALLEL_SORT_FWD_H
#define CPP11SORT_IMPL_PARALLEL_SORT_FWD_H

#include <cstdlib>
#include <iterator>

namespace cpp11sort
{
namespace impl
{

/**
* Parallel sort divide step.
* Divides an array to the two parts, the first one contains elements smaller than
* pivot, the second one contains elements greater or equal to the pivot.
* @tparam RAIter To iterate elements of the array.
* @tparam Compare Type of comparing function.
* @param first The beginning of the data to be sorted.
* @param last The ending of the data to be sorted.
* @param comp Comparing function to determine order of two elements in the sorted
* array.
* @param pivot_rank Rank of the pivot in the array of pivot samples.
* @param num_threads Number of threads in the pool to execute the tasks.
* @return Position of the pivot after partitioning.
*/
template<typename RAIter, typename Compare, typename ThreadManager>
typename std::iterator_traits<RAIter>::difference_type
ParallelSortQsDivide(RAIter first, RAIter last, Compare comp,
                     typename std::iterator_traits<RAIter>::difference_type difference_type,
                     std::size_t num_threads,
                     ThreadManager &thread_manager);

/**
* Parallel quicksort conquere step.
* @tparam RAIter To iterate elements of the array.
* @tparam Compare Type of comparing function.
* @param first The beginning of the data to be sorted.
* @param last The ending of the data to be sorted.
* @param comp Comparing function to determine order of two elements in the sorted
* array.
* @param num_threads Number of threads in the pool to execute the tasks.
* @param sequential_threshold Arrays smaller than this size will be executed
* sequentially.
* @param thread_manager Specifies the thread managing.
*/
template<typename RAIter, typename Compare, typename ThreadManager>
void
ParallelSortQsConquer(RAIter first, RAIter last, Compare comp,
                      std::size_t num_threads, std::size_t sequential_threshold,
                      ThreadManager &thread_manager);
/**
* Entering function for sorting with parallel c++11 sort. Creates thread pools to be
* used during the sorting and spawns the first task for the pool.
* @tparam RAIter To iterate elements of the array.
* @tparam Compare Type of comparing function.
* @param first The beginning of the data to be sorted.
* @param last The ending of the data to be sorted.
* @param comp Comparing function to determine order of two elements in the sorted
* array.
* @param num_threads Number of threads in the pool to execute the tasks.
* @param sequential_threshold When the parallel partitioning should switch to the
* sequential one.
* @param thread_manager Specifies the thread managing.
*/
template<typename RAIter, typename Compare, typename ThreadManager>
void ParallelSortQs(RAIter first, RAIter last, Compare comp,
                    std::size_t num_threads, std::size_t sequential_threshold,
                    ThreadManager &thread_manager);

}

}

#endif
