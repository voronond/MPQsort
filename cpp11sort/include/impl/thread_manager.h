#ifndef CPP11SORT_IMPL_THREAD_MANAGER_H
#define CPP11SORT_IMPL_THREAD_MANAGER_H

#include <atomic>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <thread>

#include "barrier.h"
#include "parallel_partition_fwd.h"
#include "parallel_sort_fwd.h"

namespace cpp11sort
{
namespace impl
{

/**
 * Determines the thread management policy for parallel sorting and partitioning. This
 * is a default thread manager that creates many threads and do not reuse them.
 * Individual methods defines how the sort should behave in a thread managing specific
 * parts.
 */
class ThreadManager {
  public:
    /**
     * Starts recursive quicksort by calling the first conquer step.
     * @tparam RAIter To iterate elements of the array.
     * @tparam Compare Compare Type of comparing function.
     * @tparam ThreadManager Class managing threads.
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
    void StartRecursiveQs(RAIter first, RAIter last, Compare comp,
                          std::size_t num_threads, std::size_t sequential_threshold,
                          ThreadManager &thread_manager) {
        ParallelSortQsConquer(first, last, comp, num_threads,
                              sequential_threshold, thread_manager);
    }
    /**
     * Runs conquer step on the left and on the right partition of the partitioned array..
     * @tparam RAIter To iterate elements of the array.
     * @tparam Compare Compare Type of comparing function.
     * @tparam ThreadManager Class managing threads.
     * @param first The beginning of the data to be sorted.
     * @param last The ending of the data to be sorted.
     * @param comp Comparing function to determine order of two elements in the sorted
     * array.
     * @param split Position of division of the array partitions.
     * @param num_threads Number of threads in the pool to execute the tasks.
     * @param num_threads_left Number of threads to work on the left partition.
     * @param sequential_threshold Arrays smaller than this size will be executed
     * sequentially.
     * @param thread_manager
     */
    template<typename RAIter, typename Compare, typename ThreadManager>
    void RunConquer(RAIter first, RAIter last, Compare comp,
                    typename std::iterator_traits<RAIter>::difference_type split,
                    std::size_t num_threads, std::size_t num_threads_left,
                    std::size_t sequential_threshold, ThreadManager &thread_manager) {
        std::thread t1(&ParallelSortQsConquer<RAIter, Compare, ThreadManager>, first,
                       first + split, comp, num_threads_left,
                       sequential_threshold, std::ref(thread_manager));
        ParallelSortQsConquer(first + split, last, comp,
                              num_threads - num_threads_left,
                              sequential_threshold, thread_manager);

        t1.join();
    }

    /**
     * Determines how many threads are used for the blocks creating.
     * @param num_threads Number of threads left for the job.
     * @return How many threads are used for the blocks building.
     */
    std::size_t NumberOfBlockThreads(std::size_t num_threads) {
        num_threads = std::min(std::thread::hardware_concurrency(), (unsigned int)
        num_threads);
        return num_threads;
    }
    /**
     * Starts running making blocks in parallel.
     * @tparam RAIter To iterate elements of the array.
     * @tparam Compare Compare Type of comparing function.
     * @param first The beginning of the data to be sorted.
     * @param left Position of the left most not already acquired element.
     * @param right Position of the rightmost not already acquired element.
     * @param dist Number of the elements that haven't been acquired so far.
     * @param leftover_left Number of unfinished left blocks.
     * @param leftover_right Number of unfinished right blocks.
     * @param pivot Divides elements to the ones smaller or bigger than it.
     * @param n Number of elements to be partitioned.
     * @param reserved_left Determines if block in the area of unfinished left blocks is on
     * the right position.
     * @param reserved_right Determines if block in the area of unfinished right blocks
     * is on the right position.
     * @param comp Comparing function to determine order of two elements in the sorted
     * array.
     * @param block_threads Number of threads used for the building blocks.
     */
    template<typename RAIter, typename Compare>
    void RunMakingBlocks(RAIter first,
                         std::atomic<typename std::iterator_traits<RAIter>::difference_type> &left,
                         std::atomic<typename std::iterator_traits<RAIter>::difference_type> &right,
                         std::atomic<typename std::iterator_traits<RAIter>::difference_type> &dist,
                         std::atomic<typename std::iterator_traits<RAIter>::difference_type> &leftover_left,
                         std::atomic<typename std::iterator_traits<RAIter>::difference_type> &leftover_right,
                         typename std::iterator_traits<RAIter>::value_type pivot,
                         typename std::iterator_traits<RAIter>::difference_type n,
                         std::atomic<int> *&reserved_left,
                         std::atomic<int> *&reserved_right,
                         Compare comp, std::size_t block_threads) {
        auto *t = new std::thread[block_threads];
        Barrier barrier(block_threads);
        for (std::size_t i = 0; i < block_threads; ++i) {
            t[i] = std::thread(&MakeBlocks<RAIter, Compare>, first, std::ref(left),
                               std::ref(right), std::ref(dist), std::ref(leftover_left),
                               std::ref(leftover_right), pivot, n,
                               std::ref(reserved_left), std::ref(reserved_right), comp,
                               std::ref(barrier));

        }
        for (std::size_t i = 0; i < block_threads; ++i) {
            t[i].join();
        }
        delete[] t;
    }

};

}

}

#endif
