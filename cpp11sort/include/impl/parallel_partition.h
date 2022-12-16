#ifndef CPP11SORT_IMPL_PARALLEL_PARTITION_H
#define CPP11SORT_IMPL_PARALLEL_PARTITION_H

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <utility>

#include "barrier.h"
#include "thread_manager.h"

namespace cpp11sort
{
namespace impl
{

/**
* Executes parallel partitioning of an array. Calls making blocks in parallel, waits
* for them and after that does the cleaning of the leftovers.
* @tparam RAIter To iterate elements of the array.
* @tparam Compare Type of comparing function.
* @param first The beginning of the data to be sorted.
* @param last The ending of the data to be sorted.
* @param pivot Value used to divide array to values smaller than this or greater and
 * equal to this.
* @param comp Comparing function to determine order of two elements in the sorted
* array.
* @param num_threads Number of threads used for the partitioning.
* @param threads_manager Defines thread managing policy.
* @return Position of the pivot after partitioning.
*/
template<typename RAIter, typename Compare, typename ThreadManager>
std::size_t ParallelPartition(RAIter first, RAIter last,
                              typename std::iterator_traits<RAIter>::value_type pivot,
                              Compare comp, std::size_t num_threads,
                              ThreadManager &thread_manager);
template<typename RAIter, typename Compare>
std::size_t ParallelPartition(RAIter first, RAIter last,
                              typename std::iterator_traits<RAIter>::value_type pivot,
                              Compare comp, std::size_t num_threads);

template<typename RAIter, typename Compare>
std::size_t ParallelPartition(RAIter first, RAIter last,
                              typename std::iterator_traits<RAIter>::value_type pivot,
                              Compare comp);

template<typename RAIter, typename Compare>
std::size_t SequentialPartition(RAIter first, RAIter last,
                                typename std::iterator_traits<RAIter>::value_type pivot,
                                Compare comp);

template<typename RAIter, typename Compare, typename TraitsType, typename ValueType,
    typename DifferenceType>
void MakeBlocks(RAIter first, std::atomic<DifferenceType> &left,
                std::atomic<DifferenceType> &right, std::atomic<DifferenceType> &dist,
                std::atomic<DifferenceType> &leftover_left,
                std::atomic<DifferenceType> &leftover_right, ValueType pivot,
                DifferenceType n, std::atomic<int> *&reserved_left,
                std::atomic<int> *&reserved_right, Compare comp, Barrier &barrier) {
    // Private variables for single thread.
    DifferenceType thread_left, thread_left_border, thread_right,
        thread_right_border;
    DifferenceType block_size = CPP11SORT_PARTITION_BLOCK_SIZE;
    thread_left = left + 1;
    thread_left_border = thread_left - 1;
    thread_right = n - 1;
    thread_right_border = thread_right + 1;

    // Acquire blocks from left and right until the whole array is partitioned.
    while (true) {
        if (thread_left > thread_left_border) {
            DifferenceType former_dist = std::atomic_fetch_add(&dist, -block_size);
            if (former_dist < block_size) {
                dist += block_size;
                break;
            } else {
                thread_left = std::atomic_fetch_add(&left, block_size);
                thread_left_border = thread_left + (block_size - 1);
            }
        }
        if (thread_right < thread_right_border) {
            DifferenceType former_dist = std::atomic_fetch_add(&dist,
                                                               -block_size);
            if (former_dist < block_size) {
                dist += block_size;
                break;
            } else {
                thread_right = std::atomic_fetch_add(&right,
                                                     -block_size);
                thread_right_border = thread_right - (block_size - 1);
            }
        }

        // Swap elements in the right and left blocks.
        while (thread_left < thread_right) {
            while (comp(first[thread_left], pivot) &&
                thread_left <= thread_left_border) {
                ++thread_left;
            }
            while (!comp(first[thread_right], pivot) &&
                thread_right >= thread_right_border) {
                --thread_right;
            }

            if (thread_left > thread_left_border ||
                thread_right < thread_right_border) {
                // Fetch new block.
                break;
            }

            std::iter_swap(first + thread_left, first + thread_right);
            ++thread_left;
            --thread_right;
        }
    }

    // Swap the leftover blocks to the right places.
    if (thread_left <= thread_left_border) {
        leftover_left += 1;
    }
    if (thread_right >= thread_right_border) {
        leftover_right += 1;
    }
    barrier.Wait();

    DifferenceType left_old = left, left_new =
        left -
            leftover_left *
                block_size, right_old = right, right_new =
        right + leftover_right * block_size;

    if (thread_left <= thread_left_border &&
        thread_left_border >= left_new) {
        // Block already in place, reserve spot.
        reserved_left[(left - (thread_left_border + 1)) /
            block_size] = 1;
    }

    if (thread_right >= thread_right_border &&
        thread_right_border <= right_new) {
        // Block already in place, reserve spot.
        reserved_right[((thread_right_border - 1) - right) /
            block_size] = 1;
    }
    barrier.Wait();

    if (thread_left <= thread_left_border &&
        thread_left_border < left_new) {
        // Find spot and swap.
        DifferenceType swap_start = -1;
        for (int r = 0; r < leftover_left; ++r) {
            int zero_constant = 0, one_constant = 1;
            if (reserved_left[r] == 0 &&
                reserved_left[r].compare_exchange_strong(zero_constant,
                                                         one_constant)) {
                swap_start = left_old - (r + 1) * block_size;
                break;
            }
        }
        std::swap_ranges(first + thread_left_border - (block_size - 1),
                         first + thread_left_border + 1,
                         first + swap_start);
    }

    if (thread_right >= thread_right_border &&
        thread_right_border > right_new) {
        // Find spot and swap.
        DifferenceType swap_start = -1;
        for (int r = 0; r < leftover_right; ++r) {
            int zero_constant = 0, one_constant = 1;
            if (reserved_right[r] == 0 &&
                reserved_right[r].compare_exchange_strong(zero_constant,
                                                          one_constant)) {
                swap_start = right_old + r * block_size + 1;
                break;
            }
        }
        std::swap_ranges(first + thread_right_border,
                         first + thread_right_border + block_size,
                         first + swap_start);
    }

    left.store(left_new);
    right.store(right_new);
    dist.store(right.load() - left.load() + 1);
}

template<typename RAIter, typename Compare, typename ThreadManager>
std::size_t ParallelPartition(RAIter first, RAIter last,
                              typename std::iterator_traits<RAIter>::value_type pivot,
                              Compare comp, std::size_t num_threads,
                              ThreadManager &thread_manager) {
    typedef std::iterator_traits<RAIter> TraitsType;
    typedef typename TraitsType::difference_type DifferenceType;
    DifferenceType n = last - first;
    std::atomic<DifferenceType> left{0}, right{n - 1},
        dist{n}, leftover_left, leftover_right;
    // Values in the arrays can be 0 or 1, but it is int type to allow atomic operation.
    std::atomic<int> *reserved_left{nullptr}, *reserved_right{nullptr};
    static const DifferenceType block_size = CPP11SORT_PARTITION_BLOCK_SIZE;
    std::size_t block_threads = thread_manager.NumberOfBlockThreads(num_threads);
    // At least two chunks per thread.
    if (dist >= (long signed int) (2 * block_threads * block_size)) {
        reserved_left = new std::atomic<int>[block_threads];
        reserved_right = new std::atomic<int>[block_threads];
        // Make blocks in parallel.
        while (dist >= (long signed int) (2 * block_threads * block_size)) {
            for (std::size_t r = 0;
                 r < block_threads; ++r) {
                reserved_left[r] = 0; //false
                reserved_right[r] = 0; //false
            }
            leftover_left = 0;
            leftover_right = 0;
            thread_manager.RunMakingBlocks(first, std::ref(left), std::ref(right),
                                           std::ref(dist), std::ref(leftover_left),
                                           std::ref(leftover_right), pivot, n,
                                           std::ref(reserved_left),
                                           std::ref(reserved_right), comp, block_threads);
        }
    } // End of parallel section.
    DifferenceType final_left = left, final_right = right;
    while (final_left < final_right) {
        // Go right until key is geq than pivot.
        while (comp(first[final_left], pivot) &&
            final_left < final_right) {
            ++final_left;
        }

        // Go left until key is less han pivot
        while (!comp(first[final_right], pivot) &&
            final_left < final_right) {
            --final_right;
        }

        if (final_left == final_right) {
            break;
        }
        std::iter_swap(first + final_left, first + final_right);
        ++final_left;
        --final_right;
    }


    // All elements on the left side are < pivot, all elements on the right are >= pivot
    delete[] reserved_left;
    delete[] reserved_right;

    // Elements "between" final_left and final_right might not have been regarded yet
    if (final_left < n && !
        comp(first[final_left], pivot
        )) {
        // Really swapped
        return final_left;
    } else {
        return final_left + 1;
    }
}

template<typename RAIter, typename Compare>
std::size_t ParallelPartition(RAIter first, RAIter last,
                              typename std::iterator_traits<RAIter>::value_type pivot,
                              Compare comp, std::size_t num_threads) {
    ThreadManager thread_manager;
    return ParallelPartition(first, last, pivot, comp, num_threads, thread_manager);
}

template<typename RAIter, typename Compare>
std::size_t ParallelPartition(RAIter first, RAIter last,
                              typename std::iterator_traits<RAIter>::value_type pivot,
                              Compare comp) {
    std::size_t num_threads = std::thread::hardware_concurrency();
    return ParallelPartition(first, last, pivot, comp, num_threads);
}

template<typename RAIter, typename Compare>
std::size_t SequentialPartition(RAIter first, RAIter last,
                                typename std::iterator_traits<RAIter>::value_type pivot,
                                Compare comp) {
    auto it1 = first - 1, it2 = last;

    while (true) {
        do {
            ++it1;
        } while (*it1 < pivot);

        do {
            --it2;
        } while (*it2 > pivot);

        if (it1 >= it2) {
            return it2 - first + 1;
        }
        std::swap(*it1, *it2);
    }
}

}

}

#endif
