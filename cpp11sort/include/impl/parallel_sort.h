#ifndef CPP11SORT_IMPL_PARALLEL_SORT_H
#define CPP11SORT_IMPL_PARALLEL_SORT_H

#include <algorithm>
#include <cstdlib>
#include <iterator>

#include "barrier.h"
#include "parallel_partition.h"
#include "thread_pool.h"

namespace cpp11sort
{
namespace impl
{

template<typename RAIter, typename Compare, typename ThreadManager>
typename std::iterator_traits<RAIter>::difference_type
ParallelSortQsDivide(RAIter first, RAIter last, Compare comp,
                     typename std::iterator_traits<RAIter>::difference_type pivot_rank,
                     std::size_t num_threads, ThreadManager &thread_manager) {
    typedef std::iterator_traits<RAIter> TraitsType;
    typedef typename TraitsType::value_type ValueType;
    typedef typename TraitsType::difference_type DifferenceType;

    DifferenceType n = last - first;
    DifferenceType
        num_samples = std::min((DifferenceType) CPP11SORT_PIVOT_SAMPLES, n);

    auto *samples = static_cast<ValueType *> (::operator new(
        num_samples * sizeof(ValueType)));

    for (DifferenceType s = 0; s < num_samples; ++s) {
        const unsigned long long index = (unsigned long long) (((s) * n) / num_samples);
        ::new(&(samples[s])) ValueType(first[index]);
    }

    std::sort(samples, samples + num_samples, comp);
    ValueType &pivot = samples[pivot_rank * num_samples / n];

    DifferenceType split;
    if (num_threads <= 1) {
        split = SequentialPartition(first, last, pivot, comp);
    } else {
        split = ParallelPartition(first, last, pivot, comp,
                                  num_threads, thread_manager);
    }
    for (DifferenceType s = 0; s < num_samples; ++s) {
        samples[s].~ValueType();
    }
    ::operator delete(samples);
    return split;
}

template<typename RAIter, typename Compare, typename ThreadManager>
void
ParallelSortQsConquer(RAIter first, RAIter last, Compare comp,
                      std::size_t num_threads, std::size_t sequential_threshold,
                      ThreadManager &thread_manager) {
    typedef std::iterator_traits<RAIter> TraitsType;
    typedef typename TraitsType::difference_type DifferenceType;

    DifferenceType n = last - first, pivot_rank = n / 2;
    num_threads = std::max<long>(num_threads, 1lu);
    if (n <= (int) sequential_threshold) {
        std::sort(first, last, comp);
        return;
    }

    DifferenceType
        split = ParallelSortQsDivide(first, last, comp, pivot_rank, num_threads,
                                     thread_manager);
    std::size_t num_threads_left = std::max<long>((num_threads * split) / n, 1lu);

    thread_manager.RunConquer(first, last, comp, split, num_threads, num_threads_left,
                              sequential_threshold, thread_manager);
}

template<typename RAIter, typename Compare, typename ThreadManager>
void ParallelSortQs(RAIter first, RAIter last, Compare comp, std::size_t num_threads,
                    std::size_t sequential_threshold, ThreadManager &thread_manager) {
    typedef std::iterator_traits<RAIter> _TraitsType;
    typedef typename _TraitsType::difference_type _DifferenceType;

    // If sequential threshold hasn't been set, set it to the default value
    if (sequential_threshold == 0) 
        sequential_threshold = CPP11SORT_SEQUENTIAL_THRESHOLD;
    
    thread_manager.StartRecursiveQs(first, last, comp, num_threads, sequential_threshold,
                                    thread_manager);
}

}

}

#endif
