#ifndef CPP11SORT_IMPL_PARALLEL_PARTITION_FWD_H
#define CPP11SORT_IMPL_PARALLEL_PARTITION_FWD_H

#include <atomic>
#include <iterator>

namespace cpp11sort
{
namespace impl
{

/**
* Acquires blocks of the array to be partitioned. This function is meant to be executed
* by multiple threads simultaneously. Each thread takes one block from the left and one from
* the right. Rearranges elements so that in the left block there are elements smaller
* than pivot and in the right block are elements greater than pivot.
* @tparam RAIter To iterate elements of the array.
* @tparam Compare Type of comparing function.
* @tparam TraitsType
* @tparam ValueType Type of the element in the array.
* @tparam DifferenceType Type of difference of the two positions in the array.
* @param first Beginning of the array to be partitioned.
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
* @param comp Comparing function to tell the ordering of the elements.
* @param barrier Barrier to wait for all the threads in the critical sections where
* all threads must meet.
*/
template<typename RAIter, typename Compare, typename TraitsType = std::iterator_traits<
    RAIter>, typename ValueType = typename TraitsType::value_type,
    typename DifferenceType = typename TraitsType::difference_type>
void MakeBlocks(RAIter first, std::atomic<DifferenceType> &left,
                std::atomic<DifferenceType> &right, std::atomic<DifferenceType> &dist,
                std::atomic<DifferenceType> &leftover_left,
                std::atomic<DifferenceType> &leftover_right, ValueType pivot,
                DifferenceType n, std::atomic<int> *&reserved_left,
                std::atomic<int> *&reserved_right, Compare comp, Barrier &barrier);

}

}

#endif
