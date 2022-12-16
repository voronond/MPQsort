#ifndef AQSORT_IMPL_SEQUENTIAL_SORT_H
#define AQSORT_IMPL_SEQUENTIAL_SORT_H

#include <cassert>
#include <cmath>
#include <cstddef>

#include "select_pivot.h"
#include "sequential_partition.h"

namespace aqsort
{
    namespace impl
    {
        template<typename Comp, typename Swap>
        inline void insertion_sort(std::size_t start, std::size_t n, Comp* const comp, Swap* const swap)
        {
            // avoid unnecessary calling performance penalty
            assert (n > 1);
/*
            // find minimal element
            std::size_t temp = 0;
            for (std::size_t i = 1; i < n; ++i)
                if ((*comp)(start + i, start + temp)) 
                    temp = i;

            (*swap)(start, start + temp);

            temp = 0;
            while (++temp < n)
                for (std::size_t j = temp; (*comp)(start + j, start + j - 1); j--) 
                    (*swap)(start + j - 1, start + j);
*/
            for (std::size_t i = 1; i < n; i++) {
                std::size_t j = i;
                while ((j > 0) && ((*comp)(start + j, start + j - 1))) {
                    (*swap)(start + j - 1, start + j);
                    j--;
                }
            }
        }

        template<typename Comp, typename Swap>
        void heap_sort(std::size_t start, std::size_t n, Comp* const comp, Swap* const swap);

        template<typename Comp, typename Swap>
        void quick_sort(std::size_t start, std::size_t n, Comp* const comp, Swap* const swap, std::size_t level)
        {
            // avoid unnecessary calling performance penalty
            assert (n > 0);

            // quick sort need to be invoked with at least 2 elements
            assert (AQSORT_INSERTION_SORT_THRESHOLD > 1);

            // tail-recursion removal loop
            while (true) {
                // check quicksort worst case
                if (level == 0) {
                    heap_sort(start, n, comp, swap);
                    return;
                }
                level--;

                // choose median-of-3 pivot and put it at the end (position n - 1)
                std::size_t pivot = select_pivot_mo3(start, n, comp);

                // swap to the end
                (*swap)(pivot, start + n - 1);
                pivot = start + n - 1;

                // partitioning
                std::size_t less_than = sequential_partition(start, n, pivot, comp, swap);

                assert(less_than >= 0);
                assert(less_than <= n);

                // swap pivot to its final position
                (*swap)(start + less_than, pivot);

                // do not process equal-to-pivot elements
                std::size_t greater_than = n - less_than - 1;
                while ((greater_than > 0) &&
                        ((*comp)(start + less_than, start + n - greater_than) == false) &&
                        ((*comp)(start + n - greater_than, start + less_than) == false))
                    greater_than--;

                if (less_than > greater_than) {
                    if (greater_than > AQSORT_INSERTION_SORT_THRESHOLD)
                        quick_sort(start + n - greater_than, greater_than, comp, swap, level);
                    if (less_than <= AQSORT_INSERTION_SORT_THRESHOLD)
                        break;
                    n = less_than;
                }
                else {
                    if (less_than > AQSORT_INSERTION_SORT_THRESHOLD)
                        quick_sort(start,                 less_than,       comp, swap, level);
                    if (greater_than <= AQSORT_INSERTION_SORT_THRESHOLD)
                        break;
                    start += n - greater_than;
                    n = greater_than;
                }
            }
        }

        template<typename Comp, typename Swap>
        inline void sequential_sort(std::size_t start, std::size_t n,
                Comp* const comp, Swap* const swap, std::size_t max_level)
        {
            // nothing-to-sort check
            if (n <= 1)
                return;

            if (n > AQSORT_INSERTION_SORT_THRESHOLD) {
                quick_sort(start, n, comp, swap, max_level);
            }
            insertion_sort(start, n, comp, swap);
        }

        template<typename Comp, typename Swap>
        inline void sequential_sort(std::size_t n, Comp* const comp, Swap* const swap)
        {
            std::size_t max_level = (std::size_t)(2.0 * floor(log2(double(n))));
            sequential_sort(0, n, comp, swap, max_level);
        }

        template<typename Comp, typename Swap>
        inline void sift_down(std::size_t start, std::size_t start_, std::size_t end,
                Comp* const comp, Swap* const swap)
        {
            std::size_t root = start_;
            std::size_t child = (root << 1) + 1;
            while (child <= end) {
                std::size_t swap_ = root;

                if ((*comp)(start + swap_, start + child))
                    swap_ = child;
                if ((child + 1 <= end) && ((*comp)(start + swap_, start + child + 1)))
                    swap_ = child + 1;
                if (swap_ == root)
                    return;
                (*swap)(start + root, start + swap_);
                root = swap_;
                child = (root << 1) + 1;
            }
        }
/*
        template<typename Comp, typename Swap>
        void heapify(std::size_t start, std::size_t count, Comp* const comp, Swap* const swap)
        {
            std::size_t start_ = (count - 2) >> 1;
            if (start_ >= 0)
                while(true) {
                    sift_down(start, start_, count - 1, comp, swap);
                    if (start_ == 0)
                        break;
                    start_--;
                }
        }
*/
        template<typename Comp, typename Swap>
        void heap_sort(std::size_t start, std::size_t n, Comp* const comp, Swap* const swap)
        {
         // heapify(start, n, comp, swap);
            std::size_t start_ = (n - 2) >> 1;
            if (start_ >= 0)
                while(true) {
                    sift_down(start, start_, n - 1, comp, swap);
                    if (start_ == 0)
                        break;
                    start_--;
                }
            
            std::size_t end = n - 1;
            while (end > 0) {
                (*swap)(start + end, start);
                end--;
                sift_down(start, 0, end, comp, swap);
            }
        }
    }
}

#endif
