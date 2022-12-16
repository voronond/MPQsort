#ifndef AQSORT_IMPL_SEQUENTIAL_PARTITION_H
#define AQSORT_IMPL_SEQUENTIAL_PARTITION_H

#include <cassert>
#include <cstddef>

namespace aqsort
{
    namespace impl
    {
        template<typename Comp, typename Swap>
        inline std::size_t sequential_partition(
                std::size_t start, std::size_t n, std::size_t pivot, Comp* const comp, Swap* const swap)
        {
            // avoid unnecessary calling performance penalty
            assert (n > 0);

            std::size_t left  = start;
            std::size_t right = start + n - 1;

            while (true) {
                while ((left < right) && ((*comp)(left, pivot)))
                    left++;
                while ((left < right) && (!(*comp)(right, pivot)))
                    right--;
                if (left < right) {
                    (*swap)(left, right);
                    left++;
                 // if (left < right) // cannot lead to underflow since left >= 1 and right >= left here
                    right--;
                }
                else 
                    break;
            }
            if ((*comp)(left, pivot))
                left++;

            // return number of elements less than pivot
            return left - start;
/*
            while (true) {
                while ((*comp)(left, pivot)) {
                    left++;
                    if (left > right)
                        return left - start;
                }

                while (!(*comp)(right, pivot)) {
                    if (right == left)
                        return left - start;
                    right--;
                }

                if (right == left)
                    return left - start;
                
                (*swap)(left, right);
                left++;
             // right--;
            }
*/
        }
    }
}

#endif
