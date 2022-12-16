#ifndef AQSORT_IMPL_SELECT_PIVOT_H
#define AQSORT_IMPL_SELECT_PIVOT_H

#include <cassert>
#include <cstddef>

namespace aqsort
{
    namespace impl
    {
        template<typename Comp>
        inline std::size_t median(std::size_t i1, std::size_t i2, std::size_t i3, Comp* const comp)
        {
            if ((*comp)(i1, i3)) {
                if ((*comp)(i2, i1)) 
                    return i1;
                else if ((*comp)(i2, i3))
                    return i2;
                return i3;
            }

            if ((*comp)(i1, i2))
                return i1;
            else if ((*comp)(i3, i2))
                return i2;
            return i3;
        }

        // median-of-3 pivot
        template<typename Comp>
        inline std::size_t select_pivot_mo3(std::size_t start, std::size_t n, Comp* const comp)
        {
            // avoid unnecessary calling performance penalty
            assert (n > 0);

            std::size_t left   = start;
            std::size_t middle = start + (n >> 1);
            std::size_t right  = start + n - 1;

            return median(left, middle, right, comp);
        }

        // median-of-medians pivot
        template<typename Comp>
        inline std::size_t select_pivot_mom(std::size_t start, std::size_t n, Comp* const comp)
        {
            // avoid unnecessary calling performance penalty
            assert (n > 0);

            std::size_t third1 =     n / 3;
            std::size_t third2 = 2 * n / 3;

            std::size_t i1 = start;
            std::size_t i2 = start + (third1 >> 1);
            std::size_t i3 = start + third1;
            std::size_t median1 = median(i1, i2, i3, comp);

            i1 = start + third1 + 1;
            i2 = start + third1 + ((third2 - third1) >> 1);
            i3 = start + third2;
            std::size_t median2 = median(i1, i2, i3, comp);

            i1 = start + third2 + 1;
            i2 = start + third2 + ((n - third2) >> 1);
            i3 = start + n - 1;
            std::size_t median3 = median(i1, i2, i3, comp);
            
            return median(median1, median2, median3, comp);
        }
    }
}

#endif
