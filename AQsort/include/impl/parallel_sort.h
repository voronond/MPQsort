#ifndef AQSORT_IMPL_PARALLEL_SORT_H
#define AQSORT_IMPL_PARALLEL_SORT_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cassert>
#include <cmath>
#include <cstddef>

#include "parallel_partition.h"
#include "sequential_sort.h"

namespace aqsort
{
    namespace impl
    {
        template<typename Comp, typename Swap>
        void parallel_quick_sort(std::size_t total_n, std::size_t total_P, std::size_t start, std::size_t n,
                Comp* const comp, Swap* const swap, std::size_t level)
        {
            // avoid unnecessary calling performance penalty
            assert (n > 1); 

            // tail-recursion removal loop
            while (true) {
                // check quicksort worst case 
                if (level == 0) {
                    sequential_sort(start, n, comp, swap, level);
                    return;
                }
                level--;

                // actual number of threads for this chunk
                std::size_t P = total_P * n / total_n;
                if (P < 1)
                    P = 1;

                if (P < 2) {
                    // sequential sort
//#pragma omp task untied firstprivate(comp, swap)
//#pragma omp task firstprivate(comp, swap)
                    sequential_sort(start, n, comp, swap, level);
                    return;
                }

                // choose median-of-medians pivot and put it at the end (position n - 1)
                std::size_t pivot = select_pivot_mom(start, n, comp);

                // swap to the end
                (*swap)(pivot, start + n - 1);
                pivot = start + n - 1;

                // parallel partitioning using P threads
                std::size_t less_than = parallel_partition(start, n, pivot, comp, swap, P);

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
//#pragma omp task untied firstprivate(comp, swap)
#pragma omp task firstprivate(comp, swap)
                    parallel_quick_sort(total_n, total_P, start + n - greater_than, greater_than, comp, swap, level);
                    n = less_than;
                }
                else {
//#pragma omp task untied firstprivate(comp, swap)
#pragma omp task firstprivate(comp, swap)
                    parallel_quick_sort(total_n, total_P, start, less_than, comp, swap, level);
                    start += n - greater_than;
                    n = greater_than;
                }

/*
                if (less_than > 1)
#pragma omp task untied firstprivate(comp, swap)
                    parallel_quick_sort(total_n, total_P, start,                    less_than,    comp, swap);
                if (greater_than > 1)
                    parallel_quick_sort(total_n, total_P, start + n - greater_than, greater_than, comp, swap);
*/
            }
        }

        template<typename Comp, typename Swap>
        void parallel_sort(std::size_t n, Comp* const comp, Swap* const swap)
        {
            std::size_t max_level = (std::size_t)(2.0 * floor(log2(double(n))));

            int nested_saved = omp_get_nested();
            omp_set_nested(1);

            // nested is a must for parallel algorithm
            if (omp_get_nested() == 0) {
                sequential_sort(0, n, comp, swap, max_level);
                omp_set_nested(nested_saved);
                return;
            }

            std::size_t P = omp_get_num_threads();

            if (P > 1) {
                // already in parallel region
#pragma omp master
                parallel_quick_sort(n, P, 0, n, comp, swap, max_level);
#pragma omp barrier
            }
            else {
#pragma omp parallel
                {
                    std::size_t P = omp_get_num_threads();
                    // at least 2 threads is a must for parallel sorting
                    if (P < 2) 
                        sequential_sort(0, n, comp, swap, max_level);
                    else {
#pragma omp master
                        parallel_quick_sort(n, P, 0, n, comp, swap, max_level);
                    }
                } // implied barrier at the end of parallel construct
            }

            omp_set_nested(nested_saved);
        }
    }
}

#endif
