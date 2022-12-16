#ifndef AQSORT_IMPL_PARALLEL_PARTITION_H
#define AQSORT_IMPL_PARALLEL_PARTITION_H

#include <omp.h>

#include <cassert>
#include <cstddef>

#include "sequential_partition.h"

namespace aqsort
{
    namespace impl
    {
        // alias for block size
        static const std::size_t B = AQSORT_PARALLEL_PARTITION_BLOCK_SIZE;

        template<typename Comp, typename Swap>
        std::size_t parallel_partition(const std::size_t start, const std::size_t n, const std::size_t pivot,
                Comp* const comp, Swap* const swap, const std::size_t P)
        {
            // avoid unnecessary calling performance penalty
            assert (n > 0);

            // do not call for single thread only
            assert (P > 1);

            // number of blocks
            const std::size_t m = n / B;
            // must not be called for small number of elements
            assert (m >= P);

            // auxiliary data
            std::size_t  threads_temp[2 * P + 1];

            std::size_t* const threads_left  = threads_temp;
            std::size_t* const threads_start = threads_temp + P;
            std::size_t* const threads_end   = threads_temp + P + 1;

            for (std::size_t p = 0; p < P; p++) 
                threads_start[p] = start + B * (p * m / P);
            threads_start[P] = start + B * m;

            // auxiliary data
            std::size_t less_than;
            std::size_t left_thread, right_thread;
            std::size_t global_left, global_right;

// must be called by single thread only
#pragma omp parallel num_threads(P) default(none) \
    firstprivate(comp, swap, threads_left, threads_start, threads_end, pivot, P, start, n) \
    shared(less_than, left_thread, right_thread, global_left, global_right)
            {
                const std::size_t p = omp_get_thread_num();

// partitioning thread parts
                std::size_t left  = threads_start[p];
                std::size_t right = threads_start[p + 1] - 1;
                assert(right >= left);

                std::size_t seq_n = right - left + 1;

                threads_left[p] = left;
                if (seq_n > 0)
                    threads_left[p] += sequential_partition(left, seq_n, pivot, comp, swap);

// neutralize grey blocks
#pragma omp barrier

#pragma omp single
                {
                    std::size_t i = 0;
                    std::size_t j = P - 1;

                    while (i < j) {
                        // check if no grey blocks
                        std::size_t mod_i = (threads_left[i] - start) % B;
                        std::size_t mod_j = (threads_left[j] - start) % B;

                        if (mod_i == 0) {
                            i++;
                            continue;
                        }
                        if (mod_j == 0) {
                            j--;
                            continue;
                        }

                        std::size_t i_last  = threads_left[i] - mod_i + B - 1;
                        std::size_t j_first = threads_left[j] - mod_j;

                        while ((threads_left[i] <= i_last) && (threads_left[j] - 1 >= j_first)) {
                            (*swap)(threads_left[i], threads_left[j] - 1);
                            threads_left[i]++;
                            threads_left[j]--;
                        }
                    }

                    // find pivot position
                    less_than = 0;
                    for (std::size_t k = 0; k < P; k++)
                        less_than += threads_left[k] - threads_start[k];

                    // place the remaining grey block to its final position if exists and not already in place
                    std::size_t temp = (threads_left[i] - start) % B;
                    if ((temp != 0) &&
                            ((start + less_than <  threads_left[i] - temp) ||
                             (start + less_than >= threads_left[i] - temp + B))) {
                        if ((*comp)(start + less_than, pivot)) {
                            std::size_t count = B - temp;
                            while (count > 0) {
                                (*swap)(threads_left[i], start + less_than + count - 1);
                                threads_left[i]++;
                                count--;
                            }
                        }
                        else {
                            std::size_t count = temp;
                            while (count > 0) {
                                (*swap)(threads_left[i] - 1, start + less_than - count);
                                threads_left[i]--;
                                count--;
                            }
                        }
                    }

                    left_thread  = 0;
                    right_thread = P - 1;
                    global_left  = threads_left[0];
                    global_right = threads_left[P - 1];
                }
                // implicit barrier

// swap misplaced blocks in parallel
                bool done = false;
                while (true) {
                    std::size_t my_left, my_right;

                    // find next blocks to be swapped
#pragma omp critical
                    {
                        if ((global_left >= global_right) || (global_left      >= (start + less_than))
                                                          || (global_right - B <  (start + less_than))) {
                            done = true;
                        }
                        else {
                            while (global_left >= threads_end[left_thread]) {
                                if (++left_thread >= P) {
                                    done = true;
                                    break;
                                }
                                global_left = threads_left[left_thread];
                            }
                            my_left = global_left;
                            global_left += B;

                            // check underflow
                            while ((global_right > B) && (global_right - B < threads_start[right_thread])) {
                                if (right_thread == 0) {
                                    done = true;
                                    break;
                                }
                                global_right = threads_left[--right_thread];
                            }
                            if (global_right <= B) {
                                done = true;
                            }
                            else {
                                my_right = global_right - B;
                                global_right -= B;
                            }
                        }
                    }

                    if (done)
                        break;

                    // swap blocks
                    if (my_left < my_right) 
                        for (std::size_t k = 0; k < B; k++)
                            (*swap)(my_left + k, my_right + k);
                }
#pragma omp barrier

// process the remainder
#pragma omp single
                {
                    for (std::size_t k = threads_end[P - 1]; k < start + n; k++) {
                        if ((*comp)(k, pivot)) {
                            (*swap)(k, start + less_than);
                            less_than++;
                        }
                    }
                }
            } // end of parallel region

            return less_than;
        }
    }
}

#endif
