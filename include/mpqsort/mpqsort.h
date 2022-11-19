#pragma once

#include <omp.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <execution>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <type_traits>

// TODO: remove after all methods implemented
#define UNUSED(x) (void)(x)
// TODO: Remove when done
#define DEBUG

#ifdef DEBUG
#    define PRINT_ITERS(base, lp, rp, msg) mpqsort::helpers::print(base, lp, rp, msg)
#    define PRINT_VECTOR(vec, msg)                   \
        std::cout << msg << ":\n";                   \
        for (auto& el : vec) std::cout << el << " "; \
        std::cout << std::endl
#    define PRINT_TABLE(vec, msg)                                                  \
        std::cout << msg << ":\n";                                                 \
        for (size_t i = 0; i < vec[0].size(); ++i) {                               \
            for (size_t j = 0; j < vec.size(); ++j) std::cout << vec[j][i] << " "; \
            std::cout << std::endl;                                                \
        }
#else
#    define PRINT_ITERS(base, lp, rp, msg)
#    define PRINT_VECTOR(vec, msg)
#    define PRINT_TABLE(vec, msg)
#endif

#define INDEX(a, b, size) (a) * (size) + (b)

/**
 * @brief Contains execution policy definitions for sort algorithms, helper functions and type
 * traits.
 */
namespace mpqsort::execution {
    // Execution policy class declaration

    class sequenced_policy_three_way {};
    class sequenced_policy_four_way {};
    class sequenced_policy_max_way {};
    class parallel_policy_max_way {};

    // TODO: maybe add unseq if possible (standard and self defined)
    /**
     * @brief Sequenced execution policy from std
     */
    using std::execution::seq;

    /**
     * @brief Parallel execution policy from std
     */
    using std::execution::par;

    // Execution policy objects

    /**
     * @brief Sequenced execution policy using two pivots
     */
    inline constexpr sequenced_policy_three_way seq_three_way{};

    /**
     * @brief Sequenced execution policy using three pivots
     */
    inline constexpr sequenced_policy_four_way seq_four_way{};

    /**
     * @brief Sequenced execution policy using multiple pivots. The number of pivots will be
     * determined by the algorithm
     */
    inline constexpr sequenced_policy_max_way seq_max_way{};

    // Define own type trait to determine if we got policy type

    template <typename T> struct _is_execution_policy_helper : std::false_type {};

    /**
     * @brief Checks whether T is a execution policy type defined in MPQsort library. Automatically
     * applies decay.
     * @tparam T - a type to check
     */
    template <typename T> using is_execution_policy = _is_execution_policy_helper<std::decay_t<T>>;

    template <typename T> struct _is_parallel_execution_policy_helper : std::false_type {};

    /**
     * @brief Checks whether T is a parallel execution policy type defined in MPQsort library.
     * @tparam T - a type to check
     */
    template <typename T> using is_parallel_execution_policy
        = _is_parallel_execution_policy_helper<std::decay_t<T>>;

    /**
     * @brief Helper variable template to check the execution policy. Automatically applies decay.
     * type trait.
     * @tparam T - a type to check
     */
    template <typename T> inline constexpr bool is_execution_policy_v
        = _is_execution_policy_helper<std::decay_t<T>>::value;

    /**
     * @brief Helper variable template to check the execution policy. Automatically applies decay
     * type trait.
     * @tparam T - a type to check
     */
    template <typename T> inline constexpr bool is_parallel_execution_policy_v
        = _is_parallel_execution_policy_helper<std::decay_t<T>>::value;

    // Allowed execution policies
    template <> struct _is_execution_policy_helper<sequenced_policy_three_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<sequenced_policy_four_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<sequenced_policy_max_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<parallel_policy_max_way> : std::true_type {};

    template <> struct _is_parallel_execution_policy_helper<parallel_policy_max_way>
        : std::true_type {};

    // Allowed STD execution policies
    template <> struct _is_execution_policy_helper<std::execution::sequenced_policy>
        : std::true_type {};
    template <> struct _is_execution_policy_helper<std::execution::parallel_policy>
        : std::true_type {};

    template <> struct _is_parallel_execution_policy_helper<std::execution::parallel_policy>
        : std::true_type {};
}  // namespace mpqsort::execution

/**
 * @brief This namespace contains parameters for algorithm tuning
 */
namespace mpqsort::parameters {
    /**
     * @brief Size of block for parallel multiway partitioning
     * Threads process an array in blocks to prevent false sharing and reduce number of
     * atomic/critical sections.
     */
    static long BLOCK_SIZE = 64;
    static long SEQ_THRESHOLD = 1 << 17;
    static long NO_RECURSION_THRESHOLD = 128;  // based on benchmarks
    /**
     * @brief Maximum supported number of pivots by MPQsort
     * This is the max number of pivots that can MPQsort use during computation. This can't be
     * changed by the user.
     */
    const static long MAX_NUMBER_OF_PIVOTS = 3;  // can't be changed during runtime

    /**
     * @brief Number of pivots used in multiway parallel partitioning
     * Number for pivots used to divide array in PAR_PARTITION_NUM_PIVOTS + 1 segments. This is done only once and parallel multiway qsort is called afterwards on each such segment.
     */
    const static long PAR_PARTITION_NUM_PIVOTS = 255;

    /**
     * @brief Number of elements to chose a pivot from
     * Determined hos many elements from an input array we want to consider to chose one pivot
     */
    const static long ONE_PIVOT_SAMPLE_SIZE = 20;
}  // namespace mpqsort::parameters

namespace mpqsort::helpers {
    // Extension of is_same applying decay_t automatically
    template <typename T, typename U> inline constexpr bool _is_same_decay_v
        = std::is_same<std::decay_t<T>, std::decay_t<U>>::value;

    template <typename RandomBaseIt>
    void print(RandomBaseIt base, size_t lp, size_t rp, const std::string name = "") {
        if (!name.empty()) std::cout << name << ": ";

        std::cout << "{ ";

        if (lp > rp) {
            std::cout << "}" << std::endl;
            return;
        }

        std::cout << base[lp];
        ++lp;
        while (lp <= rp) {
            std::cout << ", " << base[lp];
            ++lp;
        }

        std::cout << " }" << std::endl;
    }

    template <typename RandomBaseIt, typename Index>
    inline void _cyclic_shift_left(RandomBaseIt base, Index first, Index second, Index third) {
        auto tmp = base[first];
        base[first] = base[second];
        base[second] = base[third];
        base[third] = tmp;
    }

    template <typename RandomBaseIt, typename Index>
    inline void _cyclic_shift_right(RandomBaseIt base, Index first, Index second, Index third) {
        auto tmp = base[third];
        base[third] = base[second];
        base[second] = base[first];
        base[first] = tmp;
    }

    template <typename RandomBaseIt, typename Index>
    inline void _cyclic_shift_right(RandomBaseIt base, Index first, Index second, Index third,
                                    Index forth) {
        auto tmp = base[forth];
        base[forth] = base[third];
        base[third] = base[second];
        base[second] = base[first];
        base[first] = tmp;
    }

    template <typename RandomBaseIt, typename Compare>
    inline void _push_heap(RandomBaseIt base, long size, Compare& comp) {
        for (auto i = size - 1;;) {
            auto parent = (i - 1) / 2;
            auto ai = base[i], ap = base[parent];
            if (comp(ap, ai)) break;
            base[parent] = ai;
            base[i] = ap;
            if (parent == 0) break;
            i = parent;
        }
    }

    template <typename RandomBaseIt, typename Compare>
    inline void _make_heap(RandomBaseIt base, long size, Compare& comp) {
        auto first_parent = (size - 3) / 2;
        auto first_right_kid = (first_parent + 1) * 2;

        for (;; --first_parent, first_right_kid -= 2) {
            auto drop_down = base[first_parent];
            auto parent = first_parent;
            auto right_kid = first_right_kid;
            for (;;) {
                auto tmp = right_kid - comp(base[right_kid - 1], base[right_kid]);
                auto crt = base[tmp];
                if (comp(drop_down, crt)) break;
                base[parent] = crt;
                parent = tmp;
                right_kid = (tmp + 1) * 2;
                if (right_kid >= size) goto write;
            }
            if (parent != first_parent)
            write:  // Ugly but efficient
                base[parent] = drop_down;
            if (first_parent == 0) break;
        }

        if (size & 1) return;

        // Fix only child
        _push_heap(base, size, comp);
    }

    template <typename RandomBaseIt, typename Compare>
    inline void _unguarded_insertion_sort(RandomBaseIt base, long lp, long rp, Compare& comp) {
        // Index of first sorted element
        for (auto i = lp + 1;;) {
            auto el = base[i];
            auto j = i - 1;

            while (j >= lp && comp(el, base[j])) {
                base[j + 1] = base[j];
                --j;
            }

            base[j + 1] = el;

            // Save one increment
            if (i >= rp) return;
            ++i;
        }
    }

    template <typename RandomBaseIt, typename Compare>
    inline void _heap_insertion_sort(RandomBaseIt base, long lp, long rp, Compare& comp) {
        using std::swap;

        long size = rp - lp;
        // size == 1 if sort 2 elements, otherwise swap element with itself
        if (size < 2) {
            if (!comp(base[lp], base[lp + (size == 1)])) swap(base[lp], base[lp + (size == 1)]);
            return;
        }

        _make_heap(base + lp, size + 1, comp);
        _unguarded_insertion_sort(base, lp + 1, rp, comp);
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivots_indexes_two(RandomBaseIt base, long lp, long rp, Comparator& comp) {
        auto size = rp - lp + 1;
        const auto sample_size = parameters::ONE_PIVOT_SAMPLE_SIZE * 2;
        std::vector<std::pair<typename std::iterator_traits<RandomBaseIt>::value_type, long>>
            samples;
        samples.reserve(sample_size);

        // If not enough elements for sampling
        if (size < sample_size) {
            return std::tuple{size * 1 / 3 + lp, size * 2 / 3 + lp};
        } else {
            // Get sample elements
            for (long i = 0; i < sample_size; ++i) {
                auto index = size * i / sample_size + lp;
                samples[i] = std::make_pair(base[index], index);
            }

            // Sort samples based on provided comp
            std::sort(samples.begin(), samples.end(),
                      [&](auto& a, auto& b) { return comp(a.first, b.first); });

            return std::tuple{samples[sample_size * 1 / 3].second,
                              samples[sample_size * 2 / 3].second};
        }
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivot_indexes_three(RandomBaseIt base, long lp, long rp, Comparator& comp) {
        auto size = rp - lp + 1;
        const auto sample_size = parameters::ONE_PIVOT_SAMPLE_SIZE * 3;
        std::vector<std::pair<typename std::iterator_traits<RandomBaseIt>::value_type, long>>
            samples;
        samples.reserve(sample_size);

        // If not enough elements for sampling
        if (size < sample_size) {
            return std::tuple{size * 1 / 4 + lp, size * 2 / 4 + lp, size * 3 / 4 + lp};
        } else {
            // Get sample elements
            for (long i = 0; i < sample_size; ++i) {
                auto index = size * i / sample_size + lp;
                samples[i] = std::make_pair(base[index], index);
            }

            // Sort samples based on provided comp
            std::sort(samples.begin(), samples.end(),
                      [&](auto& a, auto& b) { return comp(a.first, b.first); });

            return std::tuple{samples[sample_size * 1 / 4].second,
                              samples[sample_size * 2 / 4].second,
                              samples[sample_size * 3 / 4].second};
        }
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivots(RandomBaseIt base, long lp, long rp, long num_pivots,
                            Comparator& comp) {
        auto size = rp - lp + 1;
        const auto sample_size = parameters::ONE_PIVOT_SAMPLE_SIZE * num_pivots;

        assert("num_pivots must be 2^k - 1" && (num_pivots + 1) > 0
               && ((num_pivots + 1) & (num_pivots)) == 0);
        assert("Array smaller than num_pivots" && size >= num_pivots);

        using ValueType = typename std::iterator_traits<RandomBaseIt>::value_type;
        std::vector<ValueType> samples(sample_size);
        std::vector<ValueType> pivots(num_pivots);

        // If not enough elements for sampling
        if (size < sample_size) {
            for (long i = 1; i < num_pivots + 1; ++i) {
                pivots[i - 1] = base[size * i / (num_pivots + 1) + lp];
            }

            std::sort(pivots.begin(), pivots.end(), comp);

            return pivots;
        } else {
            // Get sample elements
            for (long i = 0; i < sample_size; ++i) {
                auto index = size * i / sample_size + lp;
                samples[i] = base[index];
            }

            std::sort(samples.begin(), samples.end(), comp);

            for (long i = 1; i < num_pivots + 1; ++i) {
                pivots[i - 1] = samples[sample_size * i / (num_pivots + 1)];
            }

            return pivots;
        }
    }

    // Number of pivots needs to be 2^k - 1
    template <typename RandomBaseIt, typename Comparator, typename Element>
    inline auto _find_element_segment_id(int num_of_comparisons, RandomBaseIt base, int size,
                                         Element& el, Comparator& comp) {
        int lp = 0, rp = size - 1;
        int idx;

        for (int i = 0; i < num_of_comparisons; ++i) {
            idx = (lp + rp) / 2;

            if (comp(el, base[idx])) {
                rp = idx - 1;
            } else {
                lp = idx + 1;
            }
        }

        idx = (lp + rp) / 2;

        return comp(el, base[idx]) ? idx : idx + 1;
    }

}  // namespace mpqsort::helpers

/**
 * @brief Implementation of multiway quicksort and multiway parallel quicksort
 */
namespace mpqsort::impl {
    // ----- Main implementation START -----
    // OpenMP mergeable possible if shared variables (code gets executed as separate task or not)
    // SEQ

    template <typename RandomBaseIt, typename Compare>
    inline auto _seq_multiway_partition_two_pivots(RandomBaseIt base, long lp, long rp,
                                                   Compare& comp) {
        // Use optimal swap method
        using std::swap;

        // Get pivots
        auto [idx1, idx2] = helpers::_get_pivots_indexes_two(base, lp + 1, rp - 1, comp);

        if (comp(base[idx1], base[idx2])) {
            swap(base[lp], base[idx1]);
            swap(base[rp], base[idx2]);
        } else {
            swap(base[lp], base[idx2]);
            swap(base[rp], base[idx1]);
        }

        const auto p1 = base[lp], p2 = base[rp];

        // Set boundaries
        auto k2 = lp + 1, k = k2, g = rp - 1;

        while (k <= g) {
            if (comp(base[k], p1)) {
                swap(base[k2], base[k]);
                ++k2;
            } else {
                if (!comp(base[k], p2)) {
                    while (k < g && comp(p2, base[g])) --g;

                    if (!comp(base[g], p1)) {
                        swap(base[k], base[g]);
                    } else {
                        helpers::_cyclic_shift_left(base, k, k2, g);
                        ++k2;
                    }
                    --g;
                }
            }
            ++k;
        }

        // Move pivots at right place
        swap(base[rp], base[g + 1]);
        swap(base[lp], base[k2 - 1]);

        return std::tuple{k2 - 1, g + 1};
    }

    template <typename NumPivot, typename RandomBaseIt, typename Compare>
    void _seq_multiway_qsort_inner_ybb(NumPivot pivot_num, RandomBaseIt base, long lp, long rp,
                                       Compare& comp, long depth) {
        while (rp - lp > parameters::NO_RECURSION_THRESHOLD && depth > 0) {
            auto [index_p1, index_p2] = _seq_multiway_partition_two_pivots(base, lp, rp, comp);
            _seq_multiway_qsort_inner_ybb(pivot_num, base, lp, index_p1 - 1, comp, depth - 1);
            _seq_multiway_qsort_inner_ybb(pivot_num, base, index_p1 + 1, index_p2 - 1, comp,
                                          depth - 1);
            lp = index_p2 + 1;
        }

        helpers::_heap_insertion_sort(base, lp, rp, comp);
    }

    template <typename RandomBaseIt, typename Compare>
    inline auto _seq_multiway_partition_three_pivots(RandomBaseIt base, const long lp,
                                                     const long rp, Compare& comp) {
        // Use optimal swap method
        using std::swap;

        // Get pivots
        auto [idx1, idx2, idx3] = helpers::_get_pivot_indexes_three(base, lp + 2, rp - 1, comp);

        // Sort pivots
        if (comp(base[idx2], base[idx1])) swap(base[idx1], base[idx2]);
        if (comp(base[idx3], base[idx1])) swap(base[idx1], base[idx3]);
        if (comp(base[idx3], base[idx2])) swap(base[idx2], base[idx3]);

        // Move pivots to array edges
        swap(base[idx1], base[lp]);
        swap(base[idx2], base[lp + 1]);
        swap(base[idx3], base[rp]);

        // Indexes
        long k2, k, g, g2;
        k2 = k = lp + 2;
        g = g2 = rp - 1;

        const auto p1 = base[lp], p2 = base[lp + 1], p3 = base[rp];

        while (k <= g) {
            while (k <= g && comp(base[k], p2)) {
                if (comp(base[k], p1)) {
                    swap(base[k2], base[k]);
                    ++k2;
                }
                ++k;
            }
            while (k <= g && comp(p2, base[g])) {
                if (comp(p3, base[g])) {
                    swap(base[g], base[g2]);
                    --g2;
                }
                --g;
            }
            if (k <= g) {
                if (comp(p3, base[k])) {
                    if (comp(base[g], p1)) {
                        helpers::_cyclic_shift_right(base, k2, k, g2, g);
                        ++k2;
                    } else {
                        helpers::_cyclic_shift_left(base, k, g, g2);
                    }
                    ++k;
                    --g;
                    --g2;
                } else {
                    if (comp(base[g], p1)) {
                        helpers::_cyclic_shift_right(base, k2, k, g);
                        ++k2;
                    } else {
                        swap(base[k], base[g]);
                    }
                    ++k;
                    --g;
                }
            }
        }

        --k2;
        --k;
        ++g;
        ++g2;

        swap(base[lp + 1], base[k2]);
        swap(base[k2], base[k]);
        --k2;
        swap(base[lp], base[k2]);
        swap(base[rp], base[g2]);

        return std::tuple{k2, k, g2};
    }

    template <typename NumPivot, typename RandomBaseIt, typename Compare>
    void _seq_multiway_qsort_inner_waterloo(NumPivot pivot_num, RandomBaseIt base, long lp, long rp,
                                            Compare& comp, long depth) {
        while (rp - lp > parameters::NO_RECURSION_THRESHOLD && depth > 0) {
            auto [index_p1, index_p2, index_p3]
                = _seq_multiway_partition_three_pivots(base, lp, rp, comp);
            _seq_multiway_qsort_inner_waterloo(pivot_num, base, lp, index_p1 - 1, comp, depth - 1);
            _seq_multiway_qsort_inner_waterloo(pivot_num, base, index_p1 + 1, index_p2 - 1, comp,
                                               depth - 1);
            _seq_multiway_qsort_inner_waterloo(pivot_num, base, index_p2 + 1, index_p3 - 1, comp,
                                               depth - 1);
            lp = index_p3 + 1;
        }

        helpers::_heap_insertion_sort(base, lp, rp, comp);
    }

    template <typename NumPivot, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _seq_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        if (last - first <= 1) return;

        // TODO: Remove one pivot algorithm
        if (pivot_num == 2 || pivot_num == 1) {
            _seq_multiway_qsort_inner_ybb(pivot_num, first, 0, last - first - 1, comp,
                                          1.5 * std::log(last - first) / std::log(3));
        } else if (pivot_num == 3) {
            _seq_multiway_qsort_inner_waterloo(pivot_num, first, 0, last - first - 1, comp,
                                               1.5 * std::log(last - first) / std::log(4));
        } else {
            throw std::invalid_argument("Unknown number of pivots, this should never happen");
        }
    }

    // PAR
    template <typename RandomBaseIt, typename Cores, typename Compare>
    auto _par_multiway_partition(long num_pivots, Cores cores, RandomBaseIt base, long lp, long rp,
                                 Compare& comp) {
        using std::swap;
        using ValueType = typename std::iterator_traits<RandomBaseIt>::value_type;

        const int num_segments = num_pivots + 1;
        int num_segments_left = num_segments;
        // Precompute number of comparisons
        int num_element_comparisons = std::log2(num_pivots);

        // Get pivots
        auto pivots = helpers::_get_pivots(base, lp, rp, num_pivots, comp);

        // Indexes in blocks
        std::vector<long> segment_idx(num_pivots + 1);
        auto idx_ptr = segment_idx.data();

// Find boundaries of segments
#pragma omp parallel for reduction(+ : idx_ptr[:num_segments])
        for (long i = lp; i <= rp; ++i) {
            auto segment_id = helpers::_find_element_segment_id(num_element_comparisons, pivots,
                                                                pivots.size(), base[i], comp);
            ++idx_ptr[segment_id];
        }

        // Compute beginning of each segment
        std::exclusive_scan(segment_idx.begin(), segment_idx.end(), segment_idx.begin(), 0);

        // Create boundaries for each segment
        std::vector<long> segment_boundary(segment_idx.begin() + 1, segment_idx.end());
        segment_boundary.emplace_back(rp + 1);

        // If index already >= boundary => segment clean and no elements belong in it
        for (size_t i = 0; i < segment_idx.size(); ++i) {
            if (segment_idx[i] == segment_boundary[i]) --num_segments_left;
        }

        // Returns if element belongs to given segment
        auto element_in_segment = [&](auto& el, size_t segment) {
            // Compound logical comparison to prevent redundant branching
            return (segment == 0 || !comp(el, pivots[segment - 1]))
                   && (segment == segment_idx.size() || comp(el, pivots[segment]));
        };

        auto find_unprocessed_segment = [&](auto& segment) -> bool {
            for (size_t i = 0; i < segment_idx.size(); ++i) {
                ++segment %= segment_idx.size();
                if (segment_idx[segment] < segment_boundary[segment]) return true;
            }
            // All segments processed
            return false;
        };

        // This is the max number of elements that can be inserted in one table "segment"
        int table_height = (cores - 1) * parameters::BLOCK_SIZE;

        // Table of indexes where is the emtpy space to insert an element
        std::vector<long> elements_insertions(num_segments * table_height);
        std::vector<long> elements_insertions_index(num_segments, 0);

        // Table of elements belonging to segment
        std::vector<ValueType> elements_table(num_segments * table_height);
        std::vector<long> elements_table_index(num_segments, 0);

        // Find first unclean segment
        int current_segment = 0;
        find_unprocessed_segment(current_segment);

#pragma omp parallel shared(elements_insertions, elements_insertions_index, elements_table,    \
                            elements_table_index, segment_idx, segment_boundary, base, pivots) \
    firstprivate(num_segments, num_segments_left, current_segment) num_threads(cores)
        {
            // Private start and end of block for given segment
            std::vector<long> block_start(pivots.size() + 1, 0);
            std::vector<long> block_end(pivots.size() + 1, 0);
            ValueType tmp_el;
            int tmp_el_index = -1;
            int tmp_el_segment = -1;
            int dept_segment = -1;

            auto find_unprocessed_block = [&](auto& segment) -> bool {
                for (size_t i = 0; i < block_start.size() - 1; ++i) {
                    ++segment %= block_start.size();
                    if (block_start[segment] < block_end[segment]) return true;
                }
                // All blocks processed
                return false;
            };

            while (num_segments_left > 0) {
                // Block already marked as a "clean", so only move element in a global table
                // This can happen only if this thread as an element belonging in this segment but
                // is owned by another thread
                if (block_start[current_segment] == -1) {
                    // Do not know the index of this element in other block, it will be provided by
                    // another thread after blocks clean
                    if (tmp_el_index != -1) {  // Only if tmp_el set
                        int index;
#pragma omp atomic capture
                        index = elements_table_index[current_segment]++;
                        elements_table[INDEX(current_segment, index, table_height)] = tmp_el;
                        // tmp_el was from dept segment, segment no longer in dept
                        if (dept_segment != -1) {
                            // Insert element index in global table so that other threads can copy
                            // element to that place
#pragma omp atomic capture
                            index = elements_insertions_index[tmp_el_segment]++;
                            elements_insertions[INDEX(tmp_el_segment, index, table_height)]
                                = tmp_el_index;
                            ++block_start[dept_segment];
                            dept_segment = -1;
                        }

                        tmp_el_index = -1;
                    }

                    // Find next block/segment
                    if (!find_unprocessed_block(current_segment)) {
                        if (!find_unprocessed_segment(current_segment)) {
                            num_segments_left = 0;
                            continue;
                        }
                    }
                } else if (block_start[current_segment] >= block_end[current_segment]) {
// Get private block from global segment
#pragma omp atomic capture
                    {
                        block_start[current_segment] = segment_idx[current_segment];
                        segment_idx[current_segment] += parameters::BLOCK_SIZE;
                    }
                    // Block end can't be greater than segment_boundary
                    block_end[current_segment]
                        = std::min(block_start[current_segment] + parameters::BLOCK_SIZE,
                                   segment_boundary[current_segment]);

                    // All blocks taken, segment "clean" for this thread
                    if (block_start[current_segment] >= block_end[current_segment]) {
                        --num_segments_left;
                        // Set to know that this segment was already cleaned once so only insert
                        // element in a table
                        block_start[current_segment] = -1;
                        block_end[current_segment] = -1;
                        continue;
                    }
                }

                // Segment was in dept, save tmp_el
                if (dept_segment == current_segment) {
                    base[block_start[current_segment]] = tmp_el;
                    ++block_start[current_segment];
                    dept_segment = -1;
                    tmp_el_index = -1;
                }

                // While elements in current segment belong to it
                while (block_start[current_segment] < block_end[current_segment]
                       && element_in_segment(base[block_start[current_segment]], current_segment)) {
                    ++block_start[current_segment];
                }

                // All elements of this segment processed, continue with next one
                if (block_start[current_segment] >= block_end[current_segment]) {
                    continue;
                }

                // If started partitioning, no element to swap with
                if (dept_segment == -1) {
                    dept_segment = current_segment;
                    tmp_el_index = block_start[current_segment];
                    tmp_el_segment = current_segment;
                    tmp_el = base[block_start[current_segment]];
                } else {
                    // Found segment for tmp_el, swap them
                    swap(tmp_el, base[block_start[current_segment]]);
                    ++block_start[current_segment];
                }

                current_segment = helpers::_find_element_segment_id(
                    num_element_comparisons, pivots.begin(), pivots.size(), tmp_el, comp);
            }

            // All segments clean from this thread perspective, but some blocks have still
            // unprocessed elements. Insert those indexes in a global table so that other threads
            // know where they can place their elements from a global table
            for (size_t i = 0; i < block_start.size(); ++i) {
                if (block_start[i] >= block_end[i]) continue;

                int num_of_indexes = block_end[i] - block_start[i];
                int start;

#pragma omp atomic capture
                {
                    start = elements_insertions_index[i];
                    elements_insertions_index[i] += num_of_indexes;
                }

                auto s = block_start[i], e = block_end[i];
                while (s < e) {
                    elements_insertions[INDEX(i, start, table_height)] = s;

                    ++start;
                    ++s;
                }
            }
        }

// Insert elements from tables in an array
#pragma omp parallel for
        for (size_t i = 0; i < elements_insertions_index.size(); ++i) {
            for (int j = 0; j < elements_insertions_index[i]; ++j) {
                base[elements_insertions[INDEX(i, j, table_height)]]
                    = elements_table[INDEX(i, j, table_height)];
            }
        }

        return segment_boundary;
    }

    template <typename NumPivot, typename Cores, typename RandomBaseIt, typename Compare>
    void _par_multiway_qsort_inner(NumPivot pivot_num, Cores cores, RandomBaseIt base, size_t lp,
                                   size_t rp, Compare& comp) {
        // Parallel multiway partition using 2^k - 1 pivots
        auto boundaries = _par_multiway_partition(pivot_num, cores, base, lp, rp, comp);

        // Create ranges {start, end} for segments
        std::vector<std::pair<long, long>> segment_ranges(boundaries.size());
        for (size_t i = 1; i < boundaries.size(); ++i)
            segment_ranges.emplace_back(std::make_pair(boundaries[i - 1], boundaries[i]));

        // Sort segments based on their length from the longest to shortest
        std::sort(segment_ranges.begin(), segment_ranges.end(), [](const auto& p1, const auto& p2) {
            return (p1.second - p1.first) > (p2.second - p2.first);
        });

        // Parallel multiway sort of each segment
// TODO: Decide if use 2 or 3 pivots sort
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < segment_ranges.size(); i++) {
            _seq_multiway_qsort_inner_waterloo(3, base, segment_ranges[i].first,
                                               segment_ranges[i].second - 1, comp,
                                               1.5 * std::log(rp - lp) / std::log(4));
        }
    }

    // ----- Main implementation END -----

    template <typename NumPivot, typename Cores, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(NumPivot pivot_num, Cores cores, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        if (last - first <= 1) return;

        // Allow MAX nested parallelism
        omp_set_max_active_levels(std::numeric_limits<int>::max());

        // Should not normally happen, but we set SEQ_THRESHOLD to 0 during testing
        if (last - first > std::max(parameters::SEQ_THRESHOLD, parameters::PAR_PARTITION_NUM_PIVOTS)) {
            _par_multiway_qsort_inner(pivot_num, cores, first, 0, last - first - 1, comp);
        } else {
            _seq_multiway_qsort(3, first, last, comp);
        }
    }

    // Wrapper for different arguments
    template <typename NumPivot, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        _par_multiway_qsort(pivot_num, omp_get_max_threads(), first, last, comp);
    }

    // Call sort based on policy type
    template <typename ExecutionPolicy, typename... T>
    constexpr void _call_sort(ExecutionPolicy&& policy, T... args) {
        static_assert(
            execution::is_execution_policy_v<decltype(policy)>,
            "Provided ExecutionPolicy is not valid. Use predefined policies from namespace "
            "mpqsort::execution.");

        if constexpr (helpers::_is_same_decay_v<ExecutionPolicy, decltype(execution::seq)>) {
            // Call with one pivot
            _seq_multiway_qsort(1, std::forward<T>(args)...);
        } else if constexpr (helpers::_is_same_decay_v<ExecutionPolicy,
                                                       decltype(execution::seq_three_way)>) {
            // Call with two pivots
            _seq_multiway_qsort(2, std::forward<T>(args)...);
        } else if constexpr (helpers::_is_same_decay_v<ExecutionPolicy,
                                                       decltype(execution::seq_four_way)>) {
            // Call with three pivots
            _seq_multiway_qsort(3, std::forward<T>(args)...);
        } else if constexpr (helpers::_is_same_decay_v<ExecutionPolicy,
                                                       decltype(execution::seq_max_way)>) {
            // Let algorithm decide how many pivots to use
            _seq_multiway_qsort(parameters::MAX_NUMBER_OF_PIVOTS, std::forward<T>(args)...);
        } else if constexpr (helpers::_is_same_decay_v<ExecutionPolicy, decltype(execution::par)>) {
            // Call with one pivot
            _par_multiway_qsort(parameters::PAR_PARTITION_NUM_PIVOTS, std::forward<T>(args)...);
        } else {
            throw std::invalid_argument(
                "Unknown policy. This should never happen as we test policy type at a beginning of "
                "a function.");
        }
    }
}  // namespace mpqsort::impl

/**
 * @brief Contains public functions of the library
 */
namespace mpqsort {
    /**
     * @brief Sorts data with their default comparator
     * Runs alway sequentially
     * @tparam RandomIt
     * @param first First element of a container
     * @param last Last element of a container
     */
    template <typename RandomIt> void sort(RandomIt first, RandomIt last) {
        impl::_seq_multiway_qsort(2, first, last);
    }

    /**
     * @brief  Sorts data with their default comparator
     * Can run in parallel or sequentially. If run in parallel, uses all available cores.
     * @tparam ExecutionPolicy
     * @tparam RandomIt
     * @param policy Determines if sort runs in parallel or sequentially
     * @param first First element of a container
     * @param last Last element of a container
     */
    template <typename ExecutionPolicy, typename RandomIt>
    void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last) {
        impl::_call_sort(std::forward<ExecutionPolicy>(policy), first, last);
    }

    /**
     * @brief  Sorts data with their default comparator
     * Can run in parallel or sequentially. If run in parallel, uses all available cores.
     * @tparam ExecutionPolicy
     * @tparam Cores
     * @tparam RandomIt
     * @param policy Determines if sort runs in parallel or sequentially
     * @param cores Number of used cores for parallel version. Sequential version ignores this
     * argument.
     * @param first First element of a container
     * @param last Last element of a container
     */
    template <typename ExecutionPolicy, typename Cores, typename RandomIt>
    void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last) {
        static_assert(std::is_integral_v<Cores>, "Integral required for cores");
        impl::_call_sort(std::forward<ExecutionPolicy>(policy), cores, first, last);
    }

    /**
     * @brief Sorts data using provided comparator
     * Sorts data sequentially using provided comparator
     * @tparam RandomIt
     * @tparam Compare
     * @param first First element of a container
     * @param last Last element of a container
     * @param comp Provided comparator to compare container elements
     */
    template <typename RandomIt, typename Compare>
    void sort(RandomIt first, RandomIt last, Compare comp) {
        impl::_seq_multiway_qsort(2, first, last, comp);
    }

    /**
     * @brief Sorts data using provided comparator
     * Can run in parallel or sequentially
     * @tparam ExecutionPolicy
     * @tparam RandomIt
     * @tparam Compare
     * @param policy Determines if sort runs in parallel or sequentially
     * @param first First element of a container
     * @param last Last element of a container
     * @param comp Provided comparator to compare container elements
     */
    template <typename ExecutionPolicy, typename RandomIt, typename Compare>
    void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last, Compare comp) {
        impl::_call_sort(std::forward<ExecutionPolicy>(policy), first, last, comp);
    }

    /**
     * @brief Sorts data using provided comparator
     * Can run in parallel or sequentially
     * @tparam ExecutionPolicy
     * @tparam Cores
     * @tparam RandomIt
     * @tparam Compare
     * @param policy Determines if sort runs in parallel or sequentially
     * @param cores Number of used cores for parallel version. Sequential version ignores this
     * argument.
     * @param first First element of a container
     * @param last Last element of a container
     * @param comp Provided comparator to compare container elements
     */
    template <typename ExecutionPolicy, typename Cores, typename RandomIt, typename Compare>
    void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last, Compare comp) {
        static_assert(std::is_integral_v<Cores>, "Integral required for cores");
        impl::_call_sort(std::forward<ExecutionPolicy>(policy), cores, first, last, comp);
    }
}  // namespace mpqsort