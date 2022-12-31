#pragma once

#include <omp.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <execution>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <type_traits>

//#define DEBUG
//#define MEASURE
//#define TIME_MEASURE

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

#ifdef MEASURE
size_t NUM_SWAP = 0;
size_t NUM_COMP = 0;
#    define MEASURE_INIT() \
        NUM_SWAP = 0;      \
        NUM_COMP = 0
#    define MEASURE_SWAP() ++NUM_SWAP
#    define MEASURE_SWAP_N(n) NUM_SWAP += n
#    define MEASURE_COMP() ++NUM_COMP
#    define MEASURE_COMP_N(n) NUM_COMP += n
#    define MEASURE_RESULTS(msg)  \
        std::cout << msg << ": "; \
        std::cout << "SWAP: " << NUM_SWAP << " COMP: " << NUM_COMP << "\n"
#else
#    define MEASURE_INIT()
#    define MEASURE_SWAP()
#    define MEASURE_SWAP_N(n)
#    define MEASURE_COMP()
#    define MEASURE_COMP_N(n)
#    define MEASURE_RESULTS(msg)
#endif

#ifdef TIME_MEASURE
#    define TIME_MEASURE_START(name) \
        std::chrono::steady_clock::time_point name = std::chrono::steady_clock::now();
#    define TIME_MEASURE_END(name_start, name_end, msg)                                           \
        std::chrono::steady_clock::time_point name_end = std::chrono::steady_clock::now();        \
        std::cout << msg << " "                                                                   \
                  << std::chrono::duration_cast<std::chrono::milliseconds>(name_end - name_start) \
                         .count()                                                                 \
                  << std::endl
#else
#    define TIME_MEASURE_START(name)
#    define TIME_MEASURE_END(name_start, name_end, msg)
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
    static long BLOCK_SIZE = 1 << 10;

    /**
     * @brief Length of one cache line
     * Set to one cache line size to maximize performance (number of bytes);
     */
    // static long CACHE_LINE_SIZE = 64;

    /**
     * @brief When to switch to sequential algorithm from parallel
     * Threshold for switching from parallel version to sequential algorithm.
     */
    static long SEQ_THRESHOLD = 1 << 18;

    /**
     * @brief Switch to non-recursive algorithm
     * Threshold for switching to non-recursive algorithm (heapsort + insertion sort).
     */

    static long NO_RECURSION_THRESHOLD = 128;
    /**
     * @brief Maximum supported number of pivots by MPQsort
     * This is the max number of pivots that can MPQsort use during computation. This can't be
     * changed by the user.
     */
    const static long MAX_NUMBER_OF_PIVOTS = 3;  // can't be changed during runtime

    /**
     * @brief Number of pivots used in multiway parallel partitioning
     * Number for pivots used to divide array in PAR_PARTITION_NUM_PIVOTS + 1 segments. This is done
     * only once and parallel multiway qsort is called afterwards on each such segment.
     */
    static int PAR_PARTITION_NUM_PIVOTS = 0;

    /**
     * @brief Number of elements to chose a pivot from in a first parallel multiway partitioning
     * Determined hos many elements from an input array we want to consider to chose one pivot when
     * we do the first parallel partitioning
     */
    static long ONE_PIVOT_PAR_MULT_PARTITIONING_SAMPLE_SIZE = 100;

    /**
     * @brief Number of elements to chose a pivot from an array in parallel multiway qsorts
     * Array was already splitted up in num_pivots + 1 segments so parallel sorts does not need to
     * be as precise.
     */
    static long ONE_PIVOT_PAR_SORT_SAMPLE_SIZE = 3;
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

        if (lp >= rp) return;

        long size = rp - lp;
        // size == 1 if sort 2 elements, otherwise swap element with itself
        if (size < 2) {
            if (!comp(base[lp], base[lp + (size == 1)])) swap(base[lp], base[lp + (size == 1)]);
            return;
        }

        _make_heap(base + lp, size + 1, comp);
        _unguarded_insertion_sort(base, lp + 1, rp, comp);
    }

    template <typename Pair, typename Comparator>
    inline auto _get_median_from_three(Pair&& el1, Pair&& el2, Pair&& el3, Comparator& comp) {
        using std::swap;

        MEASURE_COMP_N(3);
        if (comp(el2.first, el1.first)) {
            MEASURE_SWAP();
            swap(el1.second, el2.second);
        }
        if (comp(el3.first, el2.first)) {
            MEASURE_SWAP();
            swap(el2.second, el3.second);
        }
        if (comp(el3.first, el1.first)) {
            MEASURE_SWAP();
            swap(el1.second, el3.second);
        }

        return el2.second;
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivot_index_median(RandomBaseIt base, long lp, long rp, Comparator& comp) {
        using ValueType = typename std::iterator_traits<RandomBaseIt>::value_type;
        long sample_size = 3;

        auto size = rp - lp + 1;

        if (size < sample_size) return size / 2 + lp;

        std::pair<ValueType, long> el1, el2, el3;
        el1 = std::make_pair(base[lp], lp);
        el2 = std::make_pair(base[size / 2], size / 2);
        el3 = std::make_pair(base[rp], rp);

        return _get_median_from_three(el1, el2, el3, comp);
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivot_index(RandomBaseIt base, long lp, long rp, Comparator& comp) {
        auto size = rp - lp + 1;
        const auto sample_size = parameters::ONE_PIVOT_PAR_SORT_SAMPLE_SIZE;
        std::vector<std::pair<typename std::iterator_traits<RandomBaseIt>::value_type, long>>
            samples;
        samples.reserve(sample_size);

        // If not enough elements for sampling
        if (size < sample_size) {
            return size * 1 / 2 + lp;
        } else {
            // Get sample elements
            for (long i = 0; i < sample_size; ++i) {
                auto index = size * i / sample_size + lp;
                samples.emplace_back(std::make_pair(base[index], index));
            }

            // Sort samples based on provided comp
            std::sort(samples.begin(), samples.end(),
                      [&](auto& a, auto& b) { return comp(a.first, b.first); });

            MEASURE_COMP_N(sample_size * std::log(sample_size));
            MEASURE_SWAP_N(sample_size * std::log(sample_size));

            return samples[sample_size * 1 / 2].second;
        }
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivots_indexes_two_medians(RandomBaseIt base, long lp, long rp,
                                                Comparator& comp) {
        using std::swap;
        // using ValueType = typename std::iterator_traits<RandomBaseIt>::value_type;
        long sample_size = 3;

        auto size = rp - lp + 1;

        if (size < sample_size) {
            auto idx1 = size * 1 / 3 + lp, idx2 = size * 2 / 3 + lp;

            MEASURE_COMP();
            if (comp(base[idx1], base[idx2])) {
                return std::tuple{idx1, idx2};
            } else {
                return std::tuple{idx2, idx1};
            }
        }

        long i1, i2, i3;

        i1 = size / 4 + lp;
        i2 = size / 4 * 2 + lp;
        i3 = size / 4 * 3 + lp;

        MEASURE_COMP_N(3);
        if (comp(base[i2], base[i1])) {
            MEASURE_SWAP();
            swap(i1, i2);
        }
        if (comp(base[i3], base[i2])) {
            MEASURE_SWAP();
            swap(i2, i3);
        }
        if (comp(base[i3], base[i1])) {
            MEASURE_SWAP();
            swap(i1, i3);
        }

        return std::tuple{i1, i3};
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivots_indexes_two(RandomBaseIt base, long lp, long rp, Comparator& comp) {
        using std::swap;

        auto size = rp - lp + 1;
        const auto sample_size = parameters::ONE_PIVOT_PAR_SORT_SAMPLE_SIZE * 2;
        std::vector<std::pair<typename std::iterator_traits<RandomBaseIt>::value_type, long>>
            samples;
        samples.reserve(sample_size);

        // If not enough elements for sampling
        if (size < sample_size) {
            auto idx1 = size * 1 / 3 + lp, idx2 = size * 2 / 3 + lp;

            MEASURE_COMP();
            if (comp(base[idx1], base[idx2])) {
                return std::tuple{idx1, idx2};
            } else {
                return std::tuple{idx2, idx1};
            }
        } else {
            // Get sample elements
            for (long i = 0; i < sample_size; ++i) {
                auto index = size * i / sample_size + lp;
                samples.emplace_back(std::make_pair(base[index], index));
            }

            // Sort samples based on provided comp
            std::sort(samples.begin(), samples.end(),
                      [&](auto& a, auto& b) { return comp(a.first, b.first); });

            MEASURE_COMP_N(sample_size * std::log(sample_size));
            MEASURE_SWAP_N(sample_size * std::log(sample_size));

            return std::tuple{samples[sample_size * 1 / 3].second,
                              samples[sample_size * 2 / 3].second};
        }
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivot_indexes_three_medians(RandomBaseIt base, long lp, long rp,
                                                 Comparator& comp) {
        using std::swap;
        // using ValueType = typename std::iterator_traits<RandomBaseIt>::value_type;
        long sample_size = 3;

        auto size = rp - lp + 1;

        if (size < sample_size) {
            auto idx1 = size * 1 / 4 + lp, idx2 = size * 2 / 4 + lp, idx3 = size * 3 / 4 + lp;

            // Sort pivots
            MEASURE_COMP_N(3);
            if (comp(base[idx2], base[idx1])) {
                swap(idx1, idx2);
            }
            if (comp(base[idx3], base[idx1])) {
                swap(idx1, idx3);
            }
            if (comp(base[idx3], base[idx2])) {
                swap(idx2, idx3);
            }

            return std::tuple{idx1, idx2, idx3};
        }

        long i1, i2, i3;

        i1 = size * 1 / 4 + lp;
        i2 = size * 2 / 4 + lp;
        i3 = size * 3 / 4 + lp;

        MEASURE_COMP_N(3);
        if (comp(base[i2], base[i1])) {
            MEASURE_SWAP();
            swap(i1, i2);
        }
        if (comp(base[i3], base[i1])) {
            MEASURE_SWAP();
            swap(i1, i3);
        }
        if (comp(base[i3], base[i2])) {
            MEASURE_SWAP();
            swap(i2, i3);
        }

        return std::tuple{i1, i2, i3};
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivot_indexes_three(RandomBaseIt base, long lp, long rp, Comparator& comp) {
        using std::swap;

        auto size = rp - lp + 1;
        const auto sample_size = parameters::ONE_PIVOT_PAR_SORT_SAMPLE_SIZE * 3;
        std::vector<std::pair<typename std::iterator_traits<RandomBaseIt>::value_type, long>>
            samples;
        samples.reserve(sample_size);

        // If not enough elements for sampling
        if (size < sample_size) {
            auto idx1 = size * 1 / 4 + lp, idx2 = size * 2 / 4 + lp, idx3 = size * 3 / 4 + lp;

            // Sort pivots
            MEASURE_COMP_N(3);
            if (comp(base[idx2], base[idx1])) {
                swap(idx1, idx2);
            }
            if (comp(base[idx3], base[idx1])) {
                swap(idx1, idx3);
            }
            if (comp(base[idx3], base[idx2])) {
                swap(idx2, idx3);
            }

            return std::tuple{idx1, idx2, idx3};
        } else {
            // Get sample elements
            for (long i = 0; i < sample_size; ++i) {
                auto index = size * i / sample_size + lp;
                samples.emplace_back(std::make_pair(base[index], index));
            }

            // Sort samples based on provided comp
            std::sort(samples.begin(), samples.end(),
                      [&](auto& a, auto& b) { return comp(a.first, b.first); });

            MEASURE_COMP_N(sample_size * std::log(sample_size));
            MEASURE_SWAP_N(sample_size * std::log(sample_size));

            return std::tuple{samples[sample_size * 1 / 4].second,
                              samples[sample_size * 2 / 4].second,
                              samples[sample_size * 3 / 4].second};
        }
    }

    template <typename RandomBaseIt, typename Comparator>
    inline auto _get_pivots(RandomBaseIt base, long lp, long rp, long num_pivots,
                            Comparator& comp) {
        auto size = rp - lp + 1;
        const auto sample_size
            = parameters::ONE_PIVOT_PAR_MULT_PARTITIONING_SAMPLE_SIZE * num_pivots;

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

    template <bool Par = false, typename RandomBaseIt, typename Compare>
    inline auto _seq_partition_one_pivot(RandomBaseIt base, long lp, long rp, Compare& comp) {
        // Use optimal swap method
        using std::swap;
        long idx;

#ifdef MEASURE
        idx = helpers::_get_pivot_index_median(base, lp, rp, comp);
#else
        if constexpr (Par) {
            idx = helpers::_get_pivot_index(base, lp, rp, comp);
        } else {
            idx = helpers::_get_pivot_index_median(base, lp, rp, comp);
        }
#endif

        auto p = base[idx];

        auto i = lp - 1, j = rp + 1;

        while (true) {
            do {
                MEASURE_COMP();
                ++i;
            } while (comp(base[i], p));

            do {
                MEASURE_COMP();
                --j;
            } while (comp(p, base[j]));

            if (i >= j) {
                return j;
            }

            MEASURE_SWAP();
            swap(base[i], base[j]);
        }
    }

    template <typename RandomBaseIt, typename Compare>
    void _seq_qsort_inner_hoare(RandomBaseIt base, long lp, long rp, Compare& comp, long depth) {
        while (rp - lp > parameters::NO_RECURSION_THRESHOLD && depth > 0) {
            auto index_p = _seq_partition_one_pivot(base, lp, rp, comp);
            _seq_qsort_inner_hoare(base, lp, index_p, comp, depth - 1);
            lp = index_p + 1;
        }

        helpers::_heap_insertion_sort(base, lp, rp, comp);
    }

    template <bool Par = false, typename RandomBaseIt, typename Compare>
    inline auto _seq_multiway_partition_two_pivots(RandomBaseIt base, long lp, long rp,
                                                   Compare& comp) {
        // Use optimal swap method
        using std::swap;
        long idx1, idx2;

// Get pivots
#ifdef MEASURE
        tie(idx1, idx2) = helpers::_get_pivots_indexes_two_medians(base, lp + 1, rp - 1, comp);
#else
        if constexpr (Par) {
            std::tie(idx1, idx2) = helpers::_get_pivots_indexes_two(base, lp + 1, rp - 1, comp);
        } else {
            std::tie(idx1, idx2)
                = helpers::_get_pivots_indexes_two_medians(base, lp + 1, rp - 1, comp);
        }
#endif

        MEASURE_SWAP_N(2);
        swap(base[lp], base[idx1]);
        swap(base[rp], base[idx2]);

        const auto p1 = base[lp], p2 = base[rp];

        // Set boundaries
        auto k2 = lp + 1, k = k2, g = rp - 1;

        while (k <= g) {
            MEASURE_COMP();
            if (comp(base[k], p1)) {
                MEASURE_SWAP();
                swap(base[k2], base[k]);
                ++k2;
            } else {
                MEASURE_COMP();
                if (!comp(base[k], p2)) {
                    MEASURE_COMP();
                    while (k < g && comp(p2, base[g])) {
                        MEASURE_COMP();
                        --g;
                    }

                    MEASURE_COMP();
                    if (!comp(base[g], p1)) {
                        MEASURE_SWAP();
                        swap(base[k], base[g]);
                    } else {
                        MEASURE_SWAP_N(2);
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
        MEASURE_SWAP_N(2);

        return std::tuple{k2 - 1, g + 1};
    }

    template <typename RandomBaseIt, typename Compare>
    void _seq_multiway_qsort_inner_ybb(RandomBaseIt base, long lp, long rp, Compare& comp,
                                       long depth) {
        while (rp - lp > parameters::NO_RECURSION_THRESHOLD && depth > 0) {
            auto [index_p1, index_p2] = _seq_multiway_partition_two_pivots(base, lp, rp, comp);
            _seq_multiway_qsort_inner_ybb(base, lp, index_p1 - 1, comp, depth - 1);
            _seq_multiway_qsort_inner_ybb(base, index_p1 + 1, index_p2 - 1, comp, depth - 1);
            lp = index_p2 + 1;
        }

        helpers::_heap_insertion_sort(base, lp, rp, comp);
    }

    template <bool Par = false, typename RandomBaseIt, typename Compare>
    inline auto _seq_multiway_partition_three_pivots(RandomBaseIt base, const long lp,
                                                     const long rp, Compare& comp) {
        // Use optimal swap method
        using std::swap;
        long idx1, idx2, idx3;

// Get pivots
#ifdef MEASURE
        std::tie(idx1, idx2, idx3)
            = helpers::_get_pivot_indexes_three_medians(base, lp + 2, rp - 1, comp);
#else
        if constexpr (Par) {
            std::tie(idx1, idx2, idx3)
                = helpers::_get_pivot_indexes_three(base, lp + 2, rp - 1, comp);
        } else {
            std::tie(idx1, idx2, idx3)
                = helpers::_get_pivot_indexes_three_medians(base, lp + 2, rp - 1, comp);
        }
#endif

        // Move pivots to array edges
        MEASURE_SWAP_N(3);
        swap(base[idx1], base[lp]);
        swap(base[idx2], base[lp + 1]);
        swap(base[idx3], base[rp]);

        // Indexes
        long k2, k, g, g2;
        k2 = k = lp + 2;
        g = g2 = rp - 1;

        const auto p1 = base[lp], p2 = base[lp + 1], p3 = base[rp];

        while (k <= g) {
            MEASURE_COMP();
            while (k <= g && comp(base[k], p2)) {
                MEASURE_COMP();
                MEASURE_COMP();
                if (comp(base[k], p1)) {
                    MEASURE_SWAP();
                    swap(base[k2], base[k]);
                    ++k2;
                }
                ++k;
            }
            MEASURE_COMP();
            while (k <= g && comp(p2, base[g])) {
                MEASURE_COMP();
                MEASURE_COMP();
                if (comp(p3, base[g])) {
                    MEASURE_SWAP();
                    swap(base[g], base[g2]);
                    --g2;
                }
                --g;
            }
            if (k <= g) {
                MEASURE_COMP();
                if (comp(p3, base[k])) {
                    MEASURE_COMP();
                    if (comp(base[g], p1)) {
                        MEASURE_SWAP_N(3);
                        helpers::_cyclic_shift_right(base, k2, k, g2, g);
                        ++k2;
                    } else {
                        MEASURE_SWAP_N(2);
                        helpers::_cyclic_shift_left(base, k, g, g2);
                    }
                    ++k;
                    --g;
                    --g2;
                } else {
                    MEASURE_COMP();
                    if (comp(base[g], p1)) {
                        MEASURE_SWAP_N(2);
                        helpers::_cyclic_shift_right(base, k2, k, g);
                        ++k2;
                    } else {
                        MEASURE_SWAP();
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

        MEASURE_SWAP_N(4);

        return std::tuple{k2, k, g2};
    }

    template <typename RandomBaseIt, typename Compare>
    void _seq_multiway_qsort_inner_waterloo(RandomBaseIt base, long lp, long rp, Compare& comp,
                                            long depth) {
        while (rp - lp > parameters::NO_RECURSION_THRESHOLD && depth > 0) {
            auto indexes = _seq_multiway_partition_three_pivots(base, lp, rp, comp);
            _seq_multiway_qsort_inner_waterloo(base, lp, std::get<0>(indexes) - 1, comp, depth - 1);
            _seq_multiway_qsort_inner_waterloo(base, std::get<0>(indexes) + 1,
                                               std::get<1>(indexes) - 1, comp, depth - 1);
            _seq_multiway_qsort_inner_waterloo(base, std::get<1>(indexes) + 1,
                                               std::get<2>(indexes) - 1, comp, depth - 1);
            lp = std::get<2>(indexes) + 1;
        }

        helpers::_heap_insertion_sort(base, lp, rp, comp);
    }

    template <typename NumPivot, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _seq_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        if (last - first <= 1) return;

        // Init measures of swaps and compares
        MEASURE_INIT();
        if (pivot_num == 1) {
            _seq_qsort_inner_hoare(first, 0, last - first - 1, comp,
                                   1.5 * std::log(last - first) / std::log(2));
            MEASURE_RESULTS("Num pivots: 1");
        } else if (pivot_num == 2) {
            _seq_multiway_qsort_inner_ybb(first, 0, last - first - 1, comp,
                                          1.5 * std::log(last - first) / std::log(3));
            MEASURE_RESULTS("Num pivots: 2");
        } else if (pivot_num == 3) {
            _seq_multiway_qsort_inner_waterloo(first, 0, last - first - 1, comp,
                                               1.5 * std::log(last - first) / std::log(4));
            MEASURE_RESULTS("Num pivots: 3");
        } else {
            throw std::invalid_argument("Unknown number of pivots, this should never happen");
        }
    }

    // PAR
    template <typename RandomBaseIt, typename Compare>
    void _par_multiway_qsort_inner_waterloo(RandomBaseIt base, long lp, long rp, Compare& comp,
                                            long depth) {
        // Only recursive calling is in parallel, partitioning does not need to be thanks to
        // multiway_parallel_partitioning at the beginning
        while (rp - lp > parameters::NO_RECURSION_THRESHOLD && depth > 0) {
            auto indexes = _seq_multiway_partition_three_pivots<true>(base, lp, rp, comp);

            if (std::get<0>(indexes) - lp > parameters::SEQ_THRESHOLD) {
#pragma omp task
                _par_multiway_qsort_inner_waterloo(base, lp, std::get<0>(indexes) - 1, comp,
                                                   depth - 1);
            } else {
                _seq_multiway_qsort_inner_waterloo(base, lp, std::get<0>(indexes) - 1, comp,
                                                   depth - 1);
            }

            if (std::get<1>(indexes) - std::get<0>(indexes) > parameters::SEQ_THRESHOLD) {
#pragma omp task
                _par_multiway_qsort_inner_waterloo(base, std::get<0>(indexes) + 1,
                                                   std::get<1>(indexes) - 1, comp, depth - 1);
            } else {
                _seq_multiway_qsort_inner_waterloo(base, std::get<0>(indexes) + 1,
                                                   std::get<1>(indexes) - 1, comp, depth - 1);
            }

            if (std::get<2>(indexes) - std::get<1>(indexes) > parameters::SEQ_THRESHOLD) {
#pragma omp task
                _par_multiway_qsort_inner_waterloo(base, std::get<1>(indexes) + 1,
                                                   std::get<2>(indexes) - 1, comp, depth - 1);
            } else {
                _seq_multiway_qsort_inner_waterloo(base, std::get<1>(indexes) + 1,
                                                   std::get<2>(indexes) - 1, comp, depth - 1);
            }

            lp = std::get<2>(indexes) + 1;
        }

        helpers::_heap_insertion_sort(base, lp, rp, comp);
    }

    template <typename RandomBaseIt, typename Cores, typename Compare>
    auto _par_multiway_partition(long num_pivots, Cores cores, RandomBaseIt base, long lp, long rp,
                                 Compare& comp) {
        using std::swap;
        using ValueType = typename std::iterator_traits<RandomBaseIt>::value_type;

        TIME_MEASURE_START(parallel_partition);
        TIME_MEASURE_START(par_preparation);

        const int num_segments = num_pivots + 1;
        int num_segments_left = num_segments;
        // Precompute number of comparisons
        int num_element_comparisons = std::log2(num_pivots);

        // Get pivots
        auto pivots = helpers::_get_pivots(base, lp, rp, num_pivots, comp);

        // Indexes in blocks
        std::vector<long> segment_idx(num_segments);
        std::vector<long> segment_boundary(num_segments);

        // If number of pivots is greater than cores, distribute pivots uniformly between threads
        if (pivots.size() > (size_t)cores) {
            // Calculate upper-bound num pivots for each thread
            int num_pivot_for_thread = std::ceil((double)pivots.size() / cores);
            size_t total_pivots = num_pivot_for_thread * cores;

            std::vector<int> t_num_pivots(cores, num_pivot_for_thread);

            // Remove extra pivots uniformly
            for (size_t i = 0; i < t_num_pivots.size() && total_pivots > pivots.size();
                 ++i, --total_pivots) {
                --t_num_pivots[i];
            }

            std::vector<int> t_first_pivots(cores);
            // Get beginning of each thread pivots
            std::exclusive_scan(t_num_pivots.begin(), t_num_pivots.end(), t_first_pivots.begin(),
                                0);

#pragma omp parallel num_threads(cores)
            {
                const int tid = omp_get_thread_num();
                const int t_num_pivot = t_num_pivots[tid];
                const int t_first_pivot = t_first_pivots[tid];
                std::vector<long> num_els_smaller(t_num_pivot, 0);

                for (long i = lp; i <= rp; ++i) {
                    for (int j = 0, first = t_first_pivot; j < t_num_pivot; ++j, ++first) {
                        if (comp(base[i], pivots[first])) {
                            for (; j < t_num_pivot; ++j) ++num_els_smaller[j];

                            break;
                        }
                    }
                }

                for (int i = 0, first = t_first_pivot; i < t_num_pivot; ++i, ++first) {
                    segment_boundary[first] = num_els_smaller[i];
                }
            }
        } else {
            // Find boundaries of segments
#pragma omp parallel num_threads(cores)
            {
                const int tid = omp_get_thread_num();

                if ((size_t)tid < pivots.size()) {
                    long num_el_smaller = 0;
                    for (long i = lp; i <= rp; ++i) {
                        if (comp(base[i], pivots[tid])) ++num_el_smaller;
                    }

                    segment_boundary[tid] = num_el_smaller;
                }
            }
        }

        segment_boundary.back() = rp + 1;
        segment_idx[0] = 0;
        for (size_t i = 0; i < segment_boundary.size() - 1; ++i) {
            segment_idx[i + 1] = segment_boundary[i];
        }

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
        std::vector<long> elements_table_insertions(num_segments * table_height);
        std::vector<long> elements_table_insertions_index(num_segments, 0);

        // Table of elements belonging to segment
        std::vector<ValueType> elements_table(num_segments * table_height);
        std::vector<long> elements_table_index(num_segments, 0);

        // Find first unclean segment
        int current_segment = 0;
        find_unprocessed_segment(current_segment);

        TIME_MEASURE_END(par_preparation, par_preparation_end, "Parallel preparation");
        TIME_MEASURE_START(par_main);

#pragma omp parallel shared(elements_table_insertions, elements_table_insertions_index,          \
                            elements_table, elements_table_index, segment_idx, segment_boundary, \
                            base, pivots)                                                        \
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
                            index = elements_table_insertions_index[tmp_el_segment]++;
                            elements_table_insertions[INDEX(tmp_el_segment, index, table_height)]
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
        }
        TIME_MEASURE_END(par_main, par_main_end, "Parallel main");
        TIME_MEASURE_START(par_cleanup);

// Insert elements from tables in an array
#pragma omp parallel for num_threads(cores)
        for (size_t i = 0; i < elements_table_insertions_index.size(); ++i) {
            for (int j = 0; j < elements_table_insertions_index[i]; ++j) {
                base[elements_table_insertions[INDEX(i, j, table_height)]]
                    = elements_table[INDEX(i, j, table_height)];
            }
        }

        TIME_MEASURE_END(par_cleanup, par_cleanup_end, "Parallel cleanup");
        TIME_MEASURE_END(parallel_partition, parallel_partition_end, "Parallel partitioning");

        return segment_boundary;
    }

    template <typename NumPivot, typename Cores, typename RandomBaseIt, typename Compare>
    void _par_multiway_qsort_inner(NumPivot pivot_num, Cores cores, RandomBaseIt base, size_t lp,
                                   size_t rp, Compare& comp) {
        auto boundaries = _par_multiway_partition(pivot_num, cores, base, lp, rp, comp);

        // Create ranges {start, end} for segments
        std::vector<std::pair<long, long>> segment_ranges;
        segment_ranges.reserve(boundaries.size());

        long start = 0;
        for (auto b : boundaries) {
            if (start == boundaries.back()) break;
            segment_ranges.emplace_back(std::make_pair(start, b - 1));
            start = b;
        }

        // Sort segments based on their length from the longest to shortest
        std::sort(segment_ranges.begin(), segment_ranges.end(), [](const auto& p1, const auto& p2) {
            return (p1.second - p1.first) > (p2.second - p2.first);
        });

        // Parallel multiway sort of each segment
#pragma omp parallel for schedule(dynamic) num_threads(cores)
        for (size_t i = 0; i < segment_ranges.size(); i++) {
            _par_multiway_qsort_inner_waterloo(base, segment_ranges[i].first,
                                               segment_ranges[i].second, comp,
                                               1.5 * std::log(rp - lp) / std::log(4));
        }
    }

    // ----- Main implementation END -----

    template <typename Cores, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(Cores cores, RandomIt first, RandomIt last, Compare comp = Compare()) {
        if (last - first <= 1) return;

        // Set to best number of pivots, always greater than num threads but power of 2 - 1
        long pivot_num
            = std::max(std::max((1 << (long)std::log2(cores)) - 1, parameters::PAR_PARTITION_NUM_PIVOTS), 1);

        // Allow MAX nested parallelism
        omp_set_max_active_levels(std::numeric_limits<int>::max());

        // Should not normally happen, but we set SEQ_THRESHOLD to 0 during testing
        if (last - first > std::max(parameters::SEQ_THRESHOLD, pivot_num)) {
            _par_multiway_qsort_inner(pivot_num, cores, first, 0, last - first - 1, comp);
        } else {
            _seq_multiway_qsort(3, first, last, comp);
        }
    }

    // Wrapper for different arguments
    template <typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(RandomIt first, RandomIt last, Compare comp = Compare()) {
        _par_multiway_qsort(omp_get_max_threads(), first, last, comp);
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
            _par_multiway_qsort(std::forward<T>(args)...);
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

