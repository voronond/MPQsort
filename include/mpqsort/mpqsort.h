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
//#define DEBUG

#ifdef DEBUG
#    define PRINT_ITERS(base, lp, rp, msg) mpqsort::helpers::print(base, lp, rp, msg)
#else
#    define PRINT_ITERS(base, lp, rp, msg)
#endif

/**
 * @brief Contains execution policy definitions for sort algorithms, helper functions and type
 * traits.
 */
namespace mpqsort::execution {
    // Execution policy class declaration

    class sequenced_policy_three_way {};
    class sequenced_policy_max_way {};
    class parallel_policy_three_way {};
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
     * @brief Sequenced execution policy using multiple pivots. The number of pivots will be
     * determined by the algorithm
     */
    inline constexpr sequenced_policy_max_way seq_max_way{};

    /**
     * @brief Parallel execution policy using two pivots
     */
    inline constexpr parallel_policy_three_way par_three_way{};

    /**
     * @brief Parallel execution policy using multiple pivots. The number of pivots will be
     * determined by the algorithm
     */
    inline constexpr parallel_policy_max_way par_max_way{};

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
    template <> struct _is_execution_policy_helper<sequenced_policy_max_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<parallel_policy_three_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<parallel_policy_max_way> : std::true_type {};

    template <> struct _is_parallel_execution_policy_helper<parallel_policy_three_way>
        : std::true_type {};
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
     * @brief Cacheline of your system
     * Threads process the array in blocks and those blocks are the multiple of cacheline size to
     * prevent false sharing between threads. The value is in bytes.
     */
    // static size_t CACHELINE_SIZE = 64;
    static size_t SEQ_THRESHOLD = 1 << 17;     // based on benchmarks
    static long NO_RECURSION_THRESHOLD = 128;  // based on benchmarks
    /**
     * @brief Maximum supported number of pivots by MPQsort
     * This is the max number of pivots that can MPQsort use during computation. This can't be
     * changed by the user.
     */
    const static long MAX_NUMBER_OF_PIVOTS = 3;  // can't be changed during runtime
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

    // Implementation based on
    // https://github.com/CppCon/CppCon2019/blob/master/Presentations/speed_is_found_in_the_minds_of_people/speed_is_found_in_the_minds_of_people__andrei_alexandrescu__cppcon_2019.pdf
    template <typename RandomBaseIt, typename Compare>
    inline void _heap_insertion_sort(RandomBaseIt base, long lp, long rp, Compare& comp) {
        using std::swap;

        long size = rp - lp;
        // size == 1 if sort 2 elements, otherwise swap element with itself
        if (size < 2) {
            if (!comp(base[lp], base[lp + (size == 1)])) swap(base[lp], base[lp + (size == 1)]);
            return;
        }

        // needs number of elements not index of the last one
        // Partially sort array for effective unbound inser sort
        _make_heap(base + lp, size + 1, comp);
        _unguarded_insertion_sort(base, lp + 1, rp, comp);
    }

    template <long NumPivot> inline auto _get_pivots_indexes(long lp, long rp) {
        using std::swap;
        static_assert(NumPivot <= parameters::MAX_NUMBER_OF_PIVOTS);

        static std::mt19937 en(0);
        std::uniform_int_distribution<long> dist(lp, rp);

        if constexpr (NumPivot == 1) {
            return dist(en);
        } else if constexpr (NumPivot == 2) {
            auto idx1 = dist(en);
            auto idx2 = dist(en);

            while (idx1 == idx2) idx2 = dist(en);  // Generate distinct pivot indexes

            return std::tuple{idx1, idx2};
        }
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
        auto [idx1, idx2] = helpers::_get_pivots_indexes<2>(lp + 1, rp - 1);

        if (comp(base[idx1], base[idx2])) {
            swap(base[lp], base[idx1]);
            swap(base[rp], base[idx2]);
        } else {
            swap(base[lp], base[idx2]);
            swap(base[rp], base[idx1]);
        }

        auto p1 = base[lp], p2 = base[rp];

        // Set boundaries
        auto k2 = lp + 1, k = k2, g = rp - 1;

        while (k <= g) {
            if (comp(base[k], p1)) {
                swap(base[k2], base[k]);
                ++k2;
            } else {
                if (!comp(base[k], p2)) {
                    while (k < g
                           && comp(p2, base[g]))  // Not the same comparison! Should be k <= g?
                        --g;

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
    void _seq_multiway_qsort_inner(NumPivot pivot_num, RandomBaseIt base, long lp, long rp,
                                   Compare& comp, long depth) {
        while (rp - lp > parameters::NO_RECURSION_THRESHOLD && depth > 0) {
            auto [index_p1, index_p2] = _seq_multiway_partition_two_pivots(base, lp, rp, comp);
            _seq_multiway_qsort_inner(pivot_num, base, lp, index_p1 - 1, comp, depth - 1);
            _seq_multiway_qsort_inner(pivot_num, base, index_p1 + 1, index_p2 - 1, comp, depth - 1);
            lp = index_p2 + 1;
        }

        helpers::_heap_insertion_sort(base, lp, rp, comp);
    }

    template <typename NumPivot, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _seq_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        if (last - first <= 1) return;

        _seq_multiway_qsort_inner(pivot_num, first, 0, last - first - 1, comp,
                                  1.5 * std::log(last - first) / std::log(3));
    }

    // PAR
    template <typename NumPivot, typename RandomBaseIt, typename Compare>
    auto _par_multiway_partition(NumPivot pivot_num, RandomBaseIt base, size_t lp, size_t rp,
                                 Compare& comp) {
        // TODO: when multiple pivots supported
        UNUSED(pivot_num);
        using std::swap;

        // Last element as pivot
        // TODO: Select random pivots
        auto pivot = base[rp];  // Access last element value
        // Indexes
        auto i = lp;
        auto j = rp - 1;

        while (i < j) {
            if (!comp(base[i], pivot) && comp(base[j], pivot)) {
                std::swap(base[i], base[j]);
                ++i;
                --j;
            } else {
                if (comp(base[i], pivot)) i++;
                if (!comp(base[j], pivot)) j--;
            }
        }

        if (comp(base[j], pivot)) ++j;
        std::swap(base[j], base[rp]);

        return j;
    }

    template <typename NumPivot, typename RandomBaseIt, typename Compare>
    void _par_multiway_qsort_inner(NumPivot pivot_num, RandomBaseIt base, size_t lp, size_t rp,
                                   Compare& comp) {
        // TODO: when multiple pivots supported
        UNUSED(pivot_num);
        while (lp < rp) {
            if ((rp - lp) <= parameters::SEQ_THRESHOLD) {
                std::sort(base + lp, base + rp + 1, comp);
                return;
            }

            size_t r = _par_multiway_partition(pivot_num, base, lp, rp, comp);
            PRINT_ITERS(base, lp, rp, "After partitioning");
            // Sort if more that 1 elements
            if ((r - lp) > 1) {
#pragma omp task
                _par_multiway_qsort_inner(pivot_num, base, lp, r - 1, comp);
            }
            lp = r + 1;
        }
    }

    // ----- Main implementation END -----

    template <typename NumPivot, typename Cores, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(NumPivot pivot_num, Cores cores, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        if (last - first <= 1) return;

        PRINT_ITERS(first, 0, last - first - 1, "START");
        UNUSED(pivot_num);
        UNUSED(cores);

        // Set OpenMP parameters
        omp_set_max_active_levels(std::numeric_limits<int>::max());  // Allow nested parallelism
                                                                     // omp_set_nested(1);
#pragma omp parallel firstprivate(pivot_num, first, last, comp) num_threads(cores)
        {
#pragma omp single
            // TODO: change when multiple pivots supported!!!
            // Convert iterators for more general format
            _par_multiway_qsort_inner(1, first, 0, last - first - 1, comp);
        }
        PRINT_ITERS(first, 0, last - first - 1, "END");
    }

    // Wrapper for different arguments
    template <typename NumPivot, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        // Use all available cores on a machine
        // TODO: Change when multiple pivots supported!!!
        UNUSED(pivot_num);
        _par_multiway_qsort(1, omp_get_num_procs(), first, last, comp);
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
                                                       decltype(execution::seq_max_way)>) {
            // Let algorithm decide how many pivots to use
            _seq_multiway_qsort(parameters::MAX_NUMBER_OF_PIVOTS, std::forward<T>(args)...);
        } else if constexpr (helpers::_is_same_decay_v<ExecutionPolicy, decltype(execution::par)>) {
            // Call with one pivot
            _par_multiway_qsort(1, std::forward<T>(args)...);
        } else if constexpr (helpers::_is_same_decay_v<ExecutionPolicy,
                                                       decltype(execution::par_three_way)>) {
            // Call with two pivots
            _par_multiway_qsort(2, std::forward<T>(args)...);
        } else if constexpr (helpers::_is_same_decay_v<ExecutionPolicy,
                                                       decltype(execution::par_max_way)>) {
            // Let algorithm decide how many pivots to use
            _par_multiway_qsort(parameters::MAX_NUMBER_OF_PIVOTS, std::forward<T>(args)...);
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