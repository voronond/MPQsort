#pragma once

#include <omp.h>

#include <execution>
#include <iterator>
#include <limits>
#include <type_traits>

// TODO remove after all methods implemented
#define UNUSED(x) (void)(x)

/**
 * @brief Contains execution policy definitions for sort algorithms, helper functions and type
 * traits.
 */
namespace mpqsort::execution {
    // Execution policy class declaration

    class sequenced_policy_two_way {};
    class sequenced_policy_multi_way {};
    class parallel_policy_two_way {};
    class parallel_policy_multi_way {};

    // TODO maybe add unseq if possible (standard and self defined)
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
    inline constexpr sequenced_policy_two_way seq_two_way{};

    /**
     * @brief Sequenced execution policy using multiple pivots. The number of pivots will be
     * determined by the algorithm
     */
    inline constexpr sequenced_policy_multi_way seq_multi_way{};

    /**
     * @brief Parallel execution policy using two pivots
     */
    inline constexpr parallel_policy_two_way par_two_way{};

    /**
     * @brief Parallel execution policy using multiple pivots. The number of pivots will be
     * determined by the algorithm
     */
    inline constexpr parallel_policy_multi_way par_multi_way{};

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
    template <> struct _is_execution_policy_helper<sequenced_policy_two_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<sequenced_policy_multi_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<parallel_policy_two_way> : std::true_type {};
    template <> struct _is_execution_policy_helper<parallel_policy_multi_way> : std::true_type {};

    template <> struct _is_parallel_execution_policy_helper<parallel_policy_two_way>
        : std::true_type {};
    template <> struct _is_parallel_execution_policy_helper<parallel_policy_multi_way>
        : std::true_type {};

    // Allowed STD execution policies
    template <> struct _is_execution_policy_helper<std::execution::sequenced_policy>
        : std::true_type {};
    template <> struct _is_execution_policy_helper<std::execution::parallel_policy>
        : std::true_type {};

    template <> struct _is_parallel_execution_policy_helper<std::execution::parallel_policy>
        : std::true_type {};
}  // namespace mpqsort::execution

namespace mpqsort::helpers {
    // Extension of is_same applying decay_t automatically
    template <typename T, typename U> inline constexpr bool _is_same_decay_v
        = std::is_same<std::decay_t<T>, std::decay_t<U>>::value;

    // Helper functions
    template <typename NumPivot, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _seq_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        // TODO implement
        UNUSED(pivot_num);
        std::sort(first, last, comp);
    }

    template <typename NumPivot, typename Cores, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(NumPivot pivot_num, Cores cores, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        // TODO implement
        UNUSED(pivot_num);
        UNUSED(cores);
        std::sort(first, last, comp);
    }

    template <typename NumPivot, typename RandomIt,
              typename Compare = std::less<typename std::iterator_traits<RandomIt>::value_type>>
    void _par_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last,
                             Compare comp = Compare()) {
        // Use all available cores on a machine
        _par_multiway_qsort(pivot_num, omp_get_num_procs(), first, last, comp);
    }

    // Call sort based on policy type
    template <typename ExecutionPolicy, typename... T>
    constexpr void _call_sort(ExecutionPolicy&& policy, T... args) {
        static_assert(
            execution::is_execution_policy_v<decltype(policy)>,
            "Provided ExecutionPolicy is not valid. Use predefined policies from namespace "
            "mpqsort::execution.");

        if constexpr (_is_same_decay_v<ExecutionPolicy, decltype(execution::seq)>) {
            // Call with one pivot
            _seq_multiway_qsort(1, std::forward<T>(args)...);
        } else if constexpr (_is_same_decay_v<ExecutionPolicy, decltype(execution::seq_two_way)>) {
            // Call with two pivots
            _seq_multiway_qsort(2, std::forward<T>(args)...);
        } else if constexpr (_is_same_decay_v<ExecutionPolicy,
                                              decltype(execution::seq_multi_way)>) {
            // Let algorithm decide how many pivots to use
            _seq_multiway_qsort(std::numeric_limits<int>::max, std::forward<T>(args)...);
        } else if constexpr (_is_same_decay_v<ExecutionPolicy, decltype(execution::par)>) {
            // Call with one pivot
            _par_multiway_qsort(1, std::forward<T>(args)...);
        } else if constexpr (_is_same_decay_v<ExecutionPolicy, decltype(execution::par_two_way)>) {
            // Call with two pivots
            _par_multiway_qsort(2, std::forward<T>(args)...);
        } else if constexpr (_is_same_decay_v<ExecutionPolicy,
                                              decltype(execution::par_multi_way)>) {
            // Let algorithm decide how many pivots to use
            _par_multiway_qsort(std::numeric_limits<int>::max, std::forward<T>(args)...);
        } else {
            throw std::invalid_argument(
                "Unknown policy. This should never happen as we test policy type at a beginning of "
                "a function.");
        }
    }
}  // namespace mpqsort::helpers

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
        helpers::_seq_multiway_qsort(1, first, last);
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
        helpers::_call_sort(std::forward<ExecutionPolicy>(policy), first, last);
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
        helpers::_call_sort(std::forward<ExecutionPolicy>(policy), cores, first, last);
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
        helpers::_seq_multiway_qsort(1, first, last, comp);
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
        helpers::_call_sort(std::forward<ExecutionPolicy>(policy), first, last, comp);
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
        helpers::_call_sort(std::forward<ExecutionPolicy>(policy), cores, first, last, comp);
    }
}  // namespace mpqsort