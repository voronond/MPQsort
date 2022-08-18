#pragma once

#include <execution>
#include <type_traits>

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
    template <typename RandomIt> void sort(RandomIt first, RandomIt last);

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
    void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last);

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
    void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last);

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
    void sort(RandomIt first, RandomIt last, Compare comp);

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
    void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last, Compare comp);

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
    void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last, Compare comp);
}  // namespace mpqsort

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

    /**
     * @brief Checks whether T is a execution policy type defined in MPQsort library.
     * @tparam T - a type to check
     */
    template <typename T> struct is_execution_policy : std::false_type {};

    /**
     * @brief Helper variable template to check the execution policy
     * @tparam T - a type to check
     */
    template <typename T> inline constexpr bool is_execution_policy_v
        = is_execution_policy<T>::value;

    // Allowed execution policies
    template <> struct is_execution_policy<sequenced_policy_two_way> : std::true_type {};

    template <> struct is_execution_policy<sequenced_policy_multi_way> : std::true_type {};

    template <> struct is_execution_policy<parallel_policy_two_way> : std::true_type {};

    template <> struct is_execution_policy<parallel_policy_multi_way> : std::true_type {};

    // Allowed STD execution policies
    template <> struct is_execution_policy<std::execution::sequenced_policy> : std::true_type {};

    template <> struct is_execution_policy<std::execution::parallel_policy> : std::true_type {};
}  // namespace mpqsort::execution
