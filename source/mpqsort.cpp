#include <mpqsort.h>
#include <omp.h>

#include <exception>
#include <limits>
#include <vector>

using namespace mpqsort;

// TODO remove after all methods implemented
#define UNUSED(x) (void)(x)

// Helper functions
// Test if received policy is known
template <typename T> inline constexpr void _is_execution_policy_assert() {
    static_assert(execution::is_execution_policy_v<T>,
                  "Provided ExecutionPolicy is not valid. Use predefined policies from namespace "
                  "mpqsort::execution.");
}

// Call sort based on policy type
template <typename ExecutionPolicy, typename... T> constexpr void _call_sort(T... args) {
    _is_execution_policy_assert<ExecutionPolicy>();

    if constexpr (std::is_same_v<ExecutionPolicy, decltype(execution::seq)>) {
        // Call with one pivot
        seq_multiway_qsort(1, std::forward<T...>(args...));
    } else if constexpr (std::is_same_v<ExecutionPolicy, decltype(execution::seq_two_way)>) {
        // Call with two pivots
        seq_multiway_qsort(2, std::forward<T...>(args...));
    } else if constexpr (std::is_same_v<ExecutionPolicy, decltype(execution::seq_multi_way)>) {
        // Let algorithm decide how many pivots to use
        seq_multiway_qsort(std::numeric_limits<int>::max, std::forward<T...>(args...));
    } else if constexpr (std::is_same_v<ExecutionPolicy, decltype(execution::par)>) {
        // Call with one pivot
        par_multi_way_qsort(1, std::forward<T...>(args...));
    } else if constexpr (std::is_same_v<ExecutionPolicy, decltype(execution::par_two_way)>) {
        // Call with two pivots
        par_multi_way_qsort(2, std::forward<T...>(args...));
    } else if constexpr (std::is_same_v<ExecutionPolicy, decltype(execution::par_multi_way)>) {
        // Let algorithm decide how many pivots to use
        par_multi_way_qsort(std::numeric_limits<int>::max, std::forward<T...>(args...));
    }
}

// Implementation of sort functions

template <typename RandomIt> void sort(RandomIt first, RandomIt last) {
    seq_multiway_qsort(1, first, last);
}

// TODO use enable_if if possible
template <typename ExecutionPolicy, typename RandomIt>
void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last) {
    _call_sort(std::forward(policy), first, last);
}

template <typename ExecutionPolicy, typename Cores, typename RandomIt>
void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last) {
    _call_sort(std::forward(policy), cores, first, last);
}

template <typename RandomIt, typename Compare>
void sort(RandomIt first, RandomIt last, Compare comp) {
    seq_multiway_qsort(1, first, last, comp);
}

template <typename ExecutionPolicy, typename RandomIt, typename Compare>
void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last, Compare comp) {
    _call_sort(std::forward(policy), first, last, comp);
}

template <typename ExecutionPolicy, typename Cores, typename RandomIt, typename Compare>
void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last, Compare comp) {
    _call_sort(std::forward(policy), cores, first, last, comp);
}

template <typename NumPivot, typename RandomIt, typename Compare = std::less<RandomIt>>
void seq_multiway_qsort(NumPivot pivot_num, RandomIt first, RandomIt last, Compare comp) {
    // TODO implement
    UNUSED(pivot_num);
    std::sort(first, last, comp);
}

template <typename NumPivot, typename Cores, typename RandomIt,
          typename Compare = std::less<RandomIt>>
void par_multiway_qsort(NumPivot pivot_num, Cores cores, RandomIt first, RandomIt last,
                        Compare comp) {
    // TODO implement
    UNUSED(pivot_num);
    UNUSED(cores);
    std::sort(first, last, comp);
}