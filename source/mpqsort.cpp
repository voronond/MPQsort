#include <mpqsort.h>
#include <omp.h>

#include <exception>
#include <vector>

// TODO remove after all methods implemented
#define UNUSED(x) (void)(x)

using namespace mpqsort;

template <typename RandomIt> void sort(RandomIt first, RandomIt last) {
    UNUSED(first);
    UNUSED(last);
    throw std::logic_error("Not implemented");
}

template <typename ExecutionPolicy, typename RandomIt>
void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last) {
    UNUSED(policy);
    UNUSED(first);
    UNUSED(last);
    throw std::logic_error("Not implemented");
}

template <typename ExecutionPolicy, typename Cores, typename RandomIt>
void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last) {
    UNUSED(policy);
    UNUSED(first);
    UNUSED(last);
    UNUSED(cores);
    throw std::logic_error("Not implemented");
}

template <typename RandomIt, typename Compare>
void sort(RandomIt first, RandomIt last, Compare comp) {
    UNUSED(first);
    UNUSED(last);
    UNUSED(comp);
    throw std::logic_error("Not implemented");
}

template <typename ExecutionPolicy, typename RandomIt, typename Compare>
void sort(ExecutionPolicy&& policy, RandomIt first, RandomIt last, Compare comp) {
    UNUSED(policy);
    UNUSED(first);
    UNUSED(last);
    UNUSED(comp);
    throw std::logic_error("Not implemented");
}

template <typename ExecutionPolicy, typename Cores, typename RandomIt, typename Compare>
void sort(ExecutionPolicy&& policy, Cores cores, RandomIt first, RandomIt last, Compare comp) {
    UNUSED(policy);
    UNUSED(cores);
    UNUSED(first);
    UNUSED(last);
    UNUSED(comp);
    throw std::logic_error("Not implemented");
}
