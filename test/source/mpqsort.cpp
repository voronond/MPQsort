#include <mpqsort/mpqsort.h>
#include <mpqsort/version.h>

#include <algorithm>
#include <any>
#include <catch2/catch.hpp>
#include <iostream>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <vector>

using namespace mpqsort;

// Declaration of helpers and variables

struct _TAG {
    const std::string IS_EXECUTION_POLICY = "[is_execution_policy]";
    const std::string IS_PARALLEL_EXECUTION_POLICY = "[is_parallel_execution_policy]";
    const std::string SEQ = "[seq_sort]";
    const std::string SEQ_three_way = "[seq_three_way_sort]";
    const std::string SEQ_MULTIWAY = "[seq_multiway_sort]";
    const std::string SEQ_ALL = SEQ + SEQ_three_way + SEQ_MULTIWAY;
    const std::string PAR = "[par_sort]";
    const std::string PAR_MULTIWAY = "[par_multiway_sort]";
    const std::string PAR_ALL = PAR + PAR_MULTIWAY;
    const std::string SORT_ALL = SEQ_ALL + PAR_ALL;
} TAG;

// Combine tags and create alias for catch2 (supports only C strings)
CATCH_REGISTER_TAG_ALIAS("[@seq_all]", TAG.SEQ_ALL.c_str())
CATCH_REGISTER_TAG_ALIAS("[@par_all]", TAG.PAR_ALL.c_str())
CATCH_REGISTER_TAG_ALIAS("[@sort_all]", TAG.SORT_ALL.c_str())

// Declare allowed execution policies
using ALLOWED_EXECUTION_POLICIES
    = std::tuple<decltype(execution::par), decltype(execution::seq_three_way),
                 decltype(execution::seq_four_way), decltype(execution::seq_max_way),
                 decltype(execution::seq)>;

// Declare parallel execution policies
using PARALLEL_EXECUTION_POLICIES = std::tuple<decltype(execution::par)>;

// Declare sequential execution policies
using SEQUENTIAL_EXECUTION_POLICIES
    = std::tuple<decltype(execution::seq_three_way), decltype(execution::seq_four_way),
                 decltype(execution::seq_max_way), decltype(execution::seq)>;

// Some disallowed execution policies
using DISALLOWED_EXECUTION_POLICIES
    = std::tuple<int, bool, char, float, std::string, decltype(std::execution::par_unseq)>;

// Setup

TEST_CASE("Not a test case, just setup!") {
    // Turn off optimization for small arrays
    mpqsort::parameters::SEQ_THRESHOLD = 0;
    mpqsort::parameters::NO_RECURSION_THRESHOLD = 6;
}
// TEST CASES

TEST_CASE("Sort version") {
    static_assert(std::string_view(MPQSORT_VERSION) == std::string_view("1.0"));
    CHECK(std::string(MPQSORT_VERSION) == std::string("1.0"));
}

// Type trait should return true for all allowed policies
TEMPLATE_LIST_TEST_CASE("Execution policy type trait allowed policies", TAG.IS_EXECUTION_POLICY,
                        ALLOWED_EXECUTION_POLICIES) {
    // Test type trait helper
    REQUIRE(execution::is_execution_policy_v<TestType>);
    REQUIRE(execution::is_execution_policy_v<const TestType>);
    REQUIRE(execution::is_execution_policy_v<const TestType&>);
    REQUIRE(execution::is_execution_policy_v<const TestType&&>);
    REQUIRE(execution::is_execution_policy_v<TestType&&>);

    // Test type trait
    REQUIRE(execution::is_execution_policy<TestType>::value);
    REQUIRE(execution::is_execution_policy<const TestType>::value);
    REQUIRE(execution::is_execution_policy<const TestType&>::value);
    REQUIRE(execution::is_execution_policy<const TestType&&>::value);
    REQUIRE(execution::is_execution_policy<TestType&&>::value);
}

// Type trait should return false for all disallowed policies
TEMPLATE_LIST_TEST_CASE("Execution policy type trait disallowed policies", TAG.IS_EXECUTION_POLICY,
                        DISALLOWED_EXECUTION_POLICIES) {
    // Test type trait helper
    REQUIRE_FALSE(execution::is_execution_policy_v<TestType>);
    REQUIRE_FALSE(execution::is_execution_policy_v<const TestType>);
    REQUIRE_FALSE(execution::is_execution_policy_v<const TestType&>);
    REQUIRE_FALSE(execution::is_execution_policy_v<const TestType&&>);
    REQUIRE_FALSE(execution::is_execution_policy_v<TestType&&>);

    // Test type trait
    REQUIRE_FALSE(execution::is_execution_policy<TestType>::value);
    REQUIRE_FALSE(execution::is_execution_policy<const TestType>::value);
    REQUIRE_FALSE(execution::is_execution_policy<const TestType&>::value);
    REQUIRE_FALSE(execution::is_execution_policy<const TestType&&>::value);
    REQUIRE_FALSE(execution::is_execution_policy<TestType&&>::value);
}

// Type trait should return true for parallel policies
TEMPLATE_LIST_TEST_CASE("Execution policy type trait parallel policies for parallel policies",
                        TAG.IS_PARALLEL_EXECUTION_POLICY, PARALLEL_EXECUTION_POLICIES) {
    // Test type trait helper
    REQUIRE(execution::is_parallel_execution_policy_v<TestType>);
    REQUIRE(execution::is_parallel_execution_policy_v<const TestType>);
    REQUIRE(execution::is_parallel_execution_policy_v<const TestType&>);
    REQUIRE(execution::is_parallel_execution_policy_v<const TestType&&>);
    REQUIRE(execution::is_parallel_execution_policy_v<TestType&&>);

    // Test type trait
    REQUIRE(execution::is_parallel_execution_policy<TestType>::value);
    REQUIRE(execution::is_parallel_execution_policy<const TestType>::value);
    REQUIRE(execution::is_parallel_execution_policy<const TestType&>::value);
    REQUIRE(execution::is_parallel_execution_policy<const TestType&&>::value);
    REQUIRE(execution::is_parallel_execution_policy<TestType&&>::value);
}

// Type trait should return true for parallel policies
TEMPLATE_LIST_TEST_CASE("Execution policy type trait parallel policies for sequential policies",
                        TAG.IS_PARALLEL_EXECUTION_POLICY, SEQUENTIAL_EXECUTION_POLICIES) {
    // Test type trait helper
    REQUIRE_FALSE(execution::is_parallel_execution_policy_v<TestType>);
    REQUIRE_FALSE(execution::is_parallel_execution_policy_v<const TestType>);
    REQUIRE_FALSE(execution::is_parallel_execution_policy_v<const TestType&>);
    REQUIRE_FALSE(execution::is_parallel_execution_policy_v<const TestType&&>);
    REQUIRE_FALSE(execution::is_parallel_execution_policy_v<TestType&&>);

    // Test type trait
    REQUIRE_FALSE(execution::is_parallel_execution_policy<TestType>::value);
    REQUIRE_FALSE(execution::is_parallel_execution_policy<const TestType>::value);
    REQUIRE_FALSE(execution::is_parallel_execution_policy<const TestType&>::value);
    REQUIRE_FALSE(execution::is_parallel_execution_policy<const TestType&&>::value);
    REQUIRE_FALSE(execution::is_parallel_execution_policy<TestType&&>::value);
}

TEST_CASE("Find a segment ID for an element") {
    // Num of pivots needs to be 2^k - 1
    SECTION("One pivot") {
        std::vector<int> pivots{3};
        auto comparator = std::less<int>();
        int num_comparisons = std::log2(pivots.size());

        int element = 1;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 0);
        element = 3;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 1);
        element = 4;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 1);
    }

    SECTION("Three pivots") {
        std::vector<int> pivots{3, 6, 9};
        auto comparator = std::less<int>();
        int num_comparisons = std::log2(pivots.size());

        int element = 1;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 0);
        element = 3;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 1);
        element = 7;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 2);
        element = 10;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 3);
    }

    SECTION("Three pivots: two are the same") {
        std::vector<int> pivots{3, 6, 6};
        auto comparator = std::less<int>();
        int num_comparisons = std::log2(pivots.size());

        int element = 1;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 0);
        element = 3;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 1);
        element = 7;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 3);
        // If two or more consequent pivots are the same, there is not any segment in between them
        // In this case, segment id 2 is never used
        element = 6;
        REQUIRE(helpers::_find_element_segment_id(num_comparisons, pivots.begin(), pivots.size(),
                                                  element, comparator)
                == 3);
    }
}

TEST_CASE("Test parallel partitioning") {
    std::vector<int> test_vector{10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    auto comp = std::less<int>();
    int num_pivots = 3;

    auto check_result = [](auto boundaries, auto test_vector, auto pivots) {
#pragma omp parallel for
        for (size_t i = 0; i < boundaries.size() - 1; ++i) {
            long start = 0;

            // Elements on left from pivot are <
            while (start < boundaries[i]) {
                REQUIRE(test_vector[start] < pivots[i]);
                ++start;
            }

            // Elements on right from pivot are >=
            while (start < boundaries.back()) {
                REQUIRE(test_vector[start] >= pivots[i]);
                ++start;
            }
        }
    };

    SECTION("Reverse order small") { test_vector = {10, 9, 8, 7, 6, 5, 4}; }

    SECTION("Reverse order small 2") { test_vector = {10, 9, 8, 7}; }

    SECTION("Reverse order") { test_vector = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1}; }

    SECTION("Already sorted") { test_vector = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; }

    SECTION("Three elements") { test_vector = {3, 2, 1}; }

    SECTION("Same elements") { test_vector = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; }

    SECTION("Three distinct values") { test_vector = {4, 4, 4, 1, 1, 1, 6, 6, 6}; }

    SECTION("Pivots same values") { test_vector = {4, 4, 4, 4, 1, 1, 4, 1, 6, 4, 6, 6}; }

    SECTION("Two pivots same values small") { test_vector = {4, 4, 4, 1, 4, 1, 6, 1, 6}; }

    SECTION("Two pivots same values") { test_vector = {4, 4, 4, 4, 1, 1, 4, 1, 6, 1, 6, 6}; }

    SECTION("Random numbers and sizes") {
        // Random length from 100 to 10000
        auto vector_length = GENERATE(take(100, random(3, 1000)));
        // Generate vector with random numbers
        test_vector = GENERATE(chunk(1000, take(1000, random(0, 1000))));
        test_vector.resize(vector_length);
    }

    SECTION("Small set of number values") {
        // Random length from 100 to 10000
        auto vector_length = GENERATE(take(100, random(3, 1000)));
        // Generate vector with random numbers
        test_vector = GENERATE(chunk(1000, take(1000, random(0, 10))));
        test_vector.resize(vector_length);
    }

    SECTION("7 pivots") {
        num_pivots = 7;
        // Random length from 100 to 10000
        auto vector_length = GENERATE(take(100, random(7, 1000)));
        // Generate vector with random numbers
        test_vector = GENERATE(chunk(1000, take(1000, random(0, 1000))));
        test_vector.resize(vector_length);
    }

    SECTION("15 pivots") {
        num_pivots = 15;
        // Random length from 100 to 10000
        auto vector_length = GENERATE(take(100, random(15, 1000)));
        // Generate vector with random numbers
        test_vector = GENERATE(chunk(1000, take(1000, random(0, 1000))));
        test_vector.resize(vector_length);
    }

    SECTION("255 pivots") {
        num_pivots = 255;
        // Random length from 100 to 10000
        auto vector_length = GENERATE(take(100, random(255, 10000)));
        // Generate vector with random numbers
        test_vector = GENERATE(chunk(10000, take(10000, random(0, 1000))));
        test_vector.resize(vector_length);
    }

    // If sample small enough, pivots are chosen like this without sampling
    auto pivots
        = helpers::_get_pivots(test_vector.begin(), 0, test_vector.size() - 1, num_pivots, comp);
    auto boundaries = mpqsort::impl::_par_multiway_partition(
        num_pivots, omp_get_max_threads(), test_vector.begin(), 0, test_vector.size() - 1, comp);

    check_result(boundaries, test_vector, pivots);
}

// Test if all sort prototypes can be called and instantiated (without policies)
TEST_CASE("Instantiation of a sort overloads without policy", TAG.SORT_ALL) {
    std::vector<int> test_vector{2, 1};
    auto first = test_vector.begin();
    auto last = test_vector.end();
    std::vector<int> test_vector_res{1, 2};

    SECTION("(iter, iter) template") { mpqsort::sort(first, last); }

    SECTION("(iter, iter, comp) template") { mpqsort::sort(first, last, std::less<int>()); }

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}

// Test if all sort prototypes can be called and instantiated (with policies)
TEMPLATE_LIST_TEST_CASE("Instantiation of a sort overloads with policy", TAG.SORT_ALL,
                        ALLOWED_EXECUTION_POLICIES) {
    std::vector<int> test_vector{2, 1};
    auto first = test_vector.begin();
    auto last = test_vector.end();
    std::vector<int> test_vector_res{1, 2};

    auto policy = TestType{};
    int cores = 4;

    SECTION("(policy, iter, iter template") { mpqsort::sort(policy, first, last); }

    SECTION("(policy, iter, iter, comp) template") {
        mpqsort::sort(policy, first, last, std::less<int>());
    }

    // Run only if parallel policy. Sequential policy + cores does not make sense and we do not have
    // a prototype for that
    if constexpr (mpqsort::execution::is_parallel_execution_policy_v<TestType>) {
        SECTION("(policy, cores, iter, iter) template") {
            mpqsort::sort(policy, cores, first, last);
        }

        SECTION("(policy, cores, iter, iter, comp) template") {
            mpqsort::sort(policy, cores, first, last, std::less<int>());
        }
    }

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}

TEST_CASE("HeapInsertion sort") {
    auto test_vector = std::vector<int>{7, 6, 5, 4, 8};
    auto comp = std::less<int>();
    mpqsort::helpers::_heap_insertion_sort(test_vector.begin(), 0, test_vector.size() - 1, comp);

    REQUIRE_THAT(test_vector, Catch::Equals(std::vector<int>{4, 5, 6, 7, 8}));
}

// Test correctness of implementations
TEMPLATE_LIST_TEST_CASE("Increasing vector sizes", TAG.SORT_ALL, ALLOWED_EXECUTION_POLICIES) {
    auto test_vector = std::vector<int>();
    auto test_vector_res = test_vector;

    SECTION("Size 0") {
        test_vector.resize(0);
        test_vector_res = test_vector;
    }

    SECTION("Size 1") {
        test_vector.resize(1);
        test_vector = {9};
        test_vector_res = test_vector;
    }

    SECTION("Size 2") {
        test_vector.resize(2);
        test_vector = {9, 8};
        test_vector_res = test_vector;
    }

    SECTION("Size 3") {
        test_vector.resize(3);
        test_vector = {9, 8, 7};
        test_vector_res = test_vector;
    }

    SECTION("Size 4") {
        test_vector.resize(4);
        test_vector = {9, 8, 7, 6};
        test_vector_res = test_vector;
    }

    SECTION("Size 5") {
        test_vector.resize(5);
        test_vector = {9, 8, 7, 6, 5};
        test_vector_res = test_vector;
    }

    SECTION("Size 6") {
        test_vector.resize(6);
        test_vector = {9, 8, 7, 6, 5, 4};
        test_vector_res = test_vector;
    }

    SECTION("Size 7") {
        test_vector.resize(7);
        test_vector = {9, 8, 7, 6, 5, 4, 3};
        test_vector_res = test_vector;
    }

    SECTION("Size 8") {
        test_vector.resize(8);
        test_vector = {9, 8, 7, 6, 5, 4, 3, 2};
        test_vector_res = test_vector;
    }

    SECTION("Size 9") {
        test_vector.resize(9);
        test_vector = {9, 8, 7, 6, 5, 4, 3, 2, 1};
        test_vector_res = test_vector;
    }

    SECTION("Size 10") {
        test_vector.resize(10);
        test_vector = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
        test_vector_res = test_vector;
    }

    // Sort array
    mpqsort::sort(TestType{}, test_vector.begin(), test_vector.end());
    // Compute correct result
    std::sort(test_vector_res.begin(), test_vector_res.end());

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}

// Test correctness of implementations
TEMPLATE_LIST_TEST_CASE("Sort already sorted vector", TAG.SORT_ALL, ALLOWED_EXECUTION_POLICIES) {
    // Length 1..10
    auto vector_length = GENERATE(range(1, 11));
    // Generate sequence from 0 to 10
    auto test_vector = GENERATE(chunk(10, range(0, 10)));
    test_vector.resize(vector_length);
    // Set expected result
    auto test_vector_res = test_vector;

    CAPTURE(vector_length);

    mpqsort::sort(TestType{}, test_vector.begin(), test_vector.end());

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}

TEMPLATE_LIST_TEST_CASE("Sort vector provided in reversed order", TAG.SORT_ALL,
                        ALLOWED_EXECUTION_POLICIES) {
    // Length 1..10
    auto vector_length = GENERATE(range(1, 11));
    // Generate sequence from 0 to 10
    auto test_vector = GENERATE(chunk(10, range(0, 10)));
    test_vector.resize(vector_length);
    // Set expected result
    auto test_vector_res = test_vector;
    // Reverse input vector
    std::reverse(test_vector.begin(), test_vector.end());

    CAPTURE(vector_length);

    mpqsort::sort(TestType{}, test_vector.begin(), test_vector.end());

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}

// Do not test when DEBUG mode activated, because it prints logging and _test_struct does not have
// operator << overloaded
#ifndef DEBUG
TEMPLATE_LIST_TEST_CASE("Sort custom type with defined comparator", TAG.SORT_ALL,
                        ALLOWED_EXECUTION_POLICIES) {
    // Declare custom type and comparators
    struct _test_struct {
        _test_struct() {}
        _test_struct(std::string&& _str, int _num) : str(_str), num(_num) {}
        std::string str;
        int num;

        bool operator==(const _test_struct& o) const { return str == o.str && num == o.num; }

        bool operator!=(const _test_struct& o) const { return !(*this == o); }
    };

    auto lambda_comparator = [](const auto& a, const auto& b) {
        if (a.str < b.str)
            return true;
        else if (a.str == b.str && a.num < b.num)
            return true;
        else {
            return false;
        }
    };

    struct functor_comparator {
        bool operator()(const _test_struct& a, const _test_struct& b) {
            if (a.str < b.str)
                return true;
            else if (a.str == b.str && a.num < b.num)
                return true;
            else {
                return false;
            }
        }
    };

    std::vector<_test_struct> test_vector{_test_struct("z", 2), _test_struct("a", 2),
                                          _test_struct("g", 6), _test_struct("z", 2),
                                          _test_struct("z", 1)};

    auto test_vector_res = test_vector;

    SECTION("Use lambda comparator") {
        std::sort(test_vector_res.begin(), test_vector_res.end(), lambda_comparator);
        mpqsort::sort(TestType{}, test_vector.begin(), test_vector.end(), lambda_comparator);

        REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
    }

    SECTION("Use functor comparator") {
        std::sort(test_vector_res.begin(), test_vector_res.end(), functor_comparator());
        mpqsort::sort(TestType{}, test_vector.begin(), test_vector.end(), functor_comparator());

        REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
    }
}
#endif

TEST_CASE("Specific problematic input", TAG.SORT_ALL) {
    auto test_vector = std::vector<int>{6, 6, 7, 6, 5, 2, 2};
    auto test_vector_res = test_vector;

    SECTION("Seq three way") {
        mpqsort::sort(mpqsort::execution::seq_three_way, test_vector.begin(), test_vector.end());
    }

    SECTION("Seq four way") {
        mpqsort::sort(mpqsort::execution::seq_four_way, test_vector.begin(), test_vector.end());
    }

    std::sort(test_vector_res.begin(), test_vector_res.end());

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}

TEMPLATE_LIST_TEST_CASE("Sort vectors, small and different sizes and random numbers.", TAG.SORT_ALL,
                        ALLOWED_EXECUTION_POLICIES) {
    // Random length from 100 to 10000
    auto vector_length = GENERATE(take(100, random(0, 10)));
    // Generate vector with random numbers
    auto test_vector = GENERATE(chunk(10000, take(10000, random(0, 10))));
    test_vector.resize(vector_length);
    // Compute expected result
    auto test_vector_res = test_vector;
    std::sort(test_vector_res.begin(), test_vector_res.end());

    CAPTURE(vector_length, test_vector, Catch::rngSeed());

    mpqsort::sort(TestType{}, test_vector.begin(), test_vector.end());

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}

// Random data, random lengths. Final test trying to catch missed bugs
TEMPLATE_LIST_TEST_CASE(
    "Sort vectors, different sizes and random numbers. Last chance to catch errors.", TAG.SORT_ALL,
    ALLOWED_EXECUTION_POLICIES) {
    // Random length from 100 to 10000
    auto vector_length = GENERATE(take(100, random(100, 10000)));
    // Generate vector with random numbers
    auto test_vector = GENERATE(chunk(10000, take(10000, random(0, 1000000))));
    test_vector.resize(vector_length);
    // Compute expected result
    auto test_vector_res = test_vector;
    std::sort(test_vector_res.begin(), test_vector_res.end());

    CAPTURE(vector_length, Catch::rngSeed());

    mpqsort::sort(TestType{}, test_vector.begin(), test_vector.end());

    REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
}
