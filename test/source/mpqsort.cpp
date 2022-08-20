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
    const std::string IS_EXECUTION_POLICY = "[is_execution_policy_v]";
    const std::string SEQ = "[seq_sort]";
    const std::string SEQ_TWO_WAY = "[seq_two_way_sort]";
    const std::string SEQ_MULTIWAY = "[seq_multiway_sort]";
    const std::string SEQ_ALL = SEQ + SEQ_TWO_WAY + SEQ_MULTIWAY;
    const std::string PAR = "[par_sort]";
    const std::string PAR_TWO_WAY = "[par_two_way_sort]";
    const std::string PAR_MULTIWAY = "[par_multiway_sort]";
    const std::string PAR_ALL = PAR + PAR_TWO_WAY + PAR_MULTIWAY;
    const std::string SORT_ALL = SEQ_ALL + PAR_ALL;
} TAG;

// Combine tags and create alias for catch2 (supports only C strings)
CATCH_REGISTER_TAG_ALIAS("[@seq_all]", TAG.SEQ_ALL.c_str());
CATCH_REGISTER_TAG_ALIAS("[@par_all]", TAG.PAR_ALL.c_str());
CATCH_REGISTER_TAG_ALIAS("[@sort_all]", TAG.SORT_ALL.c_str());

// Declare allowed execution policies
using ALLOWED_EXECUTION_POLICIES
    = std::tuple<decltype(execution::par_two_way), decltype(execution::par_multi_way),
                 decltype(execution::par), decltype(execution::seq_two_way),
                 decltype(execution::seq_multi_way), decltype(execution::seq)>;

// Some disallowed execution policies
using DISALLOWED_EXECUTION_POLICIES
    = std::tuple<int, bool, char, float, std::string, decltype(std::execution::par_unseq)>;

// TEST CASES

TEST_CASE("Sort version") {
    static_assert(std::string_view(MPQSORT_VERSION) == std::string_view("1.0"));
    CHECK(std::string(MPQSORT_VERSION) == std::string("1.0"));
}

// Type trait should return true for all allowed policies
TEMPLATE_LIST_TEST_CASE("Execution policy type trait allowed policies", TAG.IS_EXECUTION_POLICY,
                        ALLOWED_EXECUTION_POLICIES) {
    REQUIRE(std::is_swappable_v<int>);
    CHECK(execution::is_execution_policy_v<TestType>);
}

// Type trait should return false for all disallowed policies
TEMPLATE_LIST_TEST_CASE("Execution policy type trait disallowed policies", TAG.IS_EXECUTION_POLICY,
                        DISALLOWED_EXECUTION_POLICIES) {
    CHECK_FALSE(execution::is_execution_policy_v<TestType>);
}

// Test if all sort prototypes can be called and instantiated (without policies)
TEST_CASE("Instantiation of a sort overloads without policy", TAG.SORT_ALL) {
    std::vector<int> test_vector{0, 0};
    auto first = test_vector.begin();
    auto last = test_vector.end();

    mpqsort::sort(first, last);
    mpqsort::sort(first, last, std::greater<int>());
}

// Test if all sort prototypes can be called and instantiated (with policies)
TEMPLATE_LIST_TEST_CASE("Instantiation of a sort overloads with policy", TAG.SORT_ALL,
                        ALLOWED_EXECUTION_POLICIES) {
    std::vector<int> test_vector{0, 0};
    auto first = test_vector.begin();
    auto last = test_vector.end();

    auto policy = TestType{};
    int cores = 4;

    mpqsort::sort(policy, first, last);
    mpqsort::sort(policy, first, last, std::greater<int>());

    // Run only if parallel policy. Sequential policy + cores does not make sense and we do not have a prototype for that
    if constexpr (mpqsort::execution::is_parallel_execution_policy_v<TestType>)
    {
        mpqsort::sort(policy, cores, first, last);
        mpqsort::sort(policy, cores, first, last, std::greater<int>());
    }
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

TEMPLATE_LIST_TEST_CASE("Sort custom type with defined comparator", TAG.SORT_ALL,
                        ALLOWED_EXECUTION_POLICIES) {
    // Declare custom type and comparators
    struct _TEST_STRUCT {
        _TEST_STRUCT(std::string&& _str, int _num) : str(_str), num(_num) {}
        std::string str;
        int num;

        bool operator==(const _TEST_STRUCT& o) const { return str == o.str && num == o.num; }

        bool operator!=(const _TEST_STRUCT& o) const { return !(*this == o); }
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
        bool operator()(const _TEST_STRUCT& a, const _TEST_STRUCT& b) {
            if (a.str < b.str)
                return true;
            else if (a.str == b.str && a.num < b.num)
                return true;
            else {
                return false;
            }
        }
    };

    std::vector<_TEST_STRUCT> test_vector{_TEST_STRUCT("z", 2), _TEST_STRUCT("a", 2),
                                          _TEST_STRUCT("g", 6), _TEST_STRUCT("z", 2),
                                          _TEST_STRUCT("z", 1)};

    auto test_vector_res = test_vector;

    SECTION("Use lambda comparator") {
        std::sort(test_vector_res.begin(), test_vector_res.end(), lambda_comparator);
        mpqsort::sort(test_vector.begin(), test_vector.end(), lambda_comparator);

        REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
    }

    SECTION("Use functor comparator") {
        std::sort(test_vector_res.begin(), test_vector_res.end(), functor_comparator());
        mpqsort::sort(test_vector.begin(), test_vector.end(), functor_comparator());

        REQUIRE_THAT(test_vector, Catch::Equals(test_vector_res));
    }
}
// Random data, random lengths. Final test trying to catch missed bugs
TEMPLATE_LIST_TEST_CASE(
    "Sort vectors, different sizes and random numbers. Last chance to catch errors.", TAG.SORT_ALL,
    ALLOWED_EXECUTION_POLICIES) {
    // Random length from 100 to 10000
    auto vector_length = GENERATE(take(10, random(100, 10000)));
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
