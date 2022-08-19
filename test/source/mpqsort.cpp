#include <mpqsort/mpqsort.h>
#include <mpqsort/version.h>

#include <catch2/catch.hpp>
#include <string>
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
TEST_CASE("Instantiation of a sort overloads without policy") {
    std::vector<int> test_vector{1, 2, 3};
    auto first = test_vector.begin();
    auto last = test_vector.end();

    mpqsort::sort(first, last);
}

TEST_CASE("Sort") {
    std::vector<int> v{1, 2, 3, 4};
    // std::vector b{};

    // REQUIRE_THROWS(sort())
}

TEST_CASE("Sort version") {
    static_assert(std::string_view(MPQSORT_VERSION) == std::string_view("1.0"));
    CHECK(std::string(MPQSORT_VERSION) == std::string("1.0"));
}
