#include <mpqsort.h>
#include <mpqsort/version.h>

#include <catch2/catch.hpp>
#include <string>
#include <vector>

TEST_CASE("Sort") {
    std::vector<int> v{1, 2, 3, 4};
    // std::vector b{};

    // REQUIRE_THROWS(sort())
}

TEST_CASE("Sort version") {
    static_assert(std::string_view(MPQSORT_VERSION) == std::string_view("1.0"));
    CHECK(std::string(MPQSORT_VERSION) == std::string("1.0"));
}
