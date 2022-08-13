#include <mpqsort/mpqsort.h>
#include <mpqsort/version.h>

#include <catch2/catch.hpp>
#include <string>

TEST_CASE("Greeter") {
    using namespace mpqsort;

    Greeter greeter("Tests");

    CHECK(greeter.greet(LanguageCode::EN) == "Hello, Tests!");
    CHECK(greeter.greet(LanguageCode::DE) == "Hallo Tests!");
    CHECK(greeter.greet(LanguageCode::ES) == "¡Hola Tests!");
    CHECK(greeter.greet(LanguageCode::FR) == "Bonjour Tests!");
}

TEST_CASE("Greeter version") {
    static_assert(std::string_view(MPQSORT_VERSION) == std::string_view("1.0"));
    CHECK(std::string(MPQSORT_VERSION) == std::string("1.0"));
}
