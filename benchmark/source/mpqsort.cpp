#include <benchmark/benchmark.h>
#include <mpqsort.h>
#include <mpqsort/version.h>

#include <string>

static void BM_StringCreation(benchmark::State& state) {
    for (auto _ : state) std::string empty_string;
}
// Register the function as a benchmark
BENCHMARK(BM_StringCreation);

// Define another benchmark
static void BM_StringCopy(benchmark::State& state) {
    using namespace mpqsort;

    Greeter greeter("Tests");

    for (auto _ : state) {
        greeter.greet(LanguageCode::EN);
        greeter.greet(LanguageCode::DE);
        greeter.greet(LanguageCode::ES);
        greeter.greet(LanguageCode::FR);
    }
}
BENCHMARK(BM_StringCopy);
