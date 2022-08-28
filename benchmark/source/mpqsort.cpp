#include <benchmark/benchmark.h>
#include <mpqsort/mpqsort.h>
#include <mpqsort/version.h>
#include <tbb/tbb.h>
#include <thrust/sort.h>

#include <algorithm>
#include <parallel/algorithm>
#include <string>
#include <vector>

// TODO just for testing purposes
std::vector<int> vec;

static void BM_std_sort(benchmark::State& state) {
    for (auto _ : state) {
        std::sort(vec.begin(), vec.end());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_std_sort);

static void BM_par_std_sort(benchmark::State& state) {
    for (auto _ : state) {
        std::sort(std::execution::par, vec.begin(), vec.end());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_par_std_sort);

static void BM_mpqsort_sort(benchmark::State& state) {
    for (auto _ : state) {
        mpqsort::sort(mpqsort::execution::par_multi_way, vec.begin(), vec.end());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_mpqsort_sort);

static void BM_gnu_qs_sort(benchmark::State& state) {
    for (auto _ : state) {
        // In-place parallel qsort
        __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::quicksort_tag());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_gnu_qs_sort);

static void BM_gnu_bqs_sort(benchmark::State& state) {
    for (auto _ : state) {
        // In-place balanced parallel qsort
        __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::balanced_quicksort_tag());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_gnu_bqs_sort);

static void BM_gnu_mwms_sort(benchmark::State& state) {
    for (auto _ : state) {
        // Out-of-place parallel mergesort
        __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::multiway_mergesort_tag());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_gnu_mwms_sort);

static void BM_tbb_sort(benchmark::State& state) {
    for (auto _ : state) {
        // Out-of-place parallel mergesort
        tbb::parallel_sort(vec.begin(), vec.end());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_tbb_sort);

static void BM_nvidia_thrust_sort(benchmark::State& state) {
    for (auto _ : state) {
        // Out-of-place parallel mergesort
        thrust::sort(vec.begin(), vec.end());
    }
}
// Register the function as a benchmark
BENCHMARK(BM_nvidia_thrust_sort);