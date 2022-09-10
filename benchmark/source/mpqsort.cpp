#include <benchmark/benchmark.h>
#include <mpqsort/mpqsort.h>
#include <mpqsort/version.h>
#include <tbb/tbb.h>
#include <thrust/sort.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <parallel/algorithm>
#include <random>
#include <string>
#include <vector>

// TODO remove when released
#define TESTING

template <typename T> struct VectorFixture : public benchmark::Fixture {
    VectorFixture() : from(std::numeric_limits<T>::min()), to(std::numeric_limits<T>::max()) {
        this->vec.resize(this->VectorSizeToFillHalfMemory());
    }

    // Calculate how many elements to generate to fill half of the system memory
    auto VectorSizeToFillHalfMemory() const {
#ifdef TESTING
        return 1000000;
#else
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        // Memory in bytes
        long available_memory = pages * page_size;

        return available_memory / 2 / sizeof(T);
#endif
    }

    template <typename Dist> void InitVector(Dist distribution) {
        for (auto& el : vec) el = distribution();
    }

    std::vector<T> vec;
    T const from, to;
};

template <typename T, int seed = 0> struct RandomVectorFixture : public VectorFixture<T> {
    void SetUp(const benchmark::State& state) override { FillVectorRandom(); }

    void FillVectorRandom() {
        // Seed with a same value so that sort same sequences across runs
        std::mt19937 en(seed);

        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> uniform_dist(this->from, this->to);
            this->InitVector([&]() { return uniform_dist(en); });
        } else {
            std::uniform_real_distribution<T> uniform_dist(this->from, this->to);
            this->InitVector([&]() { return uniform_dist(en); });
        }
    }
};

template <typename T, int seed = 0> struct SortedVectorFixture
    : public RandomVectorFixture<T, seed> {
    SortedVectorFixture() {
        // Fill and prepare only once
        this->FillVectorRandom();
        std::sort(std::execution::par, this->vec.begin(), this->vec.end());
    }

    void SetUp(const benchmark::State& state) override {
        // Not needed, vector is generated and sorted once because its already sorted and sorting
        // functions should not change that
    }
};

template <typename T, int seed = 0> struct ReverseOrderVectorFixture
    : public RandomVectorFixture<T, seed> {
    ReverseOrderVectorFixture() {
        this->FillVectorRandom();
        std::sort(std::execution::par, this->vec.begin(), this->vec.end());
    }
    void SetUp(const benchmark::State& state) override {
        // Vector is sorted at this point, so reverse it
        // SetUp is called between benchmarks and vec is sorted by benchmarked algorithm
        std::reverse(std::execution::par, this->vec.begin(), this->vec.end());
    }
};

// Registers benchmarks with different vectors (sorted, random, reverse ) and data types
#define register_bench(name)   \
    name(Random, int);         \
    name(Random, short);       \
    name(Random, double);      \
    name(Sorted, int);         \
    name(Sorted, short);       \
    name(Sorted, double);      \
    name(ReverseOrder, int);   \
    name(ReverseOrder, short); \
    name(ReverseOrder, double);

// Run std sort benchmarks
#define std_sort(dataType, type)                                                         \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_std_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                         \
        for (auto _ : state) {                                                           \
            std::sort(vec.begin(), vec.end());                                           \
        }                                                                                \
    }

register_bench(std_sort);

// Run std parallel sort benchmarks
#define std_parallel_sort(dataType, type)                                                    \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_par_std_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                             \
        for (auto _ : state) {                                                               \
            std::sort(std::execution::par, vec.begin(), vec.end());                          \
        }                                                                                    \
    }

register_bench(std_parallel_sort);

// Run mpqsort benchmarks
#define mpqsort_sort(dataType, type)                                                         \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_mpqsort_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                             \
        for (auto _ : state) {                                                               \
            mpqsort::sort(mpqsort::execution::par_multi_way, vec.begin(), vec.end());        \
        }                                                                                    \
    }

register_bench(mpqsort_sort);

// Run gnu qs sort benchmarks
// In-place parallel qsort
#define gnu_qs_sort(dataType, type)                                                         \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_gnu_qs_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                            \
        for (auto _ : state) {                                                              \
            __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::quicksort_tag());  \
        }                                                                                   \
    }

register_bench(gnu_qs_sort);

// Run gnu_bqs_sort benchmarks
// In-place balanced parallel qsort
#define gnu_bqs_sort(dataType, type)                                                         \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_gnu_bqs_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                             \
        for (auto _ : state) {                                                               \
            __gnu_parallel::sort(vec.begin(), vec.end(),                                     \
                                 __gnu_parallel::balanced_quicksort_tag());                  \
        }                                                                                    \
    }

register_bench(gnu_bqs_sort);

// Run gnu mwms sort benchmarks
// Out-of-place parallel mergesort
#define gnu_mwms_sort(dataType, type)                                                         \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_gnu_mwms_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                              \
        for (auto _ : state) {                                                                \
            __gnu_parallel::sort(vec.begin(), vec.end(),                                      \
                                 __gnu_parallel::multiway_mergesort_tag());                   \
        }                                                                                     \
    }

register_bench(gnu_mwms_sort);

// Run tbb sort benchmarks
#define tbb_sort(dataType, type)                                                         \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_tbb_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                         \
        for (auto _ : state) {                                                           \
            tbb::parallel_sort(vec.begin(), vec.end());                                  \
        }                                                                                \
    }

register_bench(tbb_sort);

// Run nvidia thrust sort benchmarks
#define nvidia_thrust_sort(dataType, type)                                                         \
    BENCHMARK_TEMPLATE_F(dataType##VectorFixture, BM_nvidia_thrust_sort_##dataType##_##type, type) \
    (benchmark::State & state) {                                                                   \
        for (auto _ : state) {                                                                     \
            thrust::sort(vec.begin(), vec.end());                                                  \
        }                                                                                          \
    }

register_bench(nvidia_thrust_sort);