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

template <typename T, long Size = -1, int From = -1, int To = -1> struct VectorFixture
    : public benchmark::Fixture {
    VectorFixture()
        : from(From == -1 ? std::numeric_limits<T>::min() : From),
          to(To == -1 ? std::numeric_limits<T>::max() : To) {}

    void TearDown(const benchmark::State& state) override { this->DeallocateVector(); }

    void AllocateVector() {
        // If size not defined, use half of the memory as a default
        this->vec.resize(Size == -1 ? this->VectorSizeToFillHalfMemory() : Size);
    }

    void DeallocateVector() {
        this->vec.clear();
        this->vec.shrink_to_fit();
    }

    // Calculate how many elements to generate to fill half of the system memory
    auto VectorSizeToFillHalfMemory() const {
#ifdef TESTING
        return 1000000;
#else
        auto pages = sysconf(_SC_PHYS_PAGES);
        auto page_size = sysconf(_SC_PAGE_SIZE);
        // Memory in bytes
        auto available_memory = pages * page_size;
        auto allocate_memory = available_memory / (sizeof(T) * 2);

        return allocate_memory;
#endif
    }

    template <typename Dist> void InitVector(Dist distribution) {
        for (auto& el : vec) el = distribution();
    }

    std::vector<T> vec;
    T const from, to;
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct RandomVectorFixture : public VectorFixture<T, Size, From, To> {
    void SetUp(const benchmark::State& state) override {
        this->AllocateVector();
        FillVectorRandom();
    }

    void FillVectorRandom() {
        // Seed with a same value so that sort same sequences across runs
        std::mt19937 en(Seed);

        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> uniform_dist(this->from, this->to);
            this->InitVector([&]() { return uniform_dist(en); });
        } else {
            std::uniform_real_distribution<T> uniform_dist(this->from, this->to);
            this->InitVector([&]() { return uniform_dist(en); });
        }
    }
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct SortedVectorFixture : public RandomVectorFixture<T, Size, From, To, Seed> {
    void SetUp(const benchmark::State& state) override {
        this->AllocateVector();
        this->FillVectorRandom();
        std::sort(std::execution::par, this->vec.begin(), this->vec.end());
    }
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct ReverseOrderVectorFixture : public RandomVectorFixture<T, Size, From, To, Seed> {
    void SetUp(const benchmark::State& state) override {
        this->AllocateVector();
        this->FillVectorRandom();
        std::sort(std::execution::par, this->vec.begin(), this->vec.end(), std::greater<T>());
    }
};

// Stringify macro
#define str(X) #X

// Registers benchmarks with different vectors (sorted, random, reverse ), data types, sizes and
// from, to
#define register_bench_template(name, bench, size, from, to) \
    name(Random, bench, int, size, from, to);                \
    name(Random, bench, short, size, from, to);              \
    name(Random, bench, double, size, from, to);             \
    name(Sorted, bench, int, size, from, to);                \
    name(Sorted, bench, short, size, from, to);              \
    name(Sorted, bench, double, size, from, to);             \
    name(ReverseOrder, bench, int, size, from, to);          \
    name(ReverseOrder, bench, short, size, from, to);        \
    name(ReverseOrder, bench, double, size, from, to);

// Registers benchmarks with default values for different types (size, from, to stays default)
#define default default
#define register_bench_default(name) register_bench_template(name, default, -1, -1, -1);

// Registers benchmark with random values of type int
#define register_bench_int_random(name, bench, size, from, to) \
    name(Random, bench, int, size, from, to);

// Registers benchmarks with small difference between from and to => generates lot of same elements
#define small_values_range small_values_range
#define register_bench_small_values_range(name)                                          \
    register_bench_int_random(name, small_values_range##_##1, -1, 0, 0);                 \
    register_bench_int_random(name, small_values_range##_##3, -1, 0, 2);                 \
    register_bench_int_random(name, small_values_range##_##10, -1, 0, 9);                \
    register_bench_int_random(name, small_values_range##_##100, -1, 0, 99);              \
    register_bench_int_random(name, small_values_range##_##1000, -1, 0, int(1e3 - 1));   \
    register_bench_int_random(name, small_values_range##_##100000, -1, 0, int(1e5 - 1)); \
    register_bench_int_random(name, small_values_range##_##100000000, -1, 0, int(1e8 - 1));

#define small_sizes small_sizes
#define register_bench_small_sizes(name)                                         \
    register_bench_int_random(name, small_sizes##_##100000, 100000, -1, -1);     \
    register_bench_int_random(name, small_sizes##_##500000, 500000, -1, -1);     \
    register_bench_int_random(name, small_sizes##_##1000000, 1000000, -1, -1);   \
    register_bench_int_random(name, small_sizes##_##20000000, 20000000, -1, -1); \
    register_bench_int_random(name, small_sizes##_##50000000, 50000000, -1, -1); \
    register_bench_int_random(name, small_sizes##_##100000000, 100000000, -1, -1);

// Run std sort benchmarks
#define std_sort(dataType, bench, type, size, from, to)                                          \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_std_sort_##dataType##_##type##_##bench, type, size, from, to) \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            std::sort(vec.begin(), vec.end());                                                   \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_std_sort_##dataType##_##type##_##bench)     \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_std_sort_##dataType##_##type##_##bench));

register_bench_default(std_sort);
register_bench_small_sizes(std_sort);

// Run std parallel sort benchmarks
#define std_parallel_sort(dataType, bench, type, size, from, to)                                 \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_par_std_sort_##dataType##_##type##_##bench, type, size, from, \
                                to)                                                              \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            std::sort(std::execution::par, vec.begin(), vec.end());                              \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_par_std_sort_##dataType##_##type##_##bench) \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_par_std_sort_##dataType##_##type##_##bench));

register_bench_default(std_parallel_sort);
register_bench_small_sizes(std_parallel_sort);

// Run mpqsort benchmarks
#define mpqsort_sort(dataType, bench, type, size, from, to)                                      \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_mpqsort_sort_##dataType##_##type##_##bench, type, size, from, \
                                to)                                                              \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            mpqsort::sort(mpqsort::execution::par_multi_way, vec.begin(), vec.end());            \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_mpqsort_sort_##dataType##_##type##_##bench) \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_mpqsort_sort_##dataType##_##type##_##bench));

register_bench_default(mpqsort_sort);
register_bench_small_sizes(mpqsort_sort);
register_bench_small_values_range(mpqsort_sort);

// Run gnu qs sort benchmarks
// In-place parallel qsort
#define gnu_qs_sort(dataType, bench, type, size, from, to)                                      \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                        \
                                BM_gnu_qs_sort_##dataType##_##type##_##bench, type, size, from, \
                                to)                                                             \
    (benchmark::State & state) {                                                                \
        for (auto _ : state) {                                                                  \
            __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::quicksort_tag());      \
        }                                                                                       \
    }                                                                                           \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_gnu_qs_sort_##dataType##_##type##_##bench) \
        ->MeasureProcessCPUTime()                                                               \
        ->UseRealTime()                                                                         \
        ->Name(str(BM_gnu_qs_sort_##dataType##_##type##_##bench));

register_bench_default(gnu_qs_sort);
register_bench_small_sizes(gnu_qs_sort);

// Run gnu_bqs_sort benchmarks
// In-place balanced parallel qsort
#define gnu_bqs_sort(dataType, bench, type, size, from, to)                                      \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_gnu_bqs_sort_##dataType##_##type##_##bench, type, size, from, \
                                to)                                                              \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            __gnu_parallel::sort(vec.begin(), vec.end(),                                         \
                                 __gnu_parallel::balanced_quicksort_tag());                      \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_gnu_bqs_sort_##dataType##_##type##_##bench) \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_gnu_bqs_sort_##dataType##_##type##_##bench));

register_bench_default(gnu_bqs_sort);
register_bench_small_sizes(gnu_bqs_sort);

// Run gnu mwms sort benchmarks
// Out-of-place parallel mergesort
#define gnu_mwms_sort(dataType, bench, type, size, from, to)                                      \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                          \
                                BM_gnu_mwms_sort_##dataType##_##type##_##bench, type, size, from, \
                                to)                                                               \
    (benchmark::State & state) {                                                                  \
        for (auto _ : state) {                                                                    \
            __gnu_parallel::sort(vec.begin(), vec.end(),                                          \
                                 __gnu_parallel::multiway_mergesort_tag());                       \
        }                                                                                         \
    }                                                                                             \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_gnu_mwms_sort_##dataType##_##type##_##bench) \
        ->MeasureProcessCPUTime()                                                                 \
        ->UseRealTime()                                                                           \
        ->Name(str(BM_gnu_mwms_sort_##dataType##_##type##_##bench));

register_bench_default(gnu_mwms_sort);
register_bench_small_sizes(gnu_mwms_sort);

// Run tbb sort benchmarks
#define tbb_sort(dataType, bench, type, size, from, to)                                          \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_tbb_sort_##dataType##_##type##_##bench, type, size, from, to) \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            tbb::parallel_sort(vec.begin(), vec.end());                                          \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_tbb_sort_##dataType##_##type##_##bench)     \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_tbb_sort_##dataType##_##type##_##bench));

register_bench_default(tbb_sort);
register_bench_small_sizes(tbb_sort);

// Run nvidia thrust sort benchmarks
#define nvidia_thrust_sort(dataType, bench, type, size, from, to)                                \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_nvidia_thrust_sort_##dataType##_##type##_##bench, type, size, \
                                from, to)                                                        \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            thrust::sort(vec.begin(), vec.end());                                                \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture,                                                \
                         BM_nvidia_thrust_sort_##dataType##_##type##_##bench)                    \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_nvidia_thrust_sort_##dataType##_##type##_##bench));

register_bench_default(nvidia_thrust_sort);
register_bench_small_sizes(nvidia_thrust_sort);