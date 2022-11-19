#include <benchmark/benchmark.h>
#include <mpqsort/mpqsort.h>
#include <mpqsort/version.h>
#include <omp.h>
#include <tbb/tbb.h>
#include <thrust/sort.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <parallel/algorithm>
#include <random>
#include <string>
#include <vector>

// TODO: remove when released
#define TESTING

// Possible vector types
using VECTOR_TYPES = std::tuple<int, short, double>;

template <typename T, long Size = -1, int From = -1, int To = -1> struct VectorFixture
    : public benchmark::Fixture {
    VectorFixture()
        : from(From == -1 ? std::numeric_limits<T>::min() : From),
          to(To == -1 ? std::numeric_limits<T>::max() : To) {
        // GNU implementations wont' set it by itself...
        omp_set_nested(1);
    }

    void AllocateVector() {
        // If size not defined, use half of the memory as a default
        this->vec.resize(Size == -1 ? this->VectorSizeToFillHalfMemory() : Size);
    }

    virtual void Prepare() {
        throw std::logic_error("Prepare function should be overriden in child class!");
    }

    void DeallocateVector() { std::vector<T>().swap(this->vec); }

    void Destroy() { this->DeallocateVector(); }

    // Calculate how many elements to generate to fill half of the system memory
    auto VectorSizeToFillHalfMemory() const {
#ifdef TESTING
        // return 100000000;
        return 10000000;
        // return 8 * 1024 * 1024 * (1024 / sizeof(double));
#else
        auto pages = sysconf(_SC_PHYS_PAGES);
        auto page_size = sysconf(_SC_PAGE_SIZE);

        // Memory in bytes
        auto available_memory = pages * page_size;
        auto allocate_memory = available_memory / 4;

        // Get the largest type and compute number of elements based on that
        // Result: occupied memory is different but number of elements the same
        std::vector<int> sizes;
        std::apply([&](auto&&... args) { ((sizes.emplace_back(sizeof(args))), ...); },
                   VECTOR_TYPES());
        auto number_of_elements = allocate_memory / *std::max_element(sizes.begin(), sizes.end());

        return number_of_elements;
#endif
    }

    std::vector<T> vec;
    T const from, to;
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct RandomVectorFixture : public VectorFixture<T, Size, From, To> {
    void Prepare() {
        this->DeallocateVector();
        this->AllocateVector();
        FillVectorRandom();
    }

    void FillVectorRandom() {
        // Seed with a same value so that sort same sequences across runs
        std::vector<std::mt19937> ens;
        auto seed = Seed;
        for (int i = 0; i < omp_get_max_threads(); ++i) ens.emplace_back(std::mt19937(seed++));

        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> uniform_dist(this->from, this->to);

#pragma omp parallel
            {
                auto tid = omp_get_thread_num();
#pragma omp for schedule(static)
                for (size_t i = 0; i < this->vec.size(); ++i) {
                    this->vec[i] = uniform_dist(ens[tid]);
                }
            }
        } else {
            std::uniform_real_distribution<T> uniform_dist(this->from, this->to);

#pragma omp parallel
            {
                auto tid = omp_get_thread_num();
#pragma omp for schedule(static)
                for (size_t i = 0; i < this->vec.size(); ++i) {
                    this->vec[i] = uniform_dist(ens[tid]);
                }
            }
        }
    }
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct SortedVectorFixture : public RandomVectorFixture<T, Size, From, To, Seed> {
    void Prepare() {
        this->DeallocateVector();
        this->AllocateVector();
        this->FillVectorRandom();
        __gnu_parallel::sort(this->vec.begin(), this->vec.end(),
                             __gnu_parallel::balanced_quicksort_tag());
    }
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct ReverseOrderVectorFixture : public RandomVectorFixture<T, Size, From, To, Seed> {
    void Prepare() {
        this->DeallocateVector();
        this->AllocateVector();
        this->FillVectorRandom();
        __gnu_parallel::sort(this->vec.begin(), this->vec.end(), std::greater<T>(),
                             __gnu_parallel::balanced_quicksort_tag());
    }
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct OrganpipeOrderVectorFixture : public RandomVectorFixture<T, Size, From, To, Seed> {
    void Prepare() {
        this->DeallocateVector();
        this->AllocateVector();
        this->FillVectorRandom();

        auto end_middle = this->vec.size() / 2 + this->vec.size() & 1;

        __gnu_parallel::sort(this->vec.begin(), this->vec.begin() + end_middle,
                             __gnu_parallel::balanced_quicksort_tag());
        __gnu_parallel::sort(this->vec.begin() + end_middle, this->vec.end(), std::greater<T>(),
                             __gnu_parallel::balanced_quicksort_tag());
    }
};

template <typename T, long Size = -1, int From = -1, int To = -1, int Seed = 0>
struct RotatedOrderVectorFixture : public RandomVectorFixture<T, Size, From, To> {
    // 1, 2, ..., n - 1, 0
    void Prepare() {
        this->DeallocateVector();
        this->AllocateVector();
        this->FillVectorRandom();
        __gnu_parallel::sort(this->vec.begin(), this->vec.end(),
                             __gnu_parallel::balanced_quicksort_tag());
        std::rotate(this->vec.begin(), this->vec.begin() + 1, this->vec.end());
    }
};

// TODO: Teeth like structure (for p chunks generate p sorted parts)

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
    name(ReverseOrder, bench, double, size, from, to);       \
    name(OrganpipeOrder, bench, int, size, from, to);        \
    name(OrganpipeOrder, bench, short, size, from, to);      \
    name(OrganpipeOrder, bench, double, size, from, to);     \
    name(RotatedOrder, bench, int, size, from, to);          \
    name(RotatedOrder, bench, short, size, from, to);        \
    name(RotatedOrder, bench, double, size, from, to);

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

// Helps to find out threshold when to switch to sequential algorithm
#define small_sizes_threshold small_size_threshold
#define register_bench_small_size_threshold(name)                                         \
    register_bench_int_random(name, small_size_threshold##_##50000, 50000, -1, -1);       \
    register_bench_int_random(name, small_size_threshold##_##100000, 100000, -1, -1);     \
    register_bench_int_random(name, small_size_threshold##_##200000, 200000, -1, -1);     \
    register_bench_int_random(name, small_size_threshold##_##300000, 300000, -1, -1);     \
    register_bench_int_random(name, small_size_threshold##_##400000, 400000, -1, -1);     \
    register_bench_int_random(name, small_size_threshold##_##500000, 500000, -1, -1);     \
    register_bench_int_random(name, small_size_threshold##_##1000000, 1000000, -1, -1);   \
    register_bench_int_random(name, small_size_threshold##_##10000000, 10000000, -1, -1); \
    register_bench_int_random(name, small_size_threshold##_##100000000, 100000000, -1, -1);

// Run std sort benchmarks
#define std_sort(dataType, bench, type, size, from, to)                                          \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_std_sort_##dataType##_##type##_##bench, type, size, from, to) \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            state.PauseTiming();                                                                 \
            Prepare();                                                                           \
            state.ResumeTiming();                                                                \
            std::sort(vec.begin(), vec.end());                                                   \
            state.PauseTiming();                                                                 \
            Destroy();                                                                           \
            state.ResumeTiming();                                                                \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_std_sort_##dataType##_##type##_##bench)     \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_std_sort_##dataType##_##type##_##bench));

register_bench_default(std_sort);
register_bench_small_sizes(std_sort);
register_bench_small_values_range(std_sort);

// Run std parallel sort benchmarks
// STL implementation is broken as it allocates all memory!!!!
/*
#define std_parallel_sort(dataType, bench, type, size, from, to)                                 \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                         \
                                BM_par_std_sort_##dataType##_##type##_##bench, type, size, from, \
                                to)                                                              \
    (benchmark::State & state) {                                                                 \
        for (auto _ : state) {                                                                   \
            state.PauseTiming();                                                                 \
            Prepare();                                                                           \
            state.ResumeTiming();                                                                \
            std::sort(std::execution::par, vec.begin(), vec.end());                              \
            state.PauseTiming();                                                                 \
            Destroy();                                                                           \
            state.ResumeTiming();                                                                \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture, BM_par_std_sort_##dataType##_##type##_##bench) \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_par_std_sort_##dataType##_##type##_##bench));

register_bench_default(std_parallel_sort);
register_bench_small_sizes(std_parallel_sort);
*/

// Run mpqsort parallel benchmarks
#define mpqsort_par_sort(dataType, bench, type, size, from, to)                                \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                       \
                                BM_mpqsort_par_sort_##dataType##_##type##_##bench, type, size, \
                                from, to)                                                      \
    (benchmark::State & state) {                                                               \
        for (auto _ : state) {                                                                 \
            state.PauseTiming();                                                               \
            Prepare();                                                                         \
            state.ResumeTiming();                                                              \
            mpqsort::sort(mpqsort::execution::par, vec.begin(), vec.end());            \
            state.PauseTiming();                                                               \
            Destroy();                                                                         \
            state.ResumeTiming();                                                              \
        }                                                                                      \
    }                                                                                          \
    BENCHMARK_REGISTER_F(dataType##VectorFixture,                                              \
                         BM_mpqsort_par_sort_##dataType##_##type##_##bench)                    \
        ->MeasureProcessCPUTime()                                                              \
        ->UseRealTime()                                                                        \
        ->Name(str(BM_mpqsort_par_sort_##dataType##_##type##_##bench));

register_bench_default(mpqsort_par_sort);
register_bench_small_sizes(mpqsort_par_sort);
register_bench_small_values_range(mpqsort_par_sort);

// Run mpqsort sequential three way benchmarks
#define mpqsort_seq_three_way_sort(dataType, bench, type, size, from, to)                          \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                           \
                                BM_mpqsort_seq_three_way_sort_##dataType##_##type##_##bench, type, \
                                size, from, to)                                                    \
    (benchmark::State & state) {                                                                   \
        for (auto _ : state) {                                                                     \
            state.PauseTiming();                                                                   \
            Prepare();                                                                             \
            state.ResumeTiming();                                                                  \
            mpqsort::sort(mpqsort::execution::seq_three_way, vec.begin(), vec.end());              \
            state.PauseTiming();                                                                   \
            Destroy();                                                                             \
            state.ResumeTiming();                                                                  \
        }                                                                                          \
    }                                                                                              \
    BENCHMARK_REGISTER_F(dataType##VectorFixture,                                                  \
                         BM_mpqsort_seq_three_way_sort_##dataType##_##type##_##bench)              \
        ->MeasureProcessCPUTime()                                                                  \
        ->UseRealTime()                                                                            \
        ->Name(str(BM_mpqsort_seq_three_way_sort_##dataType##_##type##_##bench));

register_bench_default(mpqsort_seq_three_way_sort);
register_bench_small_sizes(mpqsort_seq_three_way_sort);
register_bench_small_values_range(mpqsort_seq_three_way_sort);

// Run mpqsort sequential four way benchmarks
#define mpqsort_seq_four_way_sort(dataType, bench, type, size, from, to)                          \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                          \
                                BM_mpqsort_seq_four_way_sort_##dataType##_##type##_##bench, type, \
                                size, from, to)                                                   \
    (benchmark::State & state) {                                                                  \
        for (auto _ : state) {                                                                    \
            state.PauseTiming();                                                                  \
            Prepare();                                                                            \
            state.ResumeTiming();                                                                 \
            mpqsort::sort(mpqsort::execution::seq_four_way, vec.begin(), vec.end());              \
            state.PauseTiming();                                                                  \
            Destroy();                                                                            \
            state.ResumeTiming();                                                                 \
        }                                                                                         \
    }                                                                                             \
    BENCHMARK_REGISTER_F(dataType##VectorFixture,                                                 \
                         BM_mpqsort_seq_four_way_sort_##dataType##_##type##_##bench)              \
        ->MeasureProcessCPUTime()                                                                 \
        ->UseRealTime()                                                                           \
        ->Name(str(BM_mpqsort_seq_four_way_sort_##dataType##_##type##_##bench));

register_bench_default(mpqsort_seq_four_way_sort);
register_bench_small_sizes(mpqsort_seq_four_way_sort);
register_bench_small_values_range(mpqsort_seq_four_way_sort);

// Run gnu qs sort benchmarks
// In-place parallel qsort
#define gnu_qs_sort(dataType, bench, type, size, from, to)                                      \
    BENCHMARK_TEMPLATE_DEFINE_F(dataType##VectorFixture,                                        \
                                BM_gnu_qs_sort_##dataType##_##type##_##bench, type, size, from, \
                                to)                                                             \
    (benchmark::State & state) {                                                                \
        for (auto _ : state) {                                                                  \
            state.PauseTiming();                                                                \
            Prepare();                                                                          \
            state.ResumeTiming();                                                               \
            __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::quicksort_tag());      \
            state.PauseTiming();                                                                \
            Destroy();                                                                          \
            state.ResumeTiming();                                                               \
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
            state.PauseTiming();                                                                 \
            Prepare();                                                                           \
            state.ResumeTiming();                                                                \
            __gnu_parallel::sort(vec.begin(), vec.end(),                                         \
                                 __gnu_parallel::balanced_quicksort_tag());                      \
            state.PauseTiming();                                                                 \
            Destroy();                                                                           \
            state.ResumeTiming();                                                                \
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
            state.PauseTiming();                                                                  \
            Prepare();                                                                            \
            state.ResumeTiming();                                                                 \
            __gnu_parallel::sort(vec.begin(), vec.end(),                                          \
                                 __gnu_parallel::multiway_mergesort_tag());                       \
            state.PauseTiming();                                                                  \
            Destroy();                                                                            \
            state.ResumeTiming();                                                                 \
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
            state.PauseTiming();                                                                 \
            Prepare();                                                                           \
            state.ResumeTiming();                                                                \
            tbb::parallel_sort(vec.begin(), vec.end());                                          \
            state.PauseTiming();                                                                 \
            Destroy();                                                                           \
            state.ResumeTiming();                                                                \
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
            state.PauseTiming();                                                                 \
            Prepare();                                                                           \
            state.ResumeTiming();                                                                \
            thrust::sort(vec.begin(), vec.end());                                                \
            state.PauseTiming();                                                                 \
            Destroy();                                                                           \
            state.ResumeTiming();                                                                \
        }                                                                                        \
    }                                                                                            \
    BENCHMARK_REGISTER_F(dataType##VectorFixture,                                                \
                         BM_nvidia_thrust_sort_##dataType##_##type##_##bench)                    \
        ->MeasureProcessCPUTime()                                                                \
        ->UseRealTime()                                                                          \
        ->Name(str(BM_nvidia_thrust_sort_##dataType##_##type##_##bench));

register_bench_default(nvidia_thrust_sort);
register_bench_small_sizes(nvidia_thrust_sort);

register_bench_small_size_threshold(std_sort);
register_bench_small_size_threshold(mpqsort_seq_three_way_sort)
    register_bench_small_size_threshold(mpqsort_seq_four_way_sort)
        register_bench_small_size_threshold(mpqsort_par_sort)