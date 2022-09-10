#include <benchmark/benchmark.h>
#include <mpqsort/mpqsort.h>
#include <mpqsort/version.h>
#include <tbb/tbb.h>
#include <thrust/sort.h>

#include <algorithm>
#include <parallel/algorithm>
#include <string>
#include <vector>
#include <random>
#include <iostream>
#include <unistd.h>

// TODO remove when released
#define TESTING 1

template <typename T, int seed = 0> class RandomVectorFixture
    : public benchmark::Fixture {
  public:
    void SetUp(const ::benchmark::State& state) override {
        // Seed with a same value so that sort same sequences across runs
        std::mt19937 en(seed);
        vec.resize(VectorSizeToFillHalfMemory());
        // Generate all numbers the type can hols
        auto from = std::numeric_limits<T>::min();
        auto to = std::numeric_limits<T>::max();

        if constexpr (std::is_integral_v<T>) {
            std::uniform_int_distribution<T> uniform_dist(from, to);
            InitVector([&]() { return uniform_dist(en); });
        } else {
            std::uniform_real_distribution<T> uniform_dist(from, to);
            InitVector([&]() { return uniform_dist(en); });
        }
    }

    void TearDown(const ::benchmark::State& state) override { vec.clear(); }

    std::vector<int> vec;

  private:
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
};

// Run std sort benchmarks
#define std_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_std_sort_##type, type)(benchmark::State& state) {\
    for (auto _ : state) {\
        std::sort(vec.begin(), vec.end());\
    }\
}

std_sort(int);
std_sort(short);
std_sort(double);

// Run std parallel sort benchmarks
#define std_parallel_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_par_std_sort_##type, type)(benchmark::State& state){\
    for (auto _ : state) {\
        std::sort(std::execution::par, vec.begin(), vec.end());\
    }\
}

std_parallel_sort(int);
std_parallel_sort(short);
std_parallel_sort(double);

// Run mpqsort benchmarks
#define mpqsort_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_mpqsort_sort_##type, type)(benchmark::State& state){\
    for (auto _ : state) {\
        mpqsort::sort(mpqsort::execution::par_multi_way, vec.begin(), vec.end());\
    }\
}

mpqsort_sort(int);
mpqsort_sort(short);
mpqsort_sort(double);

// Run gnu qs sort benchmarks
// In-place parallel qsort
#define gnu_qs_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_gnu_qs_sort_##type, type)(benchmark::State& state){\
    for (auto _ : state) {\
        __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::quicksort_tag());\
    }\
}

gnu_qs_sort(int);
gnu_qs_sort(short);
gnu_qs_sort(double);

// Run gnu_bqs_sort benchmarks
// In-place balanced parallel qsort
#define gnu_bqs_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_gnu_bqs_sort_##type, type)(benchmark::State& state){\
    for (auto _ : state) {\
        __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::balanced_quicksort_tag());\
    }\
}

gnu_bqs_sort(int);
gnu_bqs_sort(short);
gnu_bqs_sort(double);

// Run gnu mwms sort benchmarks
// Out-of-place parallel mergesort
#define gnu_mwms_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_gnu_mwms_sort_##type, type)(benchmark::State& state){\
    for (auto _ : state) {\
        __gnu_parallel::sort(vec.begin(), vec.end(), __gnu_parallel::multiway_mergesort_tag());\
    }\
}

gnu_mwms_sort(int);
gnu_mwms_sort(short);
gnu_mwms_sort(double);

// Run tbb sort benchmarks
#define tbb_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_tbb_sort_##type, type)(benchmark::State& state){\
    for (auto _ : state) {\
        tbb::parallel_sort(vec.begin(), vec.end());\
    }\
}

tbb_sort(int);
tbb_sort(short);
tbb_sort(double);

// Run nvidia thrust sort benchmarks
#define nvidia_thrust_sort(type) BENCHMARK_TEMPLATE_F(RandomVectorFixture, BM_nvidia_thrust_sort_##type, type)(benchmark::State& state) {\
    for (auto _ : state) {\
        thrust::sort(vec.begin(), vec.end());\
    }\
}

nvidia_thrust_sort(int);
nvidia_thrust_sort(short);
nvidia_thrust_sort(double);