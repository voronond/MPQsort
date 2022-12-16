#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

#define CPP11SORT_SEQUENTIAL_THRESHOLD 100000
#include <cpp11sort.h>

#ifdef WITH_STDCPP_PM
#include <parallel/algorithm>
#include <omp.h>
#endif

#ifdef WITH_TBB
#include <tbb/parallel_sort.h>
//#include <tbb/task_scheduler_init.h> // deprecated
#include <tbb/global_control.h>      // new
#endif

#ifdef WITH_AQSORT
#include <aqsort.h>
#include <omp.h>
#endif

#ifdef WITH_INTEL_PSTL
#include <pstl/algorithm>
#include <pstl/execution>
#endif

#ifdef WITH_MS_PSTL
#include <execution>
#endif

// data set types:
enum 
{
   RANDOM,           // 0: random
   SORTED,           // 1: sorted
   INV_SORTED,       // 2: inversly sorted
   REPEATED,         // 3: with repeated values
   GLOB_PART_SORTED, // 4: globally partially-sorted
   LOC_PART_SORTED   // 5: locally partially-sorted
};

std::vector<int> generator(long n, int type, long seed)
{
   std::vector<int> v;
   v.reserve(n);

   std::mt19937 eng(seed);

   switch (type)
   {
      case RANDOM:
         std::cout << "Data set type: random with seed = " << seed << std::endl;
         {
            std::uniform_int_distribution<int> dist;
            for (long i = 0; i < n; i++)
               v.emplace_back(dist(eng));
         }
         break;

      case SORTED:
         std::cout << "Data set type: sorted" << std::endl;
         for (long i = 0; i < n; i++)
            v.emplace_back((int)i);
         break;

      case INV_SORTED:
         std::cout << "Data set type: inversely-sorted" << std::endl;
         for (long i = 0; i < n; i++)
            v.emplace_back((int)(n - i - 1));
         break;

      case REPEATED:
         std::cout << "Data set type: repeated random with seed = " << seed << std::endl;
         {
            std::uniform_int_distribution<int> dist(0, 99);
            for (long i = 0; i < n; i++)
               v.emplace_back(dist(eng));
         }
         break;

      case GLOB_PART_SORTED:
         std::cout << "Data set type: globally partially-sorted with seed = " << seed << std::endl;
         {
            int eps = INT_MAX / 10000;
            std::uniform_int_distribution<int> dist(0, eps - 1);
            for (long i = 0; i < n; i++)
            {
               long chunk = i / 10000;
               v.emplace_back((int)(eps * chunk + dist(eng)));
            }
         }
         break;

      case LOC_PART_SORTED:
         std::cout << "Data set type: locally partially-sorted with seed = " << seed << std::endl;
         {
            std::uniform_int_distribution<int> dist;
            for (long i = 0; i < n; i++)
               v.emplace_back(dist(eng));
            for (long k = 0; k < n / 10000; k++)
               std::sort(v.begin() + k * 10000, v.begin() + (k + 1) * 10000);
         }
         break;

      default:
         throw std::runtime_error("invalid data set type");
   }

   return v;
}

int main(int argc, char* argv[])
{
   long n = std::atol(argv[1]);
   int method = std::atoi(argv[2]);
   unsigned int num_threads = std::atoi(argv[3]); // must be unsiged int
   int type = std::atoi(argv[4]);

   long seed;
   if (argc > 5)
      seed = std::atol(argv[5]);
   else
      seed = std::random_device{}();

   auto v = generator(n, type, seed);

   enum {
      CPP11SORT,     // 0: CPP11sort 
      STDCPP_PM,     // 1: libstdc++ paralle mode quicksort
      STDCPP_PM_B,   // 2: libstdc++ paralle mode balanced quicksort
      TBB,           // 3: Intel TBB 
      INTEL_PSTL,    // 4: Intel Parallel STL
      MS_PSTL,       // 5: Microsort Parallel STL
      AQSORT,        // 6: AQsort
      STD_SORT       // 7: sequential std::sort
   };

#if defined(WITH_STDCPP_PM) || defined(WITH_AQSORT)
   omp_set_nested(1);
#endif

#if defined(WITH_TBB) || defined(WITH_INTEL_PSTL)
// tbb::task_scheduler_init init(num_threads); // deprecated
   tbb::global_control tbb_gc(tbb::global_control::max_allowed_parallelism, num_threads);
#endif

#ifdef WITH_AQSORT
   omp_set_num_threads(num_threads);
#endif

   auto t1 = std::chrono::high_resolution_clock::now();

   switch (method)
   {
      case CPP11SORT:
         std::cout << "CPP11sort...\n";
         cpp11sort::sort(v.begin(), v.end(), num_threads);
         break;

#ifdef WITH_STDCPP_PM
      case STDCPP_PM:
         std::cout << "Libstdc++ parallel mode quicksort...\n";
         __gnu_parallel::__parallel_sort_qs(v.begin(), v.end(), std::less<>(), num_threads);
         break;

      case STDCPP_PM_B:
         std::cout << "Libstdc++ parallel mode balanced quicksort...\n";
         __gnu_parallel::__parallel_sort_qsb(v.begin(), v.end(), std::less<>(), num_threads);
         break;
#endif

#ifdef WITH_TBB
      case TBB:
         std::cout << "Intel TBB...\n";
         tbb::parallel_sort(v.begin(), v.end());
         break;
#endif

#ifdef WITH_INTEL_PSTL
      case INTEL_PSTL:
         std::cout << "Intel Parallel STL...\n";
         std::sort(pstl::execution::par, v.begin(), v.end());
         break;
#endif

#ifdef WITH_MS_PSTL
      case MS_PSTL:
         std::cout << "Microsoft Parallel STL...\n";
         std::sort(std::execution::par, v.begin(), v.end());
         break;
#endif 

#ifdef WITH_AQSORT
      case AQSORT:
         {
            std::cout << "AQsort...\n";
            auto comp = [&v](size_t i, size_t j){ return v[i] < v[j]; };
            auto swap = [&v](size_t i, size_t j){ std::swap(v[i], v[j]); };
            aqsort::parallel_sort(v.size(), &comp, &swap);
         }
         break;
#endif

      case STD_SORT:
      default:
         std::cout << "sequential std::sort...\n";
         std::sort(v.begin(), v.end());
   }

   auto t2 = std::chrono::high_resolution_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

   std::cout << "Sorting time: " << duration.count() << " [ms]" << std::endl;

   // check
   for (long i = 1; i < n; i++)
      if (v[i] < v[i - 1])
         throw std::runtime_error("Output sequence is not sorted");
}
