#ifndef CPP11SORT_IMPL_BARRIER_H
#define CPP11SORT_IMPL_BARRIER_H

#include <cstdlib>
#include <condition_variable>
#include <mutex>

namespace cpp11sort
{
namespace impl
{

/**
 * Creates a barrier for std::threads so they can continue only if all of them reach the
 * specific point. Provides synchronization of the threads. To enable creating multiple
 * points of the barrier, it either counts up or down the number of threads.
 */
class Barrier
{
  public:
    explicit Barrier(std::size_t count) : count(0),
                                          number_of_threads(count),
                                          is_increasing(true) {}
    /**
     * Waits for all the threads to reach this point.
     */
    void Wait()
    {
        std::unique_lock<std::mutex> lock(mutex);
        if (is_increasing)
        {
            if (++count != number_of_threads)
            {
                cv.wait(lock, [this] { return !is_increasing; });
            }
            else
            {
                is_increasing = false;
                cv.notify_all();
            }
        }
        else
        {
            if (--count != 0)
            {
                cv.wait(lock, [this] { return is_increasing; });
            }
            else
            {
                is_increasing = true;
                cv.notify_all();
            }
        }
    }

  private:
    std::mutex mutex;
    std::condition_variable cv;
    std::size_t count;
    std::size_t number_of_threads;
    bool is_increasing;
};

}

}

#endif
