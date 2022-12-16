/**
 * Synchronized queue, to enqueue elements thread-safe. With some edits of behaviour if
 * queue is empty from
 * https://github.com/tigerjack/prog_mgmt/blob/bf00f2b4ab3fc1310c9dca9e1cbb48c0700a017b/c_c%2B%2B/AOS/multithreading_examples/13_thread_pool_good/synchronized_queue.h
 */

#ifndef CPP11SORT_IMPL_SYNC_QUEUE
#define CPP11SORT_IMPL_SYNC_QUEUE

#include <list>
#include <functional>
#include <mutex>
#include <stdexcept>

namespace cpp11sort
{
namespace impl 
{

template<typename T>
class SynchronizedQueue {
  public:
    SynchronizedQueue() = default;

    SynchronizedQueue(const SynchronizedQueue &) = delete;

    SynchronizedQueue &operator=(const SynchronizedQueue &) = delete;

    /**
     * Enqueues one element to the queue.
     * @tparam T Type of the elements in the queue.
     * @param data Element to be enqueued to the queue.
     */
    void Push(const T &data);

    /**
     * Gets one element from the queue.
     * Throws an exception if the queue is empty.
     * @tparam T Type of the element in the queue.
     * @return One element from the queue.
     */
    T Get();

    /**
     *  Returns size of the queue.
     * @return Size of the queue.
     */
    size_t Size();

  private:

    std::list<T> queue;
    std::mutex queue_mutex;
};

template<typename T>
void SynchronizedQueue<T>::Push(const T &data) {
    std::unique_lock<std::mutex> lck(queue_mutex);
    queue.push_back(data);
}

template<typename T>
T SynchronizedQueue<T>::Get() {
    std::unique_lock<std::mutex> lck(queue_mutex);
    if (queue.empty()) {
        throw (std::logic_error("Queue is empty."));
    }
    T result = queue.front();
    queue.pop_front();
    return result;
}

template<typename T>
size_t SynchronizedQueue<T>::Size() {
    std::unique_lock<std::mutex> lck(queue_mutex);
    return queue.size();
}

template
class SynchronizedQueue<std::function<void()>>;

}

}

#endif
