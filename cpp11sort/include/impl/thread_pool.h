#ifndef CPP11SORT_THREAD_POOL_H
#define CPP11SORT_THREAD_POOL_H

#include <atomic>
#include <cstdlib>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#include "synchronized_queue.h"

namespace cpp11sort
{
namespace impl
{

/**
 * Manages pool of threads that can be used repeatedly in cycles.
 */
class ThreadPool {
  public:
    ThreadPool(int threads_count, bool recursive_tasks_generating);

    ~ThreadPool();

    /**
     * Submits task to be executed by threads.
     * @param task Task to be executed.
     */
    void PushTask(const std::function<void()> &task);

    /**
     * @return Number of threads in the pool.
     */
    int GetNumberOfThreads();

    /**
     * Sets variables for a new work cycle.
     */
    void SetNewWorkCycle();

    /**
     * Notifies threads in the thread pool that the tasks for them are ready.
     */
    void SetTasksReady();

    /**
     * Declares that no more tasks will be submitted to the thread pool queue (from
     * outside, tasks generated from the inside of the threadpool can still be submited
     * to the queue).
     */
    void FinishThreadPool();

    std::mutex exclusive_access; ///< Provides exclusive access to the pool.
    std::mutex
        mutex_finished_threads; ///< Used for waiting for all threads to finish.
    std::condition_variable cv_finished_threads; ///< Triggers if all threads are
    // finished.
    std::atomic<std::size_t> finished_workers{0}; ///< Number of workers already finished with
    // their tasks


  private:
    /**
     * Working thread that executes tasks from the queue. It runs in the "infinite"
     * loop until all tasks are done and the signal for finishing threads is sent.
     */
    void WorkerThread(int index);

    void WaitForCycleStart();

    void WaitForTaskReady();

    void AttemptToGetTask();

    SynchronizedQueue<std::function<void()>>
        work_queue; ///< Queue of tasks for threads.
    std::vector<std::thread> threads; ///< Workers processing tasks from queue.
    std::size_t total_threads_count; ///< Number of threads in the pool.
    bool recursive_tasks_generating = false; ///< Tasks are generated recursively or
    // all at once.
    std::atomic<int> working_threads_count{0}; ///< Number of threads working on some
    // task at the moment.
    std::atomic<bool> tasks_ready{false}; ///< Determines if data are ready and the work
    // cycle can start.
    std::atomic<bool> finish_threads{false}; ///< Determines that threads should end
    // their loop, stop waiting for tasks and die.
    std::condition_variable
        cv_work_cycle; ///< Triggers if working cycle should start.
    std::condition_variable cv_task_ready; ///< Triggers if task for threads is ready.
    std::mutex mutex_get_task; ///< Provides exclusive getting tasks from the queue.
    std::mutex mutex_start_working; ///< Protects start of executing tasks from the
    // queue.

};

ThreadPool::ThreadPool(int threads_count, bool recursive_tasks_generating) {
    this->recursive_tasks_generating = recursive_tasks_generating;
    total_threads_count = threads_count;
    for (int i = 0; i < threads_count; ++i) {
        threads.emplace_back(std::thread(&ThreadPool::WorkerThread, this, i));
    }
}

ThreadPool::~ThreadPool() {
    finish_threads = true;
    cv_work_cycle.notify_all();
    cv_task_ready.notify_all();
    for (auto &worker_thread: threads) {
        if (worker_thread.joinable()) {
            worker_thread.join();
        }
    }
}

void ThreadPool::PushTask(const std::function<void()> &task) {
    work_queue.Push(task);
    cv_task_ready.notify_one();
}

int ThreadPool::GetNumberOfThreads() {
    return total_threads_count;
}

void ThreadPool::SetNewWorkCycle() {
    finished_workers = 0;
    finish_threads = false;
    tasks_ready = false;
    working_threads_count = 0;
}

void ThreadPool::SetTasksReady() {
    tasks_ready = true;
    cv_work_cycle.notify_all();
}

void ThreadPool::FinishThreadPool() {
    finish_threads = true;
    cv_work_cycle.notify_all();
    cv_task_ready.notify_all();
    for (auto &worker_thread: threads) {
        if (worker_thread.joinable()) {
            worker_thread.join();
        }
    }
}

void ThreadPool::WaitForCycleStart() {
    std::unique_lock<std::mutex> lock_start_working(mutex_start_working);
    cv_work_cycle.wait(lock_start_working, [&] {
      // Can continue if work is ready or the thread should finish his
      // working loop.
      return tasks_ready && (recursive_tasks_generating || work_queue.Size() > 0 ||
          (finish_threads && working_threads_count == 0));
    });
}

void ThreadPool::WaitForTaskReady() {
    std::unique_lock<std::mutex> lock_get_job(mutex_get_task);
    cv_task_ready.wait(lock_get_job, [&] {
      return work_queue.Size() > 0 || working_threads_count == 0;
    });
}

void ThreadPool::AttemptToGetTask() {
    std::function<void()> task_to_be_done;
    std::unique_lock<std::mutex> lock_get_job(mutex_get_task);
    if (work_queue.Size() > 0) {
        try {
            working_threads_count++;
            task_to_be_done = work_queue.Get();
            lock_get_job.unlock();
            task_to_be_done();
        }
        catch (std::exception &e) {
            // Task from queue wasn't acquired, try again later.
        }
        working_threads_count--;
    }
}

void ThreadPool::WorkerThread(int index) {
    while (!finish_threads || working_threads_count > 0 ||
        work_queue.Size() > 0) {
        WaitForCycleStart();
        while (working_threads_count > 0 || work_queue.Size() > 0) {
            WaitForTaskReady();
            AttemptToGetTask();
        }
        cv_task_ready.notify_all();
        finished_workers++;
        if (finished_workers == total_threads_count) {
            cv_finished_threads.notify_one();
        }
    }
    cv_work_cycle.notify_all();
}

}

}

#endif
