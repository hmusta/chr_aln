#pragma once

#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <gperftools/profiler.h>

using namespace std::chrono_literals;

std::jthread start_profiler_flusher(std::chrono::milliseconds period = 250ms) {
    return std::jthread([period](std::stop_token st) {
        std::mutex m;
        std::condition_variable_any cv;
        std::unique_lock lk(m);
        while (!cv.wait_for(lk, st, period, [&] { return st.stop_requested(); })) {
            ProfilerFlush();
        }
        ProfilerFlush(); // capture samples since the last periodic flush
    });
}