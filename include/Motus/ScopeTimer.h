#pragma once

#include <chrono>

namespace motus {

struct ScopeTimer {
    const char* name;
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    ~ScopeTimer() {
        double ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - t0).count();
        printf("%s: %.3f ms\n", name, ms);
    }
};

}
