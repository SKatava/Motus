#pragma once

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string_view>
#include <vector>
#include <limits>
#include <unordered_map>

namespace motus {

struct SectionStats {
    std::int64_t count = 0;
    double total = 0;
    double min = std::numeric_limits<double>::max();
    double max = 0;
    double last = 0;

    void add(double ms) {
        ++count;
        total += ms;
        last = ms;
        if(ms < min) min = ms;
        if(ms > max) max = ms;
    }

    double avg() const { return count ? (total / count) : 0; }
};


class Profiler {
    public:
        static Profiler& get() {
            static Profiler p;
            return p;
        }

        void record(std::string_view name, double ms) {
            auto [it, inserted] = m_stats.try_emplace(name);
            if(inserted) m_order.push_back(name);
            it->second.add(ms);
        }

        void report() const {
            std::printf("%-22s %8s %8s %8s %8s\n", "section", "avg", "min", "max", "count");
            for(auto n : m_order) {
                const auto& s = m_stats.at(n);
                std::printf("%-22.*s %8.3f %8.3f %8.3f %8lld\n", int(n.size()), n.data(), s.avg(), s.min, s.max, (long long)s.count); 
            }
        }

    private:
        std::unordered_map<std::string_view, SectionStats> m_stats;
        std::vector<std::string_view> m_order;
};

class ScopeTimer {
public:
    explicit ScopeTimer(std::string_view name) : m_name(name), m_t0(clock::now()) {}
    ~ScopeTimer() {
        double ms = std::chrono::duration<double, std::milli>(clock::now() - m_t0).count();
        Profiler::get().record(m_name, ms);
    }
private:
    using clock = std::chrono::steady_clock;
    std::string_view m_name; clock::time_point m_t0;
};

}

#ifdef MOTUS_PROFILE
  #define MOTUS_CONCAT2(a,b) a##b
  #define MOTUS_CONCAT(a,b)  MOTUS_CONCAT2(a,b)
  #define PROFILE_SCOPE(name) ::motus::ScopeTimer MOTUS_CONCAT(_pt_, __LINE__){name}
#else
  #define PROFILE_SCOPE(name) ((void)0)
#endif
