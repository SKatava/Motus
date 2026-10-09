#pragma once


#include <span>
#include <vector>
#include <algorithm>


struct BodiesView {
    std::span<double> x;
    std::span<double> y;
    std::span<double> vx;
    std::span<double> vy;
    std::span<double> m;

    [[nodiscard]] std::size_t size() const { return m.size();}
};

struct AccelView {
    std::span<double> ax;
    std::span<double> ay;
};

struct Accel {
    std::vector<double> ax;
    std::vector<double> ay;

    [[nodiscard]] std::size_t size() const { return ax.size(); }

    void resize(size_t size) { 
        ax.resize(size);
        ay.resize(size);
    }

    void zero() {
        std::ranges::fill(ax, 0.0);
        std::ranges::fill(ay, 0.0);
    }

    AccelView view() {
        return {ax, ay};
    }

};
