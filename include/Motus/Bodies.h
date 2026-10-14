#pragma once

#include <glm/vec3.hpp>

#include <span>
#include <vector>
#include <algorithm>

namespace motus {

struct BodiesView {
    std::span<glm::dvec3> pos;
    std::span<glm::dvec3> vel;
    std::span<double> m;

    [[nodiscard]] size_t size() const { return m.size();}
};

using AccelView = std::span<glm::dvec3>;
using Accel = std::vector<glm::dvec3>;

}
/*
struct AccelView {
    std::span<double> ax;
    std::span<double> ay;
    std::span<double> az;
};

struct Accel {
    std::vector<double> ax;
    std::vector<double> ay;
    std::vector<double> az;

    [[nodiscard]] std::size_t size() const { return ax.size(); }

    void resize(size_t size) { 
        ax.resize(size);
        ay.resize(size);
        az.resize(size);
    }

    void zero() {
        std::ranges::fill(ax, 0.0);
        std::ranges::fill(ay, 0.0);
        std::ranges::fill(az, 0.0);
    }

    AccelView view() {
        return {ax, ay, az};
    }

};*/
