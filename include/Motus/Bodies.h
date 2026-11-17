#pragma once

#include <glm/vec3.hpp>

#include <span>
#include <vector>

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

