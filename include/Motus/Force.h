#pragma once

#include <cstdint>
#include "Gravity.h"

namespace motus {

enum class Forces : std::uint32_t {
    None    = 0,
    Gravity = 1u << 1
};

constexpr Forces operator|(Forces a, Forces b) {
    return static_cast<Forces>( static_cast<std::uint32_t>(a) | static_cast<std::uint32_t>(b) );
}

constexpr Forces operator&(Forces a, Forces b) {
    return static_cast<Forces>( static_cast<std::uint32_t>(a) & static_cast<std::uint32_t>(b) );
}

constexpr bool has(Forces set, Forces f) {
    return (set & f) != Forces::None;
}

struct ForceSet {
    Forces enabled = Forces::Gravity;
    Gravity gravity;

    void operator()(BodiesView b, AccelView out) const {
        if(has(Forces::Gravity, enabled)) gravity.compute(b, out); }

    double potentialEnergy(BodiesView b) const {
        double U = 0.0;
        if(has(Forces::Gravity, enabled)) { U += gravity.potentialEnergy(b); }
        return U;
    }
};

}
