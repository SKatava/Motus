#pragma once

#include <functional>
#include <algorithm>

#include "Motus/Bodies.h"


namespace motus {

enum class IntegratorMethod {
    EULER,
    VERLET
};

using ForcesFn = std::function<void(BodiesView, AccelView)>;
    
class Integrator {
    public:
        //Destrucotrs
        virtual ~Integrator() = default; 

        //Core functions
        virtual void step(BodiesView b, const ForcesFn& f, double dt) = 0;
    protected:
        static void evaluate(BodiesView b, const ForcesFn& f, Accel& out) {
            out.resize(b.size());
            std::ranges::fill(out, glm::dvec3{0, 0, 0});
            f(b, {out});
        }
};
}
