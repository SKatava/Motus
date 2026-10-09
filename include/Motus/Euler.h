#pragma once

#include "Integrator.h"

namespace motus {

class Euler final : public Integrator {
    public:
        void step(BodiesView b, const ForcesFn& f, double dt) override;
    private:
        Accel m_acc;
};

}
