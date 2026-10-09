#pragma once

#include "Integrator.h"

namespace motus {

class Verlet final : public Integrator {
    public:
        Verlet();
        void step(BodiesView b, const ForcesFn& f, double dt) override;
    private:
        Accel m_acc;
        Accel m_accNew;

};

}
