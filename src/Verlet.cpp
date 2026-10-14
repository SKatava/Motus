#include "Motus/Verlet.h"

namespace motus {

Verlet::Verlet() {}

void Verlet::step(BodiesView b, const ForcesFn &f, double dt) {
    size_t n = b.size();
    
    if(m_acc.size() != n) {
        m_acc.resize(n);
        m_accNew.resize(n);
        evaluate(b, f, m_accNew);
    }

    for(size_t i = 0; i < n; ++i) {
        b.pos[i] += b.vel[i] * dt + (m_acc[i]*dt*dt) * 0.5;
    }

    evaluate(b, f ,m_accNew);

    for(size_t i = 0; i < n; ++i) {
        b.vel[i] += 0.5 * (m_acc[i] + m_accNew[i]) * dt;
    }

    std::swap(m_acc, m_accNew);
}

}
