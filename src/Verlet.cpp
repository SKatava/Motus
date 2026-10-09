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
        b.x[i] += b.vx[i] * dt + (m_acc.ax[i]*dt*dt) / (2);
        b.y[i] += b.vy[i] * dt + (m_acc.ay[i]*dt*dt) / (2);
    }

    evaluate(b, f ,m_accNew);

    for(size_t i = 0; i < n; ++i) {
        b.vx[i] += 0.5 * (m_acc.ax[i] + m_accNew.ax[i]) * dt;
        b.vy[i] += 0.5 * (m_acc.ay[i] + m_accNew.ay[i]) * dt;
    }

    std::swap(m_acc, m_accNew);
}

}
