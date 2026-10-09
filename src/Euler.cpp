#include "Motus/Euler.h"

namespace motus {

void Euler::step(BodiesView b, const ForcesFn& f, double dt) {
    evaluate(b, f, m_acc);
    for(size_t i = 0; i < b.size(); ++i) {
        b.x[i] += b.vx[i] * dt;
        b.y[i] += b.vy[i] * dt;
        b.vx[i] += m_acc.ax[i] * dt;
        b.vy[i] += m_acc.ay[i] * dt;
    } 
}

}
