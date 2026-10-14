#include "Motus/Euler.h"

namespace motus {

void Euler::step(BodiesView b, const ForcesFn& f, double dt) {
    evaluate(b, f, m_acc);
    for(size_t i = 0; i < b.size(); ++i) {
        b.pos[i] += b.vel[i] * dt;
        b.vel[i] += m_acc[i] * dt;
    } 
}

}
