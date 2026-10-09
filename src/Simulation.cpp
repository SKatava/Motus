#include "Motus/Simulation.h"
#include "Motus/Euler.h"
#include "Motus/Verlet.h"

#include <cassert>

std::unique_ptr<motus::Integrator> makeIntegrator(motus::IntegratorMethod m) {
    switch(m) {
        case motus::IntegratorMethod::EULER: return std::make_unique<motus::Euler>();
        case motus::IntegratorMethod::VERLET: return std::make_unique<motus::Verlet>();
    }
    assert(false && "MOTUS: unknown integrator method");
    std::abort();
}

namespace motus {

Simulation::Simulation(IntegratorMethod m) : m_forceSet(), m_integrator(makeIntegrator(m)){}

void Simulation::step(BodiesView b, double dt) {
    m_integrator->step(b, m_forceSet, dt); 
    m_time += dt;
}


double Simulation::kineticEnergy(BodiesView b) {
    double K = 0.0;
    for(size_t i = 0; i < b.size(); ++i) {
        K += 0.5 * b.m[i] * (b.vx[i] * b.vx[i] + b.vy[i]*b.vy[i]);
    }
    return K;

}

double Simulation::totalEnergy(BodiesView b) {
    return kineticEnergy(b) + m_forceSet.potentialEnergy(b); 
}

}
