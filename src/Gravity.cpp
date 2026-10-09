#include "Motus/Gravity.h"
#include <cmath>

namespace motus {

Gravity::Gravity(double G, double eps) : m_G(G), m_eps(eps) {}

void Gravity::compute(BodiesView b, AccelView out) const {
    for(size_t i = 0; i < b.size(); ++i) {
        for(size_t j = i + 1; j < b.size(); ++j) {
            double dx = b.x[j] - b.x[i];
            double dy = b.y[j] - b.y[i];
            double r2  = dx*dx + dy*dy + m_eps*m_eps;

            double inv_r3 = 1 / (r2 * std::sqrt(r2));

            out.ax[i] += m_G * b.m[j] * dx * inv_r3; 
            out.ay[i] += m_G * b.m[j] * dy * inv_r3; 

            out.ax[j] += m_G * b.m[i] * (-dx) * inv_r3; 
            out.ay[j] += m_G * b.m[i] * (-dy) * inv_r3; 
        }
    }
}

double Gravity::potentialEnergy(BodiesView b) const {
    double U = 0.0;
    for(size_t i = 0; i < b.size(); ++i) {
        for(size_t j = i+1; j < b.size(); ++j) {
            double dx = b.x[j] - b.x[i];
            double dy = b.y[j] - b.y[i];
            double r  = std::sqrt(dx*dx + dy*dy + m_eps*m_eps);
            U -= m_G * b.m[i] * b.m[j] / r;
        }
    }
    return U;
}

 

void Gravity::setG(double newG) { m_G = newG; }
void Gravity::resetG() { m_G = s_defaultG; }
void Gravity::setEps(double newEps) { m_eps = newEps; }
void Gravity::resetEps() {m_eps = s_defaultEps; }

double Gravity::G() const noexcept { return m_G; }
double Gravity::eps() const noexcept { return m_eps; }

}

