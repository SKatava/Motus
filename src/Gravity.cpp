#include "Motus/Gravity.h"
#include <cmath>

#include "Motus/BarnesHutTree.h"
#include "Motus/Profiler.h"

namespace motus {

Gravity::Gravity(double G, double eps) : m_G(G), m_eps(eps) {}

void Gravity::compute(BodiesView b, AccelView out) const {
    if(m_algo == Algorithm::DIRECT_NEWTON) directNewton(b, out);    
    else if(m_algo == Algorithm::BARNES_HUT) barnesHut(b, out);    
}

double Gravity::potentialEnergy(BodiesView b) const {
    double U = 0.0;
    for(size_t i = 0; i < b.size(); ++i) {
        for(size_t j = i+1; j < b.size(); ++j) {
            double dx = b.pos[j].x - b.pos[i].x;
            double dy = b.pos[j].y - b.pos[i].y;
            double dz = b.pos[j].z - b.pos[i].z;
            double r  = std::sqrt(dx*dx + dy*dy + dz*dz + m_eps*m_eps);
            U -= m_G * b.m[i] * b.m[j] / r;
        }
    }
    return U;
}

 void Gravity::directNewton(BodiesView b, AccelView out) const {
    for(size_t i = 0; i < b.size(); ++i) {
        for(size_t j = i + 1; j < b.size(); ++j) {
            double dx = b.pos[j].x - b.pos[i].x;
            double dy = b.pos[j].y - b.pos[i].y;
            double dz = b.pos[j].z - b.pos[i].z;
            double r2  = dx*dx + dy*dy + dz*dz + m_eps*m_eps;

            double inv_r3 = 1 / (r2 * std::sqrt(r2));

            out[i].x += m_G * b.m[j] * dx * inv_r3; 
            out[i].y += m_G * b.m[j] * dy * inv_r3; 
            out[i].z += m_G * b.m[j] * dz * inv_r3; 

            out[j].x += m_G * b.m[i] * (-dx) * inv_r3; 
            out[j].y += m_G * b.m[i] * (-dy) * inv_r3; 
            out[j].z += m_G * b.m[i] * (-dz) * inv_r3; 
        }
    }
 }


void Gravity::barnesHut(BodiesView b, AccelView out) const {
    BarnesHutTree tree(b);
    tree.build();
    const std::vector<int>& indices = tree.indices();
    {
        PROFILE_SCOPE("bh/traverse");
        #pragma omp parallel for schedule(dynamic, 64)
        for(size_t i = 0; i < b.size(); ++i) {
            out[indices[i]] += tree.traverseNode(indices[i], 0, m_theta, m_G, m_eps) / b.m[indices[i]];
        }
    }
}

void Gravity::setG(double newG) { m_G = newG; }
void Gravity::resetG() { m_G = s_defaultG; }
void Gravity::setEps(double newEps) { m_eps = newEps; }
void Gravity::resetEps() {m_eps = s_defaultEps; }
void Gravity::setAlgorithm(Algorithm newAlgo) { m_algo = newAlgo; }
void Gravity::setTheta(double newTheta) { m_theta = newTheta; }
void Gravity::resetTheta() { m_theta = s_defaultTheta;}

double Gravity::G() const noexcept { return m_G; }
double Gravity::eps() const noexcept { return m_eps; }
double Gravity::theta() const noexcept { return m_theta; }

}

