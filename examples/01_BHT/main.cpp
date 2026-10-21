#include "Motus/Simulation.h"

#include <iostream>
#include <random>
#include <vector>

struct Bodies {
    std::vector<glm::dvec3> pos;
    std::vector<glm::dvec3> vel;
    std::vector<double> m;
};


void generateBodies(size_t n, Bodies& b) {
    std::random_device rd;
    std::mt19937 gen(21);
    std::uniform_real_distribution<double> real(250, 750.0);
    for(int i {}; i < n; ++i) {
        b.pos.push_back({real(gen), real(gen), real(gen)});
        b.vel.push_back({real(gen), real(gen), real(gen)});
        b.m.push_back(1);
    }
    }

int main(void) {
    Bodies b;
    generateBodies(10'000, b);

    motus::Simulation sim;
    sim.gravity().setAlgorithm(motus::Gravity::Algorithm::BARNES_HUT);
    size_t step = 1;
    

    for(;;) {
        {PROFILE_SCOPE("sim"); sim.step({b.pos, b.vel, b.m}, 0.1);}
        if(step++ % 100 == 0) {
            motus::Profiler::get().report();
        }
    }

    return 0;
}
