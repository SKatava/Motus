#include "Motus/Simulation.h"

#include <cstdio>
#include <iostream>
#include <vector>


int main(void) {

    std::vector<double> x = {250.0, -250.0};
    std::vector<double> y = {250.0, -250.0};
    std::vector<double> z = {0.0, 0.0};
    std::vector<double> vx = {0.0, 0.0};
    std::vector<double> vy = {0.0, 0.0};
    std::vector<double> vz = {0.0, 0.0};
    std::vector<double> m = {1.0, 1.0};

    std::vector<glm::dvec3> pos = {{x[0], y[0], z[0]}, {x[1], y[1], z[1]}};
    std::vector<glm::dvec3> vel = {{vx[0], vy[0], vz[0]}, {vx[1], vy[1], vz[1]}};

    motus::Simulation sim;

    size_t step = 0;

    for(;;) {
        sim.step({pos, vel, m}, 0.1);
        if(step++ % 1000000 == 0) {
            std::cout << "BODY ONE: pos=(" << pos[0].x << ", " << pos[0].y << ") vel=(" << vel[0].x << ", " << vel[0].y << ")\n";
            std::cout << "BODY TWO: pos=(" << pos[1].x << ", " << pos[1].y << ") vel=(" << vel[1].x << ", " << vel[1].y << ")\n";
            std::cout << "E: " << sim.totalEnergy({pos, vel, m}) << "\n";
        }
    }

    return 0;
}
