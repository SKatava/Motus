#include "Motus/Bodies.h"
#include "Motus/Gravity.h"
#include "Motus/Integrator.h"
#include "Motus/Simulation.h"

#include <cmath>
#include <cstdio>
#include <numbers>
#include <iostream>
#include <vector>


int main(void) {

    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {0.0, 0.0};
    std::vector<double> vx = {0.0, 0.0};
    std::vector<double> vy = {0.0, 0.0};
    std::vector<double> m = {1.0, 1.0};

    motus::Simulation sim;

    size_t step = 0;

    for(;;) {
        sim.step({x, y, vx, vy, m}, 0.1);
        
        if(step++ % 1000000 == 0) {
            std::cout << "BODY ONE: pos=(" << x[0] << ", " << y[0] << " vel=(" << vx[0] << ", " << vy[0] << "\n";
            std::cout << "BODY TWO: pos=(" << x[1] << ", " << y[1] << " vel=(" << vx[1] << ", " << vy[1] << "\n";
            std::cout << "E: " << sim.totalEnergy({x,y,vx,vy,m}) << "\n";
        }
    }

    return 0;
}
