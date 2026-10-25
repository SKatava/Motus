#include <gtest/gtest.h>

#include "Motus/Simulation.h"

#include "TestUtils.h"

using namespace motus;
using motus::test::Bodies;

namespace {

Simulation makeSim(IntegratorMethod m, double G, double eps,
                   Gravity::Algorithm algo = Gravity::Algorithm::DIRECT_NEWTON) {
    Simulation sim(m);
    sim.gravity().setG(G);
    sim.gravity().setEps(eps);
    sim.gravity().setAlgorithm(algo);
    return sim;
}

void run(Simulation& sim, Bodies& b, double dt, int steps) {
    for(int i = 0; i < steps; ++i) sim.step(b.view(), dt);
}

//Plummer-like ball of bodies with small random velocities
Bodies randomCluster(size_t n, unsigned seed) {
    Bodies b = test::randomCube(n, seed);
    std::mt19937 gen(seed + 1);
    std::uniform_real_distribution<double> v(-0.1, 0.1);
    for(auto& vel : b.vel) vel = {v(gen), v(gen), v(gen)};
    return b;
}

}

//=======================Energy=========================

TEST(Simulation, KineticEnergyUsesAllThreeComponents) {
    Simulation sim;
    Bodies b;
    b.add({0, 0, 0}, {1.0, 2.0, 3.0}, 2.0);
    b.add({1, 0, 0}, {0.0, 0.0, 4.0}, 0.5);

    //0.5*2*(1+4+9) + 0.5*0.5*16
    EXPECT_DOUBLE_EQ(sim.kineticEnergy(b.view()), 18.0);
}

TEST(Simulation, TotalEnergyIsKineticPlusPotential) {
    Simulation sim = makeSim(IntegratorMethod::VERLET, 1.5, 0.1);
    Bodies b = randomCluster(20, 4);

    EXPECT_DOUBLE_EQ(sim.totalEnergy(b.view()),
                     sim.kineticEnergy(b.view()) + sim.gravity().potentialEnergy(b.view()));
}

TEST(Simulation, CircularBinaryHasAnalyticEnergy) {
    //E = -G m^2 / (2 r) for two equal masses on a circular orbit
    const double G = 1.0, m = 1.0, r = 1.0;
    Simulation sim = makeSim(IntegratorMethod::VERLET, G, 0.0);
    Bodies b;
    test::circularBinary(b, G, m, r);

    EXPECT_NEAR(sim.totalEnergy(b.view()), -G * m * m / (2.0 * r), 1e-14);
}

//=======================Defaults=========================

TEST(Simulation, DefaultsToVerletWithPhysicalG) {
    Simulation sim;
    EXPECT_DOUBLE_EQ(sim.gravity().G(), 6.674e-11);

    //Verlet's signature: the very first step already uses the end-of-step acceleration
    Simulation verlet = makeSim(IntegratorMethod::VERLET, 1.0, 0.0);
    Simulation def;
    def.gravity().setG(1.0);
    def.gravity().setEps(0.0);

    Bodies a, c;
    test::circularBinary(a, 1.0, 1.0, 1.0);
    test::circularBinary(c, 1.0, 1.0, 1.0);
    def.step(a.view(), 0.1);
    verlet.step(c.view(), 0.1);
    EXPECT_EQ(a.pos, c.pos);
    EXPECT_EQ(a.vel, c.vel);
}

//=======================Two-body orbit=========================

TEST(SimulationOrbit, CircularBinaryReturnsAfterOnePeriod) {
    const double G = 1.0;
    Simulation sim = makeSim(IntegratorMethod::VERLET, G, 0.0);
    Bodies b;
    const double T = test::circularBinary(b, G, 1.0, 1.0);
    const Bodies start = b;

    const int steps = 2000;
    run(sim, b, T / steps, steps);

    for(size_t i = 0; i < 2; ++i) {
        EXPECT_VEC3_NEAR(b.pos[i], start.pos[i], 1e-4);
        EXPECT_VEC3_NEAR(b.vel[i], start.vel[i], 1e-4);
    }
}

TEST(SimulationOrbit, CircularBinaryKeepsItsRadius) {
    const double G = 1.0, r = 1.0;
    Simulation sim = makeSim(IntegratorMethod::VERLET, G, 0.0);
    Bodies b;
    const double T = test::circularBinary(b, G, 1.0, r);

    const int steps = 2000;
    for(int i = 0; i < 3 * steps; ++i) {
        sim.step(b.view(), T / steps);
        ASSERT_NEAR(glm::length(b.pos[1] - b.pos[0]), r, 1e-5) << "step " << i;
    }
}

TEST(SimulationOrbit, EccentricOrbitConservesEnergyAndAngularMomentum) {
    //Bound orbit with eccentricity around 0.6, integrated over ~10 periods
    const double G = 1.0;
    Simulation sim = makeSim(IntegratorMethod::VERLET, G, 0.0);
    Bodies b;
    b.add({-0.5, 0.0, 0.0}, {0.0, -0.4, 0.0}, 1.0);
    b.add({ 0.5, 0.0, 0.0}, {0.0,  0.4, 0.0}, 1.0);

    const double     E0 = sim.totalEnergy(b.view());
    const glm::dvec3 L0 = test::totalAngularMomentum(b);
    ASSERT_LT(E0, 0.0);

    double worstE = 0.0;
    for(int i = 0; i < 40'000; ++i) {
        sim.step(b.view(), 5e-4);
        worstE = std::max(worstE, std::abs(sim.totalEnergy(b.view()) - E0));
    }

    EXPECT_LT(worstE / std::abs(E0), 1e-4);
    EXPECT_VEC3_NEAR(test::totalAngularMomentum(b), L0, 1e-12);
}

TEST(SimulationOrbit, KeplersThirdLawPeriod) {
    //Measure the time between successive crossings of the x axis for two separations
    //and check T^2 / r^3 = 4 pi^2 / (G M)
    auto measurePeriod = [](double r) {
        const double G = 1.0;
        Simulation sim = makeSim(IntegratorMethod::VERLET, G, 0.0);
        Bodies b;
        const double T = test::circularBinary(b, G, 1.0, r);
        const double dt = T / 20'000;

        double t = 0.0, prevY = b.pos[1].y;
        int crossings = 0;
        double firstCrossing = 0.0;
        while(crossings < 3) {
            sim.step(b.view(), dt);
            t += dt;
            const double y = b.pos[1].y;
            if(prevY < 0.0 && y >= 0.0) {
                //linear interpolation of the crossing time
                const double tc = t - dt * y / (y - prevY);
                if(crossings == 0) firstCrossing = tc;
                ++crossings;
                if(crossings == 3) return (tc - firstCrossing) / 2.0;
            }
            prevY = y;
        }
        return 0.0;
    };

    for(double r : {1.0, 2.0}) {
        const double T = measurePeriod(r);
        EXPECT_NEAR(T * T / (r * r * r), 4.0 * M_PI * M_PI / 2.0, 1e-4) << "r = " << r;
    }
}

//=======================Many-body invariants=========================

TEST(SimulationInvariants, DirectMomentumIsConserved) {
    Simulation sim = makeSim(IntegratorMethod::VERLET, 1.0, 0.05);
    Bodies b = randomCluster(64, 8);
    const glm::dvec3 p0 = test::totalMomentum(b);

    run(sim, b, 1e-3, 500);
    EXPECT_VEC3_NEAR(test::totalMomentum(b), p0, 1e-11);
}

TEST(SimulationInvariants, CenterOfMassMovesUniformly) {
    Simulation sim = makeSim(IntegratorMethod::VERLET, 1.0, 0.05);
    Bodies b = randomCluster(64, 9);
    double M = 0.0;
    for(double m : b.m) M += m;
    const glm::dvec3 c0 = test::centerOfMass(b);
    const glm::dvec3 vcm = test::totalMomentum(b) / M;

    const double dt = 1e-3;
    const int steps = 500;
    run(sim, b, dt, steps);
    EXPECT_VEC3_NEAR(test::centerOfMass(b), c0 + vcm * (steps * dt), 1e-11);
}

TEST(SimulationInvariants, DirectEnergyIsConserved) {
    Simulation sim = makeSim(IntegratorMethod::VERLET, 1.0, 0.05);
    Bodies b = randomCluster(64, 10);
    const double E0 = sim.totalEnergy(b.view());

    run(sim, b, 1e-4, 2000);
    EXPECT_LT(std::abs(sim.totalEnergy(b.view()) - E0) / std::abs(E0), 1e-4);
}

TEST(SimulationInvariants, BarnesHutTracksDirectSimulation) {
    //Same initial state, after a short integration both algorithms must agree closely.
    //Weak G keeps the cluster from collapsing, so the difference measures force error, not chaos.
    Bodies direct = randomCluster(256, 12);
    Bodies tree   = direct;

    Simulation simDirect = makeSim(IntegratorMethod::VERLET, 0.01, 0.05);
    Simulation simTree   = makeSim(IntegratorMethod::VERLET, 0.01, 0.05, Gravity::Algorithm::BARNES_HUT);
    simTree.gravity().setTheta(0.3);

    run(simDirect, direct, 1e-3, 100);
    run(simTree,   tree,   1e-3, 100);

    double worst = 0.0;
    for(size_t i = 0; i < direct.size(); ++i) {
        worst = std::max(worst, glm::length(direct.pos[i] - tree.pos[i]));
    }
    EXPECT_LT(worst, 1e-4);

    const double Ed = simDirect.totalEnergy(direct.view());
    const double Et = simTree.totalEnergy(tree.view());
    EXPECT_LT(std::abs(Ed - Et) / std::abs(Ed), 1e-4);
}

TEST(SimulationInvariants, BarnesHutEnergyIsApproximatelyConserved) {
    Simulation sim = makeSim(IntegratorMethod::VERLET, 0.01, 0.05, Gravity::Algorithm::BARNES_HUT);
    Bodies b = randomCluster(512, 13);
    const double E0 = sim.totalEnergy(b.view());

    run(sim, b, 1e-3, 200);
    EXPECT_LT(std::abs(sim.totalEnergy(b.view()) - E0) / std::abs(E0), 1e-3);
}

TEST(SimulationInvariants, EulerDriftsMoreThanVerlet) {
    auto energyError = [](IntegratorMethod m) {
        Simulation sim = makeSim(m, 1.0, 0.0);
        Bodies b;
        const double T = test::circularBinary(b, 1.0, 1.0, 1.0);
        const double E0 = sim.totalEnergy(b.view());
        run(sim, b, T / 500, 500);
        return std::abs(sim.totalEnergy(b.view()) - E0);
    };
    EXPECT_GT(energyError(IntegratorMethod::EULER), 100.0 * energyError(IntegratorMethod::VERLET));
}

//=======================Move semantics=========================

TEST(Simulation, MovedSimulationKeepsSettingsAndState) {
    Simulation original = makeSim(IntegratorMethod::VERLET, 1.0, 0.0);
    Bodies a, b;
    test::circularBinary(a, 1.0, 1.0, 1.0);
    test::circularBinary(b, 1.0, 1.0, 1.0);

    Simulation reference = makeSim(IntegratorMethod::VERLET, 1.0, 0.0);
    original.step(a.view(), 0.01);
    reference.step(b.view(), 0.01);

    Simulation moved = std::move(original);
    EXPECT_DOUBLE_EQ(moved.gravity().G(), 1.0);
    moved.step(a.view(), 0.01);
    reference.step(b.view(), 0.01);

    EXPECT_EQ(a.pos, b.pos);
    EXPECT_EQ(a.vel, b.vel);
}
