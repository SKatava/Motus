#include <gtest/gtest.h>

#include "Motus/Euler.h"
#include "Motus/Verlet.h"

#include "TestUtils.h"

using namespace motus;
using motus::test::Bodies;

namespace {

//a = -x per component, so x(t) = x0 cos t + v0 sin t
void harmonic(BodiesView b, AccelView out) {
    for(size_t i = 0; i < b.size(); ++i) out[i] += -b.pos[i];
}

ForcesFn constantAccel(glm::dvec3 g) {
    return [g](BodiesView b, AccelView out) {
        for(size_t i = 0; i < b.size(); ++i) out[i] += g;
    };
}

Bodies oscillator() {
    Bodies b;
    b.add({1.0, 0.0, 0.5}, {0.0, 1.0, -0.5}, 1.0);
    return b;
}

glm::dvec3 harmonicExactPos(const Bodies& start, double t) {
    return start.pos[0] * std::cos(t) + start.vel[0] * std::sin(t);
}

//Position error of the oscillator at t = 2*pi after integrating with the given step
template<typename I>
double oscillatorError(double dt) {
    I integrator;
    Bodies b = oscillator();
    const Bodies start = b;
    const int steps = static_cast<int>(std::round(2.0 * M_PI / dt));
    for(int i = 0; i < steps; ++i) integrator.step(b.view(), harmonic, dt);
    return glm::length(b.pos[0] - harmonicExactPos(start, steps * dt));
}

double oscillatorEnergy(const Bodies& b) {
    return 0.5 * glm::dot(b.vel[0], b.vel[0]) + 0.5 * glm::dot(b.pos[0], b.pos[0]);
}

}

//=======================Euler=========================

TEST(Euler, SingleStepUsesOldStateForEverything) {
    //Explicit Euler: x1 = x0 + v0 dt, v1 = v0 + a(x0) dt
    Euler euler;
    Bodies b = oscillator();
    const Bodies start = b;
    const double dt = 0.1;

    euler.step(b.view(), harmonic, dt);

    EXPECT_VEC3_DOUBLE_EQ(b.pos[0], start.pos[0] + start.vel[0] * dt);
    EXPECT_VEC3_DOUBLE_EQ(b.vel[0], start.vel[0] - start.pos[0] * dt);
}

TEST(Euler, FreeParticleMovesInStraightLine) {
    Euler euler;
    Bodies b;
    b.add({1, 2, 3}, {0.5, -1.0, 2.0}, 1.0);

    for(int i = 0; i < 100; ++i) euler.step(b.view(), constantAccel({0, 0, 0}), 0.01);

    EXPECT_VEC3_NEAR(b.pos[0], glm::dvec3(1.5, 1.0, 5.0), 1e-12);
    EXPECT_VEC3_DOUBLE_EQ(b.vel[0], glm::dvec3(0.5, -1.0, 2.0));
}

TEST(Euler, EvaluatesForcesOncePerStep) {
    Euler euler;
    Bodies b = oscillator();
    int calls = 0;
    ForcesFn f = [&](BodiesView v, AccelView out) { ++calls; harmonic(v, out); };

    for(int i = 0; i < 10; ++i) euler.step(b.view(), f, 0.01);
    EXPECT_EQ(calls, 10);
}

TEST(Euler, IsFirstOrderAccurate) {
    const double e1 = oscillatorError<Euler>(1e-3);
    const double e2 = oscillatorError<Euler>(5e-4);
    EXPECT_NEAR(e1 / e2, 2.0, 0.1);
}

//=======================Verlet=========================

TEST(Verlet, FirstStepMatchesVelocityVerletFormula) {
    //x1 = x0 + v0 dt + a0 dt^2 / 2,  v1 = v0 + (a0 + a1) dt / 2
    Verlet verlet;
    Bodies b = oscillator();
    const glm::dvec3 x0 = b.pos[0], v0 = b.vel[0];
    const double dt = 0.1;

    verlet.step(b.view(), harmonic, dt);

    const glm::dvec3 a0 = -x0;
    const glm::dvec3 x1 = x0 + v0 * dt + 0.5 * a0 * dt * dt;
    const glm::dvec3 a1 = -x1;
    const glm::dvec3 v1 = v0 + 0.5 * (a0 + a1) * dt;

    EXPECT_VEC3_DOUBLE_EQ(b.pos[0], x1);
    EXPECT_VEC3_DOUBLE_EQ(b.vel[0], v1);
}

TEST(Verlet, MatchesReferenceImplementationOverManySteps) {
    //Guards the acceleration carried between steps (the m_acc / m_accNew swap)
    Verlet verlet;
    Bodies b = oscillator();
    glm::dvec3 x = b.pos[0], v = b.vel[0], a = -x;
    const double dt = 0.05;

    for(int i = 0; i < 200; ++i) {
        verlet.step(b.view(), harmonic, dt);

        x += v * dt + 0.5 * a * dt * dt;
        const glm::dvec3 aNew = -x;
        v += 0.5 * (a + aNew) * dt;
        a = aNew;

        ASSERT_NEAR(glm::length(b.pos[0] - x), 0.0, 1e-13) << "step " << i;
        ASSERT_NEAR(glm::length(b.vel[0] - v), 0.0, 1e-13) << "step " << i;
    }
}

TEST(Verlet, ConstantAccelerationIsExact) {
    //Velocity Verlet integrates a parabola exactly: x = x0 + v0 t + g t^2 / 2
    Verlet verlet;
    const glm::dvec3 g {0.0, 0.0, -9.81};
    const glm::dvec3 x0 {0.0, 0.0, 100.0}, v0 {3.0, 1.0, 20.0};
    Bodies b;
    b.add(x0, v0, 2.0);

    const double dt = 0.01;
    const int steps = 500;
    for(int i = 0; i < steps; ++i) verlet.step(b.view(), constantAccel(g), dt);

    const double t = steps * dt;
    EXPECT_VEC3_NEAR(b.pos[0], x0 + v0 * t + 0.5 * g * t * t, 1e-9);
    EXPECT_VEC3_NEAR(b.vel[0], v0 + g * t, 1e-10);
}

TEST(Verlet, EvaluatesForcesOncePerStepAfterInit) {
    Verlet verlet;
    Bodies b = oscillator();
    int calls = 0;
    ForcesFn f = [&](BodiesView v, AccelView out) { ++calls; harmonic(v, out); };

    for(int i = 0; i < 10; ++i) verlet.step(b.view(), f, 0.01);
    EXPECT_EQ(calls, 11);
}

TEST(Verlet, IsSecondOrderAccurate) {
    const double e1 = oscillatorError<Verlet>(1e-2);
    const double e2 = oscillatorError<Verlet>(5e-3);
    EXPECT_NEAR(e1 / e2, 4.0, 0.2);
}

TEST(Verlet, MoreAccurateThanEulerAtSameStep) {
    EXPECT_LT(oscillatorError<Verlet>(1e-2) * 100.0, oscillatorError<Euler>(1e-2));
}

TEST(Verlet, OscillatorEnergyStaysBounded) {
    //Symplectic: energy error oscillates at O(dt^2) instead of drifting
    Verlet verlet;
    Bodies b = oscillator();
    const double E0 = oscillatorEnergy(b);
    const double dt = 0.05;

    double worst = 0.0;
    for(int i = 0; i < 100'000; ++i) {
        verlet.step(b.view(), harmonic, dt);
        worst = std::max(worst, std::abs(oscillatorEnergy(b) - E0) / E0);
    }
    EXPECT_LT(worst, dt * dt);
}

TEST(Verlet, ReinitialisesWhenBodyCountChanges) {
    const double dt = 0.1;
    Verlet reused;
    Bodies one = oscillator();
    reused.step(one.view(), harmonic, dt);

    Bodies two;
    two.add({1, 0, 0}, {0, 1, 0}, 1.0);
    two.add({0, 2, 0}, {1, 0, 0}, 1.0);
    Bodies twoFresh = two;

    reused.step(two.view(), harmonic, dt);
    Verlet fresh;
    fresh.step(twoFresh.view(), harmonic, dt);

    for(size_t i = 0; i < 2; ++i) {
        EXPECT_VEC3_DOUBLE_EQ(two.pos[i], twoFresh.pos[i]);
        EXPECT_VEC3_DOUBLE_EQ(two.vel[i], twoFresh.vel[i]);
    }
}

//=======================Both=========================

TEST(Integrators, ForceBufferIsZeroedBeforeEvaluation) {
    //Each evaluation must start from zero, otherwise accumulating forces would pile up across steps
    Bodies be, bv;
    be.add({0, 0, 0}, {0, 0, 0}, 1.0);
    bv.add({0, 0, 0}, {0, 0, 0}, 1.0);
    Euler euler;
    Verlet verlet;
    const glm::dvec3 g {1.0, 0.0, 0.0};

    for(int i = 0; i < 10; ++i) {
        euler.step(be.view(), constantAccel(g), 1.0);
        verlet.step(bv.view(), constantAccel(g), 1.0);
    }
    EXPECT_VEC3_DOUBLE_EQ(be.vel[0], 10.0 * g);
    EXPECT_VEC3_DOUBLE_EQ(bv.vel[0], 10.0 * g);
}
