#include <gtest/gtest.h>

#include "Motus/Force.h"
#include "Motus/Gravity.h"

#include "TestUtils.h"

using namespace motus;
using motus::test::Bodies;

namespace {

Accel computeAccel(const Gravity& g, Bodies& b) {
    Accel a(b.size(), glm::dvec3{0.0});
    g.compute(b.view(), a);
    return a;
}

}

//=======================Parameters=========================

TEST(Gravity, DefaultParameters) {
    Gravity g;
    EXPECT_DOUBLE_EQ(g.G(), 6.674e-11);
    EXPECT_DOUBLE_EQ(g.eps(), 1e-3);
    EXPECT_DOUBLE_EQ(g.theta(), 0.5);
}

TEST(Gravity, ConstructorSetsGAndEps) {
    Gravity g(2.5, 0.1);
    EXPECT_DOUBLE_EQ(g.G(), 2.5);
    EXPECT_DOUBLE_EQ(g.eps(), 0.1);
}

TEST(Gravity, SettersAndResets) {
    Gravity g;
    g.setG(3.0);
    g.setEps(0.25);
    g.setTheta(0.8);
    EXPECT_DOUBLE_EQ(g.G(), 3.0);
    EXPECT_DOUBLE_EQ(g.eps(), 0.25);
    EXPECT_DOUBLE_EQ(g.theta(), 0.8);

    g.resetG();
    g.resetEps();
    g.resetTheta();
    EXPECT_DOUBLE_EQ(g.G(), 6.674e-11);
    EXPECT_DOUBLE_EQ(g.eps(), 1e-3);
    EXPECT_DOUBLE_EQ(g.theta(), 0.5);
}

//=======================Direct Newton=========================

TEST(GravityDirect, TwoBodyAccelerationMatchesNewtonsLaw) {
    //a1 = G*m2 / (r^2 + eps^2)^(3/2) * r, pointing towards the other body
    const double G = 2.0, eps = 0.1, m1 = 3.0, m2 = 5.0, r = 4.0;
    Gravity g(G, eps);

    Bodies b;
    b.add({0.0, 0.0, 0.0}, {}, m1);
    b.add({r,   0.0, 0.0}, {}, m2);
    Accel a = computeAccel(g, b);

    const double d2 = r * r + eps * eps;
    const double k  = G / (d2 * std::sqrt(d2)) * r;

    EXPECT_VEC3_DOUBLE_EQ(a[0], glm::dvec3( k * m2, 0.0, 0.0));
    EXPECT_VEC3_DOUBLE_EQ(a[1], glm::dvec3(-k * m1, 0.0, 0.0));
}

TEST(GravityDirect, InverseSquareFallOff) {
    //Without meaningful softening doubling the distance must quarter the acceleration
    Gravity g(1.0, 0.0);

    Bodies near, far;
    near.add({0, 0, 0}, {}, 1.0); near.add({1, 0, 0}, {}, 1.0);
    far.add ({0, 0, 0}, {}, 1.0); far.add ({2, 0, 0}, {}, 1.0);

    const double aNear = computeAccel(g, near)[0].x;
    const double aFar  = computeAccel(g, far)[0].x;

    EXPECT_DOUBLE_EQ(aNear, 1.0);
    EXPECT_DOUBLE_EQ(aFar,  0.25);
}

TEST(GravityDirect, AccelerationIsIndependentOfOwnMass) {
    Gravity g(1.0, 0.01);

    Bodies light, heavy;
    light.add({0, 0, 0}, {}, 1.0);    light.add({1, 2, 3}, {}, 7.0);
    heavy.add({0, 0, 0}, {}, 1000.0); heavy.add({1, 2, 3}, {}, 7.0);

    EXPECT_VEC3_DOUBLE_EQ(computeAccel(g, light)[0], computeAccel(g, heavy)[0]);
}

TEST(GravityDirect, ScalesLinearlyWithG) {
    Bodies b = test::randomCube(20, 7);
    Accel a1 = computeAccel(Gravity(1.0, 0.01), b);
    Accel a3 = computeAccel(Gravity(3.0, 0.01), b);

    for(size_t i = 0; i < b.size(); ++i) {
        EXPECT_VEC3_NEAR(a3[i], 3.0 * a1[i], 1e-12 * glm::length(a3[i]));
    }
}

TEST(GravityDirect, MatchesIndependentReferenceForManyBodies) {
    const double G = 1.3, eps = 0.05;
    Bodies b = test::randomCube(200, 42);

    Accel a   = computeAccel(Gravity(G, eps), b);
    Accel ref = test::referenceAccel(b, G, eps);

    for(size_t i = 0; i < b.size(); ++i) {
        EXPECT_VEC3_NEAR(a[i], ref[i], 1e-10 * glm::length(ref[i]));
    }
}

TEST(GravityDirect, NewtonsThirdLawNetForceIsZero) {
    Bodies b = test::randomCube(100, 3);
    Accel a = computeAccel(Gravity(1.0, 0.01), b);

    glm::dvec3 net {0.0};
    double scale = 0.0;
    for(size_t i = 0; i < b.size(); ++i) {
        net   += b.m[i] * a[i];
        scale += b.m[i] * glm::length(a[i]);
    }
    EXPECT_LT(glm::length(net), 1e-12 * scale);
}

TEST(GravityDirect, SymmetricConfigurationHasZeroForceAtCenter) {
    Gravity g(1.0, 0.01);
    Bodies b;
    b.add({0, 0, 0}, {}, 1.0);
    for(glm::dvec3 d : {glm::dvec3(1, 0, 0), glm::dvec3(-1, 0, 0), glm::dvec3(0, 1, 0),
                        glm::dvec3(0, -1, 0), glm::dvec3(0, 0, 1), glm::dvec3(0, 0, -1)}) {
        b.add(d, {}, 2.0);
    }
    EXPECT_VEC3_NEAR(computeAccel(g, b)[0], glm::dvec3(0.0), 1e-15);
}

TEST(GravityDirect, SingleBodyFeelsNoForce) {
    Bodies b;
    b.add({1, 2, 3}, {}, 5.0);
    EXPECT_VEC3_DOUBLE_EQ(computeAccel(Gravity(1.0, 0.01), b)[0], glm::dvec3(0.0));
}

TEST(GravityDirect, AccumulatesIntoOutput) {
    //Forces add to the buffer, so several forces can share one acceleration array
    Gravity g(1.0, 0.0);
    Bodies b;
    b.add({0, 0, 0}, {}, 1.0);
    b.add({1, 0, 0}, {}, 1.0);

    Accel a(2, glm::dvec3{10.0, 20.0, 30.0});
    g.compute(b.view(), a);

    EXPECT_VEC3_DOUBLE_EQ(a[0], glm::dvec3(11.0, 20.0, 30.0));
    EXPECT_VEC3_DOUBLE_EQ(a[1], glm::dvec3( 9.0, 20.0, 30.0));
}

TEST(GravityDirect, DoesNotModifyBodies) {
    Bodies b = test::randomCube(10, 1);
    Bodies copy = b;
    computeAccel(Gravity(1.0, 0.01), b);
    EXPECT_EQ(b.pos, copy.pos);
    EXPECT_EQ(b.vel, copy.vel);
    EXPECT_EQ(b.m,   copy.m);
}

//=======================Potential energy=========================

TEST(GravityPotential, TwoBodyMatchesFormula) {
    const double G = 2.0, eps = 0.1, m1 = 3.0, m2 = 5.0, r = 4.0;
    Bodies b;
    b.add({0, 0, 0}, {}, m1);
    b.add({0, r, 0}, {}, m2);

    EXPECT_DOUBLE_EQ(Gravity(G, eps).potentialEnergy(b.view()),
                     -G * m1 * m2 / std::sqrt(r * r + eps * eps));
}

TEST(GravityPotential, CountsEveryPairExactlyOnce) {
    //Equilateral triangle with side 1: U = -3 * G m^2
    const double h = std::sqrt(3.0) / 2.0;
    Bodies b;
    b.add({0.0, 0.0, 0.0}, {}, 2.0);
    b.add({1.0, 0.0, 0.0}, {}, 2.0);
    b.add({0.5, h,   0.0}, {}, 2.0);

    EXPECT_NEAR(Gravity(1.0, 0.0).potentialEnergy(b.view()), -12.0, 1e-12);
}

TEST(GravityPotential, IsConsistentWithAcceleration) {
    //a_i = -(1/m_i) dU/dx_i, checked with a central finite difference
    const double G = 1.0, eps = 0.05, h = 1e-6;
    Gravity g(G, eps);
    Bodies b = test::randomCube(8, 11);
    Accel a = computeAccel(g, b);

    for(size_t i = 0; i < b.size(); ++i) {
        for(int axis = 0; axis < 3; ++axis) {
            Bodies plus = b, minus = b;
            plus.pos[i][axis]  += h;
            minus.pos[i][axis] -= h;
            const double dU = (g.potentialEnergy(plus.view()) - g.potentialEnergy(minus.view())) / (2.0 * h);
            EXPECT_NEAR(a[i][axis], -dU / b.m[i], 1e-5 * (1.0 + std::abs(a[i][axis])));
        }
    }
}

TEST(GravityPotential, SingleBodyIsZero) {
    Bodies b;
    b.add({1, 1, 1}, {}, 3.0);
    EXPECT_DOUBLE_EQ(Gravity(1.0, 0.1).potentialEnergy(b.view()), 0.0);
}

//=======================ForceSet=========================

TEST(ForceSet, FlagOperators) {
    EXPECT_TRUE (has(Forces::Gravity, Forces::Gravity));
    EXPECT_FALSE(has(Forces::None,    Forces::Gravity));
    EXPECT_TRUE (has(Forces::None | Forces::Gravity, Forces::Gravity));
    EXPECT_EQ   (Forces::Gravity & Forces::None, Forces::None);
}

TEST(ForceSet, GravityEnabledByDefaultAndMatchesGravity) {
    ForceSet fs;
    fs.gravity = Gravity(1.0, 0.01);
    Bodies b = test::randomCube(10, 5);

    Accel viaSet(b.size(), glm::dvec3{0.0});
    fs(b.view(), viaSet);
    Accel direct = computeAccel(fs.gravity, b);

    for(size_t i = 0; i < b.size(); ++i) EXPECT_VEC3_DOUBLE_EQ(viaSet[i], direct[i]);
    EXPECT_DOUBLE_EQ(fs.potentialEnergy(b.view()), fs.gravity.potentialEnergy(b.view()));
}

TEST(ForceSet, DisabledGravityProducesNothing) {
    ForceSet fs;
    fs.enabled = Forces::None;
    fs.gravity = Gravity(1.0, 0.01);
    Bodies b = test::randomCube(10, 5);

    Accel a(b.size(), glm::dvec3{0.0});
    fs(b.view(), a);

    for(const auto& v : a) EXPECT_VEC3_DOUBLE_EQ(v, glm::dvec3(0.0));
    EXPECT_DOUBLE_EQ(fs.potentialEnergy(b.view()), 0.0);
}
