#include <gtest/gtest.h>

#include "Motus/BarnesHutTree.h"
#include "Motus/Gravity.h"

#include "TestUtils.h"

#include <algorithm>

using namespace motus;
using motus::test::Bodies;

namespace {

Gravity barnesHut(double G, double eps, double theta) {
    Gravity g(G, eps);
    g.setAlgorithm(Gravity::Algorithm::BARNES_HUT);
    g.setTheta(theta);
    return g;
}

Accel computeAccel(const Gravity& g, Bodies& b) {
    Accel a(b.size(), glm::dvec3{0.0});
    g.compute(b.view(), a);
    return a;
}

//Relative error of every body's acceleration against the exact O(N^2) result
std::vector<double> relativeErrors(Bodies& b, double G, double eps, double theta) {
    Accel bh  = computeAccel(barnesHut(G, eps, theta), b);
    Accel ref = test::referenceAccel(b, G, eps);

    std::vector<double> err(b.size());
    for(size_t i = 0; i < b.size(); ++i) {
        err[i] = glm::length(bh[i] - ref[i]) / glm::length(ref[i]);
    }
    return err;
}

double rms(const std::vector<double>& v) {
    double s = 0.0;
    for(double x : v) s += x * x;
    return std::sqrt(s / v.size());
}

}

TEST(BarnesHut, ThetaZeroIsExact) {
    //With theta = 0 no node is ever approximated, so the tree must reproduce direct summation
    const double G = 1.0, eps = 0.01;
    Bodies b = test::randomCube(300, 17);

    for(double e : relativeErrors(b, G, eps, 0.0)) EXPECT_LT(e, 1e-12);
}

TEST(BarnesHut, TwoBodiesMatchNewtonsLaw) {
    const double G = 2.0, eps = 0.1, m1 = 3.0, m2 = 5.0, r = 4.0;
    Bodies b;
    b.add({0, 0, 0}, {}, m1);
    b.add({r, 0, 0}, {}, m2);
    Accel a = computeAccel(barnesHut(G, eps, 0.5), b);

    const double d2 = r * r + eps * eps;
    const double k  = G / (d2 * std::sqrt(d2)) * r;
    EXPECT_VEC3_NEAR(a[0], glm::dvec3( k * m2, 0.0, 0.0), 1e-14);
    EXPECT_VEC3_NEAR(a[1], glm::dvec3(-k * m1, 0.0, 0.0), 1e-14);
}

TEST(BarnesHut, SingleBodyFeelsNoForce) {
    Bodies b;
    b.add({1, 2, 3}, {}, 5.0);
    EXPECT_VEC3_DOUBLE_EQ(computeAccel(barnesHut(1.0, 0.01, 0.5), b)[0], glm::dvec3(0.0));
}

TEST(BarnesHut, DefaultThetaIsAccurate) {
    //theta = 0.5 with monopole nodes is expected to give errors around half a percent
    Bodies b = test::randomCube(2000, 99);
    auto err = relativeErrors(b, 1.0, 0.01, 0.5);

    EXPECT_LT(rms(err), 1e-2);
    EXPECT_LT(*std::ranges::max_element(err), 5e-2);
}

TEST(BarnesHut, ErrorDecreasesWithSmallerTheta) {
    Bodies b = test::randomCube(1000, 5);
    const double coarse = rms(relativeErrors(b, 1.0, 0.01, 1.0));
    const double medium = rms(relativeErrors(b, 1.0, 0.01, 0.5));
    const double fine   = rms(relativeErrors(b, 1.0, 0.01, 0.2));

    EXPECT_LT(medium, coarse);
    EXPECT_LT(fine,   medium);
}

namespace {

//Cluster in [0, 1]^3 plus a probe at the opposite corner of the root cube. The root splits at ~500,
//so the whole cluster lands in octant 0 and the probe alone in octant 7.
struct ClusterAndProbe {
    Bodies      b;
    glm::dvec3  probe {1000.0, 1000.0, 1000.0};
    glm::dvec3  clusterCm;
    double      clusterMass = 0.0;

    ClusterAndProbe() {
        b = test::randomCube(50, 23, 0.0, 1.0);
        clusterCm = test::centerOfMass(b);
        for(double m : b.m) clusterMass += m;
        b.add(probe, {}, 1.0);
    }

    [[nodiscard]] glm::dvec3 pointMassAccel(double G) const {
        const glm::dvec3 r = clusterCm - probe;
        return G * clusterMass * r / std::pow(glm::length(r), 3.0);
    }
};

}

TEST(BarnesHut, DistantClusterActsAsPointMassAtCenterOfMass) {
    //Checks the mass and center of mass aggregation of the tree: with theta = 0.5 the root opens,
    //and the probe sees the octant holding the cluster as one node with the total mass at the true cm
    ClusterAndProbe s;
    const glm::dvec3 a = computeAccel(barnesHut(1.0, 0.0, 0.5), s.b).back();
    const glm::dvec3 expected = s.pointMassAccel(1.0);
    EXPECT_VEC3_NEAR(a, expected, 1e-12 * glm::length(expected));
}

TEST(BarnesHut, NodeContainingBodyIsNeverApproximated) {
    //A body must not interact with the aggregate of a node that contains itself (that is a self force).
    //Here the root's size/d is ~0.59, so any theta above that would approximate the root, probe included.
    ClusterAndProbe s;
    const glm::dvec3 a = computeAccel(barnesHut(1.0, 0.0, 0.7), s.b).back();
    const glm::dvec3 expected = s.pointMassAccel(1.0);
    EXPECT_VEC3_NEAR(a, expected, 1e-2 * glm::length(expected));
}

TEST(BarnesHut, HandlesClusteredDistributions) {
    //Two dense, far apart clumps force a deep, unbalanced tree
    Bodies b = test::randomCube(300, 1, -0.01, 0.01);
    Bodies other = test::randomCube(300, 2, 9.99, 10.01);
    for(size_t i = 0; i < other.size(); ++i) b.add(other.pos[i], other.vel[i], other.m[i]);

    for(double e : relativeErrors(b, 1.0, 1e-4, 0.0)) EXPECT_LT(e, 1e-10);
    EXPECT_LT(rms(relativeErrors(b, 1.0, 1e-4, 0.5)), 1e-2);
}

TEST(BarnesHut, BodiesOnNodeBoundariesAreInserted) {
    //A regular lattice puts bodies exactly on octant split planes
    Bodies b;
    for(int x = 0; x < 4; ++x)
        for(int y = 0; y < 4; ++y)
            for(int z = 0; z < 4; ++z)
                b.add({double(x), double(y), double(z)}, {}, 1.0);

    for(double e : relativeErrors(b, 1.0, 0.01, 0.0)) EXPECT_LT(e, 1e-12);
}

TEST(BarnesHut, AccumulatesIntoOutput) {
    //Same contract as direct summation: forces add into the buffer instead of overwriting it
    Bodies b;
    b.add({0, 0, 0}, {}, 1.0);
    b.add({1, 0, 0}, {}, 1.0);

    Accel a(2, glm::dvec3{10.0, 20.0, 30.0});
    barnesHut(1.0, 0.0, 0.5).compute(b.view(), a);

    EXPECT_VEC3_NEAR(a[0], glm::dvec3(11.0, 20.0, 30.0), 1e-14);
    EXPECT_VEC3_NEAR(a[1], glm::dvec3( 9.0, 20.0, 30.0), 1e-14);
}

TEST(BarnesHut, TreeTraversalReturnsForceNotAcceleration) {
    //traverseNode returns F = m_i * a_i; Gravity::barnesHut divides by the mass afterwards
    Bodies b;
    b.add({0, 0, 0}, {}, 3.0);
    b.add({2, 0, 0}, {}, 5.0);

    BarnesHutTree tree(b.view());
    tree.build();

    const double G = 1.0, eps = 0.0;
    const glm::dvec3 f = tree.traverseNode(0, 0, 0.5, G, eps);
    EXPECT_VEC3_NEAR(f, glm::dvec3(G * 3.0 * 5.0 / 4.0, 0.0, 0.0), 1e-14);
}
