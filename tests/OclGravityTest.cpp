#include <gtest/gtest.h>

#include "Motus/Gravity.h"
#include "Motus/ocl/Handle.h"

#include "TestUtils.h"

#include <algorithm>
#include <cmath>

using namespace motus;
using motus::test::Bodies;

//GPU kernels work in float on coordinates relative to a center, so they are compared against the
//double reference with float-level tolerances. Every test skips cleanly on machines without an OpenCL GPU.

namespace {

using Algo = Gravity::Algorithm;

//Returns false when no usable OpenCL GPU is present (Device/kernel creation throws)
bool makeGpu(Gravity& g, Algo algo) {
    try {
        g.setAlgorithm(algo);
        return true;
    } catch(const ocl::Error&) {
        return false;
    }
}

#define REQUIRE_GPU(g, algo)                                        \
    do {                                                            \
        if(!makeGpu((g), (algo))) GTEST_SKIP() << "no OpenCL GPU";  \
    } while(0)

Accel computeAccel(const Gravity& g, Bodies& b) {
    Accel a(b.size(), glm::dvec3{0.0});
    g.compute(b.view(), a);
    return a;
}

std::vector<double> relativeErrors(const Accel& a, const Accel& ref) {
    std::vector<double> err(a.size());
    for(size_t i = 0; i < a.size(); ++i) {
        err[i] = glm::length(a[i] - ref[i]) / glm::length(ref[i]);
    }
    return err;
}

double maxOf(const std::vector<double>& v) { return v.empty() ? 0.0 : *std::max_element(v.begin(), v.end()); }

double rms(const std::vector<double>& v) {
    double s = 0.0;
    for(double x : v) s += x * x;
    return std::sqrt(s / v.size());
}

bool allFinite(const Accel& a) {
    return std::all_of(a.begin(), a.end(), [](const glm::dvec3& v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    });
}

Bodies shifted(Bodies b, glm::dvec3 offset) {
    for(auto& p : b.pos) p += offset;
    return b;
}

//Float has ~7 significant digits; summing thousands of terms costs another order of magnitude or so
constexpr double kFloatTol = 1e-4;

}

//------------------------------------------------------------------------------------------------
// DIRECT_GPU
//------------------------------------------------------------------------------------------------

TEST(GravityDirectGpu, MatchesIndependentReference) {
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::DIRECT_GPU);

    Bodies b = test::randomCube(1000, 3);
    EXPECT_LT(maxOf(relativeErrors(computeAccel(g, b), test::referenceAccel(b, G, eps))), kFloatTol);
}

TEST(GravityDirectGpu, FarFromOriginKeepsFloatPrecision) {
    //A unit cube placed at 1e6: absolute float coordinates would lose all relative precision,
    //so this only passes if the kernel works on coordinates relative to the cluster
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::DIRECT_GPU);

    Bodies b = shifted(test::randomCube(1000, 4), {1e6, -1e6, 1e6});
    EXPECT_LT(maxOf(relativeErrors(computeAccel(g, b), test::referenceAccel(b, G, eps))), kFloatTol);
}

TEST(GravityDirectGpu, NewtonsThirdLawNetForceIsZero) {
    //Sum of m_i a_i vanishes for pairwise forces; the GPU computes each pair twice (i->j and j->i)
    //independently, so it only cancels to float precision
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::DIRECT_GPU);

    Bodies b = test::randomCube(500, 5);
    Accel a = computeAccel(g, b);

    glm::dvec3 net {0.0};
    double scale = 0.0;
    for(size_t i = 0; i < b.size(); ++i) {
        net   += b.m[i] * a[i];
        scale += b.m[i] * glm::length(a[i]);
    }
    EXPECT_LT(glm::length(net), 1e-5 * scale);
}

TEST(GravityDirectGpu, AccumulatesIntoOutput) {
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::DIRECT_GPU);

    Bodies b = test::randomCube(200, 6);
    Accel ref = computeAccel(g, b);

    Accel a(b.size(), glm::dvec3{1.0, -2.0, 3.0});
    g.compute(b.view(), a);
    for(size_t i = 0; i < b.size(); ++i) {
        EXPECT_VEC3_NEAR(a[i], ref[i] + glm::dvec3(1.0, -2.0, 3.0), 1e-12);
    }
}

TEST(GravityDirectGpu, SingleBodyFeelsNoForce) {
    Gravity g(1.0, 0.01);
    REQUIRE_GPU(g, Algo::DIRECT_GPU);

    Bodies b;
    b.add({1, 2, 3}, {}, 5.0);
    Accel a = computeAccel(g, b);
    EXPECT_VEC3_DOUBLE_EQ(a[0], glm::dvec3(0.0));
}

TEST(GravityDirectGpu, ZeroBodies) {
    Gravity g(1.0, 0.01);
    REQUIRE_GPU(g, Algo::DIRECT_GPU);

    Bodies b;
    Accel a;
    g.compute(b.view(), a);
    EXPECT_TRUE(a.empty());
}

//------------------------------------------------------------------------------------------------
// BARNES_HUT_GPU (group traversal: 64 bodies per work-group share one tree walk)
//------------------------------------------------------------------------------------------------

TEST(BarnesHutGpu, ThetaZeroMatchesReferenceForAnyGroupFill) {
    //theta = 0 never accepts a node, so only float rounding remains. The sizes cover a single partial
    //group, an exact multiple of the group size (64) and one body past it, where the last group is short
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);
    g.setTheta(0.0);

    for(size_t n : {2, 3, 63, 64, 65, 128, 129, 1000}) {
        Bodies b = test::randomCube(n, 7 + static_cast<unsigned>(n));
        EXPECT_LT(maxOf(relativeErrors(computeAccel(g, b), test::referenceAccel(b, G, eps))), kFloatTol) << "n = " << n;
    }
}

TEST(BarnesHutGpu, ThetaZeroSpanningManyListRounds) {
    //At theta = 0 every group must visit every leaf. With 20k bodies that is thousands of leaves,
    //more than the kernel's leaf list holds, so the walk has to stop, flush and resume several times.
    //A lost or repeated leaf across a flush shows up as a mismatch here.
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);
    g.setTheta(0.0);

    Bodies b = test::randomCube(20000, 8);
    std::vector<double> err = relativeErrors(computeAccel(g, b), test::referenceAccel(b, G, eps));
    EXPECT_LT(rms(err), 1e-5);
    EXPECT_LT(maxOf(err), 5e-4);
}

TEST(BarnesHutGpu, AtLeastAsAccurateAsCpuBarnesHut) {
    //Testing nodes against the whole group's bounding box is stricter than testing against one body,
    //so at the same theta the GPU error must not exceed the CPU tree's error (plus a float floor)
    const double G = 1.0, eps = 0.01;
    Bodies b = test::randomCube(4000, 9);
    Accel ref = test::referenceAccel(b, G, eps);

    for(double theta : {0.3, 0.5, 0.7}) {
        Gravity gpu(G, eps);
        REQUIRE_GPU(gpu, Algo::BARNES_HUT_GPU);
        gpu.setTheta(theta);

        Gravity cpu(G, eps);
        cpu.setAlgorithm(Algo::BARNES_HUT);
        cpu.setTheta(theta);

        const double gpuErr = rms(relativeErrors(computeAccel(gpu, b), ref));
        const double cpuErr = rms(relativeErrors(computeAccel(cpu, b), ref));
        EXPECT_LE(gpuErr, cpuErr + 1e-5) << "theta = " << theta;
        EXPECT_LT(gpuErr, 2e-2) << "theta = " << theta;
    }
}

TEST(BarnesHutGpu, ErrorDecreasesWithSmallerTheta) {
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);

    Bodies b = test::randomCube(4000, 10);
    Accel ref = test::referenceAccel(b, G, eps);

    double previous = 1e30;
    for(double theta : {0.9, 0.7, 0.5, 0.3}) {
        g.setTheta(theta);
        const double err = rms(relativeErrors(computeAccel(g, b), ref));
        EXPECT_LT(err, previous) << "theta = " << theta;
        previous = err;
    }
}

TEST(BarnesHutGpu, FarFromOriginKeepsFloatPrecision) {
    //Same cluster at the origin and at 1e6: relative coordinates make the GPU result independent of where it sits
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);
    g.setTheta(0.5);

    Bodies near = test::randomCube(3000, 11);
    Bodies far  = shifted(near, {1e6, -1e6, 1e6});

    const double nearErr = rms(relativeErrors(computeAccel(g, near), test::referenceAccel(near, G, eps)));
    const double farErr  = rms(relativeErrors(computeAccel(g, far),  test::referenceAccel(far,  G, eps)));
    EXPECT_LT(farErr, 1.1 * nearErr + 1e-5);
}

TEST(BarnesHutGpu, HandlesSeparatedClusters) {
    //Two dense clusters far apart: the Z-order curve jumps between them, so some groups mix bodies from both
    //and get a huge bounding box. Those groups must still be correct, just slower.
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);
    g.setTheta(0.5);

    Bodies b = test::randomCube(1500, 12, -1.0, 1.0);
    Bodies other = shifted(test::randomCube(1500, 13, -1.0, 1.0), {50.0, 0.0, 0.0});
    for(size_t i = 0; i < other.size(); ++i) b.add(other.pos[i], other.vel[i], other.m[i]);

    Accel a = computeAccel(g, b);
    ASSERT_TRUE(allFinite(a));
    EXPECT_LT(rms(relativeErrors(a, test::referenceAccel(b, G, eps))), 1e-2);
}

TEST(BarnesHutGpu, AccumulatesIntoOutput) {
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);

    Bodies b = test::randomCube(500, 14);
    Accel ref = computeAccel(g, b);

    Accel a(b.size(), glm::dvec3{1.0, -2.0, 3.0});
    g.compute(b.view(), a);
    for(size_t i = 0; i < b.size(); ++i) {
        EXPECT_VEC3_NEAR(a[i], ref[i] + glm::dvec3(1.0, -2.0, 3.0), 1e-12);
    }
}

TEST(BarnesHutGpu, SingleBodyFeelsNoForce) {
    Gravity g(1.0, 0.01);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);

    Bodies b;
    b.add({1, 2, 3}, {}, 5.0);
    Accel a = computeAccel(g, b);
    EXPECT_VEC3_DOUBLE_EQ(a[0], glm::dvec3(0.0));
}

TEST(BarnesHutGpu, ZeroBodies) {
    Gravity g(1.0, 0.01);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);

    Bodies b;
    Accel a;
    g.compute(b.view(), a);
    EXPECT_TRUE(a.empty());
}

TEST(BarnesHutGpu, AllBodiesAtOnePoint) {
    //Zero-size group box and a leaf deeper than the bucket size: every pair has zero separation,
    //so the softened force is exactly zero and nothing may turn into NaN
    Gravity g(1.0, 0.01);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);

    Bodies b;
    for(int i = 0; i < 100; ++i) b.add({5, -5, 5}, {}, 1.0 + i);
    Accel a = computeAccel(g, b);

    ASSERT_TRUE(allFinite(a));
    for(const auto& ai : a) EXPECT_VEC3_NEAR(ai, glm::dvec3(0.0), 1e-12);
}

TEST(BarnesHutGpu, CoincidentPairMatchesReference) {
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::BARNES_HUT_GPU);
    g.setTheta(0.0);

    Bodies b = test::randomCube(50, 31);
    b.add(b.pos[7], {}, 4.0);
    b.add(b.pos[7], {}, 0.5);

    Accel a = computeAccel(g, b);
    ASSERT_TRUE(allFinite(a));
    EXPECT_LT(maxOf(relativeErrors(a, test::referenceAccel(b, G, eps))), kFloatTol);
}

//------------------------------------------------------------------------------------------------
// Both GPU paths on one Gravity: they share the OpenCL device
//------------------------------------------------------------------------------------------------

TEST(GravityGpu, SwitchingAlgorithmsKeepsBothCorrect) {
    const double G = 1.0, eps = 0.01;
    Gravity g(G, eps);
    REQUIRE_GPU(g, Algo::DIRECT_GPU);
    g.setTheta(0.0);

    Bodies b = test::randomCube(300, 15);
    Accel ref = test::referenceAccel(b, G, eps);

    for(Algo algo : {Algo::DIRECT_GPU, Algo::BARNES_HUT_GPU, Algo::DIRECT_GPU, Algo::BARNES_HUT_GPU}) {
        g.setAlgorithm(algo);
        EXPECT_LT(maxOf(relativeErrors(computeAccel(g, b), ref)), kFloatTol);
    }
}
