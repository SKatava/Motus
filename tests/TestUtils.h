#pragma once

#include "Motus/Bodies.h"

#include <glm/glm.hpp>

#include <cmath>
#include <random>
#include <vector>

namespace motus::test {

//Owning SoA storage for tests, the library itself only works on views
struct Bodies {
    std::vector<glm::dvec3> pos;
    std::vector<glm::dvec3> vel;
    std::vector<double>     m;

    void add(glm::dvec3 p, glm::dvec3 v, double mass) {
        pos.push_back(p);
        vel.push_back(v);
        m.push_back(mass);
    }

    [[nodiscard]] size_t size() const { return m.size(); }

    BodiesView view() { return {pos, vel, m}; }
};

inline Bodies randomCube(size_t n, unsigned seed, double lo = -1.0, double hi = 1.0) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> coord(lo, hi);
    std::uniform_real_distribution<double> mass(0.5, 2.0);

    Bodies b;
    for(size_t i = 0; i < n; ++i) {
        b.add({coord(gen), coord(gen), coord(gen)}, {0.0, 0.0, 0.0}, mass(gen));
    }
    return b;
}

//Two equal masses on a circular orbit around their common center of mass at the origin.
//Returns the orbital period for gravitational constant G.
inline double circularBinary(Bodies& b, double G, double mass, double separation) {
    const double M      = 2.0 * mass;
    const double vRel   = std::sqrt(G * M / separation);
    const double half   = 0.5 * separation;

    b.add({-half, 0.0, 0.0}, {0.0, -0.5 * vRel, 0.0}, mass);
    b.add({ half, 0.0, 0.0}, {0.0,  0.5 * vRel, 0.0}, mass);

    return 2.0 * M_PI * std::sqrt(separation * separation * separation / (G * M));
}

inline glm::dvec3 totalMomentum(const Bodies& b) {
    glm::dvec3 p {0.0};
    for(size_t i = 0; i < b.size(); ++i) p += b.m[i] * b.vel[i];
    return p;
}

inline glm::dvec3 totalAngularMomentum(const Bodies& b) {
    glm::dvec3 L {0.0};
    for(size_t i = 0; i < b.size(); ++i) L += b.m[i] * glm::cross(b.pos[i], b.vel[i]);
    return L;
}

inline glm::dvec3 centerOfMass(const Bodies& b) {
    glm::dvec3 c {0.0};
    double M = 0.0;
    for(size_t i = 0; i < b.size(); ++i) {
        c += b.m[i] * b.pos[i];
        M += b.m[i];
    }
    return c / M;
}

//Reference O(N^2) softened gravity, written independently from the library on purpose
inline Accel referenceAccel(const Bodies& b, double G, double eps) {
    Accel a(b.size(), glm::dvec3{0.0});
    for(size_t i = 0; i < b.size(); ++i) {
        for(size_t j = 0; j < b.size(); ++j) {
            if(i == j) continue;
            const glm::dvec3 r = b.pos[j] - b.pos[i];
            const double d2 = glm::dot(r, r) + eps * eps;
            a[i] += G * b.m[j] * r / (d2 * std::sqrt(d2));
        }
    }
    return a;
}

}

#define EXPECT_VEC3_NEAR(a, b, tol)             \
    do {                                        \
        const glm::dvec3 _va = (a);             \
        const glm::dvec3 _vb = (b);             \
        EXPECT_NEAR(_va.x, _vb.x, (tol));       \
        EXPECT_NEAR(_va.y, _vb.y, (tol));       \
        EXPECT_NEAR(_va.z, _vb.z, (tol));       \
    } while(0)

#define EXPECT_VEC3_DOUBLE_EQ(a, b)             \
    do {                                        \
        const glm::dvec3 _va = (a);             \
        const glm::dvec3 _vb = (b);             \
        EXPECT_DOUBLE_EQ(_va.x, _vb.x);         \
        EXPECT_DOUBLE_EQ(_va.y, _vb.y);         \
        EXPECT_DOUBLE_EQ(_va.z, _vb.z);         \
    } while(0)
