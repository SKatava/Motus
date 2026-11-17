#include "Motus/BarnesHutTree.h"


#include "Motus/Profiler.h"

#include <glm/glm.hpp>
#include <cmath>
#include <parallel/algorithm>

namespace motus {

BarnesHutTree::BarnesHutTree(BodiesView b) : m_bodies(b) {
    if(m_bodies.size() == 0) return;

    m_nodes.reserve(4 * m_bodies.size() + 8);
    m_indices.reserve(m_bodies.size());
    m_codes.reserve(m_bodies.size());

    m_nodes.emplace_back();

    glm::dvec3 max (b.pos[0]);
    glm::dvec3 min (b.pos[0]);

    for(const auto& p : b.pos) {
        min = glm::min(min, p);
        max = glm::max(max, p);
    }

    glm::dvec3 diff = max - min;

    m_nodes[0].size   = glm::max(diff.x, glm::max(diff.y, diff.z)) * 1.01 +1e-9;
    m_nodes[0].center = (min + max) * 0.5;

    for(size_t i {}; i < m_bodies.size(); ++i) {
        m_indices.emplace_back(i);
    }

    {
        PROFILE_SCOPE("bh/zOrder");
        ZOrder();
    }
}

void BarnesHutTree::build() {
    if(m_bodies.size() == 0) return;

    buildNode(0, 0, m_bodies.size(), 0);
    accumulateMass(0);
}

glm::dvec3 BarnesHutTree::getCenter(int nodeIndex, int quad) {
    const Node& node = m_nodes[nodeIndex];
    glm::dvec3 c = node.center;
    const double off = node.size * 0.25;
    return glm::dvec3(
        c.x + ((quad & 4) ? off : -off),
        c.y + ((quad & 2) ? off : -off),
        c.z + ((quad & 1) ? off : -off)
    );
}

void BarnesHutTree::accumulateMass(int nodeIndex) {
    Node& node = m_nodes[nodeIndex];
    glm::dvec3 weighted {0.0};
    double M = 0;

    if(m_nodes[nodeIndex].firstChild == -1) {
        for(int i = node.first; i < node.first + node.count; ++i) {
            M += m_bodies.m[m_indices[i]];
            weighted += m_bodies.m[m_indices[i]] * m_bodies.pos[m_indices[i]];
        }
    }

    else {
        int first = node.firstChild;

        for(size_t i = 0; i < 8; ++i) {
            accumulateMass(first + i);
            M += m_nodes[first + i].mass;
            weighted += m_nodes[first + i].mass * m_nodes[first + i].cm;
        }
    }

    node.mass = M;
    node.cm = (M > 0.0) ? (weighted / M) : m_nodes[nodeIndex].center;

}

void BarnesHutTree::subdivide(int nodeIndex) {
    Node& node = m_nodes[nodeIndex];
    double size = node.size * 0.5;
    
    int first = static_cast<int>(m_nodes.size());

    for(size_t i = 0; i < 8; ++i) {
        Node child;
        child.size = size;
        child.center = getCenter(nodeIndex, i);
        m_nodes.emplace_back(child);
    }

    node.firstChild = first;
}

void BarnesHutTree::buildNode(int nodeIndex, int begin, int end, int depth) {
    Node& node = m_nodes[nodeIndex];
    node.first = begin;
    node.count = end - begin;

    if(end - begin <= s_bucketSize || depth == 21) {
        return;
    }

    subdivide(nodeIndex);
    int firstChild = node.firstChild;

    std::uint8_t shift = 3 * (20-depth);
    int i = begin;

    for(int j {}; j < 8; ++j) {
        int start = i;
        while(i < end && j == ((m_codes[i] >> shift) & 7)) {
            ++i;
        }
        buildNode(firstChild + j, start, i, depth + 1);
    }

    

}

glm::dvec3 BarnesHutTree::traverseNode(int bodyIndex, int nodeIndex, double theta, double G, double epsilon) {
    glm::dvec3 diff; 
    double d2;
    double invD3;
    glm::dvec3 force {0.0};
    Node& node = m_nodes[nodeIndex];

    if(m_nodes[nodeIndex].firstChild == -1) {

        for(int i = node.first; i < node.first + node.count; ++i) {
             
            if(m_indices[i] == bodyIndex) continue;

            diff = m_bodies.pos[m_indices[i]] - m_bodies.pos[bodyIndex]; 
            d2 = glm::dot(diff, diff) + epsilon*epsilon;
            invD3 = 1.0 / (d2 * std::sqrt(d2));

            force += G * m_bodies.m[m_indices[i]] * m_bodies.m[bodyIndex] * invD3 * diff; 
        }
    }
    else {
        diff = node.cm - m_bodies.pos[bodyIndex]; 
        d2 = glm::dot(diff, diff) + epsilon*epsilon;
        invD3 = 1.0 / (d2 * std::sqrt(d2));
        const double d = std::sqrt(d2);

        glm::dvec3 r = m_bodies.pos[bodyIndex] - node.center;
        if( (node.size / d < theta) && (glm::dot(r, r) > node.size * node.size)) {
            return G * node.mass * m_bodies.m[bodyIndex] * invD3 * diff;
        }

        for(int i = 0; i < 8; ++i) {
            force += traverseNode(bodyIndex, node.firstChild + i, theta, G, epsilon);
        }
    }
    
    return force;
}

void BarnesHutTree::ZOrder() {
    m_codes.resize(m_bodies.size());
    m_indices.resize(m_bodies.size());
    m_sortBuffer.resize(m_bodies.size());
  

    #pragma omp parallel for schedule(dynamic, 64)
    for(size_t i = 0; i < m_bodies.size(); ++i) {
        m_sortBuffer[i] = {code_zOrder(i), i};
    }

    //Needs to be Radix Sort
    __gnu_parallel::sort(m_sortBuffer.begin(), m_sortBuffer.end());

    #pragma omp parallel for schedule(dynamic, 64)
    for(size_t i = 0; i < m_bodies.size(); ++i) {
        m_codes[i] = m_sortBuffer[i].first; 
        m_indices[i] = m_sortBuffer[i].second; 
    }
}

std::uint64_t BarnesHutTree::code_zOrder(int index) {
    glm::dvec3 min = m_nodes[0].center - m_nodes[0].size / 2;
    double chunk = m_nodes[0].size / (1ULL << 21);
    std::uint64_t x = (m_bodies.pos[index].x - min.x) / chunk;
    std::uint64_t y = (m_bodies.pos[index].y - min.y) / chunk;
    std::uint64_t z = (m_bodies.pos[index].z - min.z) / chunk;
    std::uint64_t code = 0;
    for(int n {}; n < 21; ++n) {
        code |= ((x >> n) & 1ULL) << (3*n+2);
        code |= ((y >> n) & 1ULL) << (3*n+1);
        code |= ((z >> n) & 1ULL) << (3*n);
    }
    return code;
}

const std::vector<int>& BarnesHutTree::indices() const { return m_indices; }
const std::vector<BarnesHutTree::Node>& BarnesHutTree::nodes() const { return m_nodes; }
const glm::dvec3& BarnesHutTree::rootCenter() const { return m_nodes[0].center; }

}
