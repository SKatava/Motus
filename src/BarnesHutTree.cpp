#include "Motus/BarnesHutTree.h"

#include "Motus/ScopeTimer.h"

#include <glm/glm.hpp>
#include <cmath>


namespace motus {

BarnesHutTree::BarnesHutTree(BodiesView b) : m_bodies(b) {
    m_nodes.reserve(4 * m_bodies.size() + 8);
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
}

void BarnesHutTree::build() {
    {ScopeTimer t {"inserting nodes: "};
    for(size_t i = 0; i < m_bodies.size(); ++i) {
        insertToNode(i, 0);
    }
    }
    {ScopeTimer t {"Accumulate mass: "}; accumulateMass(0);}
}

void BarnesHutTree::insertToNode(int bodyIndex, int nodeIndex) {

    //Internal node
    if(m_nodes[nodeIndex].firstChild != -1) {
        int child = m_nodes[nodeIndex].firstChild + getQuadrant(bodyIndex, nodeIndex);
        insertToNode(bodyIndex, child);
        return;
    }

    //empty leaf
    if(m_nodes[nodeIndex].body == -1) {
        m_nodes[nodeIndex].body = bodyIndex;
        return;
    }

    //full leaf
    int existing = m_nodes[nodeIndex].body;
    m_nodes[nodeIndex].body = -1;
    subdivide(nodeIndex);
    insertToNode(existing, nodeIndex);
    insertToNode(bodyIndex, nodeIndex);
}

int BarnesHutTree::getQuadrant(int bodyIndex, int nodeIndex) {
    glm::dvec3& center = m_nodes[nodeIndex].center;
    glm::dvec3& pos = m_bodies.pos[bodyIndex];

    return  (pos.x > center.x) << 2 |
            (pos.y > center.y) << 1 |
            (pos.z > center.z);
}

glm::dvec3 BarnesHutTree::getCenter(int nodeIndex, int quad) {
    glm::dvec3 c = m_nodes[nodeIndex].center;
    const double off = m_nodes[nodeIndex].size * 0.25;
    return glm::dvec3(
        c.x + ((quad & 4) ? off : -off),
        c.y + ((quad & 2) ? off : -off),
        c.z + ((quad & 1) ? off : -off)
    );
}

void BarnesHutTree::subdivide(int nodeIndex) {
    double size = m_nodes[nodeIndex].size * 0.5;
    
    int first = static_cast<int>(m_nodes.size());

    for(size_t i = 0; i < 8; ++i) {
        Node child;
        child.size = size;
        child.center = getCenter(nodeIndex, i);
        m_nodes.emplace_back(child);
    }

    m_nodes[nodeIndex].firstChild = first;
}

void BarnesHutTree::accumulateMass(int nodeIndex) {
    if(m_nodes[nodeIndex].body != -1) {
        m_nodes[nodeIndex].mass = m_bodies.m[m_nodes[nodeIndex].body];
        m_nodes[nodeIndex].cm = m_bodies.pos[m_nodes[nodeIndex].body];
        return;
    } 

    if(m_nodes[nodeIndex].firstChild != -1) {
        int first = m_nodes[nodeIndex].firstChild;
        glm::dvec3 weighted {0.0};
        double M = 0;

        for(size_t i = 0; i < 8; ++i) {
            accumulateMass(first + i);
            M += m_nodes[first + i].mass;
            weighted += m_nodes[first + i].mass * m_nodes[first + i].cm;
        }

        m_nodes[nodeIndex].mass = M;
        m_nodes[nodeIndex].cm = (M > 0.0) ? (weighted / M) : m_nodes[nodeIndex].center;
    }
}

glm::dvec3 BarnesHutTree::traverseNode(int bodyIndex, int nodeIndex, double theta, double G, double epsilon) {
    const glm::dvec3 diff = m_nodes[nodeIndex].cm - m_bodies.pos[bodyIndex]; 
    const double d2 = glm::dot(diff, diff) + epsilon*epsilon;
    const double invD3 = 1.0 / (d2 * std::sqrt(d2));

    if(m_nodes[nodeIndex].body != -1) {
        if(m_nodes[nodeIndex].body == bodyIndex) return glm::dvec3 {0.0};

        return G * m_nodes[nodeIndex].mass * m_bodies.m[bodyIndex] * invD3 * diff;
    }
    else if(m_nodes[nodeIndex].firstChild != -1) {
        const double d = std::sqrt(d2);

        if(m_nodes[nodeIndex].size / d < theta) {
            return G * m_nodes[nodeIndex].mass * m_bodies.m[bodyIndex] * invD3 * diff;
        }

        glm::dvec3 force {0.0};
        for(int i = 0; i < 8; ++i) {
            force += traverseNode(bodyIndex, m_nodes[nodeIndex].firstChild + i, theta, G, epsilon);
        }
        return force;
    }

    return glm::dvec3 {0.0};
}

}

