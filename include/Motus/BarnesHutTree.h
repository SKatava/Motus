#pragma once

#include <glm/vec3.hpp>

#include <vector>
#include <cstdint>

#include "Motus/Bodies.h"

namespace motus {

class BarnesHutTree {
    public:
        //Construcotors
        BarnesHutTree(BodiesView b); 

        //Getters
        [[nodiscard]] const std::vector<int>& indices() const;

        //Core functions
        void build();
        glm::dvec3 traverseNode(int bodyIndex, int nodeIndex, double theta, double G, double epsilon);
    private:
        
        struct Node {
            double      mass {};            //8     bytes
            glm::dvec3  cm{};               //24    bytes
            glm::dvec3  center;             //24    bytes
            double      size;               //8     bytes
            int         first {-1};         //4     bytes
            int         firstChild {-1};    //4     bytes
            int         count {};           //4     bytes   =   74 bytes! => compiler padds to 80
 
        };

        static constexpr int s_bucketSize {1};

        std::vector<Node>           m_nodes;
        BodiesView                  m_bodies;
        std::vector<int>            m_indices;
        std::vector<std::uint64_t>  m_codes;

        void buildNode(int nodeIndex, int begin, int end, int depth);
        glm::dvec3 getCenter(int nodeIndex, int quad);
        void subdivide(int nodeIndex);
        void accumulateMass(int nodeIndex);
        void ZOrder();
        std::uint64_t code_zOrder(int index);
};

}
