#pragma once

#include <glm/vec3.hpp>

#include <vector>

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
            double mass {0};        //8     bytes
            glm::dvec3 cm{0.0};     //24    bytes
            glm::dvec3 center;      //24    bytes
            double size;            //8     bytes
            int body {-1};          //4     bytes
            int firstChild {-1};    //4     bytes   =   62 bytes!
        };

        std::vector<Node>   m_nodes;
        BodiesView          m_bodies;
        std::vector<int>    m_indices;

        
        //Helper functions
        glm::dvec3 getCenter(int nodeIndex, int quad);
        void insertToNode(int bodyIndex, int nodeIndex);
        int getQuadrant(int bodyIndex, int nodeIndex);
        void computeMassDistribution();
        void accumulateMass(int nodeIndex);
        void subdivide(int nodeIndex);
        void ZOrder();
        std::uint64_t code_zOrder(int index);

};

}
