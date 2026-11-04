#pragma once

#include <CL/cl.h>

#include <memory>
#include <vector>

#include "Motus/ocl/Device.h"
#include "Motus/ocl/Handle.h"
#include "Motus/Bodies.h"

namespace motus::ocl {

class BarnesHutKernel {
    public:
        //Constructor
        explicit BarnesHutKernel(std::shared_ptr<Device> device);

        //Copy
        BarnesHutKernel(const BarnesHutKernel&) = delete;
        BarnesHutKernel& operator=(const BarnesHutKernel&) = delete;

        //Move
        BarnesHutKernel(BarnesHutKernel&&) = default;
        BarnesHutKernel& operator=(BarnesHutKernel&&) = default;

        //Destructor
        ~BarnesHutKernel() = default;

        //Core
        void compute(BodiesView b, AccelView out, double G, double eps, double theta);


    private:
        struct GpuNode {
            cl_float4 cmMass;
            cl_float4 centerSize;
            cl_int4 meta;
        };

        std::shared_ptr<Device> m_device;
        Program m_program;
        Kernel m_kernel;
        Buffer m_posMass;
        Buffer m_nodes;
        Buffer m_acc;
        size_t m_capacity {};
        size_t m_capacityNode {};
        std::vector<cl_float4> m_staging;
        std::vector<GpuNode> m_stagingNodes;

        //Private helpers
        void ensureCapacity(std::size_t n);
        void ensureCapacityNodes(std::size_t n);
};

}
