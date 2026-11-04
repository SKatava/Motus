#pragma once

#include <CL/cl.h>

#include <memory>
#include <vector>

#include "Motus/ocl/Device.h"
#include "Motus/ocl/Handle.h"
#include "Motus/Bodies.h"

namespace motus::ocl {

class GravityKernel {
    public:
        //Constructor
        explicit GravityKernel(std::shared_ptr<Device> device);

        //Copy
        GravityKernel(const GravityKernel&) = delete;
        GravityKernel& operator=(const GravityKernel&) = delete;

        //Move
        GravityKernel(GravityKernel&&) = default;
        GravityKernel& operator=(GravityKernel&&) = default;

        //Destructor
        ~GravityKernel() = default;

        //Core
        void compute(BodiesView b, AccelView out, double G, double eps);


    private:
        std::shared_ptr<Device> m_device;
        Program m_program;
        Kernel m_kernel;
        Buffer m_posMass;
        Buffer m_acc;
        size_t m_capacity {};
        std::vector<cl_float4> m_staging;

        //Private helpers
        void ensureCapacity(std::size_t n);
};

}
