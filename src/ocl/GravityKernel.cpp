#include "Motus/ocl/GravityKernel.h"

#include <glm/glm.hpp>

#include <string>

#include "Motus/ocl/kernels/DirectSource.h"

namespace motus::ocl {

GravityKernel::GravityKernel(std::shared_ptr<Device> device) : m_device(std::move(device)) {
    const char* src = kernels::direct;
    cl_int err;
    cl_program rawProgram = clCreateProgramWithSource(m_device->context(), 1, &src, nullptr, &err);
    check(err, "clCreateProgramWithSource");
    m_program = Program{rawProgram};

    cl_device_id dev = m_device->device();
    cl_int buildErr = clBuildProgram(m_program.get(), 1, &dev, "", nullptr, nullptr);
    if(buildErr != CL_SUCCESS) {
        size_t size = 0;
        clGetProgramBuildInfo(m_program.get(), dev, CL_PROGRAM_BUILD_LOG, 0, nullptr, &size);
        std::string log(size, '\0');
        clGetProgramBuildInfo(m_program.get(), dev, CL_PROGRAM_BUILD_LOG, size, log.data(), nullptr);
        throw Error(buildErr, log.c_str());
    }

    cl_kernel rawKernel = clCreateKernel(m_program.get(), "directGravity", &err);
    check(err, "clCreateKernel");
    m_kernel = Kernel{rawKernel};
}

void GravityKernel::ensureCapacity(std::size_t n) {
   if(n <= m_capacity) return;

   cl_int err;
   cl_mem rawBuffer = clCreateBuffer(m_device->context(), CL_MEM_READ_ONLY, n * sizeof(cl_float4), nullptr, &err);
   check(err, "clCreateBuffer");
   m_posMass = Buffer{rawBuffer};

   rawBuffer = clCreateBuffer(m_device->context(), CL_MEM_WRITE_ONLY, n * sizeof(cl_float4), nullptr, &err);
   check(err, "clCreateBuffer");
   m_acc = Buffer{rawBuffer};

   m_capacity = n;
   m_staging.resize(n);
}

void GravityKernel::compute(BodiesView b, AccelView out, double G, double eps) {
    size_t n= b.size();
    if( n == 0) return;
    ensureCapacity(n);

    glm::dvec3 max = (b.pos[0]);
    glm::dvec3 min = (b.pos[0]);

    for(const auto& p : b.pos) {
        min = glm::min(min, p);
        max = glm::max(max, p);
    }

    glm::dvec3 center = (min + max) * 0.5;

    for(size_t i {}; i < n; ++i) {
        glm::dvec3 rel = b.pos[i] - center;
        m_staging[i].s[0] = static_cast<float>(rel.x);
        m_staging[i].s[1] = static_cast<float>(rel.y);
        m_staging[i].s[2] = static_cast<float>(rel.z);
        m_staging[i].s[3] = static_cast<float>(b.m[i]);
    }

    check(clEnqueueWriteBuffer(m_device->queue(), m_posMass.get(), CL_FALSE, 0, n * sizeof(cl_float4), m_staging.data(), 0, nullptr, nullptr), "clEnqueueWriteBuffer");

    cl_mem rawPosMass = m_posMass.get();
    cl_mem rawAcc = m_acc.get();
    cl_int rawN = static_cast<cl_int>(n);
    cl_float rawG = static_cast<cl_float>(G);
    cl_float rawEps2 = static_cast<cl_float>(eps*eps);

    check(clSetKernelArg(m_kernel.get(), 0, sizeof(cl_mem), &rawPosMass), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 1, sizeof(cl_mem), &rawAcc), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 2, sizeof(cl_int), &rawN), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 3, sizeof(cl_float), &rawG), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 4, sizeof(cl_float), &rawEps2), "clSetKernelArg");

    size_t global = n;
    check(clEnqueueNDRangeKernel(m_device->queue(), m_kernel.get(), 1, nullptr, &global, nullptr, 0, nullptr, nullptr), "clEnqueueNDRangeKernel");

    check(clEnqueueReadBuffer(m_device->queue(), m_acc.get(), CL_TRUE, 0, n * sizeof(cl_float4), m_staging.data(), 0, nullptr, nullptr), "clEnqueueReadBuffer");

    for(size_t i {}; i < n; ++i) {
        out[i] += glm::dvec3(m_staging[i].s[0], m_staging[i].s[1], m_staging[i].s[2]);
    }
}

}
