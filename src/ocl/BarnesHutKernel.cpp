#include "Motus/ocl/BarnesHutKernel.h"

#include <glm/glm.hpp>

#include <string>
#include <limits>

#include "Motus/ocl/kernels/BarnesHutSource.h"
#include "Motus/BarnesHutTree.h"
#include "Motus/Profiler.h"

namespace motus::ocl {

BarnesHutKernel::BarnesHutKernel(std::shared_ptr<Device> device) : m_device(std::move(device)) {
    const char* src = kernels::barnesHut;
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

    cl_kernel rawKernel = clCreateKernel(m_program.get(), "barnesHutGravity", &err);
    check(err, "clCreateKernel");
    m_kernel = Kernel{rawKernel};
}

void BarnesHutKernel::ensureCapacity(std::size_t n) {
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

void BarnesHutKernel::ensureCapacityNodes(std::size_t n) {
    if(n <= m_capacityNode) return;

    cl_int err;
    cl_mem rawBuffer = clCreateBuffer(m_device->context(), CL_MEM_READ_ONLY, n * sizeof(GpuNode), nullptr, &err);
    check(err, "clCreateBuffer");
    m_nodes = Buffer{rawBuffer};

    m_capacityNode = n;
    m_stagingNodes.resize(n);
}

void BarnesHutKernel::compute(BodiesView b, AccelView out, double G, double eps, double theta) {
    if(b.size() == 0) return;
    BarnesHutTree tree(b);
    {
        PROFILE_SCOPE("BHT-GPU/build");
        tree.build();
    }
    
    size_t n = b.size();
    size_t nNode = tree.nodes().size();
    if( n == 0 ) return;
    if( nNode == 0 ) return;
    ensureCapacity(n);
    ensureCapacityNodes(nNode);

    glm::dvec3 center = tree.rootCenter();
    const std::vector<int>& indices = tree.indices();


    for(size_t i {}; i < n; ++i) {
        glm::dvec3 rel = b.pos[indices[i]] - center;
        m_staging[i].s[0] = static_cast<float>(rel.x);
        m_staging[i].s[1] = static_cast<float>(rel.y);
        m_staging[i].s[2] = static_cast<float>(rel.z);
        m_staging[i].s[3] = static_cast<float>(b.m[indices[i]]);
    }


    const std::vector<BarnesHutTree::Node>& nodes = tree.nodes();
    for(size_t i {}; i < nNode; ++i) {
        glm::dvec3 relCenter =  nodes[i].center - center;
        glm::dvec3 relCm =  nodes[i].cm - center;
        m_stagingNodes[i].cmMass.s[0] = static_cast<float>(relCm.x);
        m_stagingNodes[i].cmMass.s[1] = static_cast<float>(relCm.y);
        m_stagingNodes[i].cmMass.s[2] = static_cast<float>(relCm.z);
        m_stagingNodes[i].cmMass.s[3] = static_cast<float>(nodes[i].mass);
        m_stagingNodes[i].centerSize.s[0] = static_cast<float>(relCenter.x);
        m_stagingNodes[i].centerSize.s[1] = static_cast<float>(relCenter.y);
        m_stagingNodes[i].centerSize.s[2] = static_cast<float>(relCenter.z);
        m_stagingNodes[i].centerSize.s[3] = static_cast<float>(nodes[i].size);
        m_stagingNodes[i].meta.s[0] = static_cast<cl_int>(nodes[i].first);
        m_stagingNodes[i].meta.s[1] = static_cast<cl_int>(nodes[i].firstChild);
        m_stagingNodes[i].meta.s[2] = static_cast<cl_int>(nodes[i].count);
        m_stagingNodes[i].meta.s[3] = 0;
    }

    check(clEnqueueWriteBuffer(m_device->queue(), m_posMass.get(), CL_FALSE, 0, n * sizeof(cl_float4), m_staging.data(), 0, nullptr, nullptr), "clEnqueueWriteBuffer");
    check(clEnqueueWriteBuffer(m_device->queue(), m_nodes.get(), CL_FALSE, 0, nNode * sizeof(GpuNode), m_stagingNodes.data(), 0, nullptr, nullptr), "clEnqueueWriteBuffer");

    cl_mem rawPosMass = m_posMass.get();
    cl_mem rawNodes = m_nodes.get();
    cl_mem rawAcc = m_acc.get();
    cl_int rawN = static_cast<cl_int>(n);
    cl_float rawG = static_cast<cl_float>(G);
    cl_float rawEps2 = static_cast<cl_float>(eps*eps);
    cl_float rawTheta2 = static_cast<cl_float>(theta*theta);

    check(clSetKernelArg(m_kernel.get(), 0, sizeof(cl_mem), &rawPosMass), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 1, sizeof(cl_mem), &rawAcc), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 2, sizeof(cl_mem), &rawNodes), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 3, sizeof(cl_int), &rawN), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 4, sizeof(cl_float), &rawG), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 5, sizeof(cl_float), &rawEps2), "clSetKernelArg");
    check(clSetKernelArg(m_kernel.get(), 6, sizeof(cl_float), &rawTheta2), "clSetKernelArg");

    size_t global = (n + s_groupSize - 1)/s_groupSize * s_groupSize;
    size_t local = s_groupSize;
    check(clEnqueueNDRangeKernel(m_device->queue(), m_kernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr), "clEnqueueNDRangeKernel");

    check(clEnqueueReadBuffer(m_device->queue(), m_acc.get(), CL_TRUE, 0, n * sizeof(cl_float4), m_staging.data(), 0, nullptr, nullptr), "clEnqueueReadBuffer");

    for(size_t i {}; i < n; ++i) {
        out[indices[i]] += glm::dvec3(m_staging[i].s[0], m_staging[i].s[1], m_staging[i].s[2]);
    }
}

}
