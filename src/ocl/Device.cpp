#include "Motus/ocl/Device.h"

#include <vector>
namespace motus::ocl {

Device::Device() {
    cl_uint count = 0;
    check(clGetPlatformIDs(0, nullptr, &count), "clGetPlatformIDs");
    std::vector<cl_platform_id> platforms(count);
    check(clGetPlatformIDs(count, platforms.data(), nullptr), "clGetPlatformIDs");
    for(cl_platform_id p : platforms) {
        cl_int err = clGetDeviceIDs(p, CL_DEVICE_TYPE_GPU, 1, &m_device, nullptr);
        if (err == CL_DEVICE_NOT_FOUND) continue;   
        check(err, "clGetDeviceIDs");
        m_platform = p;
        break;
    }

    if(m_device == nullptr) {
        throw Error(CL_DEVICE_NOT_FOUND, "Device: no GPU");
    }

    cl_int err;
    cl_context rawCtx = clCreateContext(nullptr, 1, &m_device, nullptr, nullptr, &err);
    check(err, "clCreateContext");
    m_context = Context{rawCtx};

    cl_command_queue rawQueue = clCreateCommandQueueWithProperties(m_context.get(), m_device, nullptr, &err);
    check(err, "clCreateCommandQueueWithProperties");
    m_queue = Queue{rawQueue};
}

cl_device_id Device::device() const noexcept { return m_device; }
cl_context Device::context() const noexcept { return m_context.get(); }
cl_command_queue Device::queue() const noexcept { return m_queue.get(); }    
std::string Device::name() const {
    size_t size = 0;
    check(clGetDeviceInfo(m_device, CL_DEVICE_NAME, 0, nullptr, &size), "clGetDeviceInfo");
    std::string s(size, '\0');
    check(clGetDeviceInfo(m_device, CL_DEVICE_NAME, size, s.data(), nullptr), "clGetDeviceInfo");
    s.pop_back();
    return s;
}

}
