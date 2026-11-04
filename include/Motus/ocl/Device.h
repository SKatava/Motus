#pragma once

#include <CL/cl.h>

#include <string>

#include "Motus/ocl/Handle.h"

namespace motus::ocl {

class Device {
    public:
        //Constructor
        Device();

        //Copy
        Device(const Device&) = delete;
        Device& operator=(const Device&) = delete;

        //Move
        Device(Device&&) = default;
        Device& operator=(Device&&) = default;

        //Getters
        [[nodiscard]] cl_device_id device() const noexcept;
        [[nodiscard]] cl_context context() const noexcept;
        [[nodiscard]] cl_command_queue queue() const noexcept;
        [[nodiscard]] std::string name() const;

    private:
        cl_platform_id  m_platform  {nullptr};
        cl_device_id    m_device    {nullptr};
        Context         m_context;
        Queue           m_queue;
};

}
