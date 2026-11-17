#pragma once

#include <CL/cl.h>

#include <utility>
#include <stdexcept>
#include <string>

namespace motus::ocl {

template <typename T, auto Release> 
class Handle {
    public:
        //Constructors
        Handle() = default;
        explicit Handle(T handle) noexcept : m_handle(handle) {}

        //Move and copy
        Handle(const Handle&) = delete;
        Handle& operator=(const Handle&) = delete;

        Handle(Handle&& other) noexcept : m_handle(std::exchange(other.m_handle, nullptr)) {}
        Handle& operator=(Handle&& other) noexcept {
            if(&other != this) {
                if(m_handle) Release(m_handle);
                m_handle = std::exchange(other.m_handle, nullptr);
            }
            return *this;
        }

        //Destructor
        ~Handle() { if(m_handle) Release(m_handle); }

        [[nodiscard]] T get() const noexcept { return m_handle; }
        explicit operator bool() const noexcept { return m_handle != nullptr; }

    private:
        T m_handle = nullptr;
};

class Error : public std::runtime_error {
    public:
        Error(cl_int code, const char* what) : 
            std::runtime_error(std::string(what) + " failed (" + std::to_string(code) + ")"), 
            m_code(code) {}

        [[nodiscard]] cl_int code() const noexcept { return m_code; }

    private:
        cl_int m_code;
}; 

inline void check(cl_int err, const char* what) {
    if(err != CL_SUCCESS) {
        throw Error(err, what);
    }
}

using Context   = Handle<cl_context, &clReleaseContext>;
using Queue     = Handle<cl_command_queue, &clReleaseCommandQueue>;
using Program   = Handle<cl_program, &clReleaseProgram>;
using Kernel    = Handle<cl_kernel, &clReleaseKernel>;
using Buffer    = Handle<cl_mem, &clReleaseMemObject>;

}
