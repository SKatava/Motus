cmake_minimum_required(VERSION 3.22...3.31)

set(NAME motus)

project(${NAME})

set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(${NAME} STATIC)

target_sources(${NAME} PRIVATE
    src/Gravity.cpp
    src/Simulation.cpp
    src/Verlet.cpp
    src/Euler.cpp
    src/BarnesHutTree.cpp

    src/ocl/Device.cpp
    src/ocl/GravityKernel.cpp
    src/ocl/BarnesHutKernel.cpp
)

add_subdirectory(external)
add_subdirectory(examples)
enable_testing()
add_subdirectory(tests)

target_link_libraries(${NAME} PUBLIC glm)

find_package(OpenMP REQUIRED)
target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_CXX)

find_package(OpenCL REQUIRED)
target_link_libraries(${NAME} PUBLIC OpenCL::OpenCL)

target_include_directories(${NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(${NAME} PRIVATE -Wall -Wextra -Wpedantic)

target_compile_definitions(${NAME} PUBLIC MOTUS_PROFILE)
target_compile_definitions(${NAME} PUBLIC CL_TARGET_OPENCL_VERSION=300)

file(READ ${CMAKE_CURRENT_SOURCE_DIR}/kernels/direct.cl DIRECT_SOURCE)
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/kernels/barnesHut.cl BARNESHUT_SOURCE)
configure_file(kernels/DirectSource.h.in ${CMAKE_CURRENT_BINARY_DIR}/generated/Motus/ocl/kernels/DirectSource.h)
configure_file(kernels/BarnesHutSource.h.in ${CMAKE_CURRENT_BINARY_DIR}/generated/Motus/ocl/kernels/BarnesHutSource.h)
target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS kernels/direct.cl)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS kernels/barnesHut.cl)



