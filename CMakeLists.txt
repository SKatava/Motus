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
)

add_subdirectory(external)
add_subdirectory(examples)

target_link_libraries(${NAME} PUBLIC glm)

find_package(OpenMP REQUIRED)
target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_CXX)

target_include_directories(${NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(${NAME} PRIVATE -Wall -Wextra -Wpedantic)



