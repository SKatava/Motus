include(GoogleTest)

add_executable(${NAME}-tests
    GravityTest.cpp
    BarnesHutTest.cpp
    IntegratorTest.cpp
    SimulationTest.cpp
    ProfilerTest.cpp
)

target_link_libraries(${NAME}-tests
    PRIVATE
        motus
        GTest::gtest_main
)

#Degenerate tree inputs can recurse forever; fail fast instead of waiting for ctest's 1500 s default
gtest_discover_tests(${NAME}-tests PROPERTIES TIMEOUT 10)
