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

gtest_discover_tests(${NAME}-tests)
