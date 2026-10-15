function(motus_example) 
    cmake_parse_arguments(EX "" "NAME" "SOURCES" ${ARGN})
    add_executable(${EX_NAME} ${EX_SOURCES})
    target_link_libraries(${EX_NAME} PRIVATE ${NAME})
endfunction()

motus_example(NAME 01_BHT SOURCES 01_BHT/main.cpp)
