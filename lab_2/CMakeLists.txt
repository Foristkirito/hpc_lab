cmake_minimum_required(VERSION 3.3)
project(proj1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    cpu_test/submit.sh
    cpu_test/test.c
    gpu_test/a.out
    gpu_test/submit.sh
    gpu_test/test.cu
    measure/measure_data.md
    p1/main.c
    p1/main.o
    p1/Makefile
    p1/out
    p1/stencil
    p1/submit.sh
    p1/wxl.a p1/proj_engine.h)

add_executable(proj1 ${SOURCE_FILES})