cmake_minimum_required(VERSION 3.3)
project(lab_3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.c)
add_executable(lab_3 ${SOURCE_FILES})