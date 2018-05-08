cmake_minimum_required(VERSION 3.10)
project(Assignment)

set(CMAKE_CXX_STANDARD 11)


add_executable(Assignment main.cpp test.h weighted_graph.h)

target_link_libraries(Assignment cxxtest)