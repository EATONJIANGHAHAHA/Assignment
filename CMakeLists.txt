cmake_minimum_required(VERSION 3.10)
project(Assignment)

set(CMAKE_CXX_STANDARD 11)

set(INC_DIR /cxxtest/cxxtest-4.4/cxxtest)
set(LINK_DIR /cxxtest/cxxtest-4.4/cxxtest)
include_directories(${INC_DIR})
link_directories(${LINK_DIR})
LINK_LIBRARIES(cxxtest)

add_executable(Assignment main.cpp test.h weighted_graph.h)

target_link_libraries(Assignment cxxtest)