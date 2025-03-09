cmake_minimum_required(VERSION 3.10)

project(test_jacobian)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED)

file(GLOB SOURCE "main.cpp")

add_executable(test_jacobian ${SOURCE})

target_link_libraries(test_jacobian benchmark)