cmake_minimum_required(VERSION 3.1)
project(sqp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
set(MODULE_DIR ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)

include_directories(
    ${MODULE_DIR}/core/include
    ${MODULE_DIR}/example/
    ${EIGEN3_INCLUDE_DIR}
)

SET(src_file
    ${MODULE_DIR}/core/src/sqp.cpp
    ${MODULE_DIR}/core/src/integratior.cpp
    ${MODULE_DIR}/main.cpp
)

ADD_EXECUTABLE(Result ${src_file})

target_link_libraries(Result
${catkin_LIBRARIES}  qpOASES
)