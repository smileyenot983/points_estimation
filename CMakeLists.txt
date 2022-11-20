cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)

# find_package(Eigen3 REQUIRED)

find_package(OpenMP)


add_executable(estimate_points src/estimate_points.cpp)
# target_link_libraries(estimate_points  Eigen3::Eigen OpenMP::OpenMP_CXX)
target_link_libraries(estimate_points OpenMP::OpenMP_CXX)