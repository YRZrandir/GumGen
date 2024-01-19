cmake_minimum_required(VERSION 3.0)
project(GumGen CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(Ceres REQUIRED)

aux_source_directory(src/ SRCS)
aux_source_directory(src/util UTIL_SRCS)
aux_source_directory(src/gumgen GUM_SRCS)
aux_source_directory(src/Bezier BEZIER_SRCS)

add_executable(GumGen
        ${SRCS} ${UTIL_SRCS} ${GUM_SRCS} ${BEZIER_SRCS})

target_compile_definitions(GumGen PRIVATE CGAL_PMP_USE_CERES_SOLVER)
target_compile_definitions(GumGen PRIVATE GLOG_NO_ABBREVIATED_SEVERITIES)
target_include_directories(GumGen PRIVATE dep/json/include)
target_include_directories(GumGen PRIVATE dep/libigl/include)

target_link_libraries(GumGen CGAL::CGAL Eigen3::Eigen Ceres::ceres)
