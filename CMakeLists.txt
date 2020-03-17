cmake_minimum_required(VERSION 2.8.3)
project(bnuk_opt)


add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/sdf_tools)

find_package(catkin REQUIRED COMPONENTS
        roscpp
        std_msgs
        nav_msgs
        visualization_msgs
        quadrotor_msgs
        tf
        image_transport
        cv_bridge
        )

find_package(Eigen3 REQUIRED)
message(Eigen3_dir: ${EIGEN3_INCLUDE_DIR} STATUS)
find_package(PCL REQUIRED)

set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

catkin_package(
        INCLUDE_DIRS include
)

link_directories(${PROJECT_SOURCE_DIR}/third_party/mosek/lib/mosek8_1)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/third_party/mosek/include)
include_directories(${PROJECT_SOURCE_DIR}/third_party/sdf_tools/include)

include_directories(
        include
        SYSTEM
        third_party
        ../../devel/include
        ${catkin_INCLUDE_DIRS}
        ${Eigen3_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
)

#set(CMAKE_CXX_FLAGS_DEBUG "-std=c++0x ${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS} -O3 -Wall") # -Wextra -Werror
add_library(bnuk_search
        src/bnuk_search.cpp
        )

add_executable(
        map_recorder
        src/map_recorder.cpp
)
target_link_libraries(map_recorder ${catkin_LIBRARIES} sdf_tools ${PCL_LIBRARIES})

add_executable(bnuk_traj_server src/bnuk_traj_server.cpp)
target_link_libraries(bnuk_traj_server
        ${catkin_LIBRARIES}
        )

add_executable(random_forest_sensing src/random_forest_sensing.cpp)
target_link_libraries(random_forest_sensing
        ${catkin_LIBRARIES}
        ${PCL_LIBRARIES}
        )

add_executable(odom_generator src/odom_generator.cpp)
target_link_libraries(odom_generator
        ${catkin_LIBRARIES}
        )

add_library(b_spline_optimizer
        src/b_spline_optimizer.cpp
        )

add_executable(
        bnuk_replan
        src/bnuk_replan.cpp

)
target_link_libraries(bnuk_replan ${catkin_LIBRARIES} sdf_tools ${PCL_LIBRARIES} bnuk_search b_spline_optimizer mosek64)


