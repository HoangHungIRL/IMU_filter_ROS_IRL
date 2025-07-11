cmake_minimum_required(VERSION 3.0.2)
project(imu_filter_ros_irl)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  eigen_conversions
  message_filters
)
find_package(Eigen3 REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS roscpp sensor_msgs eigen_conversions message_filters
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

# Build EKF node
add_executable(imu_ekf_node 
  src/imu_ekf_node.cpp
  src/butterworth.cpp
)
target_link_libraries(imu_ekf_node
  ${catkin_LIBRARIES}
  ${EIGEN3_LIBRARIES}
)
add_dependencies(imu_ekf_node
  ${catkin_EXPORTED_TARGETS}
)

# Build ESKF node
add_executable(imu_eskf_node
  src/imu_eskf_node.cpp
  src/butterworth.cpp
)
target_link_libraries(imu_eskf_node
  ${catkin_LIBRARIES}
  ${EIGEN3_LIBRARIES}
)
add_dependencies(imu_eskf_node
  ${catkin_EXPORTED_TARGETS}
)