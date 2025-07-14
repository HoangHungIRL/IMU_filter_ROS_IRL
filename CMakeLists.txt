cmake_minimum_required(VERSION 3.8)
project(imu_filter_ros_irl)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -pedantic)

# Find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(rosgraph_msgs REQUIRED)
find_package(message_filters REQUIRED)
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
)

# Add executable for EKF node
add_executable(imu_ekf_node
  src/imu_ekf_node.cpp
  src/butterworth.cpp
)

# Link libraries for EKF node
ament_target_dependencies(imu_ekf_node
  rclcpp
  sensor_msgs
  message_filters
)
target_link_libraries(imu_ekf_node
  Eigen3::Eigen
)

# Add executable for ESKF node
add_executable(imu_eskf_node
  src/imu_eskf_node.cpp
  src/butterworth.cpp
)

# Link libraries for ESKF node
ament_target_dependencies(imu_eskf_node
  rclcpp
  sensor_msgs
  message_filters
)
target_link_libraries(imu_eskf_node
  Eigen3::Eigen
)

# Install executables
install(TARGETS
  imu_ekf_node
  imu_eskf_node
  DESTINATION lib/${PROJECT_NAME}
)

# Install header files
install(DIRECTORY include/
  DESTINATION include/${PROJECT_NAME}
)

# Install launch files
install(FILES
  launch/imu_eskf_node.launch.xml
  launch/imu_ekf_node.launch.xml
  DESTINATION share/${PROJECT_NAME}/launch
)

# Install parameter files
install(DIRECTORY config
  DESTINATION share/${PROJECT_NAME}
)

# Install Python scripts
install(PROGRAMS
  scripts/euler_plotter.py
  DESTINATION lib/${PROJECT_NAME}
)

# Export package dependencies
ament_export_dependencies(
  rclcpp
  sensor_msgs
  message_filters
)
ament_export_include_directories(include)

# Package setup
ament_package()