# IMU EKF Filter with Butterworth Noise Filtering



## Overview

This ROS package implements an Extended Kalman Filter (EKF) for estimating the orientation of an Inertial Measurement Unit (IMU) using accelerometer, gyroscope, and optional magnetometer data. A second-order Butterworth low-pass filter can be applied to reduce noise in the IMU's linear acceleration and angular velocity data.

The Butterworth filter can be enabled or disabled via a ROS parameter, offering flexibility for different use cases.

## Requirements

Ubuntu 20.04.

ROS Noetic.

## Install

Use the following commands:

```
mkdir ~/imu_filter_ws/src
cd ~/imu_filter_ws/src
git clone https://github.com/HoangHungIRL/IMU_filter_ROS_IRL.git

cd ~/imu_filter_ws/src/IMU_filter_ROS_IRL/scripts
sudo chmod +x euler_plotter.py

cd ~/imu_filter_ws
catkin_make
```
## Launch

Use the following commands:

```
cd ~/imu_filter_ws
source devel/setup.bash
roslaunch imu_filter_ros_irl imu_ekf_node.launch
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

For questions or support, contact Hoang Quoc Hung via email hoanghung21301580@gmail.com or GitHub or open an issue in the repository.

## Acknowledgments


Developed as part of the IRL Agricultural Autonomous Vehicle project.


Implementation inspired by:

https://github.com/ZacharyTaylor/butter.

[AHRS ](https://ahrs.readthedocs.io/en/latest/)

