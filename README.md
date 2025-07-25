# IMU Attitude Estimation Using EKF and ESKF Filters with Butterworth Low-Pass Filtering



## Overview

This IMU filter ROS/ROS2 package implements both Extended Kalman Filter (EKF) and Error-State Kalman Filter (ESKF) for estimating the orientation of IMU using accelerometer, gyroscope, and optional magnetometer data. Second-order Butterworth low-pass filter is applied to reduce noise in the IMU's linear acceleration and angular velocity measurements.

The Butterworth filter can be enabled or disabled via a ROS parameter, offering flexibility for different use cases.

## Requirements

Ubuntu 20.04.

ROS Noetic.

## Install matplotlib (optional)
```
python -m pip install -U pip
python -m pip install -U matplotlib
```

## Install the package

Use the following commands:

```
mkdir -p ~/imu_filter_ws/src
cd ~/imu_filter_ws/src
git clone --branch ROS-Noetic https://github.com/HoangHungIRL/IMU_filter_IRL_ROS2.git

cd ~/imu_filter_ws/src/IMU_filter_IRL_ROS2/scripts
sudo chmod +x euler_plotter.py

cd ~/imu_filter_ws
catkin_make
```
## Launch

Use the following commands for EKF:

```
cd ~/imu_filter_ws
source devel/setup.bash
roslaunch imu_filter_ros_irl imu_ekf_node.launch
```
Use the following commands for ESKF:

```
cd ~/imu_filter_ws
source devel/setup.bash
roslaunch imu_filter_ros_irl imu_eskf_node.launch
```

## License

This project is licensed under the MIT License.

## Contact

For questions or support, contact Hoang Quoc Hung via email hoanghung21301580@gmail.com or GitHub or open an issue in the repository.

## Acknowledgments


Developed as part of the IRL Agricultural Autonomous Vehicle project.


Implementation inspired by:

https://github.com/ZacharyTaylor/butter.

https://github.com/alalagong/EKF

[AHRS ](https://ahrs.readthedocs.io/en/latest/)

