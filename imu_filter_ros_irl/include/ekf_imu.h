#ifndef EKF_IMU_H
#define EKF_IMU_H

#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/MagneticField.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <Eigen/Dense>
#include "butterworth.h"

class EKF_IMU {
public:
    EKF_IMU();
    Eigen::Matrix3d skew(const Eigen::Vector3d& x);
    Eigen::Matrix3d q2R(const Eigen::Vector4d& q);
    Eigen::Vector4d ecompass(const Eigen::Vector3d& acc, const Eigen::Vector3d& mag);
    Eigen::MatrixXd Omega(const Eigen::Vector3d& x);
    Eigen::Vector4d f(const Eigen::Vector4d& q, const Eigen::Vector3d& omega, double dt);
    Eigen::Matrix4d dfdx(const Eigen::Vector3d& omega, double dt);
    Eigen::VectorXd h(const Eigen::Vector4d& q);
    Eigen::MatrixXd dhdx(const Eigen::Vector4d& q);
    Eigen::Vector4d update(const sensor_msgs::Imu& imu_msg, const sensor_msgs::MagneticField* mag_msg, double dt);
    void syncedCallback(const sensor_msgs::Imu::ConstPtr& imu_msg, const sensor_msgs::MagneticField::ConstPtr& mag_msg);

private:
    ros::NodeHandle nh_;
    message_filters::Subscriber<sensor_msgs::Imu> imu_sub_;
    message_filters::Subscriber<sensor_msgs::MagneticField> mag_sub_;
    message_filters::TimeSynchronizer<sensor_msgs::Imu, sensor_msgs::MagneticField> sync_;
    ros::Publisher imu_pub_;
    ros::Publisher accel_comp_pub_;
    ros::Publisher gravity_pub_;
    Eigen::Matrix4d P_;
    Eigen::MatrixXd R_;
    Eigen::Matrix3d g_noise_;
    Eigen::Vector3d a_ref_;
    Eigen::Vector3d m_ref_;
    double mag_norm_ref_;
    Eigen::Vector3d mag_bias_;
    Eigen::Vector4d q_;
    Eigen::VectorXd z_;
    ros::Time last_time_;
    sensor_msgs::Imu last_imu_msg_;
    sensor_msgs::MagneticField last_mag_msg_;
    ros::Time last_mag_time_;
    uint32_t update_count_;
    bool initialized_;
    bool use_magnetometer_;
    bool has_recent_mag_;
    bool noise_filter_;
    bool filter_gyro_;  // Whether to filter gyroscope data when noise_filter is true
    Butter2 butter_ax_;
    Butter2 butter_ay_;
    Butter2 butter_az_;
    Butter2 butter_wx_;
    Butter2 butter_wy_;
    Butter2 butter_wz_;
};

#endif // EKF_IMU_H