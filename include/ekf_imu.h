#ifndef EKF_IMU_H
#define EKF_IMU_H

#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/magnetic_field.hpp>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <Eigen/Dense>
#include "butterworth.h"

class EKF_IMU : public rclcpp::Node {
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
    Eigen::Vector4d update(const sensor_msgs::msg::Imu& imu_msg, const sensor_msgs::msg::MagneticField* mag_msg, double dt);
    void syncedCallback(const sensor_msgs::msg::Imu::SharedPtr imu_msg, const sensor_msgs::msg::MagneticField::SharedPtr mag_msg);

private:
    rclcpp::Clock::SharedPtr clock_;
    std::shared_ptr<message_filters::Subscriber<sensor_msgs::msg::Imu>> imu_sub_;
    std::shared_ptr<message_filters::Subscriber<sensor_msgs::msg::MagneticField>> mag_sub_;
    std::shared_ptr<message_filters::TimeSynchronizer<sensor_msgs::msg::Imu, sensor_msgs::msg::MagneticField>> sync_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr imu_pub_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr accel_comp_pub_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr gravity_pub_;
    Eigen::Matrix4d P_;
    Eigen::MatrixXd R_;
    Eigen::Matrix3d g_noise_;
    Eigen::Vector3d a_ref_;
    Eigen::Vector3d m_ref_;
    double mag_norm_ref_;
    Eigen::Vector3d mag_bias_;
    Eigen::Vector4d q_;
    Eigen::VectorXd z_;
    rclcpp::Time last_time_;
    sensor_msgs::msg::Imu last_imu_msg_;
    sensor_msgs::msg::MagneticField last_mag_msg_;
    rclcpp::Time last_mag_time_;
    uint32_t update_count_;
    bool initialized_;
    bool use_magnetometer_;
    bool has_recent_mag_;
    bool noise_filter_;
    bool filter_gyro_;
    Butter2 butter_ax_;
    Butter2 butter_ay_;
    Butter2 butter_az_;
    Butter2 butter_wx_;
    Butter2 butter_wy_;
    Butter2 butter_wz_;
};

#endif // EKF_IMU_H