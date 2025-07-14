#ifndef ESKF_IMU_H
#define ESKF_IMU_H

#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/magnetic_field.hpp>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "butterworth.h"

#define DEG2PI (M_PI / 180.0)

struct Observation {
    Eigen::Quaterniond quat;
    Eigen::Vector3d euler;
};

class ImuEskf : public rclcpp::Node {
public:
    ImuEskf();
    Eigen::Matrix3d skew(const Eigen::Vector3d& x);
    Eigen::Matrix3d q2R(const Eigen::Vector4d& q);
    Eigen::Vector4d ecompass(const Eigen::Vector3d& acc, const Eigen::Vector3d& mag);
    Eigen::MatrixXd diffQstarvqQ(const Eigen::Quaterniond& q, const Eigen::Vector3d& v);
    void predict();
    void measurement();
    void estimate();
    void error2nominal();
    Observation getData();
    double getTime() const;
    void syncedCallback(const sensor_msgs::msg::Imu::SharedPtr imu_msg, const sensor_msgs::msg::MagneticField::SharedPtr mag_msg);

private:
    rclcpp::Clock::SharedPtr clock_;
    std::shared_ptr<message_filters::Subscriber<sensor_msgs::msg::Imu>> imu_sub_;
    std::shared_ptr<message_filters::Subscriber<sensor_msgs::msg::MagneticField>> mag_sub_;
    std::shared_ptr<message_filters::TimeSynchronizer<sensor_msgs::msg::Imu, sensor_msgs::msg::MagneticField>> sync_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr imu_pub_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr accel_comp_pub_;
    rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr gravity_pub_;
    Eigen::VectorXd x_error_;
    Eigen::VectorXd x_nominal_;
    Eigen::VectorXd z_;
    Eigen::MatrixXd P_;
    Eigen::MatrixXd Fx_;
    Eigen::MatrixXd H_;
    Eigen::MatrixXd Fi_;
    Eigen::MatrixXd Q_;
    Eigen::MatrixXd R_;
    Eigen::Vector3d a_ref_;
    Eigen::Vector3d m_ref_;
    double mag_norm_ref_;
    Eigen::Vector3d mag_bias_;
    Eigen::Vector3d raw_acc_;
    Eigen::Vector3d raw_gyro_;
    Eigen::Vector3d raw_mag_;
    rclcpp::Time last_time_;
    sensor_msgs::msg::Imu last_imu_msg_;
    sensor_msgs::msg::MagneticField last_mag_msg_;
    rclcpp::Time last_mag_time_;
    uint32_t update_count_;
    bool initialized_;
    bool euler_initialized_;
    bool use_magnetometer_;
    bool has_recent_mag_;
    bool noise_filter_;
    bool filter_gyro_;
    double current_t_;
    double dt_;
    double gyro_noise_;
    double gyro_bias_noise_;
    double acc_noise_;
    double mag_noise_;
    double theta_noise_;
    double wb_noise_;
    Butter2 butter_ax_;
    Butter2 butter_ay_;
    Butter2 butter_az_;
    Butter2 butter_wx_;
    Butter2 butter_wy_;
    Butter2 butter_wz_;
    const int n_state_ = 6;
};

#endif // ESKF_IMU_H