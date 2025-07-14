#ifndef ESKF_IMU_H
#define ESKF_IMU_H

#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/MagneticField.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>
#include <Eigen/Dense>
#include "butterworth.h"

#define DEG2PI (M_PI / 180.0)

struct Observation {
  Eigen::Quaterniond quat;
  Eigen::Vector3d euler;
};

class ImuEskf {
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
    void syncedCallback(const sensor_msgs::Imu::ConstPtr& imu_msg, const sensor_msgs::MagneticField::ConstPtr& mag_msg);

private:
    ros::NodeHandle nh_;
    message_filters::Subscriber<sensor_msgs::Imu> imu_sub_;
    message_filters::Subscriber<sensor_msgs::MagneticField> mag_sub_;
    message_filters::TimeSynchronizer<sensor_msgs::Imu, sensor_msgs::MagneticField> sync_;
    ros::Publisher imu_pub_;
    ros::Publisher accel_comp_pub_;
    ros::Publisher gravity_pub_;
    Eigen::VectorXd x_error_;     // Error state: [theta_x, theta_y, theta_z, wb_x, wb_y, wb_z]
    Eigen::VectorXd x_nominal_;   // Nominal state: [qw, qx, qy, qz, wb_x, wb_y, wb_z]
    Eigen::VectorXd z_;           // Measurement vector: [ax, ay, az, magx, magy, magz]
    Eigen::MatrixXd P_;           // State covariance
    Eigen::MatrixXd Fx_;          // State transition matrix
    Eigen::MatrixXd H_;           // Measurement Jacobian
    Eigen::MatrixXd Fi_;          // Noise input matrix
    Eigen::MatrixXd Q_;           // Process noise covariance
    Eigen::MatrixXd R_;           // Measurement noise covariance
    Eigen::Vector3d a_ref_;       // Reference gravity (m/s²)
    Eigen::Vector3d m_ref_;       // Reference magnetic field (unit vector)
    double mag_norm_ref_;         // Expected magnetic field norm (uT)
    Eigen::Vector3d mag_bias_;    // Magnetometer hard iron bias (uT)
    Eigen::Vector3d raw_acc_;     // Accelerometer (unit vector)
    Eigen::Vector3d raw_gyro_;    // Gyroscope (rad/s)
    Eigen::Vector3d raw_mag_;     // Magnetometer (unit vector)
    ros::Time last_time_;
    sensor_msgs::Imu last_imu_msg_;
    sensor_msgs::MagneticField last_mag_msg_;
    ros::Time last_mag_time_;
    uint32_t update_count_;
    bool initialized_;
    bool euler_initialized_;
    bool use_magnetometer_;
    bool has_recent_mag_;
    bool noise_filter_;
    bool filter_gyro_;           // Whether to filter gyroscope data when noise_filter is true
    double current_t_;
    double dt_;
    double gyro_noise_;           // Gyroscope noise variance
    double gyro_bias_noise_;      // Gyroscope bias noise variance
    double acc_noise_;            // Accelerometer noise variance
    double mag_noise_;            // Magnetometer noise variance
    double theta_noise_;          // Orientation error covariance
    double wb_noise_;             // Gyroscope bias error covariance
    Butter2 butter_ax_;           // Butterworth filter for accelerometer x
    Butter2 butter_ay_;           // Butterworth filter for accelerometer y
    Butter2 butter_az_;           // Butterworth filter for accelerometer z
    Butter2 butter_wx_;           // Butterworth filter for gyroscope x
    Butter2 butter_wy_;           // Butterworth filter for gyroscope y
    Butter2 butter_wz_;           // Butterworth filter for gyroscope z
    const int n_state_ = 6;       // Number of error state variables
};

#endif // ESKF_IMU_H