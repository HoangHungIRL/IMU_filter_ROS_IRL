#include "ekf_imu.h"
#include <cmath>
#include <ros/ros.h>

EKF_IMU::EKF_IMU() : nh_("/imu_ekf_node"), sync_(imu_sub_, mag_sub_, 10), update_count_(0), initialized_(false), use_magnetometer_(false), has_recent_mag_(false),
                     noise_filter_(false), filter_gyro_(true),
                     butter_ax_(nh_), butter_ay_(nh_), butter_az_(nh_), butter_wx_(nh_), butter_wy_(nh_), butter_wz_(nh_)
{
    // Initialize state covariance
    P_ = Eigen::Matrix4d::Identity();

    // Initialize measurement noise covariance
    double acc_noise = 0.25; // Default σ_a^2 = 0.5^2 m^2/s^4
    double mag_noise = 0.01; // Default σ_m^2 = 0.1^2 µT^2
    bool param_loaded = false;
    param_loaded |= nh_.getParam("acc_noise", acc_noise);
    param_loaded |= nh_.getParam("mag_noise", mag_noise);
    R_ = Eigen::MatrixXd(6, 6).setZero();
    R_.block<3,3>(0,0) = acc_noise * Eigen::Matrix3d::Identity();
    R_.block<3,3>(3,3) = mag_noise * Eigen::Matrix3d::Identity();
    if (!param_loaded) {
        ROS_WARN("Failed to load acc_noise or mag_noise, using defaults: acc_noise=%.2f, mag_noise=%.2f", acc_noise, mag_noise);
    } else {
        ROS_INFO("Loaded noise parameters: acc_noise=%.2f, mag_noise=%.2f", acc_noise, mag_noise);
    }

    // Initialize process noise covariance (gyroscope)
    double gyro_noise = 0.09; // Default σ_ω^2 = 0.3^2 rad^2/s^2
    if (!nh_.getParam("gyro_noise", gyro_noise)) {
        ROS_WARN("Failed to load gyro_noise, using default: %.2f", gyro_noise);
    } else {
        ROS_INFO("Loaded gyro_noise: %.2f", gyro_noise);
    }
    g_noise_ = Eigen::Matrix3d::Zero();
    g_noise_.diagonal() << gyro_noise, gyro_noise, gyro_noise;

    // Reference vectors (NED or ENU)
    std::string coordinate_frame;
    if (!nh_.getParam("coordinate_frame", coordinate_frame)) {
        coordinate_frame = "NED";
        ROS_WARN("Failed to load coordinate_frame, using default: %s", coordinate_frame.c_str());
    }
    double mag_ref_x = 29.14, mag_ref_y = -4.46, mag_ref_z = 45.00; // Defaults for NED, Gongju-si
    param_loaded = false;
    param_loaded |= nh_.getParam("magnetic_reference_x", mag_ref_x);
    param_loaded |= nh_.getParam("magnetic_reference_y", mag_ref_y);
    param_loaded |= nh_.getParam("magnetic_reference_z", mag_ref_z);
    if (coordinate_frame == "ENU") {
        a_ref_ = Eigen::Vector3d(0.0, 0.0, -9.81); // Gravity in ENU (m/s²)
        m_ref_ = Eigen::Vector3d(mag_ref_x, mag_ref_y, -mag_ref_z); // Magnetic field in ENU
    } else {
        a_ref_ = Eigen::Vector3d(0.0, 0.0, 9.81); // Gravity in NED (m/s²)
        m_ref_ = Eigen::Vector3d(mag_ref_x, mag_ref_y, mag_ref_z); // Magnetic field in NED
    }
    mag_norm_ref_ = m_ref_.norm(); // Compute expected magnetic field norm
    m_ref_.normalize();
    if (!param_loaded) {
        ROS_WARN("Failed to load magnetic_reference parameters, using defaults: [%.2f, %.2f, %.2f] uT", mag_ref_x, mag_ref_y, mag_ref_z);
    }
    ROS_INFO("Loaded coordinate_frame: %s, magnetic_reference: [%.2f, %.2f, %.2f] uT, norm: %.2f uT",
             coordinate_frame.c_str(), mag_ref_x, mag_ref_y, mag_ref_z, mag_norm_ref_);

    // Initialize hard iron bias
    mag_bias_ = Eigen::Vector3d(0.0, 0.0, 0.0);
    param_loaded = false;
    param_loaded |= nh_.getParam("mag_bias_x", mag_bias_(0));
    param_loaded |= nh_.getParam("mag_bias_y", mag_bias_(1));
    param_loaded |= nh_.getParam("mag_bias_z", mag_bias_(2));
    if (!param_loaded) {
        ROS_WARN("Failed to load mag_bias parameters, using defaults: [%.2f, %.2f, %.2f] uT", mag_bias_(0), mag_bias_(1), mag_bias_(2));
    }

    // Initialize quaternion
    q_ = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0);

    // Check if magnetometer should be used
    if (!nh_.getParam("use_magnetometer", use_magnetometer_)) {
        use_magnetometer_ = true;
        ROS_WARN("Failed to load use_magnetometer, using default: %s", use_magnetometer_ ? "true" : "false");
    }

    // Check if Butterworth filter should be used
    if (!nh_.getParam("noise_filter", noise_filter_)) {
        noise_filter_ = false;
        ROS_WARN("Failed to load noise_filter, using default: %s", noise_filter_ ? "true" : "false");
    }
    // Check if gyroscope data should be filtered
    if (!nh_.getParam("filter_gyro", filter_gyro_)) {
        filter_gyro_ = true;
        ROS_WARN("Failed to load filter_gyro, using default: %s", filter_gyro_ ? "true" : "false");
    }
    ROS_INFO("Loaded filter parameters: noise_filter=%s, filter_gyro=%s",
             noise_filter_ ? "true" : "false", filter_gyro_ ? "true" : "false");

    // Setup ROS subscribers and publishers
    std::string imu_topic = "/imu/data";
    std::string mag_topic = "/mag";
    param_loaded = false;
    param_loaded |= nh_.getParam("imu_topic", imu_topic);
    param_loaded |= nh_.getParam("mag_topic", mag_topic);
    imu_sub_.subscribe(nh_, imu_topic, 10);
    mag_sub_.subscribe(nh_, mag_topic, 10);
    sync_.registerCallback(boost::bind(&EKF_IMU::syncedCallback, this, _1, _2));
    imu_pub_ = nh_.advertise<sensor_msgs::Imu>("/imu/data_ekf", 10);
    accel_comp_pub_ = nh_.advertise<sensor_msgs::Imu>("/imu/accel_compensated", 10);
    gravity_pub_ = nh_.advertise<sensor_msgs::Imu>("/imu/gravity_ekf", 10);
    if (!param_loaded) {
        ROS_WARN("Failed to load imu_topic or mag_topic, using defaults: imu_topic=%s, mag_topic=%s", imu_topic.c_str(), mag_topic.c_str());
    }
    ROS_INFO("Subscribed to IMU topic: %s, Magnetometer topic: %s", imu_topic.c_str(), mag_topic.c_str());

    ROS_INFO("Initializing IMU EKF node with magnetometer support and hard iron bias compensation...");
    ROS_INFO("use_magnetometer: %s, mag_bias: [%.2f, %.2f, %.2f] uT, noise_filter: %s, filter_gyro: %s",
             use_magnetometer_ ? "true" : "false", mag_bias_(0), mag_bias_(1), mag_bias_(2), noise_filter_ ? "true" : "false", filter_gyro_ ? "true" : "false");
}

Eigen::Matrix3d EKF_IMU::skew(const Eigen::Vector3d& x)
{
    Eigen::Matrix3d S;
    S << 0.0, -x(2), x(1),
         x(2), 0.0, -x(0),
         -x(1), x(0), 0.0;
    return S;
}

Eigen::Matrix3d EKF_IMU::q2R(const Eigen::Vector4d& q)
{
    Eigen::Vector4d q_norm = q / q.norm();
    double qw = q_norm(0), qx = q_norm(1), qy = q_norm(2), qz = q_norm(3);
    Eigen::Matrix3d R;
    R << qw*qw + qx*qx - qy*qy - qz*qz, 2.0 * (qx*qy - qw*qz), 2.0 * (qx*qz + qw*qy),
         2.0 * (qx*qy + qw*qz), qw*qw - qx*qx + qy*qy - qz*qz, 2.0 * (qy*qz - qw*qx),
         2.0 * (qx*qz - qw*qy), 2.0 * (qw*qx + qy*qz), qw*qw - qx*qx - qy*qy + qz*qz;
    return R;
}

Eigen::Vector4d EKF_IMU::ecompass(const Eigen::Vector3d& acc, const Eigen::Vector3d& mag)
{
    Eigen::Vector3d a = acc / acc.norm();
    Eigen::Vector3d m = mag / mag.norm();
    Eigen::Vector3d e1 = (a.cross(m)).cross(a);
    e1.normalize();
    Eigen::Vector3d e2 = a.cross(m);
    e2.normalize();
    Eigen::Vector3d e3 = a;
    e3.normalize();
    Eigen::Matrix3d C;
    C << e1, e2, e3;
    double trace = C.trace();
    Eigen::Vector4d q;
    q(0) = 0.5 * std::sqrt(trace + 1.0);
    q(1) = 0.5 * (C(2,1) - C(1,2) > 0 ? 1 : -1) * std::sqrt(std::max(0.0, C(0,0) - C(1,1) - C(2,2) + 1.0));
    q(2) = 0.5 * (C(0,2) - C(2,0) > 0 ? 1 : -1) * std::sqrt(std::max(0.0, C(1,1) - C(2,2) - C(0,0) + 1.0));
    q(3) = 0.5 * (C(1,0) - C(0,1) > 0 ? 1 : -1) * std::sqrt(std::max(0.0, C(2,2) - C(0,0) - C(1,1) + 1.0));
    q.normalize();
    return q;
}

Eigen::MatrixXd EKF_IMU::Omega(const Eigen::Vector3d& x)
{
    Eigen::MatrixXd O(4, 4);
    O << 0.0, -x(0), -x(1), -x(2),
         x(0), 0.0, x(2), -x(1),
         x(1), -x(2), 0.0, x(0),
         x(2), x(1), -x(0), 0.0;
    return O;
}

Eigen::Vector4d EKF_IMU::f(const Eigen::Vector4d& q, const Eigen::Vector3d& omega, double dt)
{
    Eigen::MatrixXd Omega_t = Omega(omega);
    return (Eigen::Matrix4d::Identity() + 0.5 * dt * Omega_t) * q;
}

Eigen::Matrix4d EKF_IMU::dfdx(const Eigen::Vector3d& omega, double dt)
{
    Eigen::Vector3d x = 0.5 * dt * omega;
    return Eigen::Matrix4d::Identity() + Omega(x);
}

Eigen::VectorXd EKF_IMU::h(const Eigen::Vector4d& q)
{
    Eigen::Matrix3d C = q2R(q).transpose();
    if (use_magnetometer_ && has_recent_mag_) {
        Eigen::VectorXd y(6);
        y.head(3) = C * a_ref_ / a_ref_.norm();
        y.tail(3) = C * m_ref_;
        return y;
    }
    return C * a_ref_ / a_ref_.norm();
}

Eigen::MatrixXd EKF_IMU::dhdx(const Eigen::Vector4d& q)
{
    double qw = q(0), qx = q(1), qy = q(2), qz = q(3);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(use_magnetometer_ && has_recent_mag_ ? 6 : 3);
    if (use_magnetometer_ && has_recent_mag_) {
        v << a_ref_ / a_ref_.norm(), m_ref_;
    } else {
        v << a_ref_ / a_ref_.norm();
    }
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(use_magnetometer_ && has_recent_mag_ ? 6 : 3, 4);
    H.block(0,0,3,4) << v(0)*qw + v(1)*qz - v(2)*qy, v(0)*qx + v(1)*qy + v(2)*qz, -v(0)*qy + v(1)*qx - v(2)*qw, -v(0)*qz + v(1)*qw + v(2)*qx,
                        -v(0)*qz + v(1)*qw + v(2)*qx, v(0)*qy - v(1)*qx + v(2)*qw, v(0)*qx + v(1)*qy + v(2)*qz, -v(0)*qw - v(1)*qz + v(2)*qy,
                        v(0)*qy - v(1)*qx + v(2)*qw, v(0)*qz - v(1)*qw - v(2)*qx, v(0)*qw + v(1)*qz - v(2)*qy, v(0)*qx + v(1)*qy + v(2)*qz;
    if (use_magnetometer_ && has_recent_mag_) {
        H.block(3,0,3,4) << v(3)*qw + v(4)*qz - v(5)*qy, v(3)*qx + v(4)*qy + v(5)*qz, -v(3)*qy + v(4)*qx - v(5)*qw, -v(3)*qz + v(4)*qw + v(5)*qx,
                           -v(3)*qz + v(4)*qw + v(5)*qx, v(3)*qy - v(4)*qx + v(5)*qw, v(3)*qx + v(4)*qy + v(5)*qz, -v(3)*qw - v(4)*qz + v(5)*qy,
                           v(3)*qy - v(4)*qx + v(5)*qw, v(3)*qz - v(4)*qw - v(5)*qx, v(3)*qw + v(4)*qz - v(5)*qy, v(3)*qx + v(4)*qy + v(5)*qz;
    }
    H *= 2.0;
    return H;
}

Eigen::Vector4d EKF_IMU::update(const sensor_msgs::Imu& imu_msg, const sensor_msgs::MagneticField* mag_msg, double dt)
{
    // Extract accelerometer data
    double ax = imu_msg.linear_acceleration.x;
    double ay = imu_msg.linear_acceleration.y;
    double az = imu_msg.linear_acceleration.z;
    double a_norm = std::sqrt(ax * ax + ay * ay + az * az);
    if (update_count_ % 10 == 0) {
        ROS_INFO("Accelerometer: [%.2f, %.2f, %.2f], Norm: %.2f", ax, ay, az, a_norm);
    }
    if (a_norm < 1e-6 || std::abs(a_norm - 9.81) > 0.3 * 9.81) {
        ROS_WARN("Invalid accelerometer reading, norm not close to 9.81 m/s^2");
        return q_;
    }
    Eigen::Vector3d a(ax / a_norm, ay / a_norm, az / a_norm);

    // Extract gyroscope data
    Eigen::Vector3d g(imu_msg.angular_velocity.x, imu_msg.angular_velocity.y, imu_msg.angular_velocity.z);
    if (update_count_ % 10 == 0) {
        ROS_INFO("Gyroscope: [%.2f, %.2f, %.2f]", g(0), g(1), g(2));
    }

    // Prepare measurement vector
    z_ = Eigen::VectorXd(use_magnetometer_ && mag_msg && has_recent_mag_ ? 6 : 3);
    z_.head(3) = a;
    if (use_magnetometer_ && mag_msg && has_recent_mag_) {
        double mx = (mag_msg->magnetic_field.x * 1e6) - mag_bias_(0); // Convert T to µT
        double my = (mag_msg->magnetic_field.y * 1e6) - mag_bias_(1);
        double mz = (mag_msg->magnetic_field.z * 1e6) - mag_bias_(2);
        double m_norm = std::sqrt(mx * mx + my * my + mz * mz);
        if (update_count_ % 10 == 0) {
            ROS_INFO("Magnetometer (bias-compensated): [%.2f, %.2f, %.2f], Norm: %.2f", mx, my, mz, m_norm);
        }
        if (m_norm < 1e-6 || std::abs(m_norm - mag_norm_ref_) > 0.3 * mag_norm_ref_) {
            ROS_WARN("Invalid magnetometer reading, norm not close to %.2f uT, skipping magnetometer update", mag_norm_ref_);
            z_.resize(3);
            has_recent_mag_ = false;
        } else {
            z_.tail(3) = Eigen::Vector3d(mx / m_norm, my / m_norm, mz / m_norm);
        }
    }

    // Update measurement noise covariance
    R_ = Eigen::MatrixXd(z_.size(), z_.size()).setZero();
    R_.block(0,0,3,3) = (0.5 * 0.5) * Eigen::Matrix3d::Identity();
    if (z_.size() == 6) {
        R_.block(3,3,3,3) = (0.1 * 0.1) * Eigen::Matrix3d::Identity();
    }

    // Prediction
    Eigen::Vector4d q_t = f(q_, g, dt);
    q_t.normalize();
    Eigen::Matrix4d F = dfdx(g, dt);
    Eigen::MatrixXd W = Eigen::MatrixXd(4, 3);
    double qw = q_(0), qx = q_(1), qy = q_(2), qz = q_(3);
    W << -qx, -qy, -qz,
         qw, -qz, qy,
         qz, qw, -qx,
         -qy, qx, qw;
    W *= 0.5 * dt;
    Eigen::Matrix4d Q_t = W * g_noise_ * W.transpose();
    P_ = F * P_ * F.transpose() + Q_t;

    // Correction
    Eigen::VectorXd y = h(q_t);
    Eigen::VectorXd v = z_ - y;
    Eigen::MatrixXd H = dhdx(q_t);
    Eigen::MatrixXd S = H * P_ * H.transpose() + R_;
    Eigen::MatrixXd K = P_ * H.transpose() * S.inverse();
    P_ = (Eigen::Matrix4d::Identity() - K * H) * P_;
    q_ = q_t + K * v;
    q_.normalize();

    if (update_count_ % 10 == 0) {
        ROS_INFO("State updated: Quaternion = [%.2f, %.2f, %.2f, %.2f]", q_(0), q_(1), q_(2), q_(3));
    }
    return q_;
}

void EKF_IMU::syncedCallback(const sensor_msgs::Imu::ConstPtr& imu_msg, const sensor_msgs::MagneticField::ConstPtr& mag_msg)
{
    // Check for valid timestamp
    if (!imu_msg->header.stamp.isValid()) {
        ROS_WARN("Received IMU message with invalid timestamp, skipping");
        return;
    }

    // Apply Butterworth filter if enabled
    sensor_msgs::Imu filtered_imu = *imu_msg;
    if (noise_filter_) {
        filtered_imu.linear_acceleration.x = butter_ax_.apply(imu_msg->linear_acceleration.x);
        filtered_imu.linear_acceleration.y = butter_ay_.apply(imu_msg->linear_acceleration.y);
        filtered_imu.linear_acceleration.z = butter_az_.apply(imu_msg->linear_acceleration.z);
        if (filter_gyro_) {
            filtered_imu.angular_velocity.x = butter_wx_.apply(imu_msg->angular_velocity.x);
            filtered_imu.angular_velocity.y = butter_wy_.apply(imu_msg->angular_velocity.y);
            filtered_imu.angular_velocity.z = butter_wz_.apply(imu_msg->angular_velocity.z);
        }
        if (update_count_ % 10 == 0) {
            ROS_INFO("Applied Butterworth filter: Accel=[%.2f, %.2f, %.2f], Gyro=[%.2f, %.2f, %.2f]",
                     filtered_imu.linear_acceleration.x, filtered_imu.linear_acceleration.y, filtered_imu.linear_acceleration.z,
                     filtered_imu.angular_velocity.x, filtered_imu.angular_velocity.y, filtered_imu.angular_velocity.z);
        }
    }

    bool valid_mag = mag_msg && mag_msg->header.stamp.isValid();
    if (use_magnetometer_ && valid_mag) {
        double time_diff = std::abs((imu_msg->header.stamp - mag_msg->header.stamp).toSec());
        if (time_diff < 0.02) { // 20ms tolerance
            has_recent_mag_ = true;
            last_mag_msg_ = *mag_msg;
            last_mag_time_ = mag_msg->header.stamp;
        } else {
            ROS_WARN("Magnetometer message too old (diff: %.2fs), processing IMU only", time_diff);
            valid_mag = false;
            has_recent_mag_ = false;
        }
    } else {
        has_recent_mag_ = false;
    }

    // Initialize on first valid message
    if (!initialized_ || !last_time_.isValid()) {
        last_time_ = imu_msg->header.stamp;
        last_imu_msg_ = filtered_imu;
        if (valid_mag) {
            last_mag_msg_ = *mag_msg;
            last_mag_time_ = mag_msg->header.stamp;
        }
        double ax = filtered_imu.linear_acceleration.x;
        double ay = filtered_imu.linear_acceleration.y;
        double az = filtered_imu.linear_acceleration.z;
        double a_norm = std::sqrt(ax * ax + ay * ay + az * az);
        ROS_INFO("Initial accelerometer: [%.2f, %.2f, %.2f], Norm: %.2f", ax, ay, az, a_norm);
        if (a_norm > 1e-6 && std::abs(a_norm - 9.81) < 0.3 * 9.81) {
            if (use_magnetometer_ && valid_mag) {
                double mx = (mag_msg->magnetic_field.x * 1e6) - mag_bias_(0); // Convert T to µT
                double my = (mag_msg->magnetic_field.y * 1e6) - mag_bias_(1);
                double mz = (mag_msg->magnetic_field.z * 1e6) - mag_bias_(2);
                double m_norm = std::sqrt(mx * mx + my * my + mz * mz);
                ROS_INFO("Initial magnetometer (bias-compensated): [%.2f, %.2f, %.2f], Norm: %.2f", mx, my, mz, m_norm);
                if (m_norm > 1e-6 && std::abs(m_norm - mag_norm_ref_) < 0.3 * mag_norm_ref_) {
                    q_ = ecompass(Eigen::Vector3d(ax, ay, az), Eigen::Vector3d(mx, my, mz));
                    initialized_ = true;
                    ROS_INFO("Initialization successful with magnetometer");
                } else {
                    ROS_WARN("Magnetometer initialization failed, norm not close to %.2f uT, falling back to accelerometer-only", mag_norm_ref_);
                    use_magnetometer_ = false;
                }
            }
            if (!use_magnetometer_) {
                ax /= a_norm;
                ay /= a_norm;
                az /= a_norm;
                double ex = std::atan2(-ay, -az); // Roll for NED
                double ey = std::atan2(ax, std::sqrt(ay * ay + az * az)); // Pitch for NED
                double cx2 = std::cos(ex / 2.0);
                double sx2 = std::sin(ex / 2.0);
                double cy2 = std::cos(ey / 2.0);
                double sy2 = std::sin(ey / 2.0);
                q_ = Eigen::Vector4d(cx2 * cy2, sx2 * cy2, cx2 * sy2, -sx2 * sy2);
                q_.normalize();
                initialized_ = true;
                ROS_INFO("Initialization successful with accelerometer only");
            }
        } else {
            ROS_WARN("Initialization failed: invalid accelerometer data");
            // Publish raw data for debugging
            sensor_msgs::Imu orientation_msg = filtered_imu;
            orientation_msg.header.frame_id = "imu_ekf";
            imu_pub_.publish(orientation_msg);
            sensor_msgs::Imu accel_comp_msg = filtered_imu;
            accel_comp_msg.header.frame_id = "imu_accel_compensated";
            accel_comp_pub_.publish(accel_comp_msg);
            sensor_msgs::Imu gravity_msg = filtered_imu;
            gravity_msg.header.frame_id = "imu_gravity_ekf";
            gravity_msg.linear_acceleration.x = ax;
            gravity_msg.linear_acceleration.y = ay;
            gravity_msg.linear_acceleration.z = az;
            gravity_pub_.publish(gravity_msg);
            return;
        }
        return;
    }

    // Compute dt
    double dt = (imu_msg->header.stamp - last_time_).toSec();
    if (dt <= 0.0 || dt > 0.2) {
        ROS_WARN("Invalid dt: %.2fs, resetting timestamp", dt);
        last_time_ = imu_msg->header.stamp;
        last_imu_msg_ = filtered_imu;
        if (valid_mag) {
            last_mag_msg_ = *mag_msg;
            last_mag_time_ = mag_msg->header.stamp;
        }
        return;
    }
    if (update_count_ % 10 == 0) {
        ROS_INFO("dt: %.2fs", dt);
    }
    last_time_ = imu_msg->header.stamp;
    last_imu_msg_ = filtered_imu;
    if (valid_mag) {
        last_mag_msg_ = *mag_msg;
        last_mag_time_ = mag_msg->header.stamp;
    }

    q_ = update(filtered_imu, use_magnetometer_ && valid_mag && has_recent_mag_ ? &(*mag_msg) : nullptr, dt);

    // Compute gravity-compensated acceleration and gravity estimate
    Eigen::Vector3d a_measured(filtered_imu.linear_acceleration.x, filtered_imu.linear_acceleration.y, filtered_imu.linear_acceleration.z);
    Eigen::Matrix3d R = q2R(q_);
    Eigen::Matrix3d R_body_to_world = R.transpose();
    Eigen::Vector3d a_compensated = a_measured - R * a_ref_;
    Eigen::Vector3d g_estimate = R_body_to_world * a_ref_;
    if (update_count_ % 10 == 0) {
        ROS_INFO("Gravity-compensated accel: [%.2f, %.2f, %.2f]", a_compensated(0), a_compensated(1), a_compensated(2));
        ROS_INFO("Gravity estimate (world): [%.2f, %.2f, %.2f]", g_estimate(0), g_estimate(1), g_estimate(2));
    }

    // Publish orientation IMU message
    sensor_msgs::Imu orientation_msg = filtered_imu;
    orientation_msg.header.frame_id = "imu_ekf";
    orientation_msg.orientation.w = q_(0);
    orientation_msg.orientation.x = q_(1);
    orientation_msg.orientation.y = q_(2);
    orientation_msg.orientation.z = q_(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            orientation_msg.orientation_covariance[i * 3 + j] = P_(i, j);
    orientation_msg.angular_velocity = filtered_imu.angular_velocity;
    orientation_msg.linear_acceleration = filtered_imu.linear_acceleration;
    imu_pub_.publish(orientation_msg);
    if (update_count_ % 10 == 0) {
        ROS_INFO("Published orientation message");
    }

    // Publish gravity-compensated acceleration IMU message
    sensor_msgs::Imu accel_comp_msg = filtered_imu;
    accel_comp_msg.header.frame_id = "imu_accel_compensated";
    accel_comp_msg.orientation.w = q_(0);
    accel_comp_msg.orientation.x = q_(1);
    accel_comp_msg.orientation.y = q_(2);
    accel_comp_msg.orientation.z = q_(3);
    accel_comp_msg.linear_acceleration.x = a_compensated(0);
    accel_comp_msg.linear_acceleration.y = a_compensated(1);
    accel_comp_msg.linear_acceleration.z = a_compensated(2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            accel_comp_msg.orientation_covariance[i * 3 + j] = P_(i, j);
    accel_comp_msg.angular_velocity = filtered_imu.angular_velocity;
    accel_comp_msg.linear_acceleration_covariance = filtered_imu.linear_acceleration_covariance;
    accel_comp_pub_.publish(accel_comp_msg);
    if (update_count_ % 10 == 0) {
        ROS_INFO("Published gravity-compensated acceleration message");
    }

    // Publish gravity estimate IMU message
    sensor_msgs::Imu gravity_msg = filtered_imu;
    gravity_msg.header.frame_id = "imu_gravity_ekf";
    gravity_msg.orientation.w = q_(0);
    gravity_msg.orientation.x = q_(1);
    gravity_msg.orientation.y = q_(2);
    gravity_msg.orientation.z = q_(3);
    gravity_msg.linear_acceleration.x = g_estimate(0);
    gravity_msg.linear_acceleration.y = g_estimate(1);
    gravity_msg.linear_acceleration.z = g_estimate(2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            gravity_msg.orientation_covariance[i * 3 + j] = P_(i, j);
    gravity_msg.angular_velocity = filtered_imu.angular_velocity;
    gravity_msg.linear_acceleration_covariance = filtered_imu.linear_acceleration_covariance;
    gravity_pub_.publish(gravity_msg);
    if (update_count_ % 10 == 0) {
        ROS_INFO("Published gravity estimate message");
    }

    update_count_++;
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "imu_ekf_node");
    ROS_INFO("Initializing IMU EKF node...");
    EKF_IMU ekf;
    ros::spin();
    return 0;
}