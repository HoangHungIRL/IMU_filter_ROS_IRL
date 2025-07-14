#include "ekf_imu.h"
#include <cmath>

EKF_IMU::EKF_IMU()
    : Node("imu_ekf_node"),
      butter_ax_(12.0, 100.0),
      butter_ay_(12.0, 100.0),
      butter_az_(12.0, 100.0),
      butter_wx_(12.0, 100.0),
      butter_wy_(12.0, 100.0),
      butter_wz_(12.0, 100.0)
{
    // Declare parameters
    this->declare_parameter("clock_type", "ros");
    this->declare_parameter("imu_topic", "/imu/data_newIMU");
    this->declare_parameter("mag_topic", "/mag_newIMU");
    this->declare_parameter("use_magnetometer", true);
    this->declare_parameter("noise_filter", true);
    this->declare_parameter("filter_gyro", false);
    this->declare_parameter("fc", 25.0);
    this->declare_parameter("fs", 100.0);
    this->declare_parameter("mag_bias_x", 0.0);
    this->declare_parameter("mag_bias_y", 0.0);
    this->declare_parameter("mag_bias_z", 0.0);
    this->declare_parameter("acc_noise", 0.25);
    this->declare_parameter("mag_noise", 0.01);
    this->declare_parameter("gyro_noise", 0.09);
    this->declare_parameter("magnetic_reference_x", 29.14);
    this->declare_parameter("magnetic_reference_y", -4.46);
    this->declare_parameter("magnetic_reference_z", 45.00);
    this->declare_parameter("coordinate_frame", "NED");

    // Get parameters
    bool use_sim_time;
    this->get_parameter_or("use_sim_time", use_sim_time, false);
    std::string clock_type = this->get_parameter("clock_type").as_string();
    std::string imu_topic = this->get_parameter("imu_topic").as_string();
    std::string mag_topic = this->get_parameter("mag_topic").as_string();
    use_magnetometer_ = this->get_parameter("use_magnetometer").as_bool();
    noise_filter_ = this->get_parameter("noise_filter").as_bool();
    filter_gyro_ = this->get_parameter("filter_gyro").as_bool();
    double fc = this->get_parameter("fc").as_double();
    double fs = this->get_parameter("fs").as_double();
    mag_bias_ = Eigen::Vector3d(
        this->get_parameter("mag_bias_x").as_double(),
        this->get_parameter("mag_bias_y").as_double(),
        this->get_parameter("mag_bias_z").as_double()
    );
    double acc_noise = this->get_parameter("acc_noise").as_double();
    double mag_noise = this->get_parameter("mag_noise").as_double();
    double gyro_noise = this->get_parameter("gyro_noise").as_double();
    double mag_ref_x = this->get_parameter("magnetic_reference_x").as_double();
    double mag_ref_y = this->get_parameter("magnetic_reference_y").as_double();
    double mag_ref_z = this->get_parameter("magnetic_reference_z").as_double();
    std::string coordinate_frame = this->get_parameter("coordinate_frame").as_string();

    // Reinitialize Butterworth filters with parameters
    butter_ax_ = Butter2(fc, fs);
    butter_ay_ = Butter2(fc, fs);
    butter_az_ = Butter2(fc, fs);
    butter_wx_ = Butter2(fc, fs);
    butter_wy_ = Butter2(fc, fs);
    butter_wz_ = Butter2(fc, fs);

    // Initialize clock
    if (clock_type == "system") {
        clock_ = std::make_shared<rclcpp::Clock>(RCL_SYSTEM_TIME);
        RCLCPP_INFO(get_logger(), "Using SystemTime clock");
    } else if (clock_type == "steady") {
        clock_ = std::make_shared<rclcpp::Clock>(RCL_STEADY_TIME);
        RCLCPP_INFO(get_logger(), "Using SteadyTime clock");
    } else {
        clock_ = std::make_shared<rclcpp::Clock>(RCL_ROS_TIME);
        RCLCPP_INFO(get_logger(), "Using ROSTime clock (use_sim_time: %s)", use_sim_time ? "true" : "false");
    }

    // Initialize state
    q_ = Eigen::Vector4d(1.0, 0.0, 0.0, 0.0);
    P_ = Eigen::Matrix4d::Identity();
    R_ = Eigen::MatrixXd(6, 6).setZero();
    R_.block<3,3>(0,0) = acc_noise * Eigen::Matrix3d::Identity();
    R_.block<3,3>(3,3) = mag_noise * Eigen::Matrix3d::Identity();
    g_noise_ = Eigen::Matrix3d::Zero();
    g_noise_.diagonal() << gyro_noise, gyro_noise, gyro_noise;

    // Reference vectors
    if (coordinate_frame == "ENU") {
        a_ref_ = Eigen::Vector3d(0.0, 0.0, -9.81);
        m_ref_ = Eigen::Vector3d(mag_ref_x, mag_ref_y, -mag_ref_z);
    } else {
        a_ref_ = Eigen::Vector3d(0.0, 0.0, 9.81);
        m_ref_ = Eigen::Vector3d(mag_ref_x, mag_ref_y, mag_ref_z);
    }
    mag_norm_ref_ = m_ref_.norm();
    m_ref_.normalize();

    // Setup subscribers and publishers
    imu_sub_ = std::make_shared<message_filters::Subscriber<sensor_msgs::msg::Imu>>(this, imu_topic);
    mag_sub_ = std::make_shared<message_filters::Subscriber<sensor_msgs::msg::MagneticField>>(this, mag_topic);
    sync_ = std::make_shared<message_filters::TimeSynchronizer<sensor_msgs::msg::Imu, sensor_msgs::msg::MagneticField>>(*imu_sub_, *mag_sub_, 10);
    sync_->registerCallback(&EKF_IMU::syncedCallback, this);
    imu_pub_ = create_publisher<sensor_msgs::msg::Imu>("/imu/data_ekf", 10);
    accel_comp_pub_ = create_publisher<sensor_msgs::msg::Imu>("/imu/accel_compensated", 10);
    gravity_pub_ = create_publisher<sensor_msgs::msg::Imu>("/imu/gravity_ekf", 10);

    // Initialize time
    last_time_ = rclcpp::Time(0, 0, clock_->get_clock_type());
    last_mag_time_ = rclcpp::Time(0, 0, clock_->get_clock_type());

    RCLCPP_INFO(get_logger(), "Initializing IMU EKF node with magnetometer support and hard iron bias compensation...");
    RCLCPP_INFO(get_logger(), "use_sim_time: %s, use_magnetometer: %s, noise_filter: %s, filter_gyro: %s, coordinate_frame: %s, mag_bias: [%.2f, %.2f, %.2f] uT",
                use_sim_time ? "true" : "false", use_magnetometer_ ? "true" : "false", noise_filter_ ? "true" : "false",
                filter_gyro_ ? "true" : "false", coordinate_frame.c_str(), mag_bias_(0), mag_bias_(1), mag_bias_(2));
}

void EKF_IMU::syncedCallback(const sensor_msgs::msg::Imu::SharedPtr imu_msg, const sensor_msgs::msg::MagneticField::SharedPtr mag_msg)
{
    if (imu_msg->header.stamp.sec == 0 && imu_msg->header.stamp.nanosec == 0) {
        RCLCPP_WARN(get_logger(), "Received IMU message with invalid timestamp, skipping");
        return;
    }

    sensor_msgs::msg::Imu filtered_imu = *imu_msg;
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
            RCLCPP_INFO(get_logger(), "Applied Butterworth filter: Accel=[%.2f, %.2f, %.2f], Gyro=[%.2f, %.2f, %.2f]",
                        filtered_imu.linear_acceleration.x, filtered_imu.linear_acceleration.y, filtered_imu.linear_acceleration.z,
                        filtered_imu.angular_velocity.x, filtered_imu.angular_velocity.y, filtered_imu.angular_velocity.z);
        }
    }

    bool valid_mag = mag_msg && (mag_msg->header.stamp.sec != 0 || mag_msg->header.stamp.nanosec != 0);
    if (use_magnetometer_ && valid_mag) {
        rclcpp::Time imu_time(imu_msg->header.stamp, clock_->get_clock_type());
        rclcpp::Time mag_time(mag_msg->header.stamp, clock_->get_clock_type());
        double time_diff = std::abs((imu_time - mag_time).seconds());
        if (time_diff < 0.02) {
            has_recent_mag_ = true;
            last_mag_msg_ = *mag_msg;
            last_mag_time_ = mag_time;
        } else {
            RCLCPP_WARN(get_logger(), "Magnetometer message too old (diff: %.2fs), processing IMU only", time_diff);
            valid_mag = false;
            has_recent_mag_ = false;
        }
    } else {
        has_recent_mag_ = false;
    }

    if (!initialized_ || last_time_ == rclcpp::Time(0, 0, clock_->get_clock_type())) {
        last_time_ = rclcpp::Time(imu_msg->header.stamp, clock_->get_clock_type());
        last_imu_msg_ = filtered_imu;
        if (valid_mag) {
            last_mag_msg_ = *mag_msg;
            last_mag_time_ = rclcpp::Time(mag_msg->header.stamp, clock_->get_clock_type());
        }
        double ax = filtered_imu.linear_acceleration.x;
        double ay = filtered_imu.linear_acceleration.y;
        double az = filtered_imu.linear_acceleration.z;
        double a_norm = std::sqrt(ax * ax + ay * ay + az * az);
        RCLCPP_INFO_STREAM(get_logger(), "Initial accelerometer: [" << ax << ", " << ay << ", " << az << "], Norm: " << a_norm);
        if (a_norm > 1e-6 && std::abs(a_norm - 9.81) < 0.3 * 9.81) {
            if (use_magnetometer_ && valid_mag) {
                double mx = (mag_msg->magnetic_field.x * 1e6) - mag_bias_(0);
                double my = (mag_msg->magnetic_field.y * 1e6) - mag_bias_(1);
                double mz = (mag_msg->magnetic_field.z * 1e6) - mag_bias_(2);
                double m_norm = std::sqrt(mx * mx + my * my + mz * mz);
                RCLCPP_INFO_STREAM(get_logger(), "Initial magnetometer (bias-compensated): [" << mx << ", " << my << ", " << mz << "], Norm: " << m_norm);
                if (m_norm > 1e-6 && std::abs(m_norm - mag_norm_ref_) < 0.3 * mag_norm_ref_) {
                    q_ = ecompass(Eigen::Vector3d(ax, ay, az), Eigen::Vector3d(mx, my, mz));
                    initialized_ = true;
                    RCLCPP_INFO(get_logger(), "Initialization successful with magnetometer");
                } else {
                    RCLCPP_WARN(get_logger(), "Magnetometer initialization failed, falling back to accelerometer-only");
                    use_magnetometer_ = false;
                }
            }
            if (!use_magnetometer_) {
                ax /= a_norm;
                ay /= a_norm;
                az /= a_norm;
                double ex = std::atan2(-ay, -az);
                double ey = std::atan2(ax, std::sqrt(ay * ay + az * az));
                double cx2 = std::cos(ex / 2.0);
                double sx2 = std::sin(ex / 2.0);
                double cy2 = std::cos(ey / 2.0);
                double sy2 = std::sin(ey / 2.0);
                q_ = Eigen::Vector4d(cx2 * cy2, sx2 * cy2, cx2 * sy2, -sx2 * sy2);
                q_.normalize();
                initialized_ = true;
                RCLCPP_INFO(get_logger(), "Initialization successful with accelerometer only");
            }
        } else {
            RCLCPP_WARN(get_logger(), "Initialization failed: invalid accelerometer data");
            auto orientation_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
            orientation_msg->header.frame_id = "imu_ekf";
            imu_pub_->publish(std::move(orientation_msg));
            auto accel_comp_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
            accel_comp_msg->header.frame_id = "imu_accel_compensated";
            accel_comp_pub_->publish(std::move(accel_comp_msg));
            auto gravity_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
            gravity_msg->header.frame_id = "imu_gravity_ekf";
            gravity_msg->linear_acceleration.x = ax;
            gravity_msg->linear_acceleration.y = ay;
            gravity_msg->linear_acceleration.z = az;
            gravity_pub_->publish(std::move(gravity_msg));
            return;
        }
        return;
    }

    double dt = (rclcpp::Time(imu_msg->header.stamp, clock_->get_clock_type()) - last_time_).seconds();
    if (dt <= 0.0 || dt > 0.2) {
        RCLCPP_WARN(get_logger(), "Invalid dt: %.2fs, resetting timestamp", dt);
        last_time_ = rclcpp::Time(imu_msg->header.stamp, clock_->get_clock_type());
        last_imu_msg_ = filtered_imu;
        if (valid_mag) {
            last_mag_msg_ = *mag_msg;
            last_mag_time_ = rclcpp::Time(mag_msg->header.stamp, clock_->get_clock_type());
        }
        return;
    }
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "dt: %.2fs", dt);
    }
    last_time_ = rclcpp::Time(imu_msg->header.stamp, clock_->get_clock_type());
    last_imu_msg_ = filtered_imu;
    if (valid_mag) {
        last_mag_msg_ = *mag_msg;
        last_mag_time_ = rclcpp::Time(mag_msg->header.stamp, clock_->get_clock_type());
    }

    q_ = update(filtered_imu, use_magnetometer_ && valid_mag && has_recent_mag_ ? &(*mag_msg) : nullptr, dt);

    Eigen::Vector3d a_measured(filtered_imu.linear_acceleration.x, filtered_imu.linear_acceleration.y, filtered_imu.linear_acceleration.z);
    Eigen::Matrix3d R = q2R(q_);
    Eigen::Matrix3d R_body_to_world = R.transpose();
    Eigen::Vector3d a_compensated = a_measured - R * a_ref_;
    Eigen::Vector3d g_estimate = R_body_to_world * a_ref_;
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "Gravity-compensated accel: [%.2f, %.2f, %.2f]", a_compensated(0), a_compensated(1), a_compensated(2));
        RCLCPP_INFO(get_logger(), "Gravity estimate (world): [%.2f, %.2f, %.2f]", g_estimate(0), g_estimate(1), g_estimate(2));
    }

    auto orientation_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
    orientation_msg->header.frame_id = "imu_ekf";
    orientation_msg->orientation.w = q_(0);
    orientation_msg->orientation.x = q_(1);
    orientation_msg->orientation.y = q_(2);
    orientation_msg->orientation.z = q_(3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            orientation_msg->orientation_covariance[i * 3 + j] = P_(i, j);
    orientation_msg->angular_velocity = filtered_imu.angular_velocity;
    orientation_msg->linear_acceleration = filtered_imu.linear_acceleration;
    imu_pub_->publish(std::move(orientation_msg));
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "Published orientation message");
    }

    auto accel_comp_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
    accel_comp_msg->header.frame_id = "imu_accel_compensated";
    accel_comp_msg->orientation.w = q_(0);
    accel_comp_msg->orientation.x = q_(1);
    accel_comp_msg->orientation.y = q_(2);
    accel_comp_msg->orientation.z = q_(3);
    accel_comp_msg->linear_acceleration.x = a_compensated(0);
    accel_comp_msg->linear_acceleration.y = a_compensated(1);
    accel_comp_msg->linear_acceleration.z = a_compensated(2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            accel_comp_msg->orientation_covariance[i * 3 + j] = P_(i, j);
    accel_comp_msg->angular_velocity = filtered_imu.angular_velocity;
    accel_comp_msg->linear_acceleration_covariance = filtered_imu.linear_acceleration_covariance;
    accel_comp_pub_->publish(std::move(accel_comp_msg));
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "Published gravity-compensated acceleration message");
    }

    auto gravity_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
    gravity_msg->header.frame_id = "imu_gravity_ekf";
    gravity_msg->orientation.w = q_(0);
    gravity_msg->orientation.x = q_(1);
    gravity_msg->orientation.y = q_(2);
    gravity_msg->orientation.z = q_(3);
    gravity_msg->linear_acceleration.x = g_estimate(0);
    gravity_msg->linear_acceleration.y = g_estimate(1);
    gravity_msg->linear_acceleration.z = g_estimate(2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            gravity_msg->orientation_covariance[i * 3 + j] = P_(i, j);
    gravity_msg->angular_velocity = filtered_imu.angular_velocity;
    gravity_msg->linear_acceleration_covariance = filtered_imu.linear_acceleration_covariance;
    gravity_pub_->publish(std::move(gravity_msg));
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "Published gravity estimate message");
    }

    update_count_++;
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

Eigen::Vector4d EKF_IMU::update(const sensor_msgs::msg::Imu& imu_msg, const sensor_msgs::msg::MagneticField* mag_msg, double dt)
{
    double ax = imu_msg.linear_acceleration.x;
    double ay = imu_msg.linear_acceleration.y;
    double az = imu_msg.linear_acceleration.z;
    double a_norm = std::sqrt(ax * ax + ay * ay + az * az);
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO_STREAM(get_logger(), "Accelerometer: [" << ax << ", " << ay << ", " << az << "], Norm: " << a_norm);
    }
    if (a_norm < 1e-6 || std::abs(a_norm - 9.81) > 0.3 * 9.81) {
        RCLCPP_WARN(get_logger(), "Invalid accelerometer reading, norm not close to 9.81 m/s^2");
        return q_;
    }
    Eigen::Vector3d a(ax / a_norm, ay / a_norm, az / a_norm);

    Eigen::Vector3d g(imu_msg.angular_velocity.x, imu_msg.angular_velocity.y, imu_msg.angular_velocity.z);
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO_STREAM(get_logger(), "Gyroscope: [" << g(0) << ", " << g(1) << ", " << g(2) << "]");
    }

    z_ = Eigen::VectorXd(use_magnetometer_ && mag_msg && has_recent_mag_ ? 6 : 3);
    z_.head(3) = a;
    if (use_magnetometer_ && mag_msg && has_recent_mag_) {
        double mx = (mag_msg->magnetic_field.x * 1e6) - mag_bias_(0);
        double my = (mag_msg->magnetic_field.y * 1e6) - mag_bias_(1);
        double mz = (mag_msg->magnetic_field.z * 1e6) - mag_bias_(2);
        double m_norm = std::sqrt(mx * mx + my * my + mz * mz);
        if (update_count_ % 10 == 0) {
            RCLCPP_INFO_STREAM(get_logger(), "Magnetometer (bias-compensated): [" << mx << ", " << my << ", " << mz << "], Norm: " << m_norm);
        }
        if (m_norm < 1e-6 || std::abs(m_norm - mag_norm_ref_) > 0.3 * mag_norm_ref_) {
            RCLCPP_WARN(get_logger(), "Invalid magnetometer reading, norm not close to %.2f uT, skipping magnetometer update", mag_norm_ref_);
            z_.resize(3);
            has_recent_mag_ = false;
        } else {
            z_.tail(3) = Eigen::Vector3d(mx / m_norm, my / m_norm, mz / m_norm);
        }
    }

    R_ = Eigen::MatrixXd(z_.size(), z_.size()).setZero();
    R_.block(0,0,3,3) = (0.5 * 0.5) * Eigen::Matrix3d::Identity();
    if (z_.size() == 6) {
        R_.block(3,3,3,3) = (0.1 * 0.1) * Eigen::Matrix3d::Identity();
    }

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

    Eigen::VectorXd y = h(q_t);
    Eigen::VectorXd v = z_ - y;
    Eigen::MatrixXd H = dhdx(q_t);
    Eigen::MatrixXd S = H * P_ * H.transpose() + R_;
    Eigen::MatrixXd K = P_ * H.transpose() * S.inverse();
    P_ = (Eigen::Matrix4d::Identity() - K * H) * P_;
    q_ = q_t + K * v;
    q_.normalize();

    if (update_count_ % 10 == 0) {
        RCLCPP_INFO_STREAM(get_logger(), "State updated: Quaternion = [" << q_(0) << ", " << q_(1) << ", " << q_(2) << ", " << q_(3) << "]");
    }
    return q_;
}

int main(int argc, char** argv)
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<EKF_IMU>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}