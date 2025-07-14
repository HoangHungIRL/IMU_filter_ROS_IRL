#include "eskf_imu.h"
#include <cmath>

ImuEskf::ImuEskf()
    : Node("imu_eskf_node"),
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
    this->declare_parameter("acc_noise", 0.0008);
    this->declare_parameter("mag_noise", 0.00000001);
    this->declare_parameter("gyro_noise", 0.00008);
    this->declare_parameter("gyro_bias_noise", 0.0000000007);
    this->declare_parameter("theta_noise", 0.00001);
    this->declare_parameter("wb_noise", 0.0000000001);
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
    acc_noise_ = this->get_parameter("acc_noise").as_double();
    mag_noise_ = this->get_parameter("mag_noise").as_double();
    gyro_noise_ = this->get_parameter("gyro_noise").as_double();
    gyro_bias_noise_ = this->get_parameter("gyro_bias_noise").as_double();
    theta_noise_ = this->get_parameter("theta_noise").as_double();
    wb_noise_ = this->get_parameter("wb_noise").as_double();
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

    // Initialize state and covariance
    x_error_ = Eigen::VectorXd::Zero(n_state_);
    x_nominal_ = Eigen::VectorXd::Zero(n_state_ + 1);
    x_nominal_(0) = 1.0;
    z_ = Eigen::VectorXd::Zero(6);
    P_ = Eigen::MatrixXd::Identity(n_state_, n_state_);
    Fx_ = Eigen::MatrixXd::Zero(n_state_, n_state_);
    H_ = Eigen::MatrixXd::Zero(6, n_state_);
    Fi_ = Eigen::MatrixXd::Identity(n_state_, n_state_);
    Q_ = Eigen::MatrixXd::Identity(n_state_, n_state_);
    R_ = Eigen::MatrixXd::Identity(6, 6);

    // Set covariance matrices
    P_.block<3, 3>(0, 0) = theta_noise_ * Eigen::Matrix3d::Identity();
    P_.block<3, 3>(3, 3) = wb_noise_ * Eigen::Matrix3d::Identity();
    R_.block<3, 3>(0, 0) = acc_noise_ * Eigen::Matrix3d::Identity();
    R_.block<3, 3>(3, 3) = mag_noise_ * Eigen::Matrix3d::Identity();
    Q_.block<3, 3>(0, 0) = gyro_noise_ * Eigen::Matrix3d::Identity();
    Q_.block<3, 3>(3, 3) = gyro_bias_noise_ * Eigen::Matrix3d::Identity();

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
    sync_->registerCallback(&ImuEskf::syncedCallback, this);
    imu_pub_ = create_publisher<sensor_msgs::msg::Imu>("/imu/data_ekf", 10);
    accel_comp_pub_ = create_publisher<sensor_msgs::msg::Imu>("/imu/accel_compensated", 10);
    gravity_pub_ = create_publisher<sensor_msgs::msg::Imu>("/imu/gravity_ekf", 10);

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

    // Initialize time
    last_time_ = rclcpp::Time(0, 0, clock_->get_clock_type());
    last_mag_time_ = rclcpp::Time(0, 0, clock_->get_clock_type());

    RCLCPP_INFO(get_logger(), "Initializing IMU ESKF node with magnetometer support, hard iron bias compensation, and Butterworth filter...");
    RCLCPP_INFO(get_logger(), "use_sim_time: %s, use_magnetometer: %s, noise_filter: %s, filter_gyro: %s, coordinate_frame: %s, mag_bias: [%.2f, %.2f, %.2f] uT",
                use_sim_time ? "true" : "false", use_magnetometer_ ? "true" : "false", noise_filter_ ? "true" : "false",
                filter_gyro_ ? "true" : "false", coordinate_frame.c_str(), mag_bias_(0), mag_bias_(1), mag_bias_(2));
}

Eigen::Matrix3d ImuEskf::skew(const Eigen::Vector3d& x) {
    Eigen::Matrix3d S;
    S << 0.0, -x(2), x(1),
         x(2), 0.0, -x(0),
         -x(1), x(0), 0.0;
    return S;
}

Eigen::Matrix3d ImuEskf::q2R(const Eigen::Vector4d& q) {
    Eigen::Vector4d q_norm = q / q.norm();
    double qw = q_norm(0), qx = q_norm(1), qy = q_norm(2), qz = q_norm(3);
    Eigen::Matrix3d R;
    R << qw*qw + qx*qx - qy*qy - qz*qz, 2.0 * (qx*qy - qw*qz), 2.0 * (qx*qz + qw*qy),
         2.0 * (qx*qy + qw*qz), qw*qw - qx*qx + qy*qy - qz*qz, 2.0 * (qy*qz - qw*qx),
         2.0 * (qx*qz - qw*qy), 2.0 * (qw*qx + qy*qz), qw*qw - qx*qx - qy*qy + qz*qz;
    return R;
}

Eigen::Vector4d ImuEskf::ecompass(const Eigen::Vector3d& acc, const Eigen::Vector3d& mag) {
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

Eigen::MatrixXd ImuEskf::diffQstarvqQ(const Eigen::Quaterniond& q, const Eigen::Vector3d& v) {
    double q0 = q.w();
    Eigen::Vector3d qv = q.vec();
    Eigen::MatrixXd D(3, 4);
    D.col(0) = 2 * (q0 * v + skew(v) * qv);
    D.block<3, 3>(0, 1) = 2 * (-v * qv.transpose() + qv * v.transpose() +
                               v.dot(qv) * Eigen::Matrix3d::Identity() + q0 * skew(v));
    return D;
}

void ImuEskf::predict() {
    Eigen::Vector3d wb_error = x_error_.segment<3>(3);
    Eigen::Vector3d wb_nominal = x_nominal_.segment<3>(4);
    Eigen::Vector3d wm = raw_gyro_;
    Eigen::Quaterniond q_nominal;
    q_nominal.w() = x_nominal_(0);
    q_nominal.vec() = x_nominal_.segment<3>(1);

    Fx_.setZero();
    Eigen::Quaterniond w_q;
    w_q.w() = 1;
    w_q.vec() = 0.5 * (wm - wb_nominal) * dt_;
    q_nominal *= w_q;
    q_nominal.normalize();
    x_nominal_(0) = q_nominal.w();
    x_nominal_.segment<3>(1) = q_nominal.vec();

    Eigen::Vector3d omega = (wm - wb_error) * dt_;
    double theta = omega.norm();
    Eigen::Vector3d n = omega / (theta + 1e-10);
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() * cos(theta) +
                        (1 - cos(theta)) * n * n.transpose() + sin(theta) * skew(n);
    Fx_.block<3, 3>(0, 0) = R.transpose();
    Fx_.block<3, 3>(0, 3) = -Eigen::Matrix3d::Identity() * dt_;
    Fx_.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
    x_error_ = Fx_ * x_error_;
    P_ = Fx_ * P_ * Fx_.transpose() + Fi_ * Q_ * dt_ * Fi_.transpose();
}

void ImuEskf::measurement() {
    Eigen::Quaterniond q_nominal;
    q_nominal.w() = x_nominal_(0);
    q_nominal.vec() = x_nominal_.segment<3>(1);
    int meas_size = (use_magnetometer_ && has_recent_mag_) ? 6 : 3;
    z_ = Eigen::VectorXd(meas_size);
    H_ = Eigen::MatrixXd::Zero(meas_size, n_state_);
    z_.segment<3>(0) = q_nominal.inverse() * (a_ref_ / a_ref_.norm());
    if (use_magnetometer_ && has_recent_mag_) {
        Eigen::Vector3d mag_global = q_nominal * raw_mag_;
        double mag_norm = mag_global.norm();
        if (mag_norm > 1e-6) {
            mag_global /= mag_norm;
        } else {
            mag_global.setZero();
        }
        z_.segment<3>(3) = q_nominal.inverse() * mag_global;
    }

    Eigen::MatrixXd h1(meas_size, 7);
    h1.setZero();
    h1.block<3, 4>(0, 0) = diffQstarvqQ(q_nominal, a_ref_ / a_ref_.norm());
    if (use_magnetometer_ && has_recent_mag_) {
        h1.block<3, 4>(3, 0) = diffQstarvqQ(q_nominal, m_ref_);
    }
    Eigen::MatrixXd h2 = Eigen::MatrixXd::Zero(7, n_state_);
    double w = q_nominal.w(), x = q_nominal.x(), y = q_nominal.y(), z = q_nominal.z();
    Eigen::Matrix<double, 4, 3> q_left;
    q_left << -x, -y, -z,
              w, -z, y,
              z, w, -x,
             -y, x, w;
    h2.block<4, 3>(0, 0) = 0.5 * q_left;
    h2.block<3, 3>(4, 3) = Eigen::Matrix3d::Identity();
    H_ = h1 * h2;

    R_ = Eigen::MatrixXd::Identity(meas_size, meas_size);
    R_.block<3, 3>(0, 0) = acc_noise_ * Eigen::Matrix3d::Identity();
    if (meas_size == 6) {
        R_.block<3, 3>(3, 3) = mag_noise_ * Eigen::Matrix3d::Identity();
    }
}

void ImuEskf::error2nominal() {
    Eigen::Quaterniond q_nominal;
    q_nominal.w() = x_nominal_(0);
    q_nominal.vec() = x_nominal_.segment<3>(1);
    Eigen::Vector3d omega = x_error_.segment<3>(0);
    double omega_norm = omega.norm();
    Eigen::Quaterniond delt_q;
    delt_q.w() = cos(omega_norm / 2);
    if (omega_norm > 1e-10) {
        delt_q.vec() = (omega / omega_norm * sin(omega_norm / 2)).eval();
    } else {
        delt_q.vec() = Eigen::Vector3d::Zero();
    }
    q_nominal = q_nominal * delt_q;
    q_nominal.normalize();
    x_nominal_(0) = q_nominal.w();
    x_nominal_.segment<3>(1) = q_nominal.vec();
    x_nominal_.segment<3>(4) += x_error_.segment<3>(3);
}

void ImuEskf::estimate() {
    int meas_size = (use_magnetometer_ && has_recent_mag_) ? 6 : 3;
    Eigen::MatrixXd K = P_ * H_.transpose() * (H_ * P_ * H_.transpose() + R_).inverse();
    Eigen::VectorXd z_meas(meas_size);
    z_meas.segment<3>(0) = raw_acc_;
    if (use_magnetometer_ && has_recent_mag_) {
        z_meas.segment<3>(3) = raw_mag_;
    }
    x_error_ = K * (z_meas - z_);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_state_, n_state_);
    P_ = (I - K * H_) * P_;
    error2nominal();
    x_error_.setZero();
    P_ = Eigen::MatrixXd::Identity(n_state_, n_state_) * P_ * Eigen::MatrixXd::Identity(n_state_, n_state_).transpose();
}

Observation ImuEskf::getData() {
    Observation data;
    data.quat.w() = x_nominal_(0);
    data.quat.vec() = x_nominal_.segment<3>(1);
    data.euler = Eigen::Vector3d::Zero();
    return data;
}

double ImuEskf::getTime() const {
    return current_t_;
}

void ImuEskf::syncedCallback(const sensor_msgs::msg::Imu::SharedPtr imu_msg, const sensor_msgs::msg::MagneticField::SharedPtr mag_msg) {
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

    double dt = (rclcpp::Time(imu_msg->header.stamp, clock_->get_clock_type()) - last_time_).seconds();
    if (last_time_ != rclcpp::Time(0, 0, clock_->get_clock_type()) && (dt <= 0.0 || dt > 0.2)) {
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
    dt_ = dt;
    current_t_ = rclcpp::Time(imu_msg->header.stamp, clock_->get_clock_type()).seconds();
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
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "Accelerometer: [%.2f, %.2f, %.2f], Norm: %.2f", ax, ay, az, a_norm);
    }
    if (a_norm < 1e-6 || std::abs(a_norm - 9.81) > 0.3 * 9.81) {
        RCLCPP_WARN(get_logger(), "Invalid accelerometer reading, norm not close to 9.81 m/s^2");
        return;
    }
    raw_acc_ << ax / a_norm, ay / a_norm, az / a_norm;
    raw_gyro_ << filtered_imu.angular_velocity.x, filtered_imu.angular_velocity.y, filtered_imu.angular_velocity.z;
    if (valid_mag) {
        double mx = (mag_msg->magnetic_field.x * 1e6) - mag_bias_(0);
        double my = (mag_msg->magnetic_field.y * 1e6) - mag_bias_(1);
        double mz = (mag_msg->magnetic_field.z * 1e6) - mag_bias_(2);
        double m_norm = std::sqrt(mx * mx + my * my + mz * mz);
        if (update_count_ % 10 == 0) {
            RCLCPP_INFO(get_logger(), "Magnetometer (bias-compensated): [%.2f, %.2f, %.2f], Norm: %.2f", mx, my, mz, m_norm);
        }
        if (m_norm < 1e-6) {
            RCLCPP_WARN(get_logger(), "Invalid magnetometer reading, norm too small, skipping magnetometer update");
            has_recent_mag_ = false;
        } else {
            raw_mag_ << mx / m_norm, my / m_norm, mz / m_norm;
        }
    }

    if (!initialized_ || !euler_initialized_) {
        if (a_norm > 1e-6 && std::abs(a_norm - 9.81) < 0.3 * 9.81) {
            if (use_magnetometer_ && valid_mag && has_recent_mag_) {
                double mx = (mag_msg->magnetic_field.x * 1e6) - mag_bias_(0);
                double my = (mag_msg->magnetic_field.y * 1e6) - mag_bias_(1);
                double mz = (mag_msg->magnetic_field.z * 1e6) - mag_bias_(2);
                double m_norm = std::sqrt(mx * mx + my * my + mz * mz);
                if (m_norm > 1e-6) {
                    Eigen::Vector4d q = ecompass(Eigen::Vector3d(ax, ay, az), Eigen::Vector3d(mx, my, mz));
                    x_nominal_(0) = q(0);
                    x_nominal_.segment<3>(1) = q.segment<3>(1);
                    euler_initialized_ = true;
                    initialized_ = true;
                    RCLCPP_INFO(get_logger(), "Initialization successful with magnetometer");
                } else {
                    RCLCPP_WARN(get_logger(), "Magnetometer initialization failed, norm too small, falling back to accelerometer-only");
                    use_magnetometer_ = false;
                }
            }
            if (!use_magnetometer_) {
                double ex = std::atan2(-ay / a_norm, -az / a_norm);
                double ey = std::atan2(ax / a_norm, std::sqrt(ay * ay + az * az) / a_norm);
                double cx2 = std::cos(ex / 2.0);
                double sx2 = std::sin(ex / 2.0);
                double cy2 = std::cos(ey / 2.0);
                double sy2 = std::sin(ey / 2.0);
                x_nominal_(0) = cx2 * cy2;
                x_nominal_(1) = sx2 * cy2;
                x_nominal_(2) = cx2 * sy2;
                x_nominal_(3) = -sx2 * sy2;
                Eigen::Quaterniond q(x_nominal_(0), x_nominal_(1), x_nominal_(2), x_nominal_(3));
                q.normalize();
                x_nominal_(0) = q.w();
                x_nominal_.segment<3>(1) = q.vec();
                euler_initialized_ = true;
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

    // Run ESKF
    predict();
    measurement();
    estimate();

    // Compute gravity-compensated acceleration and gravity estimate
    Eigen::Vector3d a_measured(filtered_imu.linear_acceleration.x, filtered_imu.linear_acceleration.y, filtered_imu.linear_acceleration.z);
    Eigen::Matrix3d R = q2R(x_nominal_.head(4));
    Eigen::Matrix3d R_body_to_world = R.transpose();
    Eigen::Vector3d a_compensated = a_measured - R * a_ref_;
    Eigen::Vector3d g_estimate = R_body_to_world * a_ref_;
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "Gravity-compensated accel: [%.2f, %.2f, %.2f]", a_compensated(0), a_compensated(1), a_compensated(2));
        RCLCPP_INFO(get_logger(), "Gravity estimate (world): [%.2f, %.2f, %.2f]", g_estimate(0), g_estimate(1), g_estimate(2));
    }

    // Publish orientation IMU message
    Observation data = getData();
    auto orientation_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
    orientation_msg->header.frame_id = "imu_ekf";
    orientation_msg->orientation.x = data.quat.x();
    orientation_msg->orientation.y = data.quat.y();
    orientation_msg->orientation.z = data.quat.z();
    orientation_msg->orientation.w = data.quat.w();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            orientation_msg->orientation_covariance[i * 3 + j] = P_(i, j);
    orientation_msg->angular_velocity = filtered_imu.angular_velocity;
    orientation_msg->linear_acceleration = filtered_imu.linear_acceleration;
    imu_pub_->publish(std::move(orientation_msg));
    if (update_count_ % 10 == 0) {
        RCLCPP_INFO(get_logger(), "Published orientation message");
    }

    // Publish gravity-compensated acceleration IMU message
    auto accel_comp_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
    accel_comp_msg->header.frame_id = "imu_accel_compensated";
    accel_comp_msg->orientation.x = data.quat.x();
    accel_comp_msg->orientation.y = data.quat.y();
    accel_comp_msg->orientation.z = data.quat.z();
    accel_comp_msg->orientation.w = data.quat.w();
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

    // Publish gravity estimate IMU message
    auto gravity_msg = std::make_unique<sensor_msgs::msg::Imu>(filtered_imu);
    gravity_msg->header.frame_id = "imu_gravity_ekf";
    gravity_msg->orientation.x = data.quat.x();
    gravity_msg->orientation.y = data.quat.y();
    gravity_msg->orientation.z = data.quat.z();
    gravity_msg->orientation.w = data.quat.w();
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

int main(int argc, char** argv) {
    rclcpp::init(argc, argv);
    auto node = std::make_shared<ImuEskf>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}