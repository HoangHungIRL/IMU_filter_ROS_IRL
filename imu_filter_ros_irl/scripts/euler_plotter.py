#!/usr/bin/env python3

import rospy
from sensor_msgs.msg import Imu
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Explicitly set TkAgg backend
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math
import csv
import os
from datetime import datetime
import threading

class EulerPlotter:
    def __init__(self):
        # Initialize ROS 1 node
        rospy.init_node('euler_plotter', anonymous=True)
        rospy.loginfo('Initializing Euler Plotter node')

        # Get parameters
        self.use_sim_time = rospy.get_param('~use_sim_time', False)
        rospy.loginfo(f'use_sim_time: {self.use_sim_time}')

        # Initialize data storage (lists for all data, no limit)
        self.times = []
        self.rolls = []
        self.pitches = []
        self.yaws = []
        self.start_time = None

        # Subscribe to /imu/data_ekf
        self.imu_sub = rospy.Subscriber(
            '/imu/data_ekf',
            Imu,
            self.imu_callback,
            queue_size=10
        )

        # Setup matplotlib figure for real-time plotting
        self.fig, (self.ax_roll, self.ax_pitch, self.ax_yaw) = plt.subplots(3, 1, figsize=(10, 8))
        self.fig.suptitle('Real-Time Euler Angles from EKF IMU (All Data)')
        self.ax_roll.set_ylabel('Roll (deg)')
        self.ax_pitch.set_ylabel('Pitch (deg)')
        self.ax_yaw.set_ylabel('Yaw (deg)')
        self.ax_yaw.set_xlabel('Time (s)')
        self.line_roll, = self.ax_roll.plot([], [], 'r-', label='Roll')
        self.line_pitch, = self.ax_pitch.plot([], [], 'g-', label='Pitch')
        self.line_yaw, = self.ax_yaw.plot([], [], 'b-', label='Yaw')
        self.ax_roll.legend()
        self.ax_pitch.legend()
        self.ax_yaw.legend()
        self.ax_roll.grid(True)
        self.ax_pitch.grid(True)
        self.ax_yaw.grid(True)

        # Initialize plot limits
        self.ax_roll.set_ylim(-180, 180)
        self.ax_pitch.set_ylim(-180, 180)
        self.ax_yaw.set_ylim(-180, 180)
        self.ax_roll.set_xlim(0, 10)  # Initial x-axis limit, will adjust dynamically
        self.ax_pitch.set_xlim(0, 10)
        self.ax_yaw.set_xlim(0, 10)

        # Register shutdown hook
        rospy.on_shutdown(self.save_data_and_plot)

        # Start animation
        self.ani = FuncAnimation(self.fig, self.update_plot, interval=100, blit=True)

    def quaternion_to_euler(self, w, x, y, z):
        """Convert quaternion (w, x, y, z) to Euler angles (roll, pitch, yaw) in NED frame."""
        # Normalize quaternion
        norm = math.sqrt(w**2 + x**2 + y**2 + z**2)
        if norm == 0:
            rospy.logwarn('Received invalid quaternion with zero norm')
            return 0.0, 0.0, 0.0
        w, x, y, z = w/norm, x/norm, y/norm, z/norm

        # Roll (x-axis rotation)
        sinr_cosp = 2.0 * (w * x + y * z)
        cosr_cosp = 1.0 - 2.0 * (x * x + y * y)
        roll = math.atan2(sinr_cosp, cosr_cosp)

        # Pitch (y-axis rotation)
        sinp = 2.0 * (w * y - z * x)
        if abs(sinp) >= 1:
            pitch = math.copysign(math.pi / 2, sinp)  # Handle singularity
        else:
            pitch = math.asin(sinp)

        # Yaw (z-axis rotation)
        siny_cosp = 2.0 * (w * z + x * y)
        cosy_cosp = 1.0 - 2.0 * (y * y + z * z)
        yaw = math.atan2(siny_cosp, cosy_cosp)

        # Convert to degrees
        roll = math.degrees(roll)
        pitch = math.degrees(pitch)
        yaw = math.degrees(yaw)

        return roll, pitch, yaw

    def imu_callback(self, msg):
        """Callback for /imu/data_ekf messages."""
        if msg.header.stamp == rospy.Time(0):
            rospy.logwarn('Received IMU message with invalid timestamp, skipping')
            return

        # Convert quaternion to Euler angles
        w = msg.orientation.w
        x = msg.orientation.x
        y = msg.orientation.y
        z = msg.orientation.z
        roll, pitch, yaw = self.quaternion_to_euler(w, x, y, z)

        # Get time
        current_time = rospy.get_rostime()
        if self.start_time is None:
            self.start_time = current_time
        time_s = (current_time - self.start_time).to_sec()

        # Store data
        self.times.append(time_s)
        self.rolls.append(roll)
        self.pitches.append(pitch)
        self.yaws.append(yaw)

        # Log Euler angles
        if len(self.times) % 10 == 0:
            rospy.loginfo(f'Roll: {roll:.2f}°, Pitch: {pitch:.2f}°, Yaw: {yaw:.2f}°')

    def update_plot(self, frame):
        """Update the real-time plot with all data."""
        # Check if data is available
        if not self.times or not self.yaws:
            return self.line_roll, self.line_pitch, self.line_yaw

        # Update data
        self.line_roll.set_data(self.times, self.rolls)
        self.line_pitch.set_data(self.times, self.pitches)
        self.line_yaw.set_data(self.times, self.yaws)

        # Adjust x-axis limits dynamically to show all data
        max_time = max(self.times) if self.times else 1.0
        self.ax_roll.set_xlim(0, max(max_time, 1))  # Ensure at least 1s for initial plot
        self.ax_pitch.set_xlim(0, max(max_time, 1))
        self.ax_yaw.set_xlim(0, max(max_time, 1))

        return self.line_roll, self.line_pitch, self.line_yaw

    def save_data_and_plot(self):
        """Save all Euler angles to CSV and generate a final plot on shutdown."""
        if not self.times:
            rospy.logwarn('No data to save.')
            return

        # Generate timestamp for filenames
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_dir = os.path.expanduser('~/euler_plotter_output')
        os.makedirs(output_dir, exist_ok=True)

        # Save data to CSV
        csv_file = os.path.join(output_dir, 'euler_angles_{}.csv'.format(timestamp))
        try:
            with open(csv_file, 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['Time (s)', 'Roll (deg)', 'Pitch (deg)', 'Yaw (deg)'])
                for t, r, p, y in zip(self.times, self.rolls, self.pitches, self.yaws):
                    writer.writerow([t, r, p, y])
            rospy.loginfo('Saved Euler angles to {}'.format(csv_file))
        except Exception as e:
            rospy.logerr('Failed to save CSV: {}'.format(e))

        # Generate and save final plot
        fig, (ax_roll, ax_pitch, ax_yaw) = plt.subplots(3, 1, figsize=(10, 8))
        fig.suptitle('All Euler Angles from EKF IMU')
        ax_roll.plot(self.times, self.rolls, 'r-', label='Roll')
        ax_pitch.plot(self.times, self.pitches, 'g-', label='Pitch')
        ax_yaw.plot(self.times, self.yaws, 'b-', label='Yaw')
        ax_roll.set_ylabel('Roll (deg)')
        ax_pitch.set_ylabel('Pitch (deg)')
        ax_yaw.set_ylabel('Yaw (deg)')
        ax_yaw.set_xlabel('Time (s)')
        ax_roll.legend()
        ax_pitch.legend()
        ax_yaw.legend()
        ax_roll.grid(True)
        ax_pitch.grid(True)
        ax_yaw.grid(True)
        ax_roll.set_ylim(-180, 180)
        ax_pitch.set_ylim(-180, 180)
        ax_yaw.set_ylim(-180, 180)
        if self.times:
            ax_roll.set_xlim(0, max(self.times))
            ax_pitch.set_xlim(0, max(self.times))
            ax_yaw.set_xlim(0, max(self.times))

        # Save plot
        plot_file = os.path.join(output_dir, 'euler_angles_plot_{}.png'.format(timestamp))
        try:
            fig.savefig(plot_file)
            rospy.loginfo('Saved final plot to {}'.format(plot_file))
            plt.close(fig)
        except Exception as e:
            rospy.logerr('Failed to save plot: {}'.format(e))

def main():
    try:
        plotter = EulerPlotter()
        # Show plot and spin node in separate threads
        def spin():
            rospy.spin()

        spin_thread = threading.Thread(target=spin)
        spin_thread.start()

        # Show plot (blocking)
        plt.show()

        # Wait for spin thread to finish
        spin_thread.join()

    except rospy.ROSInterruptException:
        pass

if __name__ == '__main__':
    main()