#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

#include <ros/ros.h>
#include <cmath>

class Butter2 {
public:
    Butter2(ros::NodeHandle& nh, double fc = 12.0, double fs = 100.0);
    double apply(double sample);
    double reset(double sample);

private:
    bool initalised;
    double a_[2];
    double gain_;
    double xs_[3];
    double ys_[3];
    const int n = 2;
    double fc_;
    double fs_;
    double tan_;
    double gamma;
    double alphak;
};

#endif // BUTTERWORTH_H