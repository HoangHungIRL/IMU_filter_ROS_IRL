#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

#include <cmath>

class Butter2 {
public:
    Butter2(double fc = 12.0, double fs = 100.0);
    double apply(double sample);
    double reset(double sample);

private:
    void update_parameters();
    bool initalised;
    double fc_;
    double fs_;
    double tan_;
    double gamma;
    double alphak;
    int n = 2;
    double a_[2];
    double gain_;
    double xs_[3];
    double ys_[3];
};

#endif // BUTTERWORTH_H