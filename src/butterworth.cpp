#include "butterworth.h"

Butter2::Butter2(double fc, double fs)
    : initalised(false), fc_(fc), fs_(fs)
{
    update_parameters();
}

void Butter2::update_parameters()
{
    tan_ = tan(M_PI * fc_ / fs_);
    gamma = 1 / tan_;
    alphak = 2 * cos(2 * M_PI * (n + 1) / (4 * n));
    gain_ = gamma * gamma - alphak * gamma + 1;
    a_[0] = -(gamma * gamma + alphak * gamma + 1) / gain_;
    a_[1] = -(2 - 2 * gamma * gamma) / gain_;
}

double Butter2::apply(double sample) {
    if (!initalised) {
        initalised = true;
        return reset(sample);
    }
    xs_[0] = xs_[1];
    xs_[1] = xs_[2];
    xs_[2] = sample / gain_;
    ys_[0] = ys_[1];
    ys_[1] = ys_[2];
    ys_[2] = (xs_[0] + xs_[2]) + 2 * xs_[1] + (a_[0] * ys_[0]) + (a_[1] * ys_[1]);
    return ys_[2];
}

double Butter2::reset(double sample) {
    xs_[0] = sample;
    xs_[1] = sample;
    xs_[2] = sample;
    ys_[0] = sample;
    ys_[1] = sample;
    ys_[2] = sample;
    return sample;
}