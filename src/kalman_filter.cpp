#include "kalman_filter.h"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

VectorXd KalmanFilter::h(const VectorXd x_prime) {
    /* to calculate y for the radar sensor, we need to convert x′
     ​​  to polar coordinates. In other words, the function h(x) maps values from Cartesian coordinates to polar coordinates. So the equation for radar becomes y=z radar−h(x′).
     */
    float px = x_prime(0);
    float py = x_prime(1);
    float vx = x_prime(2);
    float vy = x_prime(3);
    double rho = sqrt(px*px + py*py);
    double theta = std::atan2(py, px);
    double rho_dot = (px*vx + py*vy) / rho;
    
    VectorXd x_polar = VectorXd(3);
    x_polar << rho, theta, rho_dot;
    return x_polar;
}

void KalmanFilter::Predict() {
    /**
     * predict the state
     */
    x_ = F_ * x_; // + u_
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
     * update the state by using Kalman Filter equations
     */
    const VectorXd y = z - H_ * x_;
    GenericUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /*
     For radar measurements, the functions that map the x vector [px, py, vx, vy] to polar coordinates are non-linear. Instead of using H to calculate y = z - H * x', for radar measurements you'll have to use the equations that map from cartesian to polar coordinates: y = z - h(x').
     */
    
    VectorXd y = z - h(x_);
    // Calculations are essentially the same to the Update function
    GenericUpdate(y);
}

void KalmanFilter::GenericUpdate(const VectorXd &y) {
    /**
     * update the state by using Kalman Filter equations
     */
    MatrixXd PHt_ = P_ * H_.transpose();
    MatrixXd S = H_ * PHt_ + R_; // Total covariance
    MatrixXd K = PHt_ * S.inverse(); // Kalman gain
    // measurement update
    x_ += (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
