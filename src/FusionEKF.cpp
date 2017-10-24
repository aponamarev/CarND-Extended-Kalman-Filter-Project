#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    
    previous_timestamp_ = 0;
    
    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    H_laser_ <<
    1,0,0,0,
    0,1,0,0;
    Hj_ = MatrixXd(3, 4);
    
    //* Finish initializing the FusionEKF.
    // Normally state covariance matrix should be calculated based on observations
    // using the following formula: Dataset - Unity Matrix * Dataset * 1/n_observations
    // however, in the absence of the data, we can initialize state covariance with some
    // default values.
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ <<
    1,0,0,0,
    0,1,0,0,
    0,0,1000,0,
    0,0,0,1000;
    
    ekf_.F_ = Eigen::MatrixXd().Identity(4,4);
    
    // Set the process and measurement noises
    //measurement covariance matrix - laser
    R_laser_ <<
    0.0225, 0,
    0, 0.0225;
    
    //measurement covariance matrix - radar
    R_radar_ <<
    0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;
    
    
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

Eigen::MatrixXd FusionEKF::update_Q(float dt, Eigen::VectorXd a) {
    
    Eigen::MatrixXd Q = Eigen::MatrixXd().Zero(4, 4);
    
    if (dt==0 || a.size()!=2) {
        cout << "Error: incorrect inputs. dt should be greater than 0 (dt value = " << dt
        << ") and 'a' size should be equal to 2 (a size = " << a.size() << ").";
        return Q;
    }
    
    float dt2_05 = 0.5*dt*dt;
    Eigen::MatrixXd Qv = a * a.transpose();
    Qv(0,1) = 0;
    Qv(1,0) = 0;
    Eigen::MatrixXd G = Eigen::MatrixXd(4,2);
    G <<
    dt2_05,0,
    0, dt2_05,
    dt,0,
    0,dt;
    Q = G * Qv * G.transpose();
    return Q;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            // Define measurements explicitly
            float r = measurement_pack.raw_measurements_[0];
            float p = measurement_pack.raw_measurements_[1];
            
            // Convert polar measurements into cartesian
            float px = r * cos(p);
            float py = r * sin(p);
            // Assign EKF measurements: x' = f(x) + v
            ekf_.x_ << px,py,0,0;
        }
        
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            ekf_.x_ <<
            measurement_pack.raw_measurements_[0],
            measurement_pack.raw_measurements_[1], 0,0;
        }
    
        // done initializing, no need to predict or update
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }
    
    
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    // - Time is measured in seconds.
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // /=10^6 converts to miliseconds to seconds
    previous_timestamp_ = measurement_pack.timestamp_;
    // * Update the state transition matrix F according to the new elapsed time.
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    // * Update the process noise covariance matrix.
    // * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    Eigen::VectorXd a = Eigen::VectorXd(2);
    a << 9, 9;
    
    ekf_.Q_ = update_Q(dt, a);
    
    ekf_.Predict();
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    // Use the sensor type to perform the update step.
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        //Update the state and covariance matrices.
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        //Update the state and covariance matrices.
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    // print the output
    //cout << "x_ = " << ekf_.x_.transpose() << endl;
    //cout << "P_ = " << endl << ekf_.P_ << endl;
}
