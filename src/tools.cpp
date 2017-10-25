#include <iostream>
#include "tools.h"

#define eps 0.0001

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
     Returns RMSE
     */
    VectorXd RMSE(4);
    RMSE << 0,0,0,0;
    
    if (estimations.size()>0) {
        if (estimations.size()==ground_truth.size()) {
            
            for (unsigned int i=0; i<estimations.size(); i++) {
                
                VectorXd d = estimations[i] - ground_truth[i];
                d = d.array() * d.array();
                RMSE += d;
            }
            
            RMSE /= estimations.size();
            RMSE = RMSE.array().sqrt();
            
        } else {
            cout << "Error: Incompatible input size. Estimations size is 0 and it should be greater." << endl;
        }
    }
    else {
        cout << "Error: Incorrect input size. Estimations size is 0 and it should be greater." << endl;
    }
    
    return RMSE;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    
    // Jacobean matrix solution from lecture materials"
    
    MatrixXd Hj(3,4);
    //recover state parameters
    const float px = x_state(0);
    const float py = x_state(1);
    const float vx = x_state(2);
    const float vy = x_state(3);
    
    //pre-compute a set of terms to avoid repeated calculation
    const double c1 = max(float(eps), px*px + py*py);
    const float c2 = sqrt(c1);
    const float c3 = (c1*c2);
    
    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }
    float px_c2 = px/c2;
    float py_c2 = py/c2;
    //compute the Jacobian matrix
    Hj << px_c2, py_c2, 0, 0,
    -(py/c1), (px/c1), 0, 0,
    py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px_c2, py_c2;
    
    return Hj;
}
