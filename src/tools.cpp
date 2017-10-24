#include <iostream>
#include "tools.h"

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
                d = d.array().sqrt();
                RMSE += d;
            }
            
            RMSE /= estimations.size();
            
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
    
    // Calculate a Jacobian here.
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //TODO: YOUR CODE HERE
    float px2 = px*px;
    float py2 = py*py;
    float p2sum = px2+py2;
    //check division by zero
    if (p2sum<0.00001) {
        cout << "Wrong input resulting in div/0." << endl;
        Hj <<
        0,0,0,0,
        0,0,0,0,
        0,0,0,0;
    } else {
        
        float sqrt_p2sum = sqrt(p2sum);
        float vxpy = vx*py;
        float vypx = vy*px;
        
        Hj <<
        px/sqrt_p2sum,py/sqrt_p2sum,0,0,
        -py/p2sum,px/p2sum,0,0,
        py*(vxpy-vypx)/pow(p2sum,3/2),
        px*(vypx-vxpy)/pow(p2sum,3/2),px/sqrt_p2sum,py/sqrt_p2sum;
        //compute the Jacobian matrix
    };
    
    return Hj;
}
