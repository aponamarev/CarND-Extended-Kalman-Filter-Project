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
                cout << "RMSE: " << d.transpose() << endl;
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
  /**
  TODO:
    * Calculate a Jacobian here.
  */
}
