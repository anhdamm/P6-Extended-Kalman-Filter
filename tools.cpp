#include "tools.h"
#include <iostream>
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;
using std::sqrt;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  size_t est_size = estimations.size();
  size_t gt_size = ground_truth.size();

  if(est_size == 0 || est_size != gt_size) {
	  cout << "Invalid estimation or ground_truth data!" << endl;
	  return rmse;

  }

  // accumulate squared errors
  VectorXd error(4);
  error << 0, 0, 0, 0;
  for(size_t i = 0 ; i < est_size ; i ++) {
	  error = estimations[i] - ground_truth[i];
	  // coefficient-wise multiplication
	  error = error.array() * error.array();
	  rmse += error;

  } 

  // calculate the mean
  rmse = rmse / est_size;
  // calculate the squared root
  rmse = rmse.array().sqrt();
  // return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // Initialize the Jacobian Matrix Hj
  MatrixXd Hj(3, 4);
  Hj << 0, 0, 0, 0,
	      0, 0, 0, 0,
		    0, 0, 0, 0;

  // pre-compute a set of terms to avoid repeated calculation
  double denom1 = sqrt(px * px + py * py);
  double denom2 = px * px + py * py;
  double denom3 = denom1 * denom2;

  if (denom1 == 0) {
	  cout << "CalculateJacobian()-Error-The px and py cannot be zero simultaneously." << endl;
	  return Hj;

  }

  // compute the Jacobian matrix
  Hj << px / denom1, py / denom1, 0, 0,
	      -py / denom2, px / denom2, 0, 0,
		    py * (vx * py - vy * px) / denom3, px * (vy * px - vx * py) / denom3, px / denom1, py / denom1;
  // return the result
  return Hj;

}