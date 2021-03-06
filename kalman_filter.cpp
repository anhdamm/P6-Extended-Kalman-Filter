#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::sqrt;
using std::atan2;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */

   VectorXd y = z - H_ * x_;
   MatrixXd S = H_ * P_ * H_.transpose() + R_;
   MatrixXd K = P_ * H_.transpose() * S.inverse();
   x_ = x_ + K * y;
   MatrixXd I;
   I = MatrixXd::Identity(x_.size(), x_.size());
   P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

   double px = x_(0);
   double py = x_(1);
   double vx = x_(2);
   double vy = x_(3);

   double rho_predicted = sqrt(px * px + py * py);
   double phi_predicted = atan2(py, px);
   double rho_rate_predicted = (px * vx + py * vy) / rho_predicted;

   VectorXd hx(3);
   hx << rho_predicted, phi_predicted, rho_rate_predicted;
   VectorXd y = z - hx;

   //Normalizing Angles between -Pi and Pi
   double Pi = 3.14159265358979323846;

   while(y(1) > Pi) {
     y(1) = y(1) - 2 * Pi;

   }

   while(y(1) < -Pi) {
     y(1) = y(1) + 2 * Pi;

   }

   MatrixXd S = H_ * P_ * H_.transpose() + R_;
   MatrixXd K = P_ * H_.transpose() * S.inverse();

   x_ = x_ + K * y;
   MatrixXd I;
   I = MatrixXd::Identity(x_.size(), x_.size());
   P_ = (I - K * H_) * P_;
   
}