#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  
  x_ = F_ * x_; //section 8 lesson 5 
  MatrixXd Ft = F_.transpose(); //section 9 lesson 5 
  P_ = F_ * P_ * Ft + Q_;  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  //section 7 lesson 5 
  //KF measurement update step
  VectorXd z_pred = H_ * x_; //z_prediction 
  VectorXd y = z - z_pred; // z - z prediction = transition state
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  //section 14, lesson 5 h(x')
  float rho = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  float rho_dot = (px*vx + py*vy)/rho;
  
  
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  
  VectorXd y = z - z_pred; //we calculate z_pred (h) by ourselves in updateEKF
  
  // normalize y between pi and -pi
    while(y(1) > M_PI)
    {
      y(1) -= 2 * M_PI;
    }

    while(y(1) < -M_PI)
    {
      y(1) += 2 * M_PI;
    }
  

  
  //section 7 lesson 5 identical to the function before
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;  
  
 
}
