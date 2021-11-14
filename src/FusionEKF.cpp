#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1.0, 0.0, 0.0, 0.0,		//section 10, lesson 5 - H matrix projects your belief about objects current state into the measurement space of the sensor. 
              0.0, 1.0, 0.0, 0.0;		// since lidar only meaures position, we are keeping only px and py as 1 

  Hj_ << 1, 1, 0, 0,			// this is H matrix for radar 
         1, 1, 0, 0,
         1, 1, 1, 1; 

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,		//state transition matrix
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
  
  //measurement noise
//   noise_ax_ = 9.0;  //values from section 13 lesson 5 
//   noise_ay_ = 9.0;
  
  //Initializing fusion EKF
  

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) { 
  //takes the measurement package we create in main.cpp
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1; //first two values are overwritten from state vector (first measurement), tweak the last 2 elements if RMSE is very high 

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      // px, py values from RADAR
      float ro     = measurement_pack.raw_measurements_(0);
      float phi    = measurement_pack.raw_measurements_(1);
      float ro_dot = measurement_pack.raw_measurements_(2);
      ekf_.x_(0) = ro     * cos(phi); // x
      ekf_.x_(1) = ro     * sin(phi); // y 
      ekf_.x_(2) = ro_dot * cos(phi); // vx
      ekf_.x_(3) = ro_dot * sin(phi); // vy     

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state. // px, py values from LASER
      ekf_.x_(0) = measurement_pack.raw_measurements_(0); // this is x 
      ekf_.x_(1) = measurement_pack.raw_measurements_(1); // this is y 
    }
    
    ekf_.F_ << 1, 0, 0, 0,//to be modified later
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//time difference dt is in seconds. dt calculated for ekf_.F_ matrix
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated , section 8 lesson 5 
  //state transition matrix 
  ekf_.F_(0, 2) = dt; // replace dt at two positions in F matrix which was defined as unit diagnoal matrix 4x4
  ekf_.F_(1, 3) = dt;

  //set the acceleration noise components
  float noise_ax = 9;
  float noise_ay = 9;

  //set the process covariance matrix Q
  //section 9 lesson 5 
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
             0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
             dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;  
  
  //finish
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    Tools tools;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_; //calculated in the jacobian in tools.cpp
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_; //uncertainty in laser measurement 
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
