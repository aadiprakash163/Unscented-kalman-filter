#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3-n_aug_;

  n_sigma_ = 2*n_aug_ + 1;

  weights_ = VectorXd(n_sigma_);

  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

  time_us_ = 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  // Initialize matrices with first measurement
  if(!is_initialized_){
    
    cout<<"Intializing the filter........."<<endl;

    P_ << 1,0,0,0,0,
          0,1,0,0,0,
          0,0,1,0,0,
          0,0,0,1,0,
          0,0,0,0,1;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER){

      x_ << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0;

      // Deal with zero px and py
      if(fabs(x_(0)) < 0.001 and fabs(x_(1)) < 0.001){
        x_(0) = 0.001;
        x_(1) = 0.001;
      }
    } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      
      // Convert radar measurement to states
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      float px = rho*cos(phi);
      float py = rho*sin(phi);
      float vx = rho_dot*cos(phi);
      float vy = rho_dot*sin(phi);
      float v = sqrt(vx*vx + vy*vy);

      x_ << px, py, v, 0, 0;

    }

    // Set the time stamp
    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    // cout<<"Filter Initialized with following state: "<<endl<<x_<<endl;

    // Predict and Update is not done now
    return;
  }

  
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  cout<<"Predicting state............"<<endl;
  Prediction(delta_t);

  // cout<<"State is predicted, calling respective function"<<endl;

  if(use_radar_ and meas_package.sensor_type_ == MeasurementPackage::RADAR){
    cout<<"Received a RADAR measurement package"<<endl;

    UpdateRadar(meas_package);
    

  }

  if(use_laser_ and meas_package.sensor_type_ == MeasurementPackage::LASER){
    cout<<"Received a LASER measurement package"<<endl;
    UpdateLidar(meas_package);

  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Make augmented state mean and covariance matrices
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  MatrixXd Xsigma_aug = MatrixXd(n_aug_, n_sigma_);
  MatrixXd L = P_aug.llt().matrixL();

  // cout<<"L :"<<endl<<L<<endl;

  // Xsigma_aug_.fill(0);
  Xsigma_aug.col(0) = x_aug;

  for(int i=0; i<n_aug_ ;i++){

    Xsigma_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsigma_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // cout<<"Xsigma_aug: "<<endl<<Xsigma_aug<<endl;


  //predict sigma points
  for (int i = 0; i< n_sigma_; i++)
  {
    //extract values for better readability
    double p_x = Xsigma_aug(0,i);
    double p_y = Xsigma_aug(1,i);
    double v = Xsigma_aug(2,i);
    double yaw = Xsigma_aug(3,i);
    double yawd = Xsigma_aug(4,i);
    double nu_a = Xsigma_aug(5,i);
    double nu_yawdd = Xsigma_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + (v/yawd) * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + (v/yawd) * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  x_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  cout<<"Mean State Value after prediction:"<<endl<<x_<<endl;
  // x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

  }
  cout<<"Covariance matrx after prediction: "<<P_<<endl;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Lidar measure only two variales
  int n_z = 2;
  MatrixXd Zsigma = Xsig_pred_.topRows(n_z);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
      z_pred = z_pred + weights_(i) * Zsigma.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  
    //residual
    VectorXd z_diff = Zsigma.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
   
  S = S + R;

  // Update with RADAR measurement
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  

    //residual
    VectorXd z_diff = Zsigma.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  // Kalman Gain
  MatrixXd K = Tc * S.inverse();


  VectorXd z = meas_package.raw_measurements_;
  //residual
  VectorXd z_diff = z - z_pred;

  // Calculate NIS value
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  cout<<"NIS laser: "<<NIS_laser_<<endl;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
    

}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Predict Radar measurement first
  int n_z = 3;
  MatrixXd Zsigma = MatrixXd(n_z, n_sigma_);

  for(int i = 0; i<n_sigma_;i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsigma(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsigma(1,i) = atan2(p_y,p_x);                                 //phi
    Zsigma(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  // cout<<"1"<<endl;

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
      z_pred = z_pred + weights_(i) * Zsigma.col(i);
  }

  // cout<<"2"<<endl;

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsigma.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // cout<<"3"<<endl;

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  S = S + R;

  // cout<<"4"<<endl;

  // Update with RADAR measurement
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsigma.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // cout<<"5"<<endl;

  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // Calculate NIS value
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  
  cout<<"NIS radar: "<<NIS_radar_<<endl;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  // cout<<"6"<<endl;  

}

