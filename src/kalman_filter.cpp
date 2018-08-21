#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

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

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
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

  // calculate z_pred in polar coordinates
  float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
  
  float rho = sqrt(px*px + py*py);
  if (rho < 1e-7) { 
    cout << "UpdateEKF() - Error - Division by Zero" << endl;
    return;
	}
  float phi = atan2(py,px);
  float rho_dot = (px*vx + py*vy)/rho;
  
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  cout << z << "" << z_pred << endl;
  
	VectorXd y = z - z_pred;
  y(1) = NormalizeAngle(y(1));
  /**
  if (y(1) < -M_PI) {
    y(1) += 2*M_PI;
  }
  else if (y(1) > M_PI) { 
    y(1) -= 2*M_PI;
  }
  */
  cout << z.transpose() << " " << z_pred.transpose() << " " << y.transpose() << endl;

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

double KalmanFilter::NormalizeAngle(double theta) {
  
  if (theta > M_PI) {
    int n_2pi = (theta+M_PI)/(2*M_PI);
    theta -= n_2pi*2*M_PI;
  }
  if (theta < -M_PI) { 
    int n_2pi = abs(theta-M_PI)/(2*M_PI);
    theta += n_2pi*2*M_PI;
  }
  
  return theta;
}
