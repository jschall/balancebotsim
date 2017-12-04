#include <Eigen/Dense>

Eigen::Matrix<double,19,19> get_mm(Eigen::Matrix<double,19,1> states, Eigen::Matrix<double,2,1> inputs);
Eigen::Matrix<double,19,1> get_fo(Eigen::Matrix<double,19,1> states, Eigen::Matrix<double,2,1> inputs);