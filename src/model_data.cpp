#include "model_data.hpp"

// Dimensions of paths
int d;
Eigen::VectorXd n;
int N;

// Matrix of time and paths
Eigen::VectorXd t_values;
Eigen::MatrixXd t;
Eigen::MatrixXd x;

// Random variable V
Eigen::MatrixXd v;

// Number of betas
int p;