#ifndef ESTIMATIONS_HPP
#define ESTIMATIONS_HPP

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include "io_handler.hpp"
#include "model_data.hpp"
#include "regression.hpp"

double obtainMu0();
double obtainSigma0();
Eigen::VectorXd arithmeticMean();
Eigen::VectorXd geometricMean();
Eigen::VectorXd vectorSigma2();
Eigen::VectorXd valuesPolymomialRegression();
double estimateSigma();
Eigen::VectorXd estimateEtaBeta(int degree);

#endif // ESTIMATIONS_HPP
