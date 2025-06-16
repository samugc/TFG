#ifndef REGRESSION_HPP
#define REGRESSION_HPP

#include <Eigen/Dense>
#include <vector>

// Realiza regresión lineal y devuelve {pendiente, intercepto}
std::vector<double> linearRegression(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

// Construye la matriz de Vandermonde para regresión polinomial
Eigen::MatrixXd buildVandermonde(const Eigen::VectorXd& x, int degree);

// Realiza regresión polinomial y devuelve los coeficientes del polinomio
Eigen::VectorXd polynomialRegression(const Eigen::VectorXd& x, const Eigen::VectorXd& y, int degree);

#endif // REGRESSION_HPP