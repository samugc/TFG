#include "regression.hpp"

using namespace std;
using namespace Eigen;

// Regresión lineal
vector<double> linearRegression(const VectorXd& x, const VectorXd& y) {

    // Create the design matrix (X), which will have two columns: ones and x
    MatrixXd X(x.size(), 2);  // Two columns: one of 1's and one of x values
    X.col(0) = VectorXd::Ones(x.size());  // Column of 1's (for the intercept)
    X.col(1) = x;  // Column of x (independent variable)

    // Estimate the coefficients (beta) using the least squares formula
    // beta = (X^T * X)^-1 * X^T * y
    VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);

    // The first coefficient is the intercept
    double intercept = beta(0);
    // The second coefficient is the slope of the line
    double slope = beta(1);

    //cout << "The equation of the line is: y = " << intercept << " + " << slope << " * x\n";

    vector<double> result;
    result.push_back(slope);
    result.push_back(intercept);

    return result;
}

// Función que construye la matriz de Vandermonde
MatrixXd buildVandermonde(const VectorXd& x, int degree) {
    const int n = x.size();
    MatrixXd V(n, degree + 1);

    for (int i = 0; i < n; ++i) {
        V(i, 0) = 1.0;  
        for (int j = 1; j <= degree; ++j) {
            V(i, j) = V(i, j - 1) * x(i);  
        }
    }

    return V;
}

// Regresión polinomial
VectorXd polynomialRegression(const VectorXd& x, const VectorXd& y, int degree) {
    if (x.size() != y.size() || x.size() == 0) {
        throw std::invalid_argument("Vectors x and y must have the same size and must not be empty.");
    }
    const MatrixXd V = buildVandermonde(x, degree);
    return (V.transpose() * V).ldlt().solve(V.transpose() * y);
}
