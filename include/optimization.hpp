#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>
#include "model_functions.hpp"
#include "equations.hpp"

// Función de logverosimilitud
double logLikelihood(const Eigen::VectorXd& alpha, const Eigen::VectorXd& epsilon);

// Acepta
bool accept(double currentVal, double candidateVal, double T);

// Genera un vecino dentro de los límites
Eigen::VectorXd generateNeighbor(const Eigen::VectorXd& current, const Eigen::MatrixXd& bounds, double stepScale = 0.1);

// Genera un valor aleatorio entre dos valores
double randDouble(double a, double b);

// Búsqueda Local
Result localSearch(
    const Eigen::VectorXd& alpha,
    const Eigen::VectorXd& initial,
    const Eigen::MatrixXd& bounds,
    int maxIter,
    double baseStep,
    double minStep = 1e-8);

// Simulated Annealing
Result simulatedAnnealing(
    const Eigen::VectorXd& alpha,
    const Eigen::VectorXd& epsilon_init,
    const Eigen::MatrixXd& bounds,
    double T0 = 1.0,
    double gamma = 0.95,
    int L = 20,
    int maxIter = 1000,
    double T_min = 1e-4,
    double stepSize = 0.01,
    int restarts = 3);

// ILS
Result iteratedLocalSearch(
    const Eigen::VectorXd& alpha,
    const Eigen::VectorXd& initial,
    const Eigen::MatrixXd& bounds,
    int ilsMaxIter,
    double perturbationSize,
    int lsMaxIter,
    double lsStepSize);

#endif // OPTIMIZATION_HPP
