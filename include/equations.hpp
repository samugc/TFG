#ifndef EQUATIONS_HPP
#define EQUATIONS_HPP

#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <iostream>
#include <chrono>
#include "model_functions.hpp"

// Functor que define el sistema de ecuaciones
struct Equation1Functor {
    template <typename T>
    bool operator()(const T* const theta, const T* const sigma_squared, T* residuals) const {
        residuals[0] = sigma_squared[0] * (nTotal<T>() + (sigma_squared[0] / T(4)) * Z3<T>()) - Z1<T>() - A<T>(theta) + T(2) * B<T>(theta);
        for(int l=0;l<p+1;++l){
            residuals[l+1] = Y<T>(T(l),theta) + (sigma_squared[0] / T(2)) * W<T>(T(l),theta) + X<T>(T(l),theta); 
        }
        return true;
    }
};

// Función que calcula el vector de residuos del sistema
Eigen::VectorXd compute_residuals(const Eigen::VectorXd& vars);

// Función que calcula la matriz Jacobiana del sistema
Eigen::MatrixXd compute_jacobian(const Eigen::VectorXd& vars);

// Método de Newton con Restricciones usando un Paso de Reducción de Paso 
Eigen::VectorXd newton_raphson(const Eigen::VectorXd& vars_iniciales, int max_iter = 100, double tol = 1e-6, double alpha = 1.0, double beta = 0.5);

struct Result {
    Eigen::VectorXd solucion;
    double fitness;
    double time; 
};

Result newton_raphsonR(const Eigen::VectorXd& vars_iniciales, int max_iter = 100, double tol = 1e-6, double alpha = 1.0, double beta = 0.5);

#endif // EQUATIONS_HPP