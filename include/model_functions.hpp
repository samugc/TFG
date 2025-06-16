#ifndef MODEL_FUNCTIONS_HPP
#define MODEL_FUNCTIONS_HPP

#include "model_data.hpp"
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <boost/math/tools/polynomial.hpp>
#include <ceres/ceres.h>

// Funciones del modelo

template <typename T>
T evaluate_polynomial(const T* theta, T x) {
    T result = T(0);
    for (int i = 1; i < p+1; ++i) {
        result += theta[i] * ceres::pow(x, i); 
    }
    return result;
}

template <typename T>
T nTotal() {
    int nval = 0;
    for (int i = 1; i <= d; ++i) {
        int ni = n(i);  // Asegúrate de que n(i) devuelve un valor no cero
        nval += (ni - 1);
    }
    return T(nval);
}

template <typename T>
T delta(int i, int j, int k) {
    return T(t(i, j)) - T(t(i, k));
}

template <typename T>
T Z1() {
    T val = T(0);
    for (int i = 1; i <= d; ++i) {
        for (int j = 1; j <= (n(i) - 1); ++j) {
            val += T(pow(v(i, j), 2));  
        }
    }
    return val;
}

template <typename T>
T lambda(const T* theta, int i, int m, int n) {
    T eta = theta[0];
    T Qn = evaluate_polynomial(theta, T(t(i, n)));
    T Qm = evaluate_polynomial(theta, T(t(i, m)));

    // Usar ceres::exp() sobre las expresiones
    T exp_Qn = ceres::exp(-Qn);
    T exp_Qm = ceres::exp(-Qm);

    // Calcular el logaritmo con las exponenciales
    T num = eta + exp_Qn;
    T den = eta + exp_Qm;

    // Logaritmo de la fracción
    T result = ceres::log(num / den);
    return result;
}

double me(const Eigen::VectorXd& epsilon,int i,int k,int j); 
double phi(const Eigen::VectorXd& epsilon);
double varGamma(const Eigen::VectorXd& epsilon);

template <typename T>
T deltaBarra(T l){
    T val = T(1);
    if(l==T(0)){val = T(0);}
    return val;
}

template <typename T>
T lDelta(T l, const T* theta, int i, int m) {
    T eta = theta[0];
    T Qm = evaluate_polynomial(theta, T(t(i, m)));
    T den = eta + ceres::exp(-Qm);

    T term1;
    if(l==0 && t(i,m)==0){
        term1 = T(-1);
    }
    else{
        term1 = -ceres::pow(t(i, m), l);
    }

    return  (T(1) / den) * term1;
}

template <typename T>
T lD(T l, const T* theta, int i, int m, int n) {
    T Qm = evaluate_polynomial(theta, T(t(i, m)));
    T Qn = evaluate_polynomial(theta, T(t(i, n)));

    T delta = deltaBarra(l);
    T exp_m = ceres::exp(-Qm);
    T exp_n = ceres::exp(-Qn);

    T term1 = lDelta<T>(l, theta, i, m) * ceres::pow(exp_m, delta);
    T term2 = lDelta<T>(l, theta, i, n) * ceres::pow(exp_n, delta);

    return term1 - term2;
}

template <typename T>
T W(T l, const T* theta) {
    T val = T(0);
    for (int i = 1; i <= d; ++i){
        for(int j = 1;j <= (n(i)-1); ++j){
            val += lD<T>(l, theta, i, j+1, j);
        }
    }
    return val;
}

template <typename T>
T Y(T l, const T* theta) {
    T eta = theta[0];
    T val = T(0);

    for (int i = 1; i <= d; ++i) {
        for (int j = 1; j <= (n(i) - 1); ++j) {
            T term1 = T(1) / delta<T>(i, j + 1, j);

            T Qj1 = evaluate_polynomial(theta, T(t(i,j+1)));
            T Qj = evaluate_polynomial(theta, T(t(i,j)));

            T num = eta + ceres::exp(-Qj1);
            T den = eta + ceres::exp(-Qj);

            T term2 = ceres::log(num / den);
            T term3 = lD<T>(l, theta, i, j + 1, j);

            val += term1 * term2 * term3;
        }
    }
    return val;
}

template <typename T>
T X(T l, const T* theta) {
    T val = T(0);
    for (int i = 1; i <= d; ++i) {
        for (int j = 1; j <= (n(i) - 1); ++j) {
            T num = T(v(i, j));
            T den = ceres::pow(delta<T>(i, j + 1, j), (0.5));
            T term1 = num / den;
            val += term1 * lD<T>(l, theta, i, j + 1, j);
        }
    }
    return val;
}

template <typename T>
T Z3() {
    T val = T(0);
    for (int i = 1; i <= d; ++i) {
        val += delta<T>(i, n(i), 1);
    }   
    return val;
}

template <typename T>
T A(const T* theta) {
    T val = T(0);
    for (int i = 1; i <= d; ++i) {
        for (int j = 1; j <= (n(i) - 1); ++j) {
            T num = ceres::pow(lambda<T>(theta, i, j + 1, j), T(2));
            T den = delta<T>(i, j + 1, j);
            val += num/den;
        }
    }
    return val;
}

template <typename T>
T B(const T* theta) {
    T val = T(0);
    for (int i = 1; i <= d; ++i) {
        for (int j = 1; j <= (n(i) - 1); ++j) {
            T num = T(v(i, j)) * lambda<T>(theta, i, j + 1, j);  
            T den = ceres::pow(delta<T>(i, j + 1, j), 0.5);
            val += num / den;
        }
    }
    return val;
}

#endif  // MODEL_FUNCTIONS_HPP
