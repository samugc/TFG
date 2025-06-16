#include "model_functions.hpp"

using namespace std;

double me(const Eigen::VectorXd& epsilon,int i,int k,int j){

    double eta = epsilon(0);
    double beta1 = epsilon(1);
    double beta2 = epsilon(2);
    double beta3 = epsilon(3);
    double sigma2 = epsilon(4);

    vector<double> theta = {0.0, beta1, beta2, beta3};

    double Qn = evaluate_polynomial<double>(theta.data(), t(i, j));
    double Qd = evaluate_polynomial<double>(theta.data(), t(i, k));

    double num = eta + exp(-Qn);
    double den = eta + exp(-Qd);

    return log(num / den) - (sigma2/2) * delta<double>(i,k,j);
}

double phi(const Eigen::VectorXd& epsilon){
    double val = 0.0;
    for (int i = 1; i <= d; ++i){
        for(int j = 1;j <= (n(i)-1); ++j){
            double num = pow(me(epsilon,i,j+1,j),2);
            double den = delta<double>(i,j+1,j);
            val += num/den;
        }
    }
    return val;
}

double varGamma(const Eigen::VectorXd& epsilon){
    double val = 0;
    for (int i = 1; i <= d; ++i){
        for(int j = 1;j <= (n(i)-1); ++j){
            double num = v(i,j)*me(epsilon,i,j+1,j);
            double den = pow(delta<double>(i,j+1,j),0.5);
            val += num/den;
        }
    }
    return val;
}