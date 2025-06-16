#ifndef MODEL_DATA_HPP
#define MODEL_DATA_HPP

#include <Eigen/Dense>

// Dimensiones
extern int d;
extern Eigen::VectorXd n;
extern int N;

// Matriz de tiempo y caminos
extern Eigen::VectorXd t_values;
extern Eigen::MatrixXd t;
extern Eigen::MatrixXd x;

// Variable aleatoria V
extern Eigen::MatrixXd v;

// Coeficientes de beta
extern int p;

#endif // MODEL_DATA_HPP