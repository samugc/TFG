#ifndef IO_HANDLER_HPP
#define IO_HANDLER_HPP

#include <vector>
#include <string>
#include <Eigen/Dense>

// Lee un archivo CSV con separador ';' y coma decimal.
// Omite cabecera y convierte a vector de vectores.
std::vector<std::vector<double>> readCSV(const std::string& filename);

// Crea un vector de cierto tamaño con un valor
Eigen::VectorXd createEvaluatedVector(int tam, double value);

// Elimina el primer elemento de un VectorXd 
Eigen::VectorXd removeFirst(const Eigen::VectorXd& vec);

// Elimina el último elemento de un VectorXd
Eigen::VectorXd removeLast(const Eigen::VectorXd& vec);

// Convierte un vector de vectores en una MatrixXd.
Eigen::MatrixXd convertToMatrix(const std::vector<std::vector<double>>& data);

// Imprime una matriz Eigen por consola.
void printMatrix(const Eigen::MatrixXd& matrix);

// Imprime un vector Eigen por consola.
void printVector(const Eigen::VectorXd& vec);

// Imprime las filas de la matriz Eigen
void printRows(const Eigen::MatrixXd& matrix, int numRows);

// Imprime el vector de vectores
void printVectorRange(const std::vector<std::vector<double>>& data, size_t endRow, size_t endCol);

// Imprime un polinomio
void printPoly(const Eigen::VectorXd& poly);

#endif // IO_HANDLER_HPP