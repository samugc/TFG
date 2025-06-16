#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <cmath>
#include "io_handler.hpp"
#include "model_data.hpp"
#include "model_functions.hpp"
#include "visualization.hpp"
#include "regression.hpp"
#include "estimations.hpp"
#include "equations.hpp"
#include "optimization.hpp"
#include "evolutive.hpp"

using namespace std;
using namespace Eigen;

/**************************************************************************************************/
//                             OBTENEMOS EL TIEMPO Y LA VARIABLE X                                //
/**************************************************************************************************/

// Necesitamos añadir una fila y columna de 0
MatrixXd add0(const MatrixXd& matrix){
    
    int rows = matrix.rows();
    int cols = matrix.cols();

    // Crear una nueva matriz con una fila y una columna extra (inicializadas en cero)
    MatrixXd extended = MatrixXd::Zero(rows + 1, cols + 1);

    // Copiar la matriz original dentro de la nueva (desplazada una fila y una columna)
    extended.block(1, 1, rows, cols) = matrix;

    return extended;
}

// Obtenemos la matriz del tiempo
MatrixXd obtainTime(const MatrixXd& transposedMatrix){

    // Obtener dimensiones de la matriz transpuesta
    int filas = transposedMatrix.rows();

    // Extraemos los valores del tiempo
    RowVectorXd t_values = transposedMatrix.row(0);

    // Extendemos la fila de tiempos el número de trayectorias que haya
    MatrixXd t_matrix = t_values.replicate(filas, 1);
    
    // We add a row and column of 0
    MatrixXd t_extended = add0(t_matrix);

    return t_extended;
}

// Obtenemos la matriz de la variable aleatoria X
MatrixXd obtainX(const MatrixXd& transposedMatrix){

    // Eliminamos la primera fila de la matriz transpuesta
    MatrixXd x_matrix = transposedMatrix.bottomRows(transposedMatrix.rows() - 1);

    // We add a row and column of 0
    MatrixXd x_extended = add0(x_matrix);

    return x_extended;
}

/**************************************************************************************************/
//                                      VARIABLE ALEATORIA V                                      //
/**************************************************************************************************/

// Obtenemos la matriz de la variable aleatoria V
MatrixXd obtainV(){
    // Crear la matriz `v` de tamaño (d, n)
    MatrixXd v = MatrixXd::Zero(d + 1, n(1));

    // Llenar la primera fila de `v`
    for (int i = 1; i <= d; i++) {
        v(0, i) = x(i, 1);
    }

    // Llenar el resto de `v`
    for (int i = 1; i <= d; i++) {
        for (int j = 1; j <= (n(i) - 1); j++) {  
            v(i, j) = (1.0/pow(delta<double>(i, j + 1, j), 0.5)) * log(x(i, j + 1) / x(i, j));
        }
    }

    return v;
}

/**************************************************************************************************/
//                                           LÍMITES                                              //
/**************************************************************************************************/

// Calculamos los límites
void computeBoundsFromData(
    int degree,
    MatrixXd& beta_bounds,   
    Vector2d& eta_bounds,    
    Vector2d& sigma2_bounds, 
    const VectorXd& estimationEtaBeta
) {
    VectorXd m_arith = arithmeticMean();

    double m1 = m_arith(1);
    double mN = m_arith(N);
    double eta_hat = 1.0 / (mN / m1 - 1.0);

    // Construir los vectores para regresión polinómica
    VectorXd y(N-2);
    y(0)=0;
    for (int j = 1; j <= N-3; ++j) {
        double ratio = (mN / m_arith(j) - 1.0) * eta_hat;
        y(j) = -log(ratio);
    }

    VectorXd newY = removeFirst(y);

    VectorXd newT2 = removeLast(t_values);
    VectorXd newT1 = removeLast(newT2);
    VectorXd newT = removeLast(newT1);

    // Regresión polinómica
    VectorXd beta_est(degree);
    beta_est = polynomialRegression(newT, newY, degree);

    beta_bounds = MatrixXd(2, p);
    for (int j = 0; j < p; ++j) {
        double val = estimationEtaBeta[j + 1]; // salto η
        beta_bounds(0, j) = val * 0.9;
        beta_bounds(1, j) = val * 1.1;
    }

    double eta_min = std::numeric_limits<double>::infinity();
    double eta_max = 0.0;
    for (int i = 1; i <= d; ++i) {
        double ratio = x(i, N) / x(i, 1);
        if (ratio > 1.0) {
            double inv = 1.0 / (ratio - 1.0);
            eta_min = std::min(eta_min, inv);
            eta_max = std::max(eta_max, inv);
        }
    }
    eta_bounds << eta_min, eta_max;

    sigma2_bounds << 0.0, 0.01;
}

/**************************************************************************************************/
//                                           ERRORES                                              //
/**************************************************************************************************/






/**************************************************************************************************/
//                                            MAIN                                                //
/**************************************************************************************************/

int main(){

    // Leer CSV
    vector<vector<double>> data = readCSV("../data/TrayectoriasSimuladas_Multisigmoidal_Articulo.csv");
    
    // Convertir a Eigen matriz
    MatrixXd eigenMatrix = convertToMatrix(data);
    MatrixXd transposedMatrix = eigenMatrix.transpose();

    // Definir d y n
    d = transposedMatrix.rows()-1;
    int col = transposedMatrix.cols();
    n = createEvaluatedVector(d,col);
    N = n[1];
    
    // Calcular la matriz tiempo y de la variable aleatoria X
    t_values = transposedMatrix.row(0).transpose();
    t = obtainTime(transposedMatrix);
    x = obtainX(transposedMatrix);

    // Obtenemos V 
    v = obtainV();

    /// Coeficientes de beta
    p = 3;

    /// Vector de valores iniciales
    VectorXd estimationEtaBeta(p);
    estimationEtaBeta = estimateEtaBeta(p);

    double val_mu1 = obtainMu0();
    double val_sigma12 = obtainSigma0();

    VectorXd alpha(2);
    alpha << val_mu1, val_sigma12; 
    
    int tam_var = p+2;

    double val_eta = exp(estimationEtaBeta[0]);
    double val_beta1 = estimationEtaBeta[1];
    double val_beta2 = estimationEtaBeta[2];
    double val_beta3 = estimationEtaBeta[3];
    double val_sigma = 0.01;

    VectorXd epsilon(tam_var);
    epsilon << val_eta, val_beta1, val_beta2, val_beta3, val_sigma;  

    // Ejecutamos Newton-Raphson
    //VectorXd solution = newton_raphson(epsilon);
    //cout << solution << endl;

    // Calculamos los límites
    MatrixXd beta_bounds;
    Vector2d eta_bounds, sigma2_bounds;

    computeBoundsFromData(p, beta_bounds, eta_bounds, sigma2_bounds, estimationEtaBeta);

    MatrixXd bounds(2, tam_var); // 2 filas: min y max

    // Insertar límites: orden = [η, β1, ..., βp, σ²]
    bounds(0, 0) = eta_bounds(0); // mínimo η
    bounds(1, 0) = eta_bounds(1); // máximo η

    for (int j = 0; j < beta_bounds.cols(); ++j) {
        bounds(0, j + 1) = beta_bounds(0, j); // mínimos β
        bounds(1, j + 1) = beta_bounds(1, j); // máximos β
    }

    bounds(0, 1) = beta_bounds(0, 0); // mínimos β
    bounds(1, 1) = beta_bounds(1, 0); // máximos β

    bounds(0, 2) = beta_bounds(1, 1); // mínimos β
    bounds(1, 2) = beta_bounds(0, 1); // máximos β

    bounds(0, 3) = beta_bounds(0, 2); // mínimos β
    bounds(1, 3) = beta_bounds(1, 2); // máximos β

    bounds(0, tam_var - 1) = sigma2_bounds(0); // mínimo σ²
    bounds(1, tam_var - 1) = sigma2_bounds(1); // máximo σ²

    // Valores reales
    double real_eta = exp(-1);
    double real_beta1 = 0.1;
    double real_beta2 = -0.009;
    double real_beta3 = 0.0002;
    double real_sigma = 0.01*0.01;

    VectorXd epsilon_real(tam_var);
    epsilon_real << real_eta, real_beta1, real_beta2, real_beta3, real_sigma; 

    // Estimaciones de las trayectorias

    // NR
    double eta_nr   = 0.36312;
    double b1_nr    = 0.104492;
    double b2_nr    = -0.00928902;
    double b3_nr    = 0.000204749;

    // Búsqueda Local (BL)
    double eta_bl = 0.384842;
    double b1_bl = 0.107690;
    double b2_bl = -0.009599;
    double b3_bl = 0.000210860;

    // Simulated Annealing (SA)
    double eta_sa = 0.374303;
    double b1_sa = 0.109739;
    double b2_sa = -0.009670;
    double b3_sa = 0.00021118;

    // Iterated Local Search (ILS)
    double eta_ils = 0.415067;
    double b1_ils = 0.109651;
    double b2_ils = -0.009780;
    double b3_ils = 0.000215;

    // Evolutionary Algorithm sin elitismo (EA)
    double eta_ea = 0.374502;
    double b1_ea = 0.103642;
    double b2_ea = -0.009623;
    double b3_ea = 0.000200;

    // Evolutionary Algorithm con elitismo (EA+E)
    double eta_eae = 0.370183;
    double b1_eae = 0.099506;
    double b2_eae = -0.008949;
    double b3_eae = 0.000199;

    // Memético (MA)
    double eta_mae = 0.369941;
    double b1_mae = 0.099469;
    double b2_mae = -0.008973;
    double b3_mae = 0.000199;

    // Differential Evolution (DE)
    double eta_de = 0.368640;
    double b1_de = 0.098817;
    double b2_de = -0.008912;
    double b3_de = 0.000198;

    // Curvas estimadas por cada algoritmo
    std::vector<std::pair<std::string, VectorXd>> curvas = {
        {"Teórica", m_teorica},
        {"Muestral", m_sample},
        {"BL", media(eta_bl, b1_bl, b2_bl, b3_bl)},
        {"SA", media(eta_sa, b1_sa, b2_sa, b3_sa)},
        {"ILS", media(eta_ils, b1_ils, b2_ils, b3_ils)}
        {"EA", media(eta_ea, b1_ea, b2_ea, b3_ea)},
        {"EA+Elitismo", media(eta_eae, b1_eae, b2_eae, b3_eae)},
        {"MA", media(eta_mae, b1_mae, b2_mae, b3_mae)},
        {"DE", media(eta_de, b1_de, b2_de, b3_de)}
    };

    // Visualización
    //graphMultipleCurvesNR(t_values, curvas, "Comparación de medias estimadas por algoritmo");


    // Ejecución del algoritmo
    
    /*
    cout<<"NEWTON-RAPHSON"<<endl;
    
    ofstream outfile("resultados_newtonRaphson.txt");
    if (!outfile.is_open()) {
        cerr << "Error al abrir el archivo de salida." << endl;
        return 1;
    }

    for (int i = 0; i < 30; ++i) {
        Result result = newton_raphsonR(epsilon);
        outfile << "Ejecución " << i + 1 << ":\n"
                << "Solución: " << result.solucion.transpose() << "\n"
                << "Fitness: " << result.fitness << "\n"
                << "Tiempo: " << result.time << " segundos\n\n"
                << std::flush;
    }
    */
    return 0;
}