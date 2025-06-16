#include "equations.hpp"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

VectorXd compute_residuals(const VectorXd& vars) {
    const int n_residuals = 5;  // suponiendo p=3 → 1 + (p+1) = 5 ecuaciones
    double theta[4] = { vars[0], vars[1], vars[2], vars[3] };
    double sigma_squared = vars[4];
    double* parameters[] = { theta, &sigma_squared };

    // AutoDiff para obtener residuos
    ceres::AutoDiffCostFunction<Equation1Functor, 5, 4, 1> cost(
        new Equation1Functor());

    double residuals[5];
    cost.Evaluate(parameters, residuals, nullptr);  // nullptr: no calculamos Jacobiano aquí

    VectorXd res_vec(n_residuals);
    for (int i = 0; i < n_residuals; ++i) {
        res_vec[i] = residuals[i];
    }
    return res_vec;
}

MatrixXd compute_jacobian(const VectorXd& vars) {

    double theta4[4] = { vars(0), vars(1), vars(2), vars(3) };
    const double* theta = theta4;
    double sigma_squared = vars(4);

    // Crear punteros de parámetros
    double* parameters[] = { const_cast<double*>(theta), &sigma_squared };

    // Cost function con AutoDiff (5 ecuaciones, 4 parámetros + 1 sigma²)
    ceres::AutoDiffCostFunction<Equation1Functor, 5, 4, 1> auto_cost(
        new Equation1Functor());

    // Almacenes para residuos y derivadas
    double residuals_auto[5];
    double jac_auto[5 * 4 + 5 * 1];  // 5 filas × 5 columnas
    double* jac_auto_ptrs[2] = { jac_auto, jac_auto + 5 * 4 };

    // Evaluar AutoDiff
    auto_cost.Evaluate(parameters, residuals_auto, jac_auto_ptrs);

    // Construir matriz Jacobiana (5x5)
    MatrixXd J(5, 5);

    // Copiar derivadas respecto a theta
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 5; ++i) {
            J(i, j) = jac_auto[j * 5 + i];  // acceso columna-major
        }
    }

    // Copiar derivadas respecto a sigma²
    for (int i = 0; i < 5; ++i) {
        J(i, 4) = jac_auto[4 * 5 + i];
    }

    return J;
}

VectorXd newton_raphson(const VectorXd& vars_iniciales, int max_iter, double tol, double alpha, double beta) {
    
    auto start_time = high_resolution_clock::now();

    VectorXd vars = vars_iniciales;
    int n = vars.size();
    for (int iter = 0; iter < max_iter; ++iter) {
        VectorXd f = compute_residuals(vars);

        if (f.norm() < tol) {
            std::cout << "Solución encontrada en " << iter << " iteraciones.\n";
            return vars;
        }
        
        MatrixXd J = compute_jacobian(vars) ;
        VectorXd delta = J.colPivHouseholderQr().solve(f);

        // Backtracking para asegurar que las restricciones se cumplan
        double step = alpha;
        VectorXd new_vars;
        while (true) {
            new_vars = vars - step * delta;

            // Comprobar si las restricciones se cumplen
            if (new_vars[0] >= 0.0 && new_vars[n - 2] >= 0.0 && new_vars[n - 1] >= 0.0) {
                break;  
            }
            
            // Reducimos el paso
            step *= beta;

            // Si el paso es demasiado pequeño, detener iteración
            if (step < 1e-15) {
                std::cout << "No se puede avanzar sin violar restricciones.\n";
                auto end_time = high_resolution_clock::now();
                duration<double> elapsed = end_time - start_time;
                std::cout << "Tiempo total de ejecución: " << elapsed.count() << " segundos.\n";
                return vars;
            }
        }

        // Actualizar valores
        vars = new_vars;
        cout << "Iteración " << iter << endl;
        cout << vars << endl;
    }

    cout << "Número máximo de iteraciones alcanzado.\n";

    auto end_time = high_resolution_clock::now();
    duration<double> elapsed = end_time - start_time;
    std::cout << "Tiempo total de ejecución: " << elapsed.count() << " segundos.\n";

    return vars;
}

Result newton_raphsonR(const VectorXd& vars_iniciales, int max_iter, double tol, double alpha, double beta) {

    using namespace std::chrono;
    auto start_time = high_resolution_clock::now();

    VectorXd vars = vars_iniciales;
    int n = vars.size();
    double final_fitness = 0.0;

    for (int iter = 0; iter < max_iter; ++iter) {
        VectorXd f = compute_residuals(vars);
        final_fitness = -f.norm();  // fitness como -norma de residuos (mejor = más cercano a 0)

        if (f.norm() < tol) {
            std::cout << "Solución encontrada en " << iter << " iteraciones.\n";
            auto end_time = high_resolution_clock::now();
            double elapsed = duration<double>(end_time - start_time).count();
            return {vars, final_fitness, elapsed};
        }

        MatrixXd J = compute_jacobian(vars);
        VectorXd delta = J.colPivHouseholderQr().solve(f);

        // Backtracking para asegurar que las restricciones se cumplan
        double step = alpha;
        VectorXd new_vars;
        while (true) {
            new_vars = vars - step * delta;

            if (new_vars[0] >= 0.0 && new_vars[n - 2] >= 0.0 && new_vars[n - 1] >= 0.0) {
                break;
            }

            step *= beta;

            if (step < 1e-15) {
                std::cout << "No se puede avanzar sin violar restricciones.\n";
                auto end_time = high_resolution_clock::now();
                double elapsed = duration<double>(end_time - start_time).count();
                return {vars, final_fitness, elapsed};
            }
        }

        vars = new_vars;
        std::cout << "Iteración " << iter << "\n" << vars.transpose() << "\n";
    }

    std::cout << "Número máximo de iteraciones alcanzado.\n";
    auto end_time = high_resolution_clock::now();
    double elapsed = duration<double>(end_time - start_time).count();

    return {vars, final_fitness, elapsed};
}
