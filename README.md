# Estimación de parámetros en un modelo estocástico multisigmoidal

Este repositorio contiene el código fuente desarrollado para el Trabajo Fin de Grado (TFG) titulado **“Estimación de parámetros en un modelo de difusión basado en múltiples funciones sigmoides”**, presentado en la Universidad de Granada.

## Resumen

Se ha diseñado e implementado un conjunto de algoritmos de optimización para ajustar parámetros de un modelo estocástico a partir de datos simulados. El enfoque combina técnicas clásicas (Newton-Raphson) con metaheurísticas (Simulated Annealing, ILS, DE, EA y CMA-ES). Además, se proporciona un análisis estadístico de los resultados obtenidos.

## Características

- Implementación en C++ usando Eigen y Ceres.
- Módulos separados para modelo, optimización y visualización.
- Análisis de sensibilidad y comparación entre algoritmos.
- Tests no paramétricos (Kruskal-Wallis, Wilcoxon) aplicados con Python.

## Estructura del repositorio
  
## Requisitos

- C++17 o superior
- [Eigen](https://eigen.tuxfamily.org)

## Cómo usar

### Clonar el repositorio
git clone https://github.com/samugc/TFG_codigo.git \\
cd TFG_codigo

### Compilar (usando CMake, por ejemplo)
mkdir build && cd build  \\
cmake .. \\
make

### Ejecutar
./main

## 📊 Análisis estadístico

Para determinar si existen diferencias estadísticamente significativas entre los algoritmos de estimación utilizados, se ha realizado un análisis no paramétrico sobre los resultados obtenidos tras múltiples ejecuciones.

Concretamente, se han aplicado los siguientes tests:

- **Test de Kruskal-Wallis**, para comparar más de dos algoritmos en cuanto a precisión y error de estimación por parámetro.
- **Test de Wilcoxon con corrección de Holm**, para comparaciones por pares entre algoritmos.

## 👨‍🎓 Autor

**Samuel García**  
Doble Grado en Ingeniería Informática y Matemáticas  
Universidad de Granada  
Trabajo Fin de Grado — Curso 2024/2025

## 📁 Repositorio

Este repositorio contiene todo el código necesario para la reproducción de los experimentos y análisis del TFG:

👉 [https://github.com/samugc/TFG_codigo](https://github.com/samugc/TFG_codigo)
