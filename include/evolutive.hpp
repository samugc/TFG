#ifndef EVOLUTIVE_HPP
#define EVOLUTIVE_HPP

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>
#include "optimization.hpp" 

struct Individual {
    Eigen::VectorXd genes;
    double fitness;

    bool operator<(const Individual& other) const {
        return fitness > other.fitness; 
    }
};

Individual randomIndividual(const Eigen::MatrixXd& bounds, const Eigen::VectorXd& alpha);
Individual crossover(const Individual& p1, const Individual& p2, const Eigen::MatrixXd& bounds, const Eigen::VectorXd& alpha, double alpha_blx = 0.3);
void mutate(Individual& ind, const Eigen::MatrixXd& bounds, double mutationRate, double mutationStep);
Individual tournament(const std::vector<Individual>& pop);

Result evolutionaryAlgorithm(const Eigen::VectorXd& alpha,
                             const Eigen::MatrixXd& bounds,
                             int popSize,
                             int generations,
                             double crossoverRate,
                             double mutationRate,
                             double mutationStep);

Result evolutionaryAlgorithmWithElitism(const Eigen::VectorXd& alpha,
                                        const Eigen::MatrixXd& bounds,
                                        int popSize,
                                        int generations,
                                        double crossoverRate,
                                        double mutationRate,
                                        double mutationStep);

Result memeticAlgorithmWithElitism(const Eigen::VectorXd& alpha,
                                   const Eigen::MatrixXd& bounds,
                                   int popSize,
                                   int generations,
                                   double crossoverRate,
                                   double mutationRate,
                                   double mutationStep,
                                   int localSearchInterval,
                                   int localSearchSteps);

Result differentialEvolution(const Eigen::VectorXd& alpha,
                              const Eigen::MatrixXd& bounds,
                              int popSize,
                              int generations,
                              double F = 0.8,
                              double CR = 0.9);

Eigen::VectorXd adaptiveMemeticWithSmartRestarts(const Eigen::VectorXd& alpha,
                                                 const Eigen::MatrixXd& bounds,
                                                 int popSize,
                                                 int generations,
                                                 double crossoverRate,
                                                 double mutationRate,
                                                 double mutationStep,
                                                 int localSearchInterval,
                                                 int localSearchSteps,
                                                 int maxNoImprovement);

#endif // EVOLUTIVE_HPP
