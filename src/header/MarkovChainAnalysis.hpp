#ifndef MARKOVCHAINANALYZER_H
#define MARKOVCHAINANALYZER_H

#include <vector>

class MarkovChainAnalyzer {
 private:
  int numStates;
  double lambda;
  double mu;
  std::vector<std::vector<double>> transitionMatrix;

  // Build the transition matrix for M/M/1 queue
  void buildTransitionMatrix();

 public:
  // constructor
  MarkovChainAnalyzer(int states, double arrivalRate, double serviceRate);

  // Steady-state distribution using power iteration
  std::vector<double> computeSteadyState(int maxIterations = 10000, double tolerance = 1e-10);

  void printTransitionMatrix() const;

  // Getters
  const std::vector<std::vector<double>>& getTransitionMatrix() const {
    return transitionMatrix;
  }
};

#endif
