#ifndef MARKOVCHAINANALYZER_H
#define MARKOVCHAINANALYZER_H

#include <vector>

class MarkovChainAnalysis {
   private:
    int numStates;
    double lambda;
    double miu;
    std::vector<std::vector<double>> generatorMatrix;  // Q-matrix for CTMC

    // Build the generator matrix (Q-matrix) for M/M/1 queue CTMC
    void buildGeneratorMatrix();

    // Linear algebra helper functions
    void gaussianElimination(std::vector<std::vector<double>>& A, double tolerance);
    std::vector<double> backSubstitution(const std::vector<std::vector<double>>& A);
    void normalizeDistribution(std::vector<double>& pi);

   public:
    // constructor
    MarkovChainAnalysis(int states, double arrivalRate, double serviceRate);

    // Steady-state distribution by solving πQ = 0
    std::vector<double> computeSteadyState(double tolerance = 1e-10);

    // Get theoretical steady-state (closed-form for M/M/1)
    std::vector<double> getTheoreticalSteadyState() const;

    void printGeneratorMatrix() const;

    // Getters
    const std::vector<std::vector<double>>& getGeneratorMatrix() const {
        return generatorMatrix;
    }
};

#endif
