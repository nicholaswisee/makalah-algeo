#ifndef MARKOVCHAINANALYZER_H
#define MARKOVCHAINANALYZER_H

#include <vector>

class MarkovChainAnalysis {
   private:
    int numStates;
    double lambda;
    double miu;
    std::vector<std::vector<double>> generatorMatrix;

    void buildGeneratorMatrix();

    void gaussianElimination(std::vector<std::vector<double>>& A, double tolerance);
    std::vector<double> backSubstitution(const std::vector<std::vector<double>>& A);
    void normalizeDistribution(std::vector<double>& pi);

   public:
    MarkovChainAnalysis(int states, double arrivalRate, double serviceRate);

    std::vector<double> computeSteadyState(double tolerance = 1e-10);

    std::vector<double> getTheoreticalSteadyState() const;

    void printGeneratorMatrix() const;

    const std::vector<std::vector<double>>& getGeneratorMatrix() const {
        return generatorMatrix;
    }
};

#endif
