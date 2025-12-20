#ifndef EXPERIMENTRUNNER_H
#define EXPERIMENTRUNNER_H

#include <string>
#include <vector>

struct ExperimentResult {
    double lambda;
    double rho;
    double avgQueueLengthSim;
    double avgQueueLengthTheory;
    double relativeError;
};

class ExperimentRunner {
   private:
    double miu;
    int maxStates;
    double simTime;
    unsigned int seed;

    std::vector<ExperimentResult> results;

   public:
    // Constructor
    ExperimentRunner(double serviceRate, int maxStates, double simTime, unsigned int seed);

    void runVaryingLambdaExperiments();

    ExperimentResult runSingleExperiment(double lambda);

    // Output results
    void exportResults(const std::string& filename) const;
    void printResults() const;

    // Getters
    const std::vector<ExperimentResult>& getResults() const {
        return results;
    }
};

#endif
