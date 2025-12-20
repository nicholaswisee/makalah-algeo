#include "ExperimentRunner.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include "QueueSimulator.hpp"

ExperimentRunner::ExperimentRunner(double serviceRate, int maxStates, double simTime,
                                   unsigned int seed)
    : miu(serviceRate), maxStates(maxStates), simTime(simTime), seed(seed) {}

void ExperimentRunner::runVaryingLambdaExperiments() {
    std::vector<double> lambdaValues = {0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99};

    results.clear();

    for (double lambda : lambdaValues) {
        std::cout << "Running experiment with λ=" << lambda << "..." << std::endl;
        ExperimentResult result = runSingleExperiment(lambda);
        results.push_back(result);

        std::cout << "  ρ=" << result.rho << ", L_sim=" << result.avgQueueLengthSim
                  << ", L_theory=" << result.avgQueueLengthTheory
                  << ", Error=" << result.relativeError << "%" << std::endl;
    }
}

ExperimentResult ExperimentRunner::runSingleExperiment(double lambda) {
    QueueSimulator sim(lambda, miu, maxStates, seed);
    sim.simulate(simTime);

    ExperimentResult result;
    result.lambda = lambda;
    result.rho = lambda / miu;
    result.avgQueueLengthSim = sim.getAverageQueueLength();
    result.avgQueueLengthTheory = sim.getTheoreticalAverageQueueLength();

    if (result.avgQueueLengthTheory != std::numeric_limits<double>::infinity()) {
        result.relativeError = std::abs(result.avgQueueLengthSim - result.avgQueueLengthTheory) /
                               result.avgQueueLengthTheory * 100.0;
    } else {
        result.relativeError = -1.0;  // unstable system
    }

    return result;
}

void ExperimentRunner::printResults() const {
    std::cout << "\n=== Experiment Results Summary ===" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "λ      | ρ      | L_sim    | L_theory | Error (%)" << std::endl;
    std::cout << "-------|--------|----------|----------|----------" << std::endl;

    for (const auto& result : results) {
        std::cout << std::setw(6) << result.lambda << " | " << std::setw(6) << result.rho << " | "
                  << std::setw(8) << result.avgQueueLengthSim << " | " << std::setw(8)
                  << result.avgQueueLengthTheory << " | " << std::setw(8) << result.relativeError
                  << std::endl;
    }
}
